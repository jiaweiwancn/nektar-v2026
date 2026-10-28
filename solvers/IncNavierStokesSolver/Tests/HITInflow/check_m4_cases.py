#!/usr/bin/env python3
"""Accept raw M4 native-probe, solver, and field-comparison results on Linux."""

import argparse
import csv
import json
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

VARIABLES = ('u', 'v', 'w', 'p')
SAMPLE_TOLERANCE = 5e-12


def read_json(path):
    return json.loads(path.read_text())


def probe(args):
    result = read_json(args.summary)
    result['kind'] = 'probe'
    result['maximum_projection_error'] = args.projection_limit
    result['passed'] = (result['points'] > 0 and result['repeat_identical']
                        and result['raw_max_error'] <= 5e-13
                        and max(result['projection_max_error']) <= args.projection_limit)
    if args.nonzero:
        result['passed'] &= result['generator_modes'] > 0
        if args.h1d:
            result['passed'] &= min(result['nonzero_z_coeff_energy']) > 1e-12
    return result


def samples(path):
    stem = str(path)[:-len('.probe.json')]
    rows = {}
    files = sorted(path.parent.glob(Path(stem).name + '.probe-rank*.csv'))
    assert len(files) == read_json(path)['ranks'], 'Missing rank probe output'
    for file in files:
        with file.open(newline='') as stream:
            for row in csv.DictReader(stream):
                values = {key: float(value) for key, value in row.items()}
                assert all(math.isfinite(x) for x in values.values())
                coordinates = []
                for axis in 'xyz':
                    value = values[axis]
                    if axis in 'yz' and (abs(value) < 1e-10 or abs(value-2*math.pi) < 1e-10):
                        value = 0.0
                    coordinates.append(round(value, 12))
                key = (round(values['time'], 14), *coordinates)
                data = tuple(values[name] for name in
                             ('u', 'v', 'w', 'projected_u', 'projected_v', 'projected_w'))
                if key in rows:
                    assert max(abs(a-b) for a,b in zip(rows[key],data)) <= SAMPLE_TOLERANCE, 'Inconsistent repeated/shared sample'
                else:
                    rows[key] = data
    assert rows, 'Empty inlet probe'
    return rows


def compare_samples(args):
    first, second = samples(args.first), samples(args.second)
    common = first.keys() & second.keys()
    assert len(common) >= 20, 'Insufficient common inlet coordinates'
    if not args.common:
        assert first.keys() == second.keys(), 'MPI/transform changed the physical sample grid'
    raw = max(abs(first[key][i]-second[key][i]) for key in common for i in range(3))
    projected = max(abs(first[key][i]-second[key][i]) for key in common for i in range(3,6))
    a,b = read_json(args.first),read_json(args.second)
    identity = (a['identity_sha256'] == b['identity_sha256'] and
                a['modes_sha256'] == b['modes_sha256'])
    return {'kind':'samples','common_samples':len(common),
            'sample_counts':[len(first),len(second)], 'raw_max_difference':raw,
            'projected_max_difference':projected, 'identical_realization':identity,
            'projected_limit':args.projected_limit,
            'passed':identity and raw<=SAMPLE_TOLERANCE and projected<=args.projected_limit}


def fields(args):
    text = args.log.read_text()
    assert 'HIT_FIELD_COMPARE' in text, 'Missing native field-comparison marker'
    data,_ = json.JSONDecoder().raw_decode(text.split('HIT_FIELD_COMPARE',1)[1].lstrip())
    assert set(data) == set(VARIABLES)
    assert all(math.isfinite(v) and v>=0 for row in data.values() for v in row.values())
    assert data['u']['first_l2'] > 1 and data['u']['second_l2'] > 1, 'Field data was not imported'
    velocity = max(data[v]['linf'] for v in VARIABLES[:3])
    return {'kind':'fields','errors':data,'maximum_velocity_linf':velocity,
            'velocity_limit':args.velocity_limit,'pressure_limit':args.pressure_limit,
            'passed':velocity<=args.velocity_limit and data['p']['linf']<=args.pressure_limit}


def solver(args):
    root = ET.parse(args.session).getroot()
    parameters = {}
    for item in root.findall('./CONDITIONS/PARAMETERS/P'):
        name,value=item.text.split('=',1)
        parameters[name.strip().lower()]=float(value)
    text=args.log.read_text()
    assert 'Total Computation Time' in text, 'Incomplete solver run'
    identity=re.search(r'HITInflow restart identity: ([a-f0-9]{64}), start time: (\S+)',text)
    assert identity, 'Missing initialized restart identity/time'
    start=float(identity.group(2))
    statuses=re.findall(r'Steps:\s*(\d+)\s+Time:\s*(\S+)',text)
    assert statuses
    assert int(statuses[-1][0])==round(parameters['numsteps'])
    assert math.isclose(float(statuses[-1][1]),start+parameters['numsteps']*parameters['timestep'],abs_tol=1e-12)
    errors={'l2':{},'linf':{}}
    for metric,var,value in re.findall(r'L\s+(2|inf)\s+error\s*\(variable\s+([^)]+)\)\s*:\s*(\S+)',text):
        value=float(value)
        assert math.isfinite(value) and value>=0
        errors['l2' if metric=='2' else 'linf'][var]=value
    assert all(set(row)==set(VARIABLES) for row in errors.values())
    passed=True
    if args.mean_only:
        passed=max(v for row in errors.values() for v in row.values())<=1e-10
    else:
        passed=(math.pi < errors['l2']['u'] < 4*math.pi
                and min(errors['l2'][v] for v in ('v','w'))>1e-12
                and max(v for row in errors.values() for v in row.values())<20)
    return {'kind':'solver','session':str(args.session),'log':str(args.log),
            'start_time':start,'final_time':float(statuses[-1][1]),
            'identity_sha256':identity.group(1),**errors,'passed':passed}


def statistics(args):
    first,second=read_json(args.first),read_json(args.second)
    differences={v:abs(first['l2'][v]-second['l2'][v])/max(first['l2'][v],second['l2'][v],1e-30) for v in 'uvw'}
    maximum=max(differences.values())
    return {'kind':'statistics','relative_velocity_l2_differences':differences,
            'maximum_relative_velocity_l2_difference':maximum,
            'maximum_pressure_l2_difference':abs(first['l2']['p']-second['l2']['p']),
            'relative_velocity_limit':args.relative_limit,
            'passed':maximum<=args.relative_limit}


def refinement(args):
    first,second=read_json(args.first),read_json(args.second)
    coarse,fine=first[args.metric],second[args.metric]
    assert coarse>0, 'Zero coarse difference cannot establish convergence'
    return {'kind':'refinement','metric':args.metric,'coarse':coarse,'fine':fine,
            'ratio':fine/coarse,'maximum_ratio':args.ratio,'passed':fine/coarse<args.ratio}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json',type=Path)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('probe');p.add_argument('summary',type=Path)
    p.add_argument('--projection-limit',type=float,default=2e-3)
    p.add_argument('--nonzero',action='store_true');p.add_argument('--h1d',action='store_true');p.set_defaults(run=probe)
    p=sub.add_parser('samples');p.add_argument('first',type=Path);p.add_argument('second',type=Path)
    p.add_argument('--common',action='store_true');p.add_argument('--projected-limit',type=float,default=5e-12);p.set_defaults(run=compare_samples)
    p=sub.add_parser('fields');p.add_argument('log',type=Path)
    p.add_argument('--velocity-limit',type=float,default=2e-9);p.add_argument('--pressure-limit',type=float,default=2e-8);p.set_defaults(run=fields)
    p=sub.add_parser('solver');p.add_argument('session',type=Path);p.add_argument('log',type=Path)
    p.add_argument('--mean-only',action='store_true');p.set_defaults(run=solver)
    p=sub.add_parser('statistics');p.add_argument('first',type=Path);p.add_argument('second',type=Path)
    p.add_argument('--relative-limit',type=float,default=.01);p.set_defaults(run=statistics)
    p=sub.add_parser('refinement');p.add_argument('first',type=Path);p.add_argument('second',type=Path)
    p.add_argument('--metric',required=True);p.add_argument('--ratio',type=float,default=.5);p.set_defaults(run=refinement)
    args=parser.parse_args()
    try:
        result=args.run(args)
    except (AssertionError,ValueError,KeyError,OSError) as error:
        print(f'M4 check failed: {error}',file=sys.stderr)
        return 2
    if args.json:
        args.json.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return 0 if result['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
