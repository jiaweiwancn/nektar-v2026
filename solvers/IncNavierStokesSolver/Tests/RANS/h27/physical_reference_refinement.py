#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Run the bounded source19 physical reference refinement and compare owner fields."""

import argparse
import json
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path('/home/jay/build/wj-h27-coupling-r1-20261117')
NAME = 'cp2-refinement-v19'
SOURCE = '56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'
ELF = '37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b'
TABLE = 'b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'
SEEDS = {
    'base': ('physical-reference-h2-v1/reference-h2',
             '686bd303c511a504559ec646a37a92671695f824dffc25897d028ad40585b7c6'),
    'streamwise': ('physical-reference-x2-v1/reference-x2',
                   'e60ea109d0633d96f6b79a8d679c402228e3c9ab2c47a1261d052ccfdaeceb89'),
    'radial': ('cutoff-panels-v19/reference-r2',
               'f48243ee0d94b39eb394ab05eb6a00bc23d0caf91956783ab7e62699afcd072e'),
}
CASES = {
    'base': (450, 210), 'quadrature': (520, 210), 'tight': (420, 210),
    'streamwise': (1100, 355), 'targeted': (1300, 355),
}
COMPARISONS = (
    ('base', 'quadrature', 130), ('base', 'tight', 130),
    ('base', 'streamwise', 180), ('base', 'radial', 180),
    ('streamwise', 'targeted', 270), ('radial', 'targeted', 270),
)
H27 = ROOT/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27'
OUT = ROOT/'science'/NAME
TABLE_PATH = ROOT/'science/reference-tables-v1.json'


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def axes(path):
    vertices = [list(map(float, v.text.split())) for v in ET.parse(path).findall('GEOMETRY/VERTEX/V')]
    return [sorted({v[i] for v in vertices}) for i in (0, 1)]


def native(case):
    prefix = ROOT/'science'/SEEDS[case][0] if case == 'radial' else OUT/case
    return prefix.with_suffix('.result.json')


def audit(case):
    return ROOT/'science/cutoff-panels-v19/independent32.json' if case == 'radial' else OUT/(case+'-current.json')


def markers(path):
    from resources import identity

    fields = read(path)
    state = fields['state']
    require(fields['normal_main'] and fields['factory'] == 'WJNativeCoupledV2' and
            fields['steady_solver']['converged'] and state['step'] == state['time'] == 0 and
            state['reference_identity'] == TABLE, f'invalid steady field {path}')
    return {**identity(path), 'source_identity': state['source_identity'], 'state_sha256': state['state_sha256']}


def fixtures():
    from coupled_fixtures import generate, mesh
    from resources import atomic, identity

    OUT.mkdir()
    known = {}
    for case, (loc, digest) in SEEDS.items():
        prefix = ROOT/'science'/loc
        row = markers(prefix.with_suffix('.result.json'))
        require(row['sha256'] == digest, f'archived {case} seed bytes differ')
        known[case] = dict(seed=row, template=identity(prefix.with_suffix('.xml')))
    x, z = axes(Path(known['base']['template']['path']))
    xx, zx = axes(Path(known['streamwise']['template']['path']))
    xr, zz = axes(Path(known['radial']['template']['path']))
    require((len(x), len(z), len(xx), len(zz)) == (29, 101, 57, 201) and x == xr and z == zx and
            xx == sorted(x+[(a+b)/2 for a, b in zip(x[:-1], x[1:])]) and set(z) <= set(zz),
            'archived physical mesh geometry differs from the admitted axes')
    split = (1, 2, 3, *range(15, 24), 27)
    xt = sorted(x+[-9.921875, -9.6875, -9.375]+[(x[i]+x[i+1])/2 for i in split])
    zt = sorted(set(z) | {v for v in zz if v < .2375})
    require((len(xt), len(zt)) == (45, 133) and .2375 in zt, 'wrong localized material grid')
    template = OUT/'targeted-template.xml'
    mesh(template, xt, zt, quadrature=12)
    known['targeted'] = dict(seed=known['radial']['seed'], template=identity(template))
    sessions = {}
    for case, (cap, _) in CASES.items():
        reference = known.get(case, known['base'])
        path = OUT/(case+'.xml')
        source_template = reference['template']['path']
        seed = reference['seed']['path']
        seed_source = reference['seed']['source_identity']
        # Single-factor controls import the converged source19 base coefficients.
        if case in ('quadrature', 'tight'):
            seed, seed_source = str(native('base')), SOURCE
        qx, qr = (12, 16) if case == 'quadrature' else (8, 12)
        linear, nonlinear = (3e-13, 3e-12) if case == 'tight' else (1e-12, 1e-10)
        generate(path, TABLE_PATH, OUT/case, mode='Physical', dt=.05, steps=0, solve='SteadyReference',
                 quadrature_x=qx, quadrature_r=qr, linear=linear, nonlinear=nonlinear,
                 backend='FastFull3DV1', steady_seconds=cap-45, template=source_template,
                 initial_fields=seed, initial_source_identity=seed_source)
        root = ET.parse(path).getroot()
        ET.SubElement(root.find('CONDITIONS/PARAMETERS'), 'P').text = 'WJSteadyKrylovDimension = 80'
        ET.indent(root)
        path.write_text(ET.tostring(root, encoding='unicode')+'\n')
        xy, zy = axes(path)
        sessions[case] = dict(session=identity(path), options=identity(Path(str(path)+'.opt')),
                              seed_path=seed, seed_source=seed_source,
                              cells=[len(xy)-1, len(zy)-1], x=xy, physical_z=zy,
                              quadrature=[qx, qr], tolerances=[nonlinear, linear])
    atomic(OUT/'inputs.json', dict(schema='WJH27ReferenceRefinementV1', source=SOURCE, table=identity(TABLE_PATH),
                                  runtime=identity(ROOT/'science/runtime-v19.json'), archived=known,
                                  sessions=sessions, generator=identity(__file__)))
    print('SOURCE19_REFINEMENT_FIXTURES_READY', flush=True)


def compare(first, second, seconds, order):
    from physical_oracle import PhysicalOracle, Reference, load_native
    from resources import atomic, identity

    deadline = time.monotonic()+seconds
    reference = Reference(read(TABLE_PATH))
    rows, oracles = {}, []
    for case in (first, second):
        a = audit(case)
        p = native(case)
        report = read(a)
        data, metadata = load_native(p)
        row = dict(native=identity(p), audit=identity(a), state_sha256=data['state_sha256'])
        require(row['native']['sha256'] == report['native']['sha256'] and
                metadata['steady_solver']['converged'] and data['source_identity'] == SOURCE and
                report['native_source_identity'] == SOURCE and data['state_sha256'] == report['state_sha256'],
                f'current source19 audit binding differs for {case}')
        core = report['core']
        row['equations'] = report['core_steady_equations_status']
        row['boundary'] = report['boundary_checks']['status']
        row['physical_development'] = report['physical_development_status']
        row['residuals'] = core['residuals']
        row['budgets'] = {key: {k: v[k] for k in ('fixed_scale', 'mixed_physical_balance',
                         'all_owner_gradient_balance', 'status')} for key, v in core['budgets'].items()
                         if key in ('k', 'epsilon', 'mechanical')}
        rows[case] = row
        oracles.append(PhysicalOracle(data, reference, order, deadline))
    result = oracles[0].compare(oracles[1])
    result['current_equations_and_boundary_pass'] = all(v['equations'] == v['boundary'] == 'MET' for v in rows.values())
    result['normalized_budget_changes'] = {}
    for key, a in rows[first]['budgets'].items():
        b = rows[second]['budgets'][key]
        require(a['fixed_scale'] == b['fixed_scale'], 'physical budget scale differs')
        result['normalized_budget_changes'][key] = {name: (b[name]-a[name])/a['fixed_scale']
            for name in ('mixed_physical_balance', 'all_owner_gradient_balance')}
    path = OUT/f'{first}-to-{second}-order{order}.json'
    atomic(path, dict(schema='WJH27CommonOwnerRefinementV1', source_identity=SOURCE,
                      table=identity(TABLE_PATH), native_and_audit=rows, quadrature_order=order, **result))
    print(json.dumps(dict(path=str(path), norms=result['status'], rows=len(result['norms']),
                          equations=result['current_equations_and_boundary_pass'])), flush=True)


def run():
    from resources import atomic, identity, size, workers

    account = read(ROOT/'accounting.json')
    charged = sum(row['seconds'] for row in account['entries'].values() if row['status'] != 'RESERVED')
    require(not account['closed'] and not any(v['status'] == 'RESERVED' for v in account['entries'].values()) and
            charged+6500 <= account['ceiling']-account['evidence_protected'], 'the full CP2 matrix must fit the open account')
    require(size(ROOT/'science')+int(.9*1024**3) <= 8*1024**3 and
            shutil.disk_usage(ROOT).free >= int(15.9*1024**3) and workers()['workers'] <= 65,
            'insufficient room for the admitted matrix')
    require(identity(TABLE_PATH)['sha256'] == TABLE and
            identity(ROOT/'runtime-v19/IncNavierStokesSolver.elf')['sha256'] == ELF and
            read(ROOT/'source-v19.json')['expected_source_identity'] == SOURCE,
            'wrong pinned application, source or table')
    runner = [sys.executable, '-B', str(ROOT/'tools/run.py'), '--root', str(ROOT)]
    receipts = {}

    def operation(label, phase, seconds, *command):
        subprocess.run(runner+['--label', NAME+'-'+label, '--phase', phase,
                              '--seconds', str(seconds), '--', *map(str, command)], check=False)
        p = ROOT/'operations'/(NAME+'-'+label)/'receipt.json'
        receipts[label] = read(p)['status'] if p.exists() else 'NO_RECEIPT'
        if receipts[label] == 'NO_RECEIPT':
            raise RuntimeError(f'runner did not record {label}; inspect the live account')
        return receipts[label] == 'PASS'

    py = [sys.executable, '-B']
    require(operation('fixtures', 'analysis', 40, *py, __file__, 'fixtures'), 'fixture operation failed')
    states = {'radial': markers(native('radial'))}
    audited = {'radial'}
    for case, (cap, audit_cap) in CASES.items():
        if case in ('quadrature', 'tight') and 'base' not in states:
            continue
        xml = OUT/(case+'.xml')
        if not operation(case+'-native', 'science', cap, ROOT/'runtime-v19/IncNavierStokesSolver',
                         '--use-opt-file', str(xml)+'.opt', xml):
            continue
        try:
            row = markers(native(case))
            require(row['source_identity'] == SOURCE, 'wrong result source')
        except (ValueError, OSError) as error:
            print(f'STEADY_CANDIDATE_REJECTED {case} {error}', flush=True)
            continue
        states[case] = row
        if operation(case+'-current', 'analysis', audit_cap, *py, H27/'physical_oracle.py', '--native', native(case),
                     '--table', TABLE_PATH, '--out', audit(case), '--order', '32', '--seconds', str(audit_cap-15)):
            audited.add(case)
    for first, second, cap in COMPARISONS:
        if {first, second} <= audited:
            operation(first+'-to-'+second, 'analysis', cap, *py, __file__, 'compare', first, second,
                      '--order', '20', '--seconds', str(cap-10))
    atomic(OUT/'run-index.json', dict(schema='WJH27RefinementOperationsV1', status='MEASURED', source_identity=SOURCE,
                                     preflight_charged_seconds=charged, launch_cap_seconds=6500,
                                     states=states, audited=sorted(audited), receipts=receipts,
                                     reference_qualified=False, TASK3_READY=False))
    print('SOURCE19_REFINEMENT_MATRIX_FINISHED', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'fixtures', 'compare'))
    parser.add_argument('case', nargs='*')
    parser.add_argument('--order', type=int, choices=(20, 32), default=20)
    parser.add_argument('--seconds', type=float, default=0)
    args = parser.parse_args()
    require(sys.platform.startswith('linux') and ROOT.is_dir(), 'authoritative Linux H27 capsule required')
    sys.path[:0] = [str(ROOT/'tools'), str(H27)]
    if args.action == 'compare':
        require(len(args.case) == 2 and set(args.case) <= ({'radial'} | CASES.keys()) and args.seconds > 0,
                'two admitted cases and a finite comparison cap required')
        compare(*args.case, args.seconds, args.order)
    else:
        require(not args.case, 'unexpected case arguments')
        (run if args.action == 'run' else fixtures)()


if __name__ == '__main__':
    main()
