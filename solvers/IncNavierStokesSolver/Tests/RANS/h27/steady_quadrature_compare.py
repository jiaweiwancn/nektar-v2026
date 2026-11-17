#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind independent high-order equations and common-owner steady q changes."""
import argparse,json,math,os,resource,sys,time
from pathlib import Path
sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import history_recovery_check as hr
import physical_cache_history_check as history
import physical_oracle as po
import physical_transient_oracle as transient
import implicit_algorithm_check as fields
require,read = po.require,history.read


def run(manifest,folder,table,deadline,out):
    require(manifest['schema'] == 'WJH27SteadyQuadratureComparisonV1','wrong quadrature comparison manifest')
    root = hr.resolve(folder,manifest['root']);require(hr.identity(table)['sha256'] == po.TABLE_SHA,'reference table changed')
    comparison = manifest.get('comparison','OrderRefinementV1')
    require(comparison in ('OrderRefinementV1','MovingCutoffPanelsV1'),'unrecognized quadrature comparison')
    paths = {key:history.input_path(manifest[key],root)[0] for key in ('old','new')}
    packets = {key:read(path) for key,path in paths.items()};a,b = (packets[key]['state'] for key in ('old','new'))
    for packet in packets.values():
        require(packet['steady_mode'] is True and packet['steady_solver']['converged'] is True and
            packet['state']['step'] == packet['state']['time'] == 0,'comparison requires converged zero-time steady states')
    for key in fields.EXACT_KEYS:
        if key != 'quadratureX' and not (comparison == 'MovingCutoffPanelsV1' and key == 'source_identity'):
            require(a[key] == b[key],'quadrature comparison changed another input: '+key)
    require(a['quadratureX'] == manifest['old_quadratureX'] and b['quadratureX'] == manifest['new_quadratureX'] and
        a['signature'] != b['signature'],'quadrature/signature binding differs')
    if comparison == 'OrderRefinementV1':
        require(a['quadratureX'] < b['quadratureX'] and a['source_identity'] == manifest['source_identity'],'order/source binding differs')
    else:
        require(a['quadratureX'] == 16 and b['quadratureX'] == 8 and
            a['source_identity'] == manifest['old_source_identity'] != b['source_identity'] == manifest['source_identity'] and
            b['scalar_quadrature_policy'] == 'CubicCutoffFacePanelsV1','panel-repair source/policy binding differs')
    out['comparison'] = comparison
    ref = po.Reference(read(table));old = po.PhysicalOracle(a,ref,20,deadline);new = po.PhysicalOracle(b,ref,20,deadline)
    out['common_owner_fields'] = old.compare(new)
    out['coefficients'] = fields.coefficient_delta(transient.state(a['coefficients'],a['sizes']),transient.state(b['coefficients'],b['sizes']),a['sizes'],ref)
    out['independent_equations'] = {};equations = True
    require(set(manifest['independent']) == {'32','48'},'both independent32/48 audits are required')
    for label,spec in manifest['independent'].items():
        p,ident = history.input_path(spec,root);j = read(p)
        require(j['schema'] == po.SCHEMA and j['native']['sha256'] == hr.identity(paths['new'])['sha256'] and
            j['native_source_identity'] == b['source_identity'] and j['table']['sha256'] == po.TABLE_SHA and
            j['quadrature'] == [int(label),int(label)],'independent quadrature report identity differs')
        passed = all(value == 'MET' for value in [j['core_steady_equations_status'],j['boundary_checks']['status'],j['core']['accounting_identity_status']])
        equations = equations and passed
        out['independent_equations'][label] = dict(file=ident,passed=passed,residuals=j['core']['residuals'],boundaries=j['boundary_checks']['status'],
            accounting=j['core']['accounting_identity_status'],physical_development_status=j['physical_development_status'],
            physical_envelopes=j['core']['physical_envelope'])
    out['quadrature_equation_gate_pass'] = equations
    out['field_comparison_gate_pass'] = out['common_owner_fields']['status'] == 'MET'
    out['quadrature_foundation_gate_pass'] = equations and out['field_comparison_gate_pass']
    out['status'] = 'PASS_STEADY_QUADRATURE_FOUNDATION_ONLY' if out['quadrature_foundation_gate_pass'] else 'NOT_MET'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','table','out'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True);args = parser.parse_args();began = time.monotonic()
    require(not args.out.exists(),'quadrature comparison output is immutable')
    out = dict(schema='WJH27SteadyQuadratureComparisonCheckV1',status='FAIL',quadrature_foundation_gate_pass=False,full_H27_qualification=False,TASK3_READY=False,
        manifest=hr.identity(args.manifest),table=hr.identity(args.table),workers=1,
        sources={key:hr.identity(path) for key,path in dict(checker=__file__,reference_fields=po.__file__,coefficients=fields.__file__,JSON_reader=hr.__file__).items()})
    try:
        require(sys.platform.startswith('linux') and math.isfinite(args.seconds) and 0 < args.seconds <= 120,'bounded authoritative Linux assessment required')
        soft,hard = resource.getrlimit(resource.RLIMIT_AS);cap = min([8*1024**3]+[v for v in (soft,hard) if v != resource.RLIM_INFINITY])
        resource.setrlimit(resource.RLIMIT_AS,(cap,hard));run(read(args.manifest),args.manifest.resolve().parent,args.table.resolve(),began+args.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-began;out['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    with args.out.open('x') as stream:json.dump(out,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({key:out[key] for key in ('status','reason','elapsed_seconds') if key in out}))
    if 'exception' in out:raise SystemExit(2)


if __name__ == '__main__':main()
