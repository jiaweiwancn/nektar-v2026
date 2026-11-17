#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Cross-source regression for an explicitly bound native numerical repair.

Only source/signature may differ. Independent unchanged current-equation
reports, full histories, two owner samplings and fixed1e-8 gates are required.
"""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import history_recovery_check as hr
import implicit_algorithm_check as fields
import physical_cache_history_check as history
import physical_oracle as po

require,read = po.require,history.read


def run(manifest,folder,table,deadline,out):
    require(manifest['schema'] == 'WJH27HTraceRegressionV1' and set(manifest['cases']) == {'poly','hist','physical'},'incomplete H-trace regression matrix')
    scope = manifest.get('repair_scope','HTraceV1')
    require(scope in ('HTraceV1','MovingCutoffPanelsV1'),'unrecognized numerical repair scope')
    root = hr.resolve(folder,manifest['root']);require(hr.identity(table)['sha256'] == po.TABLE_SHA,'wrong unchanged reference table')
    source_ids = manifest['source_identities']
    require(set(source_ids) == {'old','new'} and source_ids['old'] != source_ids['new'] and
        all(len(value) == 64 and set(value) <= set('0123456789abcdef') for value in source_ids.values()),'explicit different source identities required')
    ref = po.Reference(read(table));results = {}
    for name,spec in manifest['cases'].items():
        fields.clock(deadline);paths = {key:history.input_path(spec[key],root)[0] for key in ('old','new')}
        old,new = (fields.load(paths[key],'CoupledJFNKV1') for key in ('old','new'));a,b = old['data'],new['data']
        require(a['source_identity'] == source_ids['old'] and b['source_identity'] == source_ids['new'],'regression source binding differs')
        for key in fields.EXACT_KEYS:
            if key != 'source_identity':require(a[key] == b[key],'geometry repair changed another configuration: '+key)
        require(old['scheme'] == new['scheme'] and a['reference_identity'] == po.TABLE_SHA and a['manufactured'] is (name != 'physical'),
            'regression scheme/reference/physical mode differs')
        owners,_ = fields.inventory(a,ref);fields.inventory(b,ref)
        result = dict(old=hr.identity(paths['old']),new=hr.identity(paths['new']),sources=source_ids,
            coefficients={},owner_fields={},dependent={},cache_logs={},scaled_limit=1e-8,
            pressure_policy='Use actual canonical pressure level; no alignment')
        for key in fields.STATE_NAMES+fields.CACHE_NAMES:
            delta = fields.coefficient_delta(old['states'][key],new['states'][key],a['sizes'],ref)
            require(all(v['scaled_MAX'] <= 1e-8 for v in delta.values()),'geometry repair changed a history coefficient beyond1e-8')
            result['coefficients'][key] = delta
        for key in fields.STATE_NAMES:
            fields.clock(deadline)
            values = [fields.physical_delta(old['states'][key],new['states'][key],a,ref,q,grid,deadline) for q,grid in ((8,9),(12,17))]
            require(all(v.get('status','MET') == 'MET' for row in values for v in row['norms'].values()),'geometry repair changed owner fields beyond1e-8')
            result['owner_fields'][key] = values
            dep = fields.dependent_delta(old['states'][key],new['states'][key],a,owners,ref,deadline)
            require(dep['gate_pass'],'geometry repair changed dependent K/Gamma/cutoff/Ew beyond1e-8');result['dependent'][key] = dep
        result['essential_lifts'] = {}
        for key,scale in [('essential_flow_lift',1.),('essential_k_lift',ref.Kin),('essential_epsilon_lift',ref.A/ref.H)]:
            x,y = np.asarray(a[key],float),np.asarray(b[key],float)
            require(x.shape == y.shape and np.all(np.isfinite(x)) and np.all(np.isfinite(y)),'essential lift inventory differs')
            error = float(np.max(abs(x-y),initial=0.)/scale);require(error <= 1e-8,'essential lift exceeds1e-8');result['essential_lifts'][key] = error
        duals = fields.cache_duals(old['states'],new['states'],a,ref,deadline)
        require(all(v <= 1e-8 for row in duals.values() for v in row.values()),'explicit history mass-dual difference exceeds1e-8')
        result['explicit_history_duals'] = duals
        for key,item in [('old',old),('new',new)]:
            result['cache_logs'][key] = fields.cache_log(paths[key],item['data'],'CurrentDataDriftV1',deadline)
        rp,rid = history.input_path(spec['current_equations'],root);report = read(rp)
        if name == 'physical':
            states,_ = history.binary_history(history.replace_suffix(paths['new'],'.result.json','.final.bin'),b,.005,deadline)
            result['independent_current'] = history.equation_report(spec['current_equations'],root,paths['new'],b,history.exported(b,states),po.TABLE_SHA)
        else:
            require(report['native_sha256'] == hr.identity(paths['new'])['sha256'] and report['table_sha256'] == po.TABLE_SHA and
                report['native_source_identity'] == source_ids['new'] and report['current_coupled_step_gate_pass'] is True,
                'independent continuous MMS current-equation gate failed')
            result['independent_current'] = dict(file=rid,current_coupled_step_gate_pass=True,status=report['status'])
        result['gate_pass'] = True;results[name] = result
    out.update(status='PASS_NATIVE_HTRACE_REGRESSION_ONLY' if scope == 'HTraceV1' else 'PASS_NATIVE_CUTOFF_PANEL_REGRESSION_ONLY',
        repair_scope=scope,foundation_gate_pass=True,cases=results,
        full_H27_qualification=False,TASK3_READY=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','table','out'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True);args = parser.parse_args();began = time.monotonic()
    require(not args.out.exists(),'H-trace regression output is immutable')
    out = dict(schema='WJH27HTraceRegressionCheckV1',status='FAIL',foundation_gate_pass=False,full_H27_qualification=False,TASK3_READY=False,
        manifest=hr.identity(args.manifest),table=hr.identity(args.table),workers=1,
        sources={key:hr.identity(path) for key,path in dict(checker=__file__,fields=fields.__file__,history=history.__file__,reference=po.__file__,JSON_reader=hr.__file__).items()})
    try:
        require(sys.platform.startswith('linux') and math.isfinite(args.seconds) and 0 < args.seconds <= 90,'bounded authoritative Linux analysis required')
        soft,hard = resource.getrlimit(resource.RLIMIT_AS);cap = min([8*1024**3]+[v for v in (soft,hard) if v != resource.RLIM_INFINITY])
        resource.setrlimit(resource.RLIMIT_AS,(cap,hard));run(read(args.manifest),args.manifest.resolve().parent,args.table.resolve(),began+args.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-began;out['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    with args.out.open('x') as stream:json.dump(out,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({key:out[key] for key in ('status','reason','elapsed_seconds') if key in out}))
    if not out['foundation_gate_pass']:raise SystemExit(2)


if __name__ == '__main__':main()
