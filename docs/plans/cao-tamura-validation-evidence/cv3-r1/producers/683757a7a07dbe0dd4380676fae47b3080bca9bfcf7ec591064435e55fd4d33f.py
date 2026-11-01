#!/usr/bin/env python3
"""Compare matched native checkpoints; this never certifies development."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode=True
from cao_tamura_experiment import (THREADS,allocation,environment,execute,observe,
    records,sha,write_csv,write_json)
from record_cao_tamura_r1 import field_bound
from validate_terrain import compare


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('baseline','candidate','baseline-field','candidate-field','build','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--time',type=float,required=True);p.add_argument('--cpus',required=True)
    a=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Linux is authoritative')
    os.environ.update(THREADS);allocation(1,a.cpus)
    os.sched_setaffinity(0,set(map(int,a.cpus.split(','))))
    build=environment(a.build);out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(__file__,out/'producer.py')
    write_csv(out/'point.csv',[dict(id='native-time-binding',x=.1,span=0.,z=.1)],['id','x','span','z'])
    result=dict(status='RUNNING',scope='Matched physical checkpoint only; no development or calibration claim',
        producer_sha256=sha(out/'producer.py'),cases={},time=a.time,
        binaries={n:sha(build/n) for n in ('IncNavierStokesSolver','RANSWallFieldAudit','RANSObservationAudit')})
    write_json(out/'evidence.json',result)
    try:
        configs=[]
        for name,folder,field in (('baseline',a.baseline,a.baseline_field),('candidate',a.candidate,a.candidate_field)):
            m=json.loads((folder/'manifest.json').read_text());r=json.loads((folder/'run-start.json').read_text())
            if (r['case_sha256']!=sha(folder/'case.xml') or m['case_sha256']!=r['case_sha256']
                    or r['solver_sha256']!=result['binaries']['IncNavierStokesSolver']
                    or r['geometry_sha256']!=sha(build/'geometry-lib/libSpatialDomains.so')):
                raise ValueError('Unbound producing case/binary: '+name)
            configs.append({k:v for k,v in m['config'].items() if k not in ('steps','checkpoint')})
            binding=field_bound(field)
            observe(folder/'case.xml',field,out/'point.csv',out/(name+'-time.csv'),build,1,a.cpus)
            times={float(row['time']) for row in records(out/(name+'-time.csv'))}
            if len(times)!=1 or abs(next(iter(times))-a.time)>1e-10:
                raise ValueError('Wrong native checkpoint time: '+str(times))
            execute(allocation(1,a.cpus)+[build/'RANSWallFieldAudit',folder/'case.xml',field,out/(name+'.csv')],out,name)
            rows=records(out/(name+'.csv'));epsilon='eps_tilde' if 'eps_tilde' in rows[0] else 'eps'
            if any(not math.isfinite(float(row[k])) or float(row[k])<=0 for row in rows for k in ('k',epsilon)):
                raise ValueError('Nonpositive/nonfinite checkpoint scalar')
            if field_bound(field)!=binding:raise ValueError('Checkpoint changed during audit')
            result['cases'][name]=dict(field=binding,case_sha256=r['case_sha256'],
                run_start_sha256=sha(folder/'run-start.json'),manifest_sha256=sha(folder/'manifest.json'),
                time=a.time,full_csv_sha256=sha(out/(name+'.csv')))
        if configs[0]!=configs[1]:raise ValueError('Different physical/mesh/time-step configurations')
        errors=compare(out/'baseline.csv',out/'candidate.csv')
        if not all(math.isfinite(v) and v<=1e-6 for v in errors.values()):
            raise ValueError('Matched-checkpoint fields differ: '+str(errors))
        result.update(status='PASS',scaled_field_errors=errors)
    except Exception as error:
        result.update(status='FAIL',error=repr(error));write_json(out/'evidence.json',result);raise
    write_json(out/'evidence.json',result);print('Matched native checkpoint PASS',errors,flush=True)


if __name__=='__main__':main()
