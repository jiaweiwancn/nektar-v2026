#!/usr/bin/env python3
"""Independent A/B audit of accepted current-solve and residual-map stages."""
import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys

sys.dont_write_bytecode=True
from cao_tamura_experiment import THREADS,allocation,environment,execute,observe,records,sha,write_csv,write_json
from cao_tamura_cv3_budgets import flat_momentum_budget
from cao_tamura_routes import coverage
from terrain_continuity import audit as strip_continuity
from validate_terrain import compare
from record_cao_tamura_r1 import field_bound


def audit(folder):
    r=json.loads((folder/'run.json').read_text());s=r['stages']
    if r['exit_code'] or not r['complete'] or not r['binaries_unchanged_during_run']:
        raise ValueError('Incomplete or changed producing run')
    if r['case_sha256']!=sha(folder/'case.xml'):raise ValueError('Producing XML changed')
    if (s['completed_stages']<3 or s['maximum_guard_activations'] or s['maximum_startup_weight']
            or min(s['minimum_k'],s['minimum_epsilon'])<=0 or
            max(s[k] for k in ('maximum_weak_divergence','maximum_projected_normal',
                 'maximum_normal_integral_defect','maximum_resolved_mass_defect'))>1e-7):
        raise ValueError('Native positivity/continuity/physical-time gate failed')
    values=[float(v['free_coefficient_l2']) for v in records(folder/'case.rans-residual.csv') if v['phase']=='1']
    if not values or not all(math.isfinite(v) and v<=1e-7 for v in values):
        raise ValueError('Independent full weak stage residual failed')
    return dict(run=r,run_sha256=sha(folder/'run.json'),maximum_free_residual=max(values),
        coverage=coverage(folder/'case.xml'),strip_continuity=strip_continuity(folder,'case'),
        strip_momentum=flat_momentum_budget(folder/'case.xml'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('baseline','candidate','build','output'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--cpus',required=True);a=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Linux is authoritative')
    os.environ.update(THREADS);os.sched_setaffinity(0,set(map(int,a.cpus.split(','))))
    build=environment(a.build);allocation(1,a.cpus)
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(__file__,out/'producer.py')
    write_csv(out/'point.csv',[dict(id='native-final-time',x=.1,span=0.,z=.1)],['id','x','span','z'])
    result=dict(version='cao-tamura-r1-map-equivalence-v2',status='RUNNING',cases={},
        producer_sha256=sha(out/'producer.py'))
    try:
        manifests=[json.loads((f/'manifest.json').read_text()) for f in (a.baseline,a.candidate)]
        if manifests[0]['config']!=manifests[1]['config']:
            raise ValueError('Different physical, mesh or time configuration')
        for name,folder in (('baseline',a.baseline),('candidate',a.candidate)):
            result['cases'][name]=audit(folder)
            binding=field_bound(folder/'case.fld')
            observe(folder/'case.xml',folder/'case.fld',out/'point.csv',out/(name+'-time.csv'),build,1,a.cpus)
            times={float(r['time']) for r in records(out/(name+'-time.csv'))}
            expected=result['cases'][name]['run']['stages']['final_time']
            if len(times)!=1 or abs(next(iter(times))-expected)>1e-10:
                raise ValueError('Native complete-field time differs from its producing run')
            execute(allocation(1,a.cpus)+[build/'RANSWallFieldAudit',folder/'case.xml',folder/'case.fld',out/(name+'.csv')],out,name)
            if field_bound(folder/'case.fld')!=binding:raise ValueError('Field changed during audit')
            result['cases'][name].update(field=binding,native_field_time=expected,
                full_csv_sha256=sha(out/(name+'.csv')),time_receipt_sha256=sha(out/(name+'-time.csv.json')))
        errors=compare(out/'baseline.csv',out/'candidate.csv')
        if max(errors.values())>1e-6:raise ValueError('Complete-field map equivalence failed: '+str(errors))
        old,new=[records(out/(name+'.csv')) for name in ('baseline','candidate')]
        relative={k:max(abs(float(x[k])-float(y[k]))/max(1e-20,abs(float(x[k]))) for x,y in zip(old,new))
                  for k in old[0] if k in ('k','eps','eps_tilde')}
        before,after=[float(re.findall(r'^Execute\s+([\d.eE+-]+)',(folder/'run.log').read_text(),re.M)[-1])
                      for folder in (a.baseline,a.candidate)]
        result.update(status='PASS',scaled_field_errors=errors,pointwise_relative_scalar_errors=relative,
            runtime_ratio_baseline_over_candidate=before/after,
            runtime_basis='native Execute wall-time region, identical measurement for both binaries',
            limitation='three unforced stages; same discrete equations, not a developed precursor or a convergence study',
            auditor_sha256=sha(build/'RANSWallFieldAudit'))
    except Exception as error:
        result.update(status='FAIL',error=repr(error));write_json(out/'evidence.json',result);raise
    write_json(out/'evidence.json',result);print('map equivalence PASS',result['scaled_field_errors'],
        'runtime ratio',result['runtime_ratio_baseline_over_candidate'],flush=True)


if __name__=='__main__':main()
