#!/usr/bin/env python3
"""CV3-R1 recovery of the documented unforced CLS EP-AWF reference controls."""
import argparse
import math
from pathlib import Path
import subprocess
import shutil
import sys

sys.dont_write_bytecode=True
from cao_tamura_experiment import ROOT, allocation, environment, execute, records, sha, write_json
from generate_m6 import generate_abl


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=('reference','perturbation','molecular'),required=True)
    p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpus',required=True);p.add_argument('--timeout',type=int,default=7200)
    a=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Use Linux for authoritative validation')
    build=environment(a.build);allocation(1,a.cpus)
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    case=out/'reference.xml';perturbed=a.case=='perturbation'
    spec=generate_abl(case,model='cls_switched' if perturbed else 'cls_fixed',dimension=2,
        degree=4,layers=12,nx=1,dt=.5,flow_times=5,backend='matrixfree',
        boundary='outflow',monitor=True,tolerance=1e-8,perturbation=.01 if perturbed else 0.,
        diffusion='Molecular' if a.case=='molecular' else 'HighRe')
    shutil.copy2(__file__,out/'producer.py')
    binaries={n:sha(build/n) for n in ('IncNavierStokesSolver','RANSFieldAudit')}
    result=dict(status='RUNNING',case=a.case,source_parent=subprocess.check_output(
        ['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),specification=spec,
        case_sha256=sha(case),binaries=binaries,script_sha256=sha(__file__),
        generator_sha256=sha(Path(__file__).parent/'generate_m6.py'))
    write_json(out/'evidence.json',result)
    try:
        result['run']=execute(allocation(1,a.cpus)+[build/'IncNavierStokesSolver',case.name],out,'solver',a.timeout)
        stages=records(out/'reference.rans.csv');wall=records(out/'reference.awf.csv')
        errors=records(out/'reference.rans-errors.csv')
        assert len(stages)==spec['steps'] and abs(float(stages[-1]['time'])-spec['physical_time'])<1e-9
        assert all(math.isfinite(float(v)) for r in stages for v in r.values())
        assert all(float(r['min_k'])>0 and float(r['min_eps'])>0 and
                   float(r['guard_activations'])==0 and float(r['weak_divergence_l2'])<1e-7 for r in stages)
        assert wall and all(float(r[k])==0 for r in wall for k in
            ('guard_activations','cutoff_fallback_points','regime_fallback_points','range_fallback_points'))
        assert all(float(r[k])==1. for r in wall for k in ('chi_min','chi_max'))
        fields=('u','k','eps_tilde','rans_nut','rans_R_xx','rans_R_yy','rans_R_zz','rans_R_xz')
        peak={f:{n:max(float(r[n]) for r in errors if r['field']==f) for n in
            ('normalized_l2','pointwise_relative_linf')} for f in fields}
        final={f:next(r for r in reversed(errors) if r['field']==f) for f in fields}
        if perturbed:
            initial=next(float(r['normalized_l2']) for r in errors if r['field']=='k')
            decay=initial/float(final['k']['normalized_l2'])
            assert decay>1000 and all(float(final[f]['normalized_l2'])<.001 for f in fields)
            result['k_error_decay_factor']=decay
        else:assert all(v<.001 for r in peak.values() for v in r.values())
        execute(allocation(1,a.cpus)+[build/'RANSFieldAudit',case,out/'reference.fld',out/'reference.fld'],out,'audit')
        assert binaries=={n:sha(build/n) for n in binaries}
        result.update(status='PASS',peak_errors=peak,final_errors=final,
            maximum_weak_divergence=max(float(r['weak_divergence_l2']) for r in stages),
            minimum_k=min(float(r['min_k']) for r in stages),minimum_epsilon=min(float(r['min_eps']) for r in stages),
            artifacts={f.name:sha(f) for f in out.iterdir() if f.suffix in ('.xml','.csv','.log')})
    except Exception as error:
        result.update(status='FAIL',error=repr(error));write_json(out/'evidence.json',result);raise
    write_json(out/'evidence.json',result);print(a.case,'PASS',flush=True)


if __name__=='__main__':main()
