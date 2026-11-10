#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CLI success/failure contracts for historical evidence supplementation."""
import argparse,json,os,subprocess,sys
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_pilot_evidence import epsilon_budget_gate,load
from rans_preservation_oracle import Reference
from rans_qualification import write_json

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--previous-root',type=Path);p.add_argument('--epsilon-witness',type=Path);p.add_argument('--budget-output',type=Path);a=p.parse_args()
 if a.epsilon_witness:
  d=epsilon_budget_gate(load(a.epsilon_witness),Reference('E94E95','Natural','Accepted'));write_json(a.budget_output,d);return 0 if d['status']=='pass' else 1
 if not a.root or not a.previous_root:p.error('root and previous-root required')
 out=a.root/'evidence-contract';out.mkdir(exist_ok=False);collector=Path(__file__).with_name('collect_rans_neutral_pilot.py');checks=[]
 def collect(dest,extra=()):
  with (out/(dest+'.log')).open('w') as f:return subprocess.run([sys.executable,'-B',str(collector),'--root',str(a.previous_root),'--output',str(out/dest),*extra],stdout=f,stderr=subprocess.STDOUT).returncode
 code=collect('supplement');s=load(out/'supplement/summary.json');bodies=s['required_evolved_space_comparisons'];assert code==0 and s['status']=='pass' and len(bodies)==12 and all(len(v['cases'])==3 and len(v['fields'])==5 for v in bodies.values())
 checks.append(dict(name='all-twelve-complete-bodies',status='pass',exit_code=code,receipt=str(out/'supplement/summary.json')))
 missing=out/'eleven-comparisons';missing.mkdir();omitted='A93-N-Accepted-space-0.016'
 for name in bodies:
  if name!=omitted:(missing/(name+'.json')).symlink_to(a.previous_root/'p5/natural-schedule'/(name+'.json'))
 code=collect('missing-required',('--comparison-directory',str(missing)));s=load(out/'missing-required/summary.json');bad=[k for k,v in s['required_evolved_space_comparisons'].items() if v['status']!='pass'];assert code!=0 and bad==[omitted] and 'FileNotFoundError' in s['required_evolved_space_comparisons'][omitted]['failure']
 checks.append(dict(name='missing-required-comparison-rejects',status='pass',exit_code=code,reason=s['required_evolved_space_comparisons'][omitted]['failure']))
 ref=Reference('E94E95','Natural','Accepted');scale=2*(ref.Z-ref.H)*(ref.A/ref.H)*ref.us/ref.H
 for kind,defect,expected in [('passing',0.,0),('failing',2e-10*scale,1)]:
  fixture=out/(kind+'-epsilon.json');receipt=out/(kind+'-epsilon-result.json');write_json(fixture,[dict(time=.016,alpha=2500.,free_constant_test_defect=defect)])
  code=subprocess.run([sys.executable,'-B',__file__,'--epsilon-witness',str(fixture),'--budget-output',str(receipt)]).returncode
  assert code==expected and load(receipt)['status']==('pass' if code==0 else 'fail')
  checks.append(dict(name=kind+'-epsilon-budget-cli',status='pass',exit_code=code,receipt=load(receipt)))
 write_json(out/'results.json',dict(status='pass',checks=checks));print('pass');return 0
if __name__=='__main__':sys.exit(main())
