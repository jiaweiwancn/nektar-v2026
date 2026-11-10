#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CLI admission test with real partial native evidence and missing-work resume."""
import argparse,os,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load
from run_cls_epawf_picard import sha

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'resume-cli-regression';out.mkdir(exist_ok=False);grid=r/'preservation-spatial';grid.mkdir(exist_ok=False);write_json(grid/'evaluation-grid.json',dict(grid=[[0.,0.,.5],[0.,0.,3.]],weights=[4.,4.]))
 os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'));script=Path(__file__).with_name('qualify_rans_preservation.py')
 base=[sys.executable,'-B',str(script),'--root',str(r),'--phase','temporal','--label','reuse-contract','--models','E94E95','--policies','Natural','--walls','Accepted','--horizon','.0001','--boundary-policy','FrozenReferenceV1'];campaign=r/'preservation-temporal-reuse-contract';checks=[]
 def run(name,extra):
  with (out/(name+'.log')).open('w') as f:return subprocess.run(base+extra,stdout=f,stderr=subprocess.STDOUT).returncode
 try:
  code=run('partial',['--new-check-limit','1']);assert code==2,code
  first=campaign/'E94E95-Natural-Accepted-dt0';saved={x.name:(sha(x),x.stat().st_mtime_ns) for x in first.iterdir() if x.is_file()};receipt=sha(campaign/'results.json')
  alias=r/'runtime-alias';alias.mkdir();[(alias/n).symlink_to(r/'runtime'/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','source-build.json','lib')]
  for name,extra in [('runtime',['--runtime-directory',str(alias)]),('quadrature',['--quadrature-index','0']),('perturbation',['--epsilon-perturbation','.01']),('horizon',['--horizon','.0002']),('boundary',['--boundary-identity','changed'])]:
   code=run('reject-'+name,['--resume']+extra);text=(out/('reject-'+name+'.log')).read_text();assert code!=0 and 'immutable request mismatch' in text,(name,code,text[-500:]);assert sha(campaign/'results.json')==receipt
   checks.append(dict(name='changed-'+name,status='pass',exit_code=code,intended_rejection='immutable request mismatch',aggregate_unchanged=True))
  code=run('resume',['--resume']);assert code in (0,2),code
  assert all((sha(first/x), (first/x).stat().st_mtime_ns)==v for x,v in saved.items()),'completed trajectory was mutated'
  complete=load(campaign/'results.json');assert len(complete['checks'])==3 and all(x['status']=='pass' for x in complete['checks'])
  checks.append(dict(name='identical-resume-only-missing-native-work',status='pass',completed_numerical_files_unchanged=True,new_native_cases=2,exit_code=code,accuracy_scope='request/reuse contract only; two-point evaluation is not pilot evidence'))
  # Runtime mismatches are separately exercised by the manifest unit test with
  # the exact runtime identity object. No frozen executable is overwritten.
  status='pass'
 except Exception as e:checks.append(dict(name='native-resume-contract',status='fail',failure=repr(e)));status='fail'
 write_json(out/'results.json',dict(status=status,checks=checks));print(status);return 0 if status=='pass' else 1
if __name__=='__main__':sys.exit(main())
