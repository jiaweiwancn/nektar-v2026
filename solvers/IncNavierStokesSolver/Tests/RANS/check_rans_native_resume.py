#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CLI admission test with real partial native evidence and missing-work resume."""
import argparse,os,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load
from run_cls_epawf_picard import sha

def verify_existing(r):
 from rans_request_manifest import RequestManifest
 campaign=r/'preservation-temporal-reuse-contract';request=load(campaign/'request-manifest.json')['request'];m=RequestManifest(campaign,request,True)
 folders=sorted(campaign.glob('E94E95-Natural-Accepted-dt*'));assert len(folders)==3 and all(m.completed(f) for f in folders)
 before={str(p):[sha(p),p.stat().st_mtime_ns] for f in folders for p in f.iterdir() if p.is_file()};original=sha(campaign/'results.json')
 command=[sys.executable,'-B',str(Path(__file__).with_name('qualify_rans_preservation.py')),'--root',str(r),'--phase','temporal','--label','reuse-contract','--models','E94E95','--policies','Natural','--walls','Accepted','--horizon','.0001','--boundary-policy','FrozenReferenceV1','--resume']
 with (r/'resume-cli-verified.log').open('w') as f:code=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT).returncode
 d=load(campaign/'results.json');failed=[name for t in d['temporal_tables'] for name,v in t['fields'].items() if v['status']=='fail']
 assert code==1 and d['status']=='fail' and failed==['Pi'],(code,d['status'],failed)
 assert original==sha(campaign/'results.json'),'original numerical failure receipt changed'
 assert all([sha(Path(p)),Path(p).stat().st_mtime_ns]==v for p,v in before.items()),'completed numerical artifacts changed'
 previous=load(r/'resume-cli-regression/results.json');guards=[c for c in previous['checks'] if c['name'].startswith('changed-')];assert len(guards)==5 and all(c['status']=='pass' for c in guards)
 log=(r/'resume-cli-regression/resume.log').read_text();assert 'Accepted-dt1' in log and 'Accepted-dt2' in log and 'Accepted-dt0' not in log
 result=dict(status='pass',execution_status='completed',reuse_status='verified',numerical_campaign_status='fail',numerical_failure_reason='two-point, T=.0001 Pi temporal order outside original band; retained, not pilot evidence',guards=guards,completed_cases=3,only_missing_native_cases_completed=[1,2],existing_case0_seal_verified=True,repeat_resume_reran_native_cases=False,artifacts_and_original_receipt_unchanged=True,expected_negative=dict(exit_code=code,required_failed_field='Pi',order=d['temporal_tables'][0]['fields']['Pi']['order']))
 write_json(r/'resume-cli-verification.json',result);print('pass');return 0

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--verify-existing',action='store_true');a=p.parse_args();r=a.root.resolve();
 if a.verify_existing:return verify_existing(r)
 out=r/'resume-cli-regression';out.mkdir(exist_ok=False);grid=r/'preservation-spatial';grid.mkdir(exist_ok=False);write_json(grid/'evaluation-grid.json',dict(grid=[[0.,0.,.5],[0.,0.,3.]],weights=[4.,4.]))
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
  code=run('resume',['--resume']);assert code in (0,1,2),code
  if code==1:
   d=load(campaign/'results.json');failed=[n for t in d['temporal_tables'] for n,v in t['fields'].items() if v['status']=='fail'];assert d['status']=='fail' and failed==['Pi'],failed
   checks.append(dict(name='required-numerical-failure-propagates',status='pass',exit_code=code,intended_failed_field='Pi',numerical_campaign_status='fail'))
  assert all((sha(first/x), (first/x).stat().st_mtime_ns)==v for x,v in saved.items()),'completed trajectory was mutated'
  complete=load(campaign/'results.json');assert len(complete['checks'])==3 and all(x['status']=='pass' for x in complete['checks'])
  checks.append(dict(name='identical-resume-only-missing-native-work',status='pass',completed_numerical_files_unchanged=True,new_native_cases=2,exit_code=code,accuracy_scope='request/reuse contract only; two-point evaluation is not pilot evidence'))
  # Runtime mismatches are separately exercised by the manifest unit test with
  # the exact runtime identity object. No frozen executable is overwritten.
  status='pass'
 except Exception as e:checks.append(dict(name='native-resume-contract',status='fail',failure=repr(e)));status='fail'
 write_json(out/'results.json',dict(status=status,checks=checks));print(status);return 0 if status=='pass' else 1
if __name__=='__main__':sys.exit(main())
