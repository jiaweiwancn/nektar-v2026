#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Selected prescribed-kappa, current-gradient cache and ordinary CLS controls."""
import argparse,os,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash,coverage
from rans_request_manifest import RequestManifest,runtime_identity
from qualify_rans_preservation import load
from continue_cls_epawf_pimple import header

def checkpoint(folder,steps):
 return next(p for p in folder.glob('checkpoint*') if p.is_dir() and (p/'State.bin').exists() and header(p)['steps']==steps)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--cls-baseline',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args();r=a.root.resolve();old=a.cls_baseline.resolve();out=r/'selected-controls';os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'));source=Path(__file__).parent
 specs=[('E94E95','EP'),('E94E95','Accepted'),('A95','EP'),('A93','Accepted')];cp2=checkpoint(old,2);cp4=checkpoint(old,4)
 req=dict(schema='SelectedPilotControlsV1',runtime={n:runtime_identity(r/'runtime'/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},prescribed=specs,dt=.0004,horizon=.0016,cls_prefix=header(cp2),cls_expected=header(cp4),cache='A95 full current-gradient p5')
 admit=RequestManifest(out,req,a.resume);result=load(out/'results.json') if a.resume and (out/'results.json').exists() else dict(status='running',checks=[])
 def control(name,fn):
  if a.resume and any(c['name']==name and c['status']=='pass' for c in result['checks']):return
  try:d=fn();row=dict(name=name,status='pass',result=d)
  except Exception as e:row=dict(name=name,status='fail',failure=repr(e))
  result['checks']=[x for x in result['checks'] if x['name']!=name]+[row];write_json(out/'results.json',result)
 def run(label,model,policy,wall,cache='ExactStepV1'):
  command=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model',model,'--policy',policy,'--wall',wall,'--dt','.0004','--cache',cache]
  if (r/'campaigns'/label).exists():command+=['--resume']
  with (out/(label+'.log')).open('w') as f:c=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT)
  m=load(r/'campaigns'/label/'T-0.0016/metrics.json');assert c.returncode==0 and m['status']=='pass',(c.returncode,m.get('failure'))
  return dict(label=label,metrics=m,request=load(r/'campaigns'/label/'request-manifest.json'))
 for model,wall in specs:
  name=model+'-Prescribed-'+wall
  control(name,lambda m=model,w=wall:run(('E' if m=='E94E95' else m)+'-P-'+w+'-h0p5',m,'Prescribed',w))
 def cache():
  new=run('A95-N-EP-h0p5-cache-off','A95','Natural','EP','Off');a=load(r/'campaigns/A95-N-EP-h0p5/T-0.0016/case.evolution.json');b=load(r/'campaigns/A95-N-EP-h0p5-cache-off/T-0.0016/case.evolution.json')
  assert [exact_vector_hash(x) for x in a['accepted_states']]==[exact_vector_hash(x) for x in b['accepted_states']]
  assert [exact_vector_hash(x) for x in a['native_history']]==[exact_vector_hash(x) for x in b['native_history']] and a['native_history_times']==b['native_history_times']
  assert a['state_audits']==b['state_audits'],'cache diagnostics changed'
  return dict(status='pass',state_history_and_audit_bitwise_equal=True,source_labels=['A95-N-EP-h0p5','A95-N-EP-h0p5-cache-off'])
 control('A95-p5-current-gradient-cache',cache)
 def cls():
  dest=out/'cls-reviewed-restart2';command=[sys.executable,'-B',str(source/'run_cls_epawf_pimple.py'),'--solver',str(r/'runtime/IncNavierStokesSolver'),'--output',str(dest),'--steps','2','--dt','.000625','--diagnostic','--wall-policy','UniformEPFlat','--checkpoint','--checkpoint-interval','1','--restart',str(cp2)]
  with (out/'cls-restart.log').open('w') as f:c=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT)
  assert c.returncode==0,c.returncode;expected=header(cp4,True);actual=header(checkpoint(dest,4),True)
  assert exact_vector_hash(expected['coupled_state'])==exact_vector_hash(actual['coupled_state'])
  assert [exact_vector_hash(x) for x in expected['history']]==[exact_vector_hash(x) for x in actual['history']] and expected['history_times']==actual['history_times']
  return dict(status='pass',bitwise_state_history_equal=True,baseline=str(old),reference_source='99112a1cbb9c5448dd3de69f50600bf91465062b',native_continuation=actual['sha256'])
 control('ordinary-reviewed-CLS-restart',cls)
 expected={m+'-Prescribed-'+w for m,w in specs}|{'A95-p5-current-gradient-cache','ordinary-reviewed-CLS-restart'};result['coverage']=coverage(result['checks'],expected,key=lambda x:x['name']);result['status']='pass' if result['coverage']['status']=='pass' and all(x['status']=='pass' for x in result['checks']) else 'fail';write_json(out/'results.json',result);print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
