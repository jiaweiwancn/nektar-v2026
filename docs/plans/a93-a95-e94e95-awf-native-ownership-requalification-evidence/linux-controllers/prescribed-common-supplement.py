import pathlib,sys,subprocess,os,time,gzip,json
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_worker_budget import spawn
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load
from rans_request_manifest import sha
from analyze_rans_duration_constitutive import analyze
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
events=[];checks=[]
def run(name,cmd):
 with (r/(name+'.log')).open('x') as f:
  proc=None
  while proc is None:
   proc=spawn(r,['python3','-B',*map(str,cmd)],env=dict(os.environ),stdout=f,stderr=subprocess.STDOUT)
   if proc is None:time.sleep(3)
  code=proc.wait()
 events.append(dict(name=name,exit_code=code,command=list(map(str,cmd))));return code
for model,wall in (('A95','Accepted'),('A93','EP')):
 tag=model+'-P-'+wall;screen=model+'-Prescribed-'+wall;specfile=r/'prescribed-meshes'/tag/'screen.json'
 write_json(specfile,dict(Balanced11V1=load(specfile.parent/'B11.json')))
 code=run('prescribed-v2-initial-'+tag,[s/'diagnose_rans_balanced_spatial.py','--root',r,'--analysis-runtime',r/'runtime','--model',model,'--policy','Prescribed','--wall',wall,'--label',screen,'--specifications',specfile])
 dest=r/'evolved-audits'/('prescribed-v2-'+tag+'-B11')
 if code==0:code=run('prescribed-v2-audit-'+tag,[s/'audit_rans_balanced_evolved.py','--root',r,'--analysis-runtime',r/'runtime','--campaign-root',r,'--label',tag+'-B11','--output-label',dest.name,'--times','.0016','--common-directory',r/'diagnosis'/screen/'common','--constitutive-samples'])
 value=dict(model=model,policy='Prescribed',wall=wall,space='Balanced11V1',time=.0016,status='pass' if code==0 else 'fail',actual_history='existing fresh corrected trajectory; sealed original request unchanged')
 if code==0:
  q=load(r/'campaigns'/(tag+'-B11')/'request-manifest.json')['request']
  with gzip.open(dest/'T-0.0016/packet.json.gz','rt') as f:packet=json.load(f)
  stress=analyze(packet,{k:q[k] for k in ('model','policy','wall','descriptor')},q['evaluation']['weights'])
  value.update(constitutive=stress,source_request_sha256=sha(r/'campaigns'/(tag+'-B11')/'request-manifest.json'),source_native_report_sha256=sha(r/'campaigns'/(tag+'-B11')/'T-0.0016/case.evolution.json'),native_audit=load(dest/'results.json'))
  if stress['status']!='pass':value['status']='fail'
 checks.append(value)
write_json(r/'selected-prescribed-results-v2.json',dict(status='pass' if len(checks)==2 and all(c['status']=='pass' for c in checks) else 'fail',checks=checks,execution=events,preserved_failed_receipt_sha256=sha(r/'selected-prescribed-results.json'),correction='explicit model/policy/wall active-break common space; original default common space was not nested and was correctly rejected',scope='two selected .0016 fine-space prescribed-kappa controls; no prescribed duration/four-member claim'))
