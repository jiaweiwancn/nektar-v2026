import pathlib,sys,subprocess,os,time,json
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_worker_budget import spawn
from rans_qualification import THREADS,write_json
from rans_preservation_oracle import Reference
from make_rans_balanced_meshes import specifications
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
while not (r/'runtime/source-build.json').exists():time.sleep(3)
events=[]
def run(name,cmd):
 with (r/(name+'.log')).open('x') as f:
  proc=None
  while proc is None:
   proc=spawn(r,cmd,env=dict(os.environ),stdout=f,stderr=subprocess.STDOUT)
   if proc is None:time.sleep(3)
  code=proc.wait()
 events.append(dict(name=name,command=cmd,exit_code=code));write_json(r/'targeted-native-executions.json',dict(execution=events))
 return code
old=r.with_name('a93-a95-e94e95-nonlinear-20261110')
run('compatibility-v1',['python3','-B',str(s/'check_rans_ownership_compatibility.py'),'--root',str(r),'--previous-root',str(old)])
regex='^(RANSNativeDiagnosticChecks|RANSPilotEvidenceContract|RANSDurationCoverageContract|RANSNonlinearEvidenceContract|RANSBudgetAccuracyContract|RANSOwnershipRequiredScopeContract|RANSQualificationFailureContract|RANSForwardInterfaceTolerance|RANSFrozenBoundaryInterfaceRegression|RANSNativeAdmissibilityRegression|RANSNewModelsIndependentOracle)$'
run('targeted-ctest',['ctest','--test-dir',str(r/'build'),'-j','1','--output-on-failure','-R',regex])
checks=[]
for model,wall in (('A95','Accepted'),('A93','EP')):
 tag=model+'-P-'+wall;ref=Reference(model,'Prescribed',wall)
 mesh=r/'prescribed-meshes'/tag;mesh.mkdir(parents=True)
 spec=specifications(ref)['Balanced11V1'];write_json(mesh/'B11.json',spec)
 code=run('prescribed-'+tag,['python3','-B',str(s/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',tag+'-B11','--model',model,'--policy','Prescribed','--wall',wall,'--mesh-file',str(mesh/'B11.json'),'--dt','.0002','--milestones','.0016','--max-time','.0016','--timeout','14400'])
 check=dict(model=model,policy='Prescribed',wall=wall,time=.0016,status='pass' if code==0 else 'fail',space='Balanced11V1',dt=.0002)
 if code==0:
  code=run('prescribed-audit-'+tag,['python3','-B',str(s/'audit_rans_balanced_evolved.py'),'--root',str(r),'--analysis-runtime',str(r/'runtime'),'--campaign-root',str(r),'--label',tag+'-B11','--output-label','prescribed-'+tag+'-B11','--times','.0016','--constitutive-samples'])
  check['native_audit_status']='pass' if code==0 else 'fail'
 checks.append(check)
write_json(r/'selected-prescribed-results.json',dict(status='pass' if all(q['status']=='pass' and q.get('native_audit_status')=='pass' for q in checks) else 'fail',checks=checks,scope='two selected short fine-space controls; no prescribed-kappa duration or four-member qualification'))
