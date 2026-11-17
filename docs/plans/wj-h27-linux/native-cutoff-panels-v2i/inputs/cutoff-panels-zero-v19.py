from pathlib import Path
import os,sys,time,json,subprocess,shutil,xml.etree.ElementTree as ET
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');t=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27';sys.path.insert(0,str(t));sys.path.insert(0,str(r/'tools'))
from resources import identity,atomic
import numpy as np
import physical_cache_history_check as history
import history_recovery_check as hr
o=r/'science/cutoff-panels-v19';old=o/'reference-r2.result.json';packet=history.read(old);a=packet['state'];checkpoint=o/'reference-r2.final.bin'
assert packet['steady_mode']is True and packet['steady_solver']['converged']is True and a['step']==a['time']==0
before=identity(checkpoint);original=ET.parse(o/'reference-r2.xml');x=original.getroot();si=x.find('.//SOLVERINFO');opts={i.get('PROPERTY'):i for i in si}
opts['WJNativeSolve'].set('VALUE','Advance');opts['WJOutputPrefix'].set('VALUE',str(o/'reference-r2-reload'))
for key in ('WJInitialFields','WJInitialSourceIdentity'):si.remove(opts[key])
ET.SubElement(si,'I',PROPERTY='WJRestart',VALUE=str(checkpoint));case=o/'reference-r2-reload.xml';assert not case.exists();original.write(case,encoding='unicode');shutil.copy2(str(o/'reference-r2.xml')+'.opt',str(case)+'.opt')
atomic(o/'radial-zero-input.json',dict(checkpoint=before,session=identity(case),options=identity(Path(str(case)+'.opt')),source=identity(r/'source-v19.json'),TASK3_READY=False))
subprocess.run([str(r/'runtime-v19/IncNavierStokesSolver'),'--use-opt-file',str(case)+'.opt',str(case)],check=True)
out=o/'reference-r2-reload.result.json';new=history.read(out);b=new['state'];deadline=time.monotonic()+8
for key in ('signature','source_identity','reference_identity','scalar_quadrature_policy','cutoff_radius_range','cutoff_log_slope_bounds','step','time','dt','sizes'):assert a[key]==b[key],key
assert b['scalar_quadrature_policy']=='CubicCutoffFacePanelsV1'and new['restart']is True and new['steady_mode']is False
assert not new['complete_step_seconds']and new['complete_advancing_seconds']==0 and not(o/'reference-r2-reload.steps.jsonl').read_text().strip()
cache=b['coupled_newton'];assert not cache['preconditioner_events']and cache['refreshes']==cache['reuses']==b['rebuilds']==0
aa,aid=history.binary_history(checkpoint,a,a['dt'],deadline,before['sha256']);bb,bid=history.binary_history(o/'reference-r2-reload.final.bin',b,b['dt'],deadline)
assert aid['sha256']==bid['sha256']and identity(checkpoint)['sha256']==before['sha256']
for name in aa:
 for key in aa[name]:assert np.array_equal(aa[name][key].view('u8'),bb[name][key].view('u8')),(name,key)
rows=json.loads((o/'operator-native.json').read_text());assert len(rows['cells'])==2 and all(c['points']>8*12 for c in rows['cells'])
atomic(o/'radial-zero-independent.json',dict(schema='WJH27CutoffPanelColdZeroV1',status='PASS_COLD_RADIAL_ZERO_HISTORY_ONLY',input=before,output=bid,
 result=identity(out),four_binary_histories_exact=True,advanced_steps=0,inverse_factorizations=0,source_identity=a['source_identity'],scalar_quadrature_policy=b['scalar_quadrature_policy'],
 derived_data_scope='Fresh constructor face levels and ReadHistory Evaluate for current/previous cutoff rules; no advancing residual or physical-time qualification claim',
 full_H27_qualification=False,TASK3_READY=False))
print('COLD_RADIAL_ZERO_HISTORY_PASSED')
