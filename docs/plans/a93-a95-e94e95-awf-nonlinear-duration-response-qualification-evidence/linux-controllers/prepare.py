from pathlib import Path
import os,subprocess,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';old=r.parent/'a93-a95-e94e95-balanced-20261110'
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),LD_LIBRARY_PATH=str(r/'runtime/lib'))
def run(name,cmd):
 with (r/(name+'.log')).open('x') as f:p=subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,env=env)
 print(name,p.returncode,flush=True)
 return p.returncode
assert run('configure',['cmake','-S',t,'-B',r/'build','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake','-DRANS_BUILD_SOLVER=ON','-DCMAKE_BUILD_TYPE=Release','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(r/'dependency/libSpatialDomains.so.5.9.0'),'-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'])==0
buildlog=(r/'build-probe.log').open('x');build=subprocess.Popen(['cmake','--build',str(r/'build'),'--parallel','2','--target','RANSPreservationProbe'],stdout=buildlog,stderr=subprocess.STDOUT,env=env)
checks=[]
checks.append(dict(name='coverage-contract',exit_code=run('coverage-contract',['python3','-B',t/'test_rans_duration_contract.py','-v'])))
for model in ('A93','A95'):
 tag=model+'-N-Accepted'
 q=run('frozen-quadrature-'+tag,['python3','-B',t/'check_rans_balanced_state_quadrature.py','--root',r,'--campaign-root',old,'--model',model,'--walls','Accepted','--output-label','frozen-quadrature-'+tag,'--time','.016'])
 checks.append(dict(name=tag+'-quadrature',exit_code=q))
 if q==0:checks.append(dict(name=tag+'-admission',exit_code=run('admit-'+tag,['python3','-B',t/'admit_rans_balanced_duration.py','--root',r,'--previous-root',old,'--model',model,'--wall','Accepted'])))
(r/'preparation-results.json').write_text(json.dumps(dict(status='pass' if all(x['exit_code']==0 for x in checks) else 'fail',checks=checks),indent=2)+'\n')
code=build.wait();buildlog.close();(r/'build-result.json').write_text(json.dumps(dict(exit_code=code,status='pass' if code==0 else 'fail'))+'\n');print('build',code,flush=True)
