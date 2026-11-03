import hashlib, json, os, pathlib, signal, subprocess, sys, time, xml.etree.ElementTree as ET
from datetime import datetime, timezone

ROOT = pathlib.Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
SOURCE, BUILD = ROOT/'verification/source', ROOT/'verification/build'
OUT = ROOT/'verification/curved-001'
OUT.mkdir(exist_ok=False)
THREADS = {key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
ENV = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1')
ENV['TMPDIR'] = str(ROOT.parent/'tmp')
ENV['LD_LIBRARY_PATH'] = str(ROOT/'build/geometry-lib')+':'+ENV.get('LD_LIBRARY_PATH','')

def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

def git(path,*args):
    return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()

revision=git(SOURCE,'rev-parse','HEAD')
assert revision=='6284787bf4c436203f4404498ba3dfb68aa144ef' and not git(SOURCE,'status','--porcelain')
primary_revision=git(ROOT/'source','rev-parse','HEAD')
primary_solver=ROOT/'build/IncNavierStokesSolver'
primary_solver_sha=sha(primary_solver)
geometry=ROOT/'build/geometry-lib/libSpatialDomains.so.5.9.0'
geometry_sha=sha(geometry)
cpp=[p for p in git(SOURCE,'diff','--name-only','aa854f4c5b05451a81a4e1dd60eaaca98eb547e4',revision).splitlines() if pathlib.Path(p).suffix in ('.cpp','.h','.hpp','.c','.cc','.cxx','.hxx','.inl')]
assert cpp==['library/UnitTests/SolverUtils/RANSCurvedWallProbe.cpp'],cpp
receipt=dict(schema=1,status='running',started_utc=datetime.now(timezone.utc).isoformat(),source_commit=revision,
             source_status='',primary_source_commit=primary_revision,primary_solver_sha256=primary_solver_sha,
             frozen_geometry=dict(path=str(geometry),sha256=geometry_sha),CPP_difference_from_native_cohort=cpp,
             environment={k:ENV.get(k) for k in (*THREADS,'TMPDIR','LD_LIBRARY_PATH')},
             compiler_workers=16,numerical_workers=1,numerical_cpu=50,commands=[],cases=[],
             scope='Selected all-mode curved operator and original controls; no genuine curved gamma MMS or robust/physical terrain qualification',
             runner_sha256=sha(__file__),CMake_input_sha256=sha(SOURCE/'solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt'))

def save():
    (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

def command(name,args,cwd,timeout):
    path=OUT/(name+'.log');start=time.monotonic()
    row=dict(name=name,command=list(map(str,args)),cwd=str(cwd),timeout_seconds=timeout)
    receipt['commands'].append(row);save()
    with path.open('w') as stream:
        process=subprocess.Popen(row['command'],cwd=cwd,env=ENV,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        row['pid']=process.pid;save();last=start
        while process.poll() is None:
            now=time.monotonic()
            if now-start>timeout:
                print(name+': hard timeout; terminating process group',flush=True)
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL);process.wait()
                row['exit_code']=124;break
            if now-last>=30:
                print(name+': alive, '+str(round(now-start,1))+'s; log bytes='+str(path.stat().st_size),flush=True);last=now
            time.sleep(1)
        else:
            row['exit_code']=process.returncode
    row.update(seconds=time.monotonic()-start,log=str(path),log_sha256=sha(path));save()
    print(name+': exit '+str(row['exit_code'])+' ('+str(round(row['seconds'],2))+'s)',flush=True)
    if row['exit_code']:
        print(path.read_text()[-2500:],flush=True)
    return row['exit_code']

save()
try:
    if command('configure',['cmake','-S',SOURCE/'solvers/IncNavierStokesSolver/Tests/RANS','-B',BUILD,
               '-DRANS_BUILD_SOLVER=OFF','-DCMAKE_BUILD_TYPE=Release',
               '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],OUT,180):
        raise RuntimeError('Configure failed')
    if command('build',['cmake','--build',BUILD,'--target','RANSCurvedWallProbe','--parallel','16'],OUT,900):
        raise RuntimeError('Isolated curved probe build failed')
    probe=BUILD/'RANSCurvedWallProbe'
    receipt['probe']=dict(path=str(probe),sha256=sha(probe))
    dependencies=subprocess.check_output(['ldd',str(probe)],env=ENV,text=True)
    (OUT/'probe-ldd.txt').write_text(dependencies)
    import re
    resolved={}
    for line in dependencies.splitlines():
        if '=>' in line:
            name,path=(s.strip() for s in line.split('=>',1));path=re.sub(r'\s+\(0x[0-9a-fA-F]+\)\s*$','',path)
            if path.startswith('/'):resolved[name]=dict(path=path,sha256=sha(path))
        elif line.strip().startswith('/'):
            path=re.sub(r'\s+\(0x[0-9a-fA-F]+\)\s*$','',line.strip());resolved['ELF interpreter']=dict(path=path,sha256=sha(path))
    assert 'not found' not in dependencies
    spatial=[v for k,v in resolved.items() if 'SpatialDomains' in k]
    assert len(spatial)==1 and pathlib.Path(spatial[0]['path']).resolve()==geometry.resolve() and spatial[0]['sha256']==geometry_sha
    receipt['probe_loader_bindings']=resolved;save()
    os.sched_setaffinity(0,{50})
    scripts=SOURCE/'solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(scripts))
    ENV['PYTHONPATH']=str(scripts)+(':'+ENV['PYTHONPATH'] if ENV.get('PYTHONPATH') else '')
    receipt['environment']['PYTHONPATH']=ENV['PYTHONPATH']
    from cls_epawf_m5_cases import mapped_session
    from cls_epawf_m6_cases import transition_session
    from cls_epawf_m8_cases import apply_variant
    from cls_epawf_m3_cases import write
    matrix=[
      ('original-variable','Original','variable',2,False,'m5'),
      ('scaled-variable','CubicScaled','variable',2,False,None),
      ('regularized-variable','Regularized','variable',2,False,None),
      ('original-transition','Original','feet',2,True,'m6'),
      ('regularized-transition','Regularized','feet',2,True,None),
      ('regularized-surface3d','Regularized','surface3d',3,False,None)]
    receipt['frozen_case_matrix']=matrix;save()
    legacy=OUT/'legacy-original-audit.py'
    legacy.write_text('''import json,pathlib,sys\nimport validate_cls_epawf_m5 as weak\nfrom cls_epawf_m7_energy_oracle import Audit\nfrom cls_epawf_m6_operator_oracle import verify\nraw=json.loads(pathlib.Path(sys.argv[1]).read_text());out=pathlib.Path(sys.argv[2]);mode=sys.argv[3];before=weak.checks;a=Audit()\nif mode=="m5":\n result=dict(weak=weak.residual_oracle(raw),interface=weak.interface_oracle(raw))\nelse:\n result=verify(raw,a.near,a.require,directional=True)\nresult.update(success=True,legacy_control=mode,weak_checks=weak.checks-before,additional_checks=a.checks)\nout.write_text(json.dumps(result,indent=2)+"\\n")\nprint(json.dumps(dict(success=True,legacy_control=mode,weak_checks=weak.checks-before,additional_checks=a.checks)))\n''')
    receipt['legacy_audit_sha256']=sha(legacy)
    for name,variant,family,dimension,upstream,control in matrix:
        directory=OUT/name;directory.mkdir();xml=directory/'case.xml'
        if upstream:
            transition_session(scripts/'EquilibriumWallEvolution.xml',xml,family=family,dimension=dimension,degree=3,
                steps=0,quadrature=8,normal_order=16)
        else:
            mapped_session(scripts/'EquilibriumWallEvolution.xml',xml,family=family,dimension=dimension,degree=3,
                steps=0,quadrature=8)
        root=ET.parse(xml).getroot()
        apply_variant(root,variant,beta=1. if variant=='Original' else .4,gamma_scale=1. if variant=='Regularized' else 0.)
        write(root,xml)
        row=dict(name=name,variant=variant,family=family,dimension=dimension,upstream=upstream,
                 session=str(xml),session_sha256=sha(xml),status='running');receipt['cases'].append(row);save()
        raw=directory/'probe.json'
        if command(name+'-probe',[probe,xml,raw],directory,900):
            row['status']='probe_failed';save();continue
        row['native_report']=str(raw);row['native_report_sha256']=sha(raw)
        if command(name+'-operator',['python3','-B',scripts/'cls_epawf_m8_curved_operator.py','--raw',raw,
                    '--output',directory/'operator-audit.json'],directory,2400):
            row['status']='operator_audit_failed';save();continue
        row['selected_operator_audit']=json.loads((directory/'operator-audit.json').read_text())
        row['status']='passed'
        if control:
            if command(name+'-legacy',['python3','-B',legacy,raw,directory/'legacy-audit.json',control],directory,2400):
                row['status']='legacy_control_failed'
            else:row['legacy_original_control']=json.loads((directory/'legacy-audit.json').read_text())
        save()
    assert sha(probe)==receipt['probe']['sha256']
    for value in resolved.values():assert sha(value['path'])==value['sha256']
    receipt['status']='passed' if all(row['status']=='passed' for row in receipt['cases']) else 'failed'
except Exception as error:
    receipt.update(status='failed',error=str(error))
finally:
    receipt['source_clean_end']=not git(SOURCE,'status','--porcelain')
    receipt['source_unchanged']=git(SOURCE,'rev-parse','HEAD')==revision
    receipt['primary_source_unchanged']=git(ROOT/'source','rev-parse','HEAD')==primary_revision
    receipt['primary_solver_unchanged']=sha(primary_solver)==primary_solver_sha
    receipt['primary_geometry_unchanged']=sha(geometry)==geometry_sha
    receipt['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    print(json.dumps(dict(status=receipt['status'],error=receipt.get('error'),cases=[(x['name'],x['status']) for x in receipt['cases']],output=str(OUT))),flush=True)
sys.exit(0 if receipt['status']=='passed' else 1)
