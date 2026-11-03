import fcntl, hashlib, json, os, pathlib, signal, subprocess, sys, time
from datetime import datetime, timezone

M8 = pathlib.Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
ROOT = M8/'initialization'
SOURCE, BUILD = ROOT/'source', ROOT/'build'
OUT = ROOT/'mms-runner-002'
OUT.mkdir(exist_ok=False)
lock=(ROOT/'validation.lock').open('a')
fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
THREADS={name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
ENV=dict(os.environ,**THREADS,PYTHONDONTWRITEBYTECODE='1',OMPI_MCA_hwloc_base_report_bindings='1')
ENV.pop('OMPI_MCA_hwloc_base_cpu_list',None)
ENV['LD_LIBRARY_PATH']=str(BUILD/'geometry-lib')+':'+ENV.get('LD_LIBRARY_PATH','')
ENV['TMPDIR']=str(ROOT/'tmp')
os.sched_setaffinity(0,{40,41,42,43})
os.environ.update(ENV)

def sha(path):return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
revision=git(SOURCE,'rev-parse','HEAD')
assert revision=='ed34052fe4dd967e963c22272a0ec0888cf35b70' and not git(SOURCE,'status','--porcelain')
main=pathlib.Path('/media/jay/新加卷/nektar-v2026')
primary=M8/'source';primary_head=git(primary,'rev-parse','HEAD');main_head=git(main,'rev-parse','HEAD')
program=BUILD/'IncNavierStokesSolver';program_sha=sha(program)
primary_program=M8/'build/IncNavierStokesSolver';primary_program_sha=sha(primary_program)
fingerprint=M8.parent/'transfers/fingerprint-runtime.py'
old_functions=M8/'mms-completion-001/mms-history'
old_neutral=M8/'native-completion-001/scaled-native/history'
old_reference=M8/'native-completion-001/scaled-native/joint-four/case.evolution.json'
old_files={str(path):sha(path) for path in (old_functions/'State.bin',old_neutral/'State.bin',old_reference)}
receipt=dict(schema=1,status='running',started_utc=datetime.now(timezone.utc).isoformat(),source_commit=revision,
    source_status='',program=dict(path=str(program),sha256=program_sha),runner_sha256=sha(__file__),
    fingerprint_tool_sha256=sha(fingerprint),old_inputs=old_files,maximum_MPI_ranks=2,
    numerical_threads=THREADS,affinity=sorted(os.sched_getaffinity(0)),cpu_binding_policy='Parent Python restricted to CPUs40-43; native driver original --bind-to core --map-by core; MPI may repin its two single-threaded ranks; actual bindings recorded',
    environment={name:ENV.get(name) for name in (*THREADS,'LD_LIBRARY_PATH','TMPDIR','OMPI_MCA_hwloc_base_report_bindings')},
    commands=[],controls=[],fresh_MMS_output=str(ROOT/'mms-002'),reuse=False)

def save():(OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
def command(name,args,cwd,timeout):
    path=OUT/(name+'.log');start=time.monotonic();row=dict(name=name,command=list(map(str,args)),cwd=str(cwd),timeout_seconds=timeout)
    receipt['commands'].append(row);save()
    with path.open('w') as log:
        process=subprocess.Popen(row['command'],cwd=cwd,env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
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
                print(name+': alive, '+str(round(now-start,1))+'s, logbytes='+str(path.stat().st_size),flush=True);last=now
            time.sleep(1)
        else:row['exit_code']=process.returncode
    row.update(seconds=time.monotonic()-start,log=str(path),log_sha256=sha(path));save()
    print(name+': exit '+str(row['exit_code'])+' ('+str(round(row['seconds'],2))+'s)',flush=True)
    if row['exit_code']:print(path.read_text()[-2500:],flush=True)
    return row['exit_code']

save()
try:
    if command('fingerprint-before',['python3','-B',fingerprint,'--root',ROOT,'--output',OUT/'runtime-before.json'],OUT,300):
        raise RuntimeError('Corrected solver pre-run full fingerprint failed')
    fp_before=json.loads((OUT/'runtime-before.json').read_text());assert fp_before['success']
    build_pin=json.loads((ROOT/'build-002/runtime-fingerprint.json').read_text());assert build_pin['success']
    for key in ('source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images'):
        assert fp_before[key]==build_pin[key],key
    receipt['build_fingerprint_equivalent_inputs_and_images']=True;save()
    scripts=SOURCE/'solvers/IncNavierStokesSolver/Tests/RANS'
    if command('fresh-mms',['python3','-B',scripts/'validate_cls_epawf_m8.py','--build',BUILD,
          '--mode','mms','--variant','Regularized','--beta','.4','--gamma-scale','1',
          '--output',ROOT/'mms-002'],OUT,3600):
        raise RuntimeError('Fresh corrected-solver MMS qualification failed')
    mms=json.loads((ROOT/'mms-002/validation.json').read_text())
    assert mms['implementation_verified'] and len(mms['cases'])==7 and not mms.get('replays')
    fresh_checkpoint=json.loads((ROOT/'mms-002/flat-checkpoint/case.evolution.json').read_text())
    fresh_restart=json.loads((ROOT/'mms-002/flat-restart/case.evolution.json').read_text())
    assert ':capped_initial_trace=1' in fresh_checkpoint['signature'] and ':capped_initial_trace=1' in fresh_restart['signature']
    receipt['fresh_MMS']=dict(validation_sha256=sha(ROOT/'mms-002/validation.json'),case_count=7,
        checks=mms['checks'],MMS_field_errors={name:row['field_errors'] for name,row in mms['regularized_mms'].items()},
        MPI_field_norms=mms['mms_mpi_field_norms'],restart_field_norms=mms['mms_restart_field_norms'],
        scoped_marker_present_in_checkpoint_and_restart=True)
    save()
    sys.path.insert(0,str(scripts))
    from cls_epawf_m8_mms import flat_session
    from cls_epawf_m8_cases import perturbation_session
    from validate_cls_epawf_m4 import state_difference,native_history_difference
    from cls_epawf_m8_native_oracle import verify as energy
    from validate_cls_epawf_m8 import physical_gate
    directory=ROOT/'compatibility-002';directory.mkdir(exist_ok=False)
    stale=directory/'old-functions-capped';stale.mkdir()
    xml=flat_session(scripts/'EquilibriumWallEvolution.xml',stale/'case.xml',variant='Regularized',
                     beta=.4,gamma_scale=1.,steps=1,dt=.005,restart=str(old_functions))
    code=command('reject-old-functions-capped',[program,xml],stale,900)
    text=(OUT/'reject-old-functions-capped.log').read_text()
    raw_path=stale/'case.evolution.json';raw=json.loads(raw_path.read_text()) if raw_path.exists() else {}
    assert code!=0 and code!=124 and 'checkpoint signature mismatch' in text
    assert raw.get('accepted_states')==[] and raw.get('initialization_complete') is False
    receipt['controls'].append(dict(name='oldAA_Functions_FixedFlat_Capped',status='expected_rejection_verified',
        XML_sha256=sha(xml),report_sha256=sha(raw_path),expected_reason='checkpoint signature mismatch',
        no_accepted_restart=True,interpretation='Scoped checkpoint versioning control, not a model/trajectory failure'))
    save()
    neutral=directory/'old-neutral-compatible';neutral.mkdir()
    xml=perturbation_session(scripts/'EquilibriumWallEvolution.xml',neutral/'case.xml',variant='CubicScaled',
        beta=.4,gamma_scale=0.,kind='joint',eta=1,dt=.0025,steps=2,restart=str(old_neutral))
    if command('resume-old-neutral',[program,xml],neutral,900):
        raise RuntimeError('Unaffected old Neutral checkpoint failed compatibility')
    raw_path=neutral/'case.evolution.json';raw=json.loads(raw_path.read_text());old=json.loads(old_reference.read_text())
    assert not raw.get('failure') and ':capped_initial_trace=1' not in raw['signature'] and abs(raw['time']-.01)<1e-12
    norms=state_difference(old,raw,'oldAA_Neutral_checkpoint_newsolver_fullfields',tolerance=1e-8)
    native_history_difference(old,raw,'oldAA_Neutral_checkpoint_newsolver_allhistories',tolerance=1e-8)
    answer=energy(raw);(neutral/'energy-audit.json').write_text(json.dumps(answer,indent=2)+'\n')
    receipt['controls'].append(dict(name='oldAA_Neutral_CubicScaled',status='compatibility_verified',
        XML_sha256=sha(xml),report_sha256=sha(raw_path),field_norms=norms,
        full_native_histories_verified=True,independent_energy_checks=answer['checks'],physical_gates=physical_gate(raw),
        checkpoint_marker_absent_as_required=True))
    receipt['status']='passed';save()
except Exception as error:
    receipt.update(status='failed',error=str(error));save()
finally:
    code=command('fingerprint-after',['python3','-B',fingerprint,'--root',ROOT,'--output',OUT/'runtime-after.json',
                 '--compare',OUT/'runtime-before.json'],OUT,300)
    after=json.loads((OUT/'runtime-after.json').read_text()) if (OUT/'runtime-after.json').exists() else {}
    receipt['full_runtime_and_source_unchanged']=code==0 and after.get('success') is True and after.get('identical_qualification_inputs_and_runtime') is True
    receipt['old_inputs_unchanged']=all(sha(path)==expected for path,expected in old_files.items())
    receipt['source_clean_end']=not git(SOURCE,'status','--porcelain')
    receipt['source_revision_unchanged']=git(SOURCE,'rev-parse','HEAD')==revision
    receipt['main_preserved']=git(main,'rev-parse','HEAD')==main_head and not git(main,'status','--porcelain')
    receipt['primary_preserved']=git(primary,'rev-parse','HEAD')==primary_head and sha(primary_program)==primary_program_sha
    receipt['program_unchanged']=sha(program)==program_sha
    if not all(receipt[key] for key in ('full_runtime_and_source_unchanged','old_inputs_unchanged','source_clean_end','source_revision_unchanged','main_preserved','primary_preserved','program_unchanged')):
        receipt['status']='failed'
    receipt['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    print(json.dumps(dict(status=receipt['status'],error=receipt.get('error'),MMS_cases=receipt.get('fresh_MMS',{}).get('case_count'),controls=[(x['name'],x['status']) for x in receipt['controls']],output=str(OUT))),flush=True)
sys.exit(0 if receipt['status']=='passed' else 1)
