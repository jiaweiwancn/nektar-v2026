import fcntl, hashlib, json, os
from pathlib import Path
import signal, subprocess, sys, time
from datetime import datetime, timezone

M8=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8')
R1=M8.parent/'m8-r1'
SOURCE=R1/'archive-source'
HIST=M8/'initialization'
BUILD=HIST/'build'
TOOLS=R1/'preflight-tools-002'
OUT=R1/'preflight-completion-001'
INITIAL=R1/'preflight-001'
OUT.mkdir(exist_ok=False)
lock=(HIST/'validation.lock').open('a')
fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
THREADS={key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
ENV=dict(os.environ,**THREADS,PYTHONDONTWRITEBYTECODE='1',OMPI_MCA_hwloc_base_report_bindings='1')
ENV.pop('OMPI_MCA_hwloc_base_cpu_list',None)
ENV['LD_LIBRARY_PATH']=str(BUILD/'geometry-lib')
ENV['PYTHONPATH']=str(SOURCE/'solvers/IncNavierStokesSolver/Tests/RANS')
ENV['TMPDIR']=str(OUT/'tmp')
Path(ENV['TMPDIR']).mkdir()
os.environ.update(ENV)
os.sched_setaffinity(0,{49})
original=M8/'flow-001/qualification/control-coarse/case.xml'
original_raw=original.with_name('case.evolution.json')
fingerprint=M8.parent/'transfers/fingerprint-runtime.py'
program=BUILD/'IncNavierStokesSolver'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def git(source,*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
assert git(SOURCE,'rev-parse','HEAD')=='ab5c3c6fee19761e8d88692e8e6c331e5889234e' and not git(SOURCE,'status','--porcelain')
tool_pins={str(path):sha(path) for path in TOOLS.iterdir() if path.is_file()}
input_pins={str(path):sha(path) for path in (original,original_raw,fingerprint)}
solver_sha=sha(program)
receipt=dict(schema=1,phase='C',status='running',started_utc=datetime.now(timezone.utc).isoformat(),
    reviewed_source=dict(path=str(SOURCE),commit=git(SOURCE,'rev-parse','HEAD'),clean=True),
    actual_build_source=dict(path=str(HIST/'source'),commit=git(HIST/'source','rev-parse','HEAD')),
    solver=dict(path=str(program),sha256=solver_sha),tool_pins=tool_pins,input_pins=input_pins,
    numerical_threads=THREADS,affinity=sorted(os.sched_getaffinity(0)),maximum_workers=1,
    environment={key:ENV.get(key) for key in (*THREADS,'LD_LIBRARY_PATH','PYTHONPATH','TMPDIR','OMPI_MCA_hwloc_base_report_bindings')},
    source_quad_levels=[[8,16],[12,24]],native_GLL_counts=[8,8],
    trajectories=0,zero_step_solver_initializations=4,commands=[],cases=[])

def save():(OUT/'runner.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
def command(name,argv,cwd,timeout):
    argv=list(map(str,argv));path=OUT/(name+'.log');start=time.monotonic()
    record=dict(name=name,command=argv,cwd=str(cwd),timeout_seconds=timeout)
    record['command_cwd_sha256']=hashlib.sha256(json.dumps(dict(command=argv,cwd=str(cwd)),sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
    receipt['commands'].append(record);save()
    with path.open('w') as log:
        process=subprocess.Popen(argv,cwd=cwd,env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        record['pid']=process.pid;save();last=start
        while process.poll() is None:
            now=time.monotonic()
            if now-start>timeout:
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(10)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                record['exit_code']=124;break
            if now-last>=30:
                print(name+': alive '+str(round(now-start,1))+' s; '+str(path.stat().st_size)+' log bytes',flush=True);last=now
            time.sleep(1)
        else:record['exit_code']=process.returncode
    record.update(seconds=time.monotonic()-start,log=str(path),log_sha256=sha(path));save()
    print(name+': exit '+str(record['exit_code'])+' ('+str(round(record['seconds'],2))+' s)',flush=True)
    if record['exit_code']:print(path.read_text()[-2500:],flush=True)
    return record['exit_code']

save()
try:
    if command('fingerprint-before',['python3','-B',fingerprint,'--root',HIST,'--output',OUT/'runtime-before.json'],OUT,180):
        raise RuntimeError('Preflight runtime/source fingerprint failed')
    before=json.loads((OUT/'runtime-before.json').read_text());assert before['success']
    cpp={name:digest for name,digest in before['source_inputs'].items()
         if Path(name).suffix.lower() in ('.cpp','.cc','.c','.cxx','.h','.hpp','.hxx','.inl')}
    assert all((SOURCE/name).is_file() and sha(SOURCE/name)==digest for name,digest in cpp.items())
    receipt['reviewed_ab5_Cxx_identical_to_actual_compiled_ed34']=dict(passed=True,count=len(cpp))
    buildpin=json.loads((HIST/'build-002/runtime-fingerprint.json').read_text())
    for key in ('source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images'):
        assert before[key]==buildpin[key],key
    receipt['build_pin_comparison']=dict(identical_fields=['source_inputs','runtime_artifacts','runtime_bindings','python_packages','tool_images'],
        runtime_environment_differences={key:dict(build=buildpin['runtime_environment'].get(key),preflight=before['runtime_environment'].get(key))
          for key in set(buildpin['runtime_environment'])|set(before['runtime_environment'])
          if buildpin['runtime_environment'].get(key)!=before['runtime_environment'].get(key)},
        full_environment_identity_not_claimed=True)
    prior=json.loads((INITIAL/'runner.json').read_text())
    assert prior['full_runtime_environment_source_images_unchanged'] and prior['tools_inputs_solver_unchanged']
    assert prior['solver']['sha256']==solver_sha and len(prior['cases'])==1 and prior['cases'][0]['exit_code']==0
    plan=INITIAL/'cases/preflight-plan.json';matrix=json.loads(plan.read_text())
    receipt['prior_attempt']=dict(path=str(INITIAL/'runner.json'),sha256=sha(INITIAL/'runner.json'),
        purpose='Preserve successful original-base initialization; rerun only its corrected independent metadata/numerical audit. Execute the three missing zero-step solver cases fresh.',
        original_solver_case_count=1,new_solver_case_count=3,total_zero_step_solver_initializations=4,
        original_native_report_sha256=prior['cases'][0]['report_sha256'])
    receipt['plan_sha256']=sha(plan);save()
    for case in matrix['cases']:
        xml=Path(case['xml']);directory=xml.parent;name=case['name']
        row=dict(name=name,xml=str(xml),xml_sha256=sha(xml),cost=case['cost'],steps=0,source_quadrature=case['source_quadrature'])
        receipt['cases'].append(row);save()
        raw_path=directory/'case.evolution.json'
        if name=='original-base':
            assert sha(xml)==prior['cases'][0]['xml_sha256'] and sha(raw_path)==prior['cases'][0]['report_sha256']
            code=0;row.update(execution_mode='existing_successful_zero_step_report_independently_reaudited',
                original_runner=str(INITIAL/'runner.json'),original_solver_seconds=next(item['seconds'] for item in prior['commands'] if item['name']=='solve-original-base'))
        else:
            assert not raw_path.exists()
            code=command('solve-'+name,[program,xml],directory,900)
            row['execution_mode']='fresh_zero_step_native_initialization'
        row['exit_code']=code;save()
        if raw_path.exists():row.update(report=str(raw_path),report_sha256=sha(raw_path))
        if code:raise RuntimeError('Zero-step initialization failed: '+name)
        raw=json.loads(raw_path.read_text())
        assert raw['time']==0 and len(raw['accepted_states'])==1 and not raw['step_audits'] and not raw.get('failure')
        assert raw['offsets'][-1]+1==case['cost']['unknowns']
        output=OUT/(name+'-audit.json')
        code=command('audit-'+name,['python3','-B',TOOLS/'cls_epawf_m8_r1_preflight.py',
            '--raw',raw_path,'--session',xml,'--source',SOURCE,'--plan',plan,'--output',output],directory,600)
        if output.exists():row.update(audit=str(output),audit_sha256=sha(output))
        if code:raise RuntimeError('Independent zero-step audit failed: '+name)
        audited=json.loads(output.read_text());answer=audited['result']
        row.update(status='implementation_audited',finite_admissibility=answer['initial_finite_admissibility'],
            outer_k_over_K0=answer['finite_ranges']['outer_k_over_K0'],
            outer_epsilon_over_epsilon0=answer['finite_ranges']['outer_epsilon_over_epsilon0'],
            outer_maximum_gradient_entry_deviation=answer['finite_ranges']['outer_maximum_gradient_entry_deviation'])
        save()
    receipt['status']='passed';save()
except Exception as error:
    receipt.update(status='failed',failure=dict(type=type(error).__name__,message=str(error)));save()
finally:
    code=command('fingerprint-after',['python3','-B',fingerprint,'--root',HIST,'--output',OUT/'runtime-after.json','--compare',OUT/'runtime-before.json'],OUT,180)
    after=json.loads((OUT/'runtime-after.json').read_text()) if (OUT/'runtime-after.json').exists() else {}
    receipt['full_runtime_environment_source_images_unchanged']=code==0 and after.get('success') is True and after.get('identical_qualification_inputs_and_runtime') is True
    receipt['reviewed_source_clean_and_unchanged']=git(SOURCE,'rev-parse','HEAD')==receipt['reviewed_source']['commit'] and not git(SOURCE,'status','--porcelain')
    receipt['tools_inputs_solver_unchanged']=sha(program)==solver_sha and all(sha(path)==value for path,value in dict(tool_pins,**input_pins).items())
    if not all(receipt[key] for key in ('full_runtime_environment_source_images_unchanged','reviewed_source_clean_and_unchanged','tools_inputs_solver_unchanged')):
        receipt['status']='failed'
    receipt['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    print(json.dumps(dict(status=receipt['status'],failure=receipt.get('failure'),cases=receipt['cases'],output=str(OUT))),flush=True)
sys.exit(0 if receipt['status']=='passed' else 1)
