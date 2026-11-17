from pathlib import Path
import json,sys,time,resource

r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');t=r/'tools-terminal-validation-v14/h27';sys.path.insert(0,str(t))
import terminal_output_check as check
from resources import atomic,identity
soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=min([8*1024**3]+[v for v in [soft,hard]if v!=resource.RLIM_INFINITY]);resource.setrlimit(resource.RLIMIT_AS,(cap,hard))
o=r/'science/coupled-terminal-v14';deadline=time.monotonic()+18
control=check.read(o/'aligned.result.json')['state'];ref=check.po.Reference(check.read(r/'science/reference-tables-v1.json'));results={}
for name in ['aligned','initial']:
    spec=dict(result=str(o/(name+'.result.json')),index=str(o/(name+'.science-index.json')),session=str(o/(name+'.xml')))
    results[name]=check.load_run(name,spec,o,control,ref,deadline)
check.require(all(results['aligned']['data'][key]==results['initial']['data'][key]for key in check.HASHES),'MMS output window changed histories')
failed=check.read(o/'failed.science-index.json');failure=check.read(o/'failed.failure.json')
check.require(failed['run_complete']is False and failed['initial_step']==failed['final_step']==0 and failure['state']['step']==0,
    'failed MMS output claims accepted progress')
check.require(len(failed['frames'])==1 and failed['frames'][0]['origin']=='fresh_initial_observation','failed MMS frame falsely accepted')
atomic(o/'mms-scheduling-independent.json',dict(schema='WJH27TerminalMMSSchedulingV1',status='PASS_TERMINAL_MMS_SCHEDULING_ONLY',
    full_H27_qualification=False,source=identity(__file__),checker=identity(t/'terminal_output_check.py'),
    sources={name:value['evidence']for name,value in results.items()},failed_index=identity(o/'failed.science-index.json'),
    failed_run=identity(o/'failed.failure.json'),same_final_history_hashes=True,workers=1,memory_cap_bytes=cap))
print('PASS_TERMINAL_MMS_SCHEDULING_ONLY')
