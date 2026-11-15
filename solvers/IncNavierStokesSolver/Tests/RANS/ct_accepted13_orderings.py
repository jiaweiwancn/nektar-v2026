#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed COLAMD/AMD/Natural assessment using the actual guarded native factor."""
import argparse,os,re,resource,shutil,subprocess,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import run,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def native_run(command,folder,log,cpu,wall):
    """Bound the native process itself, observe its actual maps, always reap it."""
    started=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);maps={};rss=vm=0;stop=None
    def limits():
        resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu+5))
    with (folder/log).open('x') as f:
        child=subprocess.Popen([str(x) for x in command],cwd=folder,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,**THREADS),preexec_fn=limits)
        try:
            while child.poll() is None:
                proc=Path('/proc')/str(child.pid)
                try:
                    for line in (proc/'maps').read_text().splitlines():
                        word=line.split()[-1];p=Path(word)
                        if word.startswith('/') and p.is_file() and word not in maps:maps[word]=sha(p)
                    for line in (proc/'status').read_text().splitlines():
                        if line.startswith('VmHWM:'):rss=max(rss,int(line.split()[1]))
                        if line.startswith('VmPeak:'):vm=max(vm,int(line.split()[1]))
                except (FileNotFoundError,ProcessLookupError):pass
                if time.monotonic()-started>wall:stop='wall_work_limit';child.terminate();break
                time.sleep(.05)
        finally:
            try:child.wait(timeout=10)
            except subprocess.TimeoutExpired:child.kill();child.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(command=list(map(str,command)),exit_code=child.returncode,wall_seconds=time.monotonic()-started,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,peak_RSS_KiB=rss,peak_address_space_KiB=vm,actual_mapped_files=maps,CPU_soft=cpu,CPU_hard=cpu+5,wall_limit=wall,stop_reason=stop)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('exports','runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--reuse-progress',type=Path)
    p.add_argument('--work-amendment',type=Path,required=True)
    p.add_argument('--independent-natural-completion',action='store_true')
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    identity=runtime_identity(a.runtime/'CT13SparseOrderingProbe');data=sealed(a.exports)
    cpu_limit=900 if a.independent_natural_completion else 300
    wall_limit=1000 if a.independent_natural_completion else 900
    if a.independent_natural_completion:
        require(__import__('json').loads(a.work_amendment.read_text())['schema']=='CT13IndependentNaturalCompletionV2','missing independent work contract')
    reused={};work_limited={}
    if a.reuse_progress:
        request=__import__('json').loads((a.reuse_progress/'request-manifest.json').read_text())['request']
        require(request['exports']==sha(a.exports/'results.json') and request['runtime']==identity,'foreign partial ordering assessment')
        reused=__import__('json').loads((a.reuse_progress/'progress.json').read_text())['results']
        closure=a.reuse_progress.parent/('resources-'+a.reuse_progress.name)/'execution.json'
        if not a.independent_natural_completion and closure.exists() and __import__('json').loads(closure.read_text())['exit_code']==124:
            witness=a.reuse_progress/'live-work-witness.json'
            for row in __import__('json').loads(witness.read_text())['owned_active_native']:
                if row['ordering']=='Natural' and row['command'][1]==data['systems']['original-initial']['A']['file']:
                    require(not (a.reuse_progress/'original-initial-Natural.json').exists(),'completed Natural result cannot be relabelled')
                    work_limited['Natural']=dict(status='unresolved_work_limited',reason='original1800-second family bound expired during first Natural control; no positive fill/residual claim or repeat',closed_execution_sha256=sha(closure),last_native_witness_sha256=sha(witness),last_native_CPU_lower_bound=row['CPU'])
    m=RequestManifest(a.output,dict(schema='CT13FixedOrderingAssessmentV3',source=sha(Path(__file__)),exports=sha(a.exports/'results.json'),contract=sha(a.contract),runtime=identity,orderings=['COLAMD','AMD','Natural'],work_amendment=sha(a.work_amendment),reused_progress=sha(a.reuse_progress/'progress.json') if a.reuse_progress else None,per_factor_CPU_soft=cpu_limit,per_factor_CPU_hard=cpu_limit+5,per_factor_wall=wall_limit))
    results={};fill={k:[] for k in ('COLAMD','AMD','Natural')};ineligible=set(work_limited)
    for label,s in data['systems'].items():
        require(sha(Path(s['A']['file']))==s['A']['sha256'],'foreign canonical matrix')
        rhs=a.output/(label+'.rhs');rhs.write_text('\n'.join(format(x,'.17g') for x in s['rhs'])+'\n');results[label]={};control=None
        for ordering in fill:
            if ordering in ineligible and ordering!='COLAMD':
                results[label][ordering]=work_limited[ordering] if label=='original-initial' and ordering in work_limited else dict(status='not_attempted_after_policy_guard',reason='same fixed alternative already ineligible on a required system')
                fill[ordering].append(None);continue
            out=a.output/(label+'-'+ordering+'.json')
            prior=reused.get(label,{}).get(ordering)
            if prior and not prior.get('result_sha256'):prior=None
            if prior:
                src=a.reuse_progress/out.name
                require(src.is_file() and sha(src)==prior['result_sha256'] and prior['matrix_sha256']==s['A']['sha256'] and sha(rhs)==sha(a.reuse_progress/rhs.name),'changed completed ordering primitive')
                shutil.copy2(src,out);e=prior['execution']
            else:e=native_run([a.runtime/'CT13SparseOrderingProbe',s['A']['file'],rhs,ordering,out],a.output,label+'-'+ordering+'.log',cpu_limit,wall_limit)
            row=__import__('json').loads(out.read_text()) if out.is_file() else dict(status='native_work_or_resource_failure',failure='no factor result; bounded native exit; inspect log and signal')
            compact={k:v for k,v in row.items() if k not in ('solution','tiny_solution')};compact.update(execution=e,result_sha256=sha(out) if out.exists() else None,matrix_sha256=s['A']['sha256'])
            if prior:compact['reused_completed_output']=str(a.reuse_progress/out.name)
            if e['exit_code']==0:
                require(row['factor_nnz']<=40000000 and row['true_residual']<=2e-12*row['rhs_maximum']+1e-25 and row['tiny_true_residual']<=2e-12*row['tiny_rhs_maximum']+1e-25,'pass label violates original factor/accuracy guards')
                compact['fill_margin']=40000000-row['factor_nnz'];fill[ordering].append(row['factor_nnz'])
                if ordering=='COLAMD':control=np.asarray(row['solution'])
                elif control is not None:
                    compact['maximum_solution_difference']=float(np.max(abs(np.asarray(row['solution'])-control),initial=0.))
                    require(compact['maximum_solution_difference']<=1e-8,'predictor solution disagreement')
            else:
                fill[ordering].append(None);match=re.search(r'factor_nnz=(\d+)',row.get('failure',''))
                if match:compact['observed_factor_nnz']=int(match[1])
                if ordering!='COLAMD':ineligible.add(ordering)
            results[label][ordering]=compact;write_json(a.output/'progress.json',dict(results=results,fill=fill))
            print(label,ordering,compact['status'],compact.get('factor_nnz',compact.get('observed_factor_nnz')),flush=True)
    eligible=[k for k in ('AMD','Natural') if len(fill[k])==len(data['systems']) and all(v is not None and v<=40000000 for v in fill[k])]
    winner=min(eligible,key=lambda k:max(fill[k])) if eligible else None
    result=dict(schema='CT13FixedOrderingAssessmentV1',status='ordering_selected' if winner else 'blocked_capacity',selected=winner,results=results,fill=fill,selection='one fixed alternative with smallest maximum fill across all prospectively declared systems; original guards/residual/refinement unchanged',governing_component='not_enabled_by_read_only_assessment')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print(result['status'],winner,flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
