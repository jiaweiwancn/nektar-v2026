#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze compact evidence and an indexed manifest beside a live campaign."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from run_cls_epawf_picard import load,write_json,require,sha,state_hash
from run_cls_epawf_pimple import identity
from continue_cls_epawf_pimple import header,process_identity,load_progress
from compare_cls_epawf_pimple import differences
from analyse_cls_epawf_constitutive_v2 import read_report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root.resolve()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    first=load(r/'analysis16/summary.json');last=load(r/'monitor/step-000032/summary.json')
    require(first['status']==last['status']=='analysis_completed','readback is not complete')
    joined=copy.deepcopy(last);joined['states']={k:v for k,v in first['states'].items() if k in ('initial','step8','step16')}
    joined['states']['step32']=last['states']['step32'];joined['windows']=[];joined['window_integrated_budgets']=[]
    points=[('initial',r/'analysis16/initial/samples.json'),('step16',r/'analysis16/step16/samples.json'),
        ('step32',r/'monitor/step-000032/step32/samples.json')]
    samples={name:load(path) for name,path in points};initial=load(r/'analysis16/initial/input.json')['state']
    for v in samples.values():v['initial']=initial
    for x,y in zip(('initial','step16'),('step16','step32')):
        start=samples[x]['time'];end=samples[y]['time'];rows=[v for v in joined['budgets'] if start+1e-13<v['time']<=end+1e-13]
        joined['windows'].append(dict(start=start,end=end,physical_L2=differences(dict(samples[x],time=end),samples[y])))
        physical=[(end-start)*(u+v)/2 for u,v in zip(joined['states'][x]['physical_production_destruction_global_near'],joined['states'][y]['physical_production_destruction_global_near'])]
        exact={k:[sum(v['exact_discrete_global_near_budget'][k][i] for v in rows) for i in range(2)] for k in ('input_load','implicit_increment','bound','actual','closure')}
        joined['window_integrated_budgets'].append(dict(start=start,end=end,
            correction_net=sum(v['signed_integral'] for v in rows),correction_near=sum(v['near_integral'] for v in rows),
            correction_positive=sum(v['positive_integral'] for v in rows),correction_negative=sum(v['negative_integral'] for v in rows),
            correction_absolute=sum(v['absolute_integral'] for v in rows),physical_production_destruction_global_near=physical,
            physical_time_quadrature='trapezoid at retained endpoints; not a closed local CG balance',exact_frozen_discrete=exact))
    legacy_report=load(r.parent/'cls-epawf-pimple-k-bound/fixed-horizon/bound-0/case.evolution.json')
    new_report=read_report(r/'campaign/trajectory/steps-000001-000016')
    same=[i for i,(old,new) in enumerate(zip(legacy_report['accepted_states'],new_report['accepted_states'])) if state_hash(old)==state_hash(new)]
    old12=load(r.parent/'cls-epawf-pimple-longtime-20261106/analysis-final/samples-later/result.json')
    new12=load(r/'analysis16/step12/samples.json');old12['initial']=new12['initial']=initial
    joined['legacy_comparison']=dict(bitwise_identical_available_states=same,
        physical_L2_at_T00075=differences(old12,new12),
        caveat='V1 accepted state twelve violates the V2 inventory; only its retained safe physical samples are compared')
    write_json(a.output/'analysis.json',joined)
    # Capture a real last accepted checkpoint immediately, with all native
    # history bytes. No signatures/history contents are edited or fabricated.
    progress_path=r/'campaign/trajectory/progress.json'
    for attempt in range(3):
        progress=load_progress(progress_path);source=Path(progress['checkpoint']);destination=a.output/f'checkpoint-{attempt}'
        try:
            shutil.copytree(source,destination);cp=header(destination,True)
            require(cp['sha256']==progress['checkpoint_sha256'] and cp['steps']==progress['accepted_steps'],'checkpoint moved during snapshot')
            break
        except (FileNotFoundError,ValueError):
            if attempt==2:raise
            time.sleep(.2)
    checkpoint=dict(original_path=str(source),snapshot=str(destination),sha256=cp['sha256'],step=cp['steps'],time=cp['time'],
        dt=cp['dt'],order=cp['order'],coupled_state_sha256=state_hash(cp['coupled_state']),
        native_history_sha256=[state_hash(v) for v in cp['history']],native_history_times=cp['history_times'])
    histories=load(r/'history/summary.json');fixtures=load(r/'fixtures/summary.json')
    text=(r/'campaign/trajectory/steps.jsonl').read_text();lines=text.splitlines()
    if text and not text.endswith('\n'):lines=lines[:-1]
    records={v['global_step']:v for v in map(json.loads,lines) if v['global_step']<=cp['steps']}
    require(len(records)==cp['steps'],'checkpoint accepted prefix is not contiguous')
    keys=('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves',
        'outer_sweeps','bound_calls','additional_governing_equation_solves','mass_inverse_applications',
        'bound_active_set_iterations','bound_mass_factor_applications','bound_kkt_solves')
    counts={k:sum(v['attempt'][k] for v in records.values()) for k in keys}
    for record in records.values():
        require(record['attempt']['governing_solves_complete'] and all(record['attempt'][k]==1 for k in keys[:6]),'accepted solve schedule changed')
    callback_rows=[]
    for segment in {v['segment'] for v in records.values()}:
        for line in (Path(segment)/'solver.log').read_text().splitlines():
            if line.startswith('EP-AWF PIMPLE '):
                row=json.loads(line[len('EP-AWF PIMPLE '):])
                if row['time']<=cp['time']+1e-13:callback_rows.append(row)
    require(len(callback_rows)==cp['steps'],'hidden callback or repeated sweep in accepted prefix')
    counts['implicit_callbacks_for_accepted_prefix']=len(callback_rows)

    working=r/'source';git=lambda *args:subprocess.check_output(['git','-C',str(working),*args],text=True).strip()
    source_files={}
    for path in subprocess.check_output(['git','-C',str(working),'diff','--name-only','0b86269739075af81622319402370967fdeffbe8','HEAD'],text=True).splitlines():
        f=working/path
        if f.is_file():source_files[path]=sha(f)
    # Index immutable or momentarily stable artifacts. Mutable progress/logs
    # are labelled snapshots; the checkpoint copy above remains recoverable.
    index=[]
    for top in ('validation-runtime','audit-probe','fixtures','history','analysis16','monitor','campaign/trajectory'):
        for path in sorted((r/top).rglob('*')):
            if not path.is_file() or path.name=='.campaign.lock':continue
            try:before=path.stat();digest=sha(path);after=path.stat()
            except FileNotFoundError:continue
            if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):continue
            index.append(dict(path=str(path.relative_to(r)),bytes=after.st_size,sha256=digest,
                snapshot_may_be_mutable=path.suffix in ('.log','.jsonl') or path.name in ('status.json','progress.json'),
                may_be_pruned=top=='campaign/trajectory'))
    write_json(a.output/'artifact-index.json',dict(recorded_unix=time.time(),root=str(r),files=index))
    compact=copy.deepcopy(joined)
    for b in compact['budgets']:
        locations=b.pop('active_constraints');b['active_constraint_count']=len(locations)
        b['active_extent']={axis:[min(v['coordinate'][i] for v in locations),max(v['coordinate'][i] for v in locations)] for i,axis in ((0,'x'),(2,'z'))} if locations else {}
    diagnostic_path=os.environ.get('LD_LIBRARY_PATH')
    os.environ['LD_LIBRARY_PATH']=str(r/'campaign/trajectory/runtime/lib')
    try:current_identity=identity(r/'campaign/trajectory/runtime/IncNavierStokesSolver')
    finally:
        if diagnostic_path is None:os.environ.pop('LD_LIBRARY_PATH',None)
        else:os.environ['LD_LIBRARY_PATH']=diagnostic_path
    pinned=load(r/'campaign/trajectory/runtime/identity.json')['frozen']
    require(current_identity['sha256']==pinned['sha256'] and current_identity['libraries']==pinned['libraries'],'pinned runtime changed')
    supervisor_status=load_progress(r/'campaign/status.json')
    supervisor_id=supervisor_status.get('supervisor_identity') or {}
    identities=dict(driver=process_identity(progress['pid']),supervisor=process_identity(supervisor_id['pid']) if supervisor_id else None,native=process_identity(progress['native_pid']),
        monitor=process_identity(load(r/'monitor-launch.json')['pid']))
    later=[]
    for path in sorted((r/'monitor').glob('step-??????/summary.json')):
        data=load(path)
        if data.get('last_verified_step',0)<=32:continue
        n=int(data['last_verified_step']);state=data['states'].get('step'+str(n),{})
        later.append(dict(step=n,time=data['last_verified_time'],states={str(k):v for k,v in state.items() if k in ('global_k_integral','near_k_integral','mean_wall_stress','spatial_block_norms','epsilon_trace_defect','spatial_rate_mass_norm')},
            cumulative_correction=data['cumulative_correction'],last_budget={k:v for k,v in data['budgets'][-1].items() if k not in ('active_constraints',)}))
    evidence=dict(schema=1,base='0b86269739075af81622319402370967fdeffbe8',branch=git('branch','--show-current'),
        source_revision=git('rev-parse','HEAD'),source_status=git('status','--porcelain'),source_files_sha256=source_files,
        root=str(r),collected_unix=time.time(),production_launch=load(r/'launch.json'),monitor_launch=load(r/'monitor-launch.json'),
        progress=progress,supervisor_status=load_progress(r/'campaign/status.json'),monitor_status=load_progress(r/'monitor/status.json'),process_identities=identities,
        runtime=load(r/'campaign/trajectory/runtime/identity.json'),runtime_identity_unchanged=True,
        production_kernel_revision=git('rev-parse','1c49c07'),probe=identity(r/'audit-probe/RANSWallResolutionProbe'),
        checkpoint=checkpoint,production_counts_at_checkpoint=counts,fixtures=fixtures,history_checks=histories,
        qp_log=(r/'qp.log').read_text(),
        durability=json.loads((r/('durability-recovery.log' if (r/'durability-recovery.log').exists() else 'durability-terminal.log')).read_text()),
        supervisor_recovery=load(r/'supervisor-recovery.json') if (r/'supervisor-recovery.json').exists() else None,later_analyses=later,simultaneous_resume=load(r/'simultaneous-resume.json'),
        frozen_repair_budgets=load(r/'repair-budget.json'),analysis=compact,
        ranks_per_native_job=1,thread_environment={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')},
        build_parallelism=8,worker_limit=72,
        build_cache=[v for v in (r/'build/CMakeCache.txt').read_text().splitlines() if any(v.startswith(k) for k in ('CMAKE_CXX_COMPILER:','CMAKE_CXX_FLAGS:','CMAKE_BUILD_TYPE:','Nektar++_DIR:','RANS_BUILD_SOLVER:'))],
        receipts={str(path.relative_to(r)):sha(path) for path in sorted(r.glob('*.log'))},
        artifact_index=dict(path=str(a.output/'artifact-index.json'),sha256=sha(a.output/'artifact-index.json'),files=len(index)))
    write_json(a.output/'native.json',evidence);print(a.output/'native.json')


if __name__=='__main__':main()
