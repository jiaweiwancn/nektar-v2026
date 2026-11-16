#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Actual selected-space segment retirement, fresh restore/replay and holdout."""
import argparse,copy,json,os,resource,shutil,time,xml.etree.ElementTree as ET
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_controls import enable
from ct_accepted15_controls import ordering
from ct_accepted15_retention import publish,retire,verify,atomic
from ct_accepted14_pair import numerical_receipts
from ct_accepted14_collect import pressure_receipts
from ct_accepted11_target import bounded_budgets
from ct_accepted15_exec import native_run
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha
from compare_cls_epawf_pimple import same_history

STAGES=('committed','post_flow','pre_k','raw','corrected','history','fixed_reference')
def referenced_files(value):
    files=[]
    if isinstance(value,dict):
        if 'file' in value and 'sha256' in value and isinstance(value['file'],str):
            path=Path(value['file']);require(sha(path)==value['sha256'],'corrupt referenced native payload');files.append(path)
        if 'matrix_file' in value and 'matrix_sha256' in value:
            path=Path(value['matrix_file']);require(sha(path)==value['matrix_sha256'],'corrupt stage matrix');files.append(path)
        for v in value.values():files.extend(referenced_files(v))
    elif isinstance(value,list):
        for v in value:
            if isinstance(v,(dict,list)):files.extend(referenced_files(v))
    return files
def snapshots(report):
    out=[]
    for row in report['pimple_attempts']:
        if not row.get('scheduled_complete'):continue
        path=Path(row['k_snapshot']);require(sha(path)==row['k_snapshot_sha256'],'changed actual stage')
        q=n.read(path);require(sha(Path(q['matrix_file']))==q['matrix_sha256'],'changed actual scalar matrix')
        flow=Path(row['flow_snapshot']);require(sha(flow)==row['flow_snapshot_sha256'],'changed actual flow stage')
        out.append(dict(time=q['time'],alpha=q['alpha'],stages={k:exact_vector_hash(q[k]) for k in STAGES},matrix_sha256=q['matrix_sha256'],flow={k:v for k,v in n.read(flow).items() if isinstance(v,list) and (not v or isinstance(v[0],(int,float,bool)))}))
    return out

def main():
    p=argparse.ArgumentParser()
    for k in ('root','runtime','pair','training','static','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1';os.environ.pop('NEKTAR_RANS_TIMING_FILE',None);n.reference=lru_cache(None)(n.reference)
    pair=n.read(a.pair/'results.json');require(pair['qualification']['status']=='pass','fresh pair prerequisite')
    m=RequestManifest(a.output,dict(schema='CT15ActualSegmentRecoveryV1',runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.runtime/'RANSPreservationProbe'),pair_sha256=sha(a.pair/'results.json'),training_sha256=sha(a.training/'results.json'),contract_sha256=sha(a.root/'source/docs/plans/cao-tamura-accepted-15-evidence/contract.json'),held_out='steps6-to8 fixed before training fit',source_sha256=sha(Path(__file__))))
    c=pair['case'];ini=n.read(a.static/'target/frozen.json');weights=n.inventory(c)[2];count=len(weights)
    full=n.read(a.pair/'0.00025/case.evolution.json');full_snaps=snapshots(full)
    zero=n.read(a.training/'steps-0/case.evolution.json');two=n.read(a.training/'steps-2/case.evolution.json')
    require(two['state']==full['accepted_states'][2],'training is a different trajectory')
    # Frozen training model: a conservative nonnegative startup/step envelope.
    training=n.read(a.training/'results.json')['cases'];cpu=[training[str(i)]['execution']['process_cpu_seconds'] for i in (0,1,2)];wall=[training[str(i)]['execution']['wall_seconds'] for i in (0,1,2)]
    def fit(v):return dict(startup=max(0.,v[0]),repeated_step=max(0.,v[1]-v[0],(v[2]-v[0])/2,v[2]-v[1]),training=v)
    model=dict(CPU=fit(cpu),wall=fit(wall),fit_scope='training zero/one/two only; heldout not used for fit')
    write_json(a.output/'frozen-model-before-heldout.json',model)
    def run(label,parent,steps=2,reject=0):
        dest=a.output/label;n.make(dest,a.runtime,c,.00025,steps,parent);enable(dest);ordering(dest)
        root=ET.parse(dest/'case.xml').getroot();set_parameter(root,'RANSWallPimpleHardAuditRejectStep',reject);write(root,dest/'case.xml');os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(dest/'sparse-receipts.jsonl')
        e=native_run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log',1600,2000);os.environ.pop('NEKTAR_RANS_SPARSE_RECEIPTS',None);r=n.read(dest/'case.evolution.json')
        return dest,r,e
    binding=dict(runtime=m.payload['request']['runtime'],probe=m.payload['request']['probe'],dt=.00025,signature=full['signature'])
    result=dict(schema='CT15ActualSegmentRecoveryV1',status='incomplete',policy='SelectiveReplayV1',cycles=[],model=model)
    previous=zero;current=two;current_folder=a.training/'steps-2'
    for end in (2,4):
        cp=Path(previous['checkpoint']);restored=a.output/('restored-'+str(end));shutil.copytree(cp,restored)
        require(all(sha(f)==sha(restored/f.relative_to(cp)) for f in cp.rglob('*') if f.is_file()),'native restart restore bytes differ')
        replay_folder,replayed,execution=run('replayed-'+str(end),restored)
        require(execution['exit_code']==0 and not replayed['failure'] and same_history(current,replayed),'fresh-process restoration/replay history')
        recovered=snapshots(replayed);actual=snapshots(current);require(recovered==actual,'actual stage/flow/scalar arrays not recovered exactly')
        require(recovered==full_snaps[end-2:end],'segmented versus continuous stage trajectory')
        metrics=n.physical_metrics(ini,replayed['state'],c,weights,count);expected=pair['members']['0.00025']['accepted_metrics'][end]
        require(metrics=={k:v for k,v in expected.items() if k!='time'},'full versus replayed scientific classification')
        nr=numerical_receipts(replay_folder/'sparse-receipts.jsonl');pr=pressure_receipts(replayed)
        proof=dict(status='pass',decisions_exact_equal=True,maximum_scaled_array_difference=0.,maximum_scaled_state_pressure_traction_difference=0.,replay_report_sha256=sha(replay_folder/'case.evolution.json'),actual_stage_witnesses=actual,independent_metrics=metrics,numerical_receipts=nr,pressure_receipts=pr)
        proof_path=a.output/('replay-proof-'+str(end)+'.json');write_json(proof_path,proof)
        # Retire only stage payloads inside the newly owned member's snapshot
        # directory. Reports, receipts and both native restart states remain.
        eligible=[]
        for row in current['pimple_attempts']:
            for key in ('flow_snapshot','k_snapshot'):
                if row.get(key):
                    path=Path(row[key]);require(path.is_relative_to(current_folder),'snapshot outside newly owned member');eligible.append(path)
                    eligible.extend(referenced_files(n.read(path)))
        eligible=sorted(set(eligible))
        required=[proof_path,current_folder/'case.evolution.json',current_folder/'case.xml',current_folder/'sparse-receipts.jsonl']
        for checkpoint in (Path(previous['checkpoint']),Path(current['checkpoint'])):required.extend(f for f in checkpoint.rglob('*') if f.is_file())
        publication=a.output/('retirement-'+str(end));manifest=publish(a.root,publication,binding,required,[Path(previous['checkpoint']),Path(current['checkpoint'])],eligible,proof)
        negative={}
        for label,mutate in [('foreign',lambda b:b.update(signature='foreign')),('stale_dt',lambda b:b.update(dt=.000125))]:
            altered=copy.deepcopy(binding);mutate(altered)
            try:retire(a.root,publication,altered)
            except (ValueError,RuntimeError):negative[label]='pass'
            else:raise ValueError('foreign retirement accepted')
            require(all(f.exists() for f in eligible),'negative control partially deleted payloads')
        interrupted=a.output/('interrupted-'+str(end));interrupted.mkdir();(interrupted/'manifest.json.pending').write_bytes((publication/'manifest.json').read_bytes())
        try:retire(a.root,interrupted,binding)
        except (ValueError,RuntimeError):negative['incomplete_seal']='pass'
        else:raise ValueError('incomplete publication accepted')
        require(all(f.exists() for f in eligible),'interruption pruned recovery inputs')
        for label,missing in [('missing_required',True),('corrupt_required',False)]:
            bad=copy.deepcopy(manifest);key=next(iter(bad['required_files']))
            if missing:bad['required_files']['absent-required-payload']=bad['required_files'].pop(key)
            else:bad['required_files'][key]['sha256']='0'*64
            test=a.output/(label+'-'+str(end));test.mkdir();atomic(test/'manifest.json',bad);atomic(test/'seal.json',dict(manifest_sha256=sha(test/'manifest.json')))
            try:retire(a.root,test,binding)
            except (ValueError,RuntimeError,FileNotFoundError):negative[label]='pass'
            else:raise ValueError('missing/corrupt recovery payload accepted')
            require(all(f.exists() for f in eligible),'bad payload caused partial prune')
        # Recover an interrupted publish from the independently sealed identical
        # manifest, then verify it before allowing any deletion.
        os.replace(interrupted/'manifest.json.pending',interrupted/'manifest.json');atomic(interrupted/'seal.json',dict(manifest_sha256=sha(interrupted/'manifest.json')));verify(a.root,interrupted,binding)
        removed=retire(a.root,publication,binding)
        result['cycles'].append(dict(end_step=end,execution=execution,retired_bytes=removed,publication_sha256=sha(publication/'seal.json'),negative_controls=negative,proof_sha256=sha(proof_path)))
        write_json(a.output/'progress.json',result)
        previous=current
        if end==2:
            current_folder,current,e=run('segment-4',Path(previous['checkpoint']));require(e['exit_code']==0 and not current['failure'],'genuine continuation after first cycle');result['segment4_execution']=e
    six,sr,se=run('segment-6',Path(current['checkpoint']));require(se['exit_code']==0 and not sr['failure'],'segment6 continuation')
    started=time.monotonic();self_cpu=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    eight,er,ee=run('held-out-6-to8',Path(sr['checkpoint']));require(ee['exit_code']==0 and not er['failure'] and same_history(full,er),'heldout later-state continuation')
    require(snapshots(er)==full_snaps[6:8],'heldout stage trajectory')
    budgets=bounded_budgets(a.runtime,eight,ini,er,c);require(budgets['gate']['status']=='pass','independent later-state audit replay')
    metrics=n.physical_metrics(ini,er['state'],c,weights,count);require(metrics=={k:v for k,v in pair['members']['0.00025']['accepted_metrics'][8].items() if k!='time'},'heldout physical decisions')
    after=resource.getrusage(resource.RUSAGE_CHILDREN);measured=dict(native=ee,whole_CPU=time.process_time()-self_cpu+after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,whole_wall=time.monotonic()-started,output_bytes=sum(f.stat().st_size for f in eight.rglob('*') if f.is_file()),budgets=budgets)
    rollback,rr,re=run('rejected-step',Path(sr['checkpoint']),2,reject=8)
    seven,sv,sve=run('rollback-predecessor7',Path(sr['checkpoint']),1)
    require(sve['exit_code']==0 and re['exit_code']!=0 and 'deliberate' in rr['failure'] and same_history(sv,rr),'selected full rejected-step rollback')
    result.update(status='pass',held_out=measured,segment6_execution=se,rollback_execution=re,rollback_predecessor_execution=sve,rollback_report_sha256=sha(rollback/'case.evolution.json'),live_bytes=sum(f.stat().st_size for f in a.root.rglob('*') if f.is_file()))
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
