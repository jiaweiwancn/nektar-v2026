#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify persisted moments and archive before rotating explicitly owned data.

This is a bounded physical foundation, not H27-D or a Task3 readiness claim.
The independent two-pass reduction is compared with the streaming accumulator.
No input is removed unless --commit-retention is explicitly supplied and every
preceding scientific/history check passes. Earlier experiments are never owned.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import tarfile
import time

sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import history_recovery_check as hr
import implicit_algorithm_check as fields
import physical_cache_history_check as history
import physical_oracle as po
import physical_transient_oracle as transient
import scientific_accumulator as accumulator
import terminal_output_check as terminal

require,read = po.require,history.read
LIMIT = 1e-8
ROUNDING = 2e-12
SCHEMA = 'WJH27ScientificRetentionCheckV1'


def clock(deadline):require(time.monotonic() < deadline,'scientific retention exhausted its cooperative cap')


def durable_json(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    fd = os.open(path.parent,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)


def load_epoch(spec,root,control,ref,start,deadline):
    path,result_id = history.input_path(spec['result'],root);packet = read(path);data = packet['state']
    index_path,index_id = history.input_path(spec['index'],root);index = read(index_path)
    _,session_id,parameters,solver = terminal.session(spec['session'],root)
    history.configuration(data,control)
    require(data['step'] == start+4 and data['time'] == (start+4)*data['dt'],'epoch endpoint differs')
    require(index['schema'] == 'WJH27NativeScientificIndexV1' and index['run_complete'] is True and
        index['reason_bits'] == terminal.BITS and index['terminal_output'] is True and index['coefficients_included'] is True,
        'incomplete epoch scientific index')
    require((index['terminal_begin'],index['terminal_end']) == terminal.WINDOW and index['initial_step'] == start and
        index['initial_time'] == start*data['dt'] and index['final_step'] == data['step'] and index['final_time'] == data['time'],
        'epoch scientific clock/window differs')
    require(all(index[key] == data[key] for key in ('source_identity','reference_identity','signature','dt')),'epoch index identity differs')
    for key,value in dict(TimeStep=.005,NumSteps=4,IO_InfoSteps=1,IO_CheckSteps=3,WJScientificSteps=1,
        WJTerminalBegin=terminal.WINDOW[0],WJTerminalEnd=terminal.WINDOW[1]).items():
        require(float(parameters[key]) == value,'epoch session differs: '+key)
    for key,value in dict(SolverType='WJNativeCoupledV2',WJNativeMode='Physical',WJImplicitSolver='CoupledJFNKV1',
        WJTerminalOutput='True',WJOutputCoefficients='True').items():require(solver[key] == value,'epoch solver differs: '+key)
    prefix = history.replace_suffix(path,'.result.json','')
    require(hr.resolve(root,solver['WJOutputPrefix']) == prefix and hr.resolve(root,packet['scientific_index']) == index_path and
        index_path == Path(str(prefix)+'.science-index.json'),'epoch output prefix differs')
    require(bool(solver.get('WJRestart')) is bool(start) and packet['restart'] is bool(start),'epoch restart mode differs')
    mesh = history.session_geometry(spec['session'],root,ref);owners = history.inventory(data,control,ref,mesh)
    log_path = Path(str(prefix)+'.steps.jsonl');logs = [hr.native_json(line) for line in log_path.read_text().splitlines()]
    require([row['step'] for row in logs] == list(range(start+1,start+5)),'epoch accepted log is incomplete')
    previous = None;event = None
    for row in logs:
        clock(deadline);history.accepted(row,control,False)
        _,event = history.cache_events(row,event,previous is None,False,'CurrentDataDriftV1')
        require(row['time'] == row['step']*.005,'epoch accepted clock differs')
        if previous is not None:
            require(row['previous_sha256'] == previous['state_sha256'] and row['previous_advection_sha256'] == previous['advection_sha256'],
                'epoch accepted history chain differs')
        previous = row
    require(all(logs[-1][key] == data[key] for key in history.HASH_NAMES),'epoch endpoint/log differs')
    expected = terminal.expected_schedule(start,4,.005,1,*terminal.WINDOW)
    require([f['step'] for f in index['frames']] == sorted(expected),'epoch scientific schedule differs')
    records = {};ids = [];by_step = {row['step']:row for row in logs}
    for frame in index['frames']:
        clock(deadline);p = hr.resolve(root,frame['path']);ident = hr.identity(p);record = read(p);step = frame['step']
        require(p == Path(str(prefix)+'.science-'+str(step)+'.json') and ident['sha256'] == frame['sha256'] and
            ident['bytes'] == frame['bytes'] and frame['reason_mask'] == expected[step] and frame['time'] == step*.005,
            'epoch scientific byte/clock/reason differs')
        require(frame['origin'] == 'current_process_accepted_step' and frame['write_seconds'] > 0,'epoch frame provenance differs')
        require(terminal.record_hash(record) == record['state_sha256'],'epoch current coefficient SHA differs')
        for key,header in [('current_advection','advection_sha256'),('previous_advection','previous_advection_sha256')]:
            require(transient.state_hash(transient.state(record[key],record['sizes'])) == record[header],'epoch explicit coefficient SHA differs')
        history.configuration(record,control)
        require(all(record[key] == data[key] for key in ('geometry','elements','free_velocity_ids','k_free_ids','epsilon_free_ids')),
            'epoch original-owner layout differs')
        require(record['step'] == step and record['time'] == frame['time'] and all(record[key] == by_step[step][key] for key in
            (*history.HASH_NAMES,'residuals',*history.GUARDS)),'epoch scientific/log binding differs')
        records[step] = record;ids.append(ident)
    require({p.resolve() for p in prefix.parent.glob(prefix.name+'.science-[0-9]*.json')} == {Path(v['path']) for v in ids},
        'epoch actual scientific inventory differs')
    times = packet['complete_step_seconds'];advancing = packet['complete_advancing_seconds'];reserve = 128*np.finfo(float).eps
    require(len(times) == 4 and all(math.isfinite(t) and t > 0 for t in times) and
        abs(math.fsum(times)-advancing) <= reserve*advancing,'epoch complete-step timing does not close')
    write = math.fsum(f['write_seconds'] for f in index['frames'])
    require(abs(write-packet['scientific_frame_seconds']) <= reserve*write and 0 < packet['scientific_index_seconds'] <= packet['final_export_seconds'],
        'epoch scientific serialization timing differs')
    states,binary = history.binary_history(Path(str(prefix)+'.final.bin'),data,.005,deadline)
    histories = history.exported(data,states)
    return dict(data=data,packet=packet,index=index,records=records,logs=logs,owners=owners,mesh=mesh,histories=histories,solver=solver,prefix=prefix,
        evidence=dict(result=result_id,index=index_id,session=session_id,log=hr.identity(log_path),frames=ids,checkpoint=binary,
            complete_advancing_seconds=advancing,complete_step_seconds=times,scientific_frame_seconds=write,
            scientific_index_seconds=packet['scientific_index_seconds']))


def own(path,root):
    require(not path.is_symlink() and path.resolve().parent.is_relative_to(root) and path.is_file(),'retention path is not an owned regular file')
    return hr.identity(path)


def rotate(root,epochs,generations,verified,deadline):
    keep = [];remove = [];all_frames = []
    for run in epochs:
        for item in run['evidence']['frames']:
            all_frames.append(item);(keep if read(Path(item['path']))['step'] >= 11 else remove).append(item)
        for p in sorted(run['prefix'].parent.glob(run['prefix'].name+'.recovery-*.bin')):
            step = int(p.name.rsplit('-',1)[1][:-4]);(keep if step >= 9 else remove).append(own(p,root))
    for i,run in enumerate(epochs):
        (keep if i >= 1 else remove).append(own(Path(str(run['prefix'])+'.final.bin'),root))
        (keep if i >= 1 else remove).append(own(Path(generations[i]['arrays']['path']),root))
    require(len(all_frames) == 12 and len(keep) == 8 and len(remove) == 14,'prospective retention inventory differs')
    archive_members = {item['path']:item for item in all_frames+remove}
    for item in keep+remove:require(own(Path(item['path']),root) == item,'owned retention input changed')
    protected = root/'protected-retention-v2';protected.mkdir()
    archive = protected/'qualification-before-retention.tar.gz'
    with tarfile.open(archive,'x:gz',compresslevel=5) as tar:
        for item in sorted(archive_members.values(),key=lambda x:x['path']):
            clock(deadline);path = Path(item['path']);tar.add(path,arcname=str(path.relative_to(root)),recursive=False)
    with archive.open('rb') as stream:os.fsync(stream.fileno())
    with tarfile.open(archive,'r:gz') as tar:
        require(set(tar.getnames()) == {str(Path(p).relative_to(root)) for p in archive_members},'protected archive inventory differs')
        for item in archive_members.values():
            clock(deadline);member = tar.extractfile(str(Path(item['path']).relative_to(root)))
            require(member is not None and hashlib.sha256(member.read()).hexdigest() == item['sha256'],'protected archive bytes differ')
    plan = dict(schema='WJH27OwnedRetentionPlanV1',verified_reduction=hr.identity(verified),archive=hr.identity(archive),
        archived_members=list(archive_members.values()),retained=keep,remove=remove,
        scope='Only new epoch1/2/3 raw science, recovery/final and V2 NPZ payloads; original indices immutable',
        production_saving_excludes_qualification_archive=True)
    plan_path = protected/'verified-plan.json';durable_json(plan_path,plan)
    for item in remove:
        clock(deadline);path = Path(item['path']);require(own(path,root) == item,'owned bytes changed immediately before removal');path.unlink()
    require(all(not Path(item['path']).exists() for item in remove) and all(own(Path(item['path']),root) == item for item in keep),
        'actual retention outcome differs')
    for parent in {Path(item['path']).parent for item in keep+remove}:
        fd = os.open(parent,os.O_RDONLY)
        try:os.fsync(fd)
        finally:os.close(fd)
    receipt = dict(schema='WJH27OwnedRetentionReceiptV1',status='PASS_ACTUAL_OWNED_ROTATION',plan=hr.identity(plan_path),
        retained=keep,removed=remove,archive=hr.identity(archive),removed_payload_bytes=sum(v['bytes'] for v in remove),
        retained_payload_bytes=sum(v['bytes'] for v in keep),immutable_indices=[run['evidence']['index'] for run in epochs],
        original_indices_describe_pre_retention_writes=True,full_H27_qualification=False,H27_D_pass=False,TASK3_READY=False)
    receipt_path = protected/'retention-receipt.json';durable_json(receipt_path,receipt)
    return dict(receipt=hr.identity(receipt_path),**receipt)


def run(manifest,folder,table,deadline,out,output,commit):
    require(manifest['schema'] == 'WJH27ScientificRetentionV1' and tuple(manifest['window']) == terminal.WINDOW,'wrong retention manifest/window')
    root = hr.resolve(folder,manifest['root']);require(hr.identity(table)['sha256'] == po.TABLE_SHA,'wrong scientific table')
    ref = po.Reference(read(table));cp,_ = history.input_path(manifest['control']['result'],root);control = read(cp)['state']
    require(control['source_identity'] == manifest['source_identity'],'retention source identity differs')
    dense = terminal.load_run('dense',manifest['control'],root,control,ref,deadline)
    require(len(manifest['epochs']) == 3,'retention requires three native epochs')
    epochs = [];generations = [];generation_ids = [];stitched = {};last = None
    for i,spec in enumerate(manifest['epochs']):
        require(spec['name'] == 'epoch'+str(i+1),'unexpected epoch owner')
        run = load_epoch(spec['native'],root,control,ref,4*i,deadline);epochs.append(run);stitched.update(run['records'])
        gp,gid = history.input_path(spec['generation'],root);g,arrays = accumulator.read_generation(gp)
        require(g['batch_number'] == i+1 and g['cumulative_frame_count'] == (i+1)*4 and g['last_step'] == (i+1)*4 and
            g['last_time'] == (i+1)*.02 and tuple(g['window']) == terminal.WINDOW and
            abs(g['weight']-min(max((i+1)*.02-terminal.WINDOW[0],0),.0364)) <= 4*np.finfo(float).eps,
            'scientific generation clock/count/weight differs')
        require(g['native_result'] == run['evidence']['result'] and g['native_index'] == run['evidence']['index'] and
            g['session'] == run['evidence']['session'] and all(g['last_native_hashes'][key] == run['data'][key] for key in history.HASH_NAMES),
            'scientific generation native identity differs')
        require(g['layout_sha256'] == accumulator.digest({key:control[key] for key in accumulator.METADATA}) and
            g['math'] == 'WeightedLinearSegmentWelfordV1' and g['sampling'] == 'IndependentSignedGauss5EndpointsV1' and
            all(hr.identity(Path(item['path'])) == item for item in g['helper_sources'].values()),'generation algorithm/source binding differs')
        require(g['previous_generation'] == (generation_ids[-1] if i else None),'scientific generation hash chain differs')
        require([f['step'] for f in g['frames']] == list(range(4*i+1,4*i+5)) and
            all(hr.identity(Path(f['path'])) == {key:f[key] for key in ('path','sha256','bytes')} for f in g['frames']),
            'generation copied/lost batch provenance')
        if last is not None:
            rp,rid = history.input_path(spec['restart'],root)
            _,binary = history.binary_history(rp,control,.005,deadline,rid['sha256'])
            require(rp == Path(str(last['prefix'])+'.final.bin') and binary['step'] == 4*i and
                all(binary['hashes'][role] == last['data'][key] for role,key in zip(history.HISTORY_NAMES,history.HASH_NAMES)),
                'cold epoch input does not bind complete predecessor')
            require(hr.resolve(root,run['solver']['WJRestart']) == rp and run['logs'][0]['previous_sha256'] == binary['hashes']['current'] and
                run['logs'][0]['previous_advection_sha256'] == binary['hashes']['advection'],'cold epoch lost imported history')
        generations.append(g);generation_ids.append(gid);last = run
    require(generations[-1]['window_complete'] is True and generations[-1]['covered_until'] == terminal.WINDOW[1] and
        generations[-1]['weight'] == .0364,'scientific final window incomplete')
    require(generations[-1]['cumulative_segment_count'] == sum(len(g['segments']) for g in generations),'scientific segment count differs')
    out['runs'] = dict(control=dense['evidence'],epochs=[e['evidence'] for e in epochs])
    out['generations'] = [dict(identity=ident,weight=g['weight'],batch_number=g['batch_number'],frames=len(g['frames']),segments=len(g['segments']),
        arrays=g['arrays'],covered_until=g['covered_until']) for g,ident in zip(generations,generation_ids)]
    out['endpoint'] = history.compare(dense['histories'],epochs[-1]['histories'],control,dense['owners'],ref,deadline)
    out['current_equations'] = history.equation_report(manifest['current_equation_report'],root,
        Path(epochs[-1]['evidence']['result']['path']),epochs[-1]['data'],epochs[-1]['histories'],po.TABLE_SHA)
    out['common_clock_coefficients'] = {}
    for step,record in stitched.items():
        a,b = transient.state(record['coefficients'],control['sizes']),transient.state(dense['records'][step]['coefficients'],control['sizes'])
        delta = fields.coefficient_delta(a,b,control['sizes'],ref)
        require(all(v['scaled_MAX'] <= LIMIT for v in delta.values()),'cold scientific frame exceeds fixed coefficient scale')
        out['common_clock_coefficients'][step] = delta
    data = {};reductions = {}
    for name,records in [('split',stitched),('control',dense['records'])]:
        data[name],reductions[name] = terminal.reduce(records,dense['owners'],dense['mesh'],ref,deadline,
            output.with_name(output.stem+'-'+name+'.npz'))
    _,final_arrays = accumulator.read_generation(Path(generation_ids[-1]['path']));stream = {}
    for key,g in final_arrays.items():
        group = data['split'][key]['group']
        require(all(np.array_equal(g[name],group[name]) for name in ('scale','x','r','weight','owner')),'streamed original-owner inventory differs')
        stream[key] = dict(mean=g['mean'],sd=np.sqrt(g['m2']/generations[-1]['weight']),group=group)
    out['reductions'] = reductions
    out['comparisons'] = dict(stream_vs_two_pass=terminal.compare(stream,data['split'],{key:ROUNDING for key in terminal.LIMITS}),
        split_vs_uninterrupted=terminal.compare(data['split'],data['control'],{key:LIMIT for key in terminal.LIMITS}))
    require(all(v['status'] == 'MET' for rows in out['comparisons'].values() for v in rows.values()),'scientific moment comparison exceeds prospective tolerance')
    stale = read(root/'epoch3-accumulator-manifest.json');stale['resume'] = generation_ids[0]
    stale_prefix = output.with_name(output.stem+'-stale-control')
    try:accumulator.run(stale,root,table,stale_prefix,deadline)
    except ValueError as error:
        require('immediate predecessor' in str(error),'stale generation failed for an unrelated reason');out['stale_generation_control'] = 'REJECTED_IMMEDIATE_PREDECESSOR'
    else:raise ValueError('stale scientific generation accepted')
    require(not Path(str(stale_prefix)+'.json').exists() and not Path(str(stale_prefix)+'.npz').exists(),'stale control published a generation')
    verified = output.with_name(output.stem+'-verified-reduction.json')
    durable_json(verified,dict(schema='WJH27VerifiedScientificReductionV1',status='PASS_SCIENTIFIC_REDUCTION_ONLY',
        manifest=out['manifest'],comparisons=out['comparisons'],generations=out['generations'],reductions=reductions,
        endpoint=out['endpoint'],current_equations=out['current_equations'],sources=out['sources'],
        full_H27_qualification=False,H27_D_pass=False,TASK3_READY=False))
    out['verified_reduction'] = hr.identity(verified)
    if commit:out['retention'] = rotate(root,epochs,generations,verified,deadline)
    out.update(status='PASS_SCIENTIFIC_REDUCTION_AND_ROTATION_ONLY' if commit else 'PASS_SCIENTIFIC_REDUCTION_ONLY',
        foundation_gate_pass=True,raw_payload_rotation_performed=commit,retained_recovery_cold_imports_still_required=True,
        production_storage_qualified=False,monitoring_escalation_qualified=False,long_pilots_qualified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','table','out'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True);parser.add_argument('--commit-retention',action='store_true');args = parser.parse_args()
    began = time.monotonic();require(not args.out.exists(),'retention assessment output is immutable')
    out = dict(schema=SCHEMA,status='FAIL',foundation_gate_pass=False,full_H27_qualification=False,H27_D_pass=False,TASK3_READY=False,
        manifest=hr.identity(args.manifest),table=hr.identity(args.table),workers=1,streaming_comparison_limit=ROUNDING,cold_comparison_limit=LIMIT,
        sources={key:hr.identity(path) for key,path in dict(checker=__file__,accumulator=accumulator.__file__,terminal=terminal.__file__,
            history=history.__file__,fields=fields.__file__,reference=po.__file__,transient=transient.__file__,JSON_reader=hr.__file__).items()})
    try:
        require(sys.platform.startswith('linux') and sys.byteorder == 'little' and math.isfinite(args.seconds) and 0 < args.seconds <= 180,
            'bounded authoritative little-endian Linux analysis required')
        soft,hard = resource.getrlimit(resource.RLIMIT_AS);cap = min([8*1024**3]+[v for v in (soft,hard) if v != resource.RLIM_INFINITY])
        resource.setrlimit(resource.RLIMIT_AS,(cap,hard));out['memory_cap_bytes'] = cap
        run(read(args.manifest),args.manifest.resolve().parent,args.table.resolve(),began+args.seconds,out,args.out.resolve(),args.commit_retention)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-began;out['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    durable_json(args.out,out);print(json.dumps({key:out[key] for key in ('status','reason','elapsed_seconds') if key in out}))
    if not out['foundation_gate_pass']:raise SystemExit(2)


if __name__ == '__main__':main()
