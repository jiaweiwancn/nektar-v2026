#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Persisted, original-owner scientific moments for bounded native batches.

Consumes WJH27ScientificBatchV1. It never deletes input data. A separate
independent assessment must authorize owned retention after verification.
All observations/scales come from the independent source14 physical sampler.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import xml.etree.ElementTree as ET
sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import history_recovery_check as hr
import physical_cache_history_check as history
import physical_oracle as po
import physical_transient_oracle as transient
import implicit_algorithm_check as fields
import terminal_output_check as terminal

require,read = po.require,history.read
SCHEMA = 'WJH27ScientificAccumulatorV1'
METADATA = ('schema','model','source_identity','reference_identity','signature','sizes','gauge','manufactured','basis',
            'geometry','elements','free_velocity_ids','k_free_ids','epsilon_free_ids','dt','quadratureX','quadratureR',
            'linear_tolerance','nonlinear_tolerance','constitutive_backend','implicit_solver','anderson_depth')
ARRAY_NAMES = ('mean','m2','last','scale','x','r','weight','owner')
HASHES = history.HASH_NAMES


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def clock(deadline):require(time.monotonic() < deadline,'scientific accumulator exhausted its cooperative cap')


def identity(path):return hr.identity(path)


def empty_groups(observation):
    result = {}
    for key,g in observation.items():
        result[key] = {name:g[name].copy() for name in ('scale','x','r','weight','owner')}
        result[key].update(region=g['region'],field=g['field'],mean=np.zeros_like(g['value']),
                           m2=np.zeros_like(g['value']),last=g['value'].copy())
    return result


def same_inventory(a,b):
    require(a.keys() == b.keys(),'scientific observation inventory changed')
    for key in a:
        for name in ('scale','x','r','weight','owner'):
            require(np.array_equal(a[key][name],b[key][name]),'scientific owner/scale changed: '+key+'/'+name)
        require(a[key]['mean'].shape == b[key]['value'].shape,'scientific observation shape changed')


def read_generation(path):
    data = read(path);require(data['schema'] == SCHEMA and data['status'] == 'PASS_SCIENTIFIC_ACCUMULATOR_BATCH_ONLY','unverified previous scientific generation')
    npz = hr.resolve(path.parent,data['arrays']['path']);require(identity(npz)['sha256'] == data['arrays']['sha256'],'previous scientific arrays changed')
    groups = {}
    with np.load(npz,allow_pickle=False) as loaded:
        require(set(loaded.files) == {key+'__'+name for key in data['groups'] for name in ARRAY_NAMES},'previous scientific array inventory differs')
        for key,g in data['groups'].items():
            groups[key] = {name:loaded[key+'__'+name].copy() for name in ARRAY_NAMES}
            groups[key].update(region=g['region'],field=g['field'])
            require(all(np.all(np.isfinite(groups[key][name])) for name in ARRAY_NAMES),'previous scientific array contains nonfinite values')
            require(np.all(groups[key]['m2'] >= 0),'previous centered moment is negative')
    return data,groups


def save_generation(prefix,metadata,groups):
    path = Path(str(prefix)+'.json');array_path = Path(str(prefix)+'.npz')
    temporary = Path(str(array_path)+'.tmp')
    require(not path.exists() and not array_path.exists() and not temporary.exists(),'scientific generation output is immutable')
    arrays = {key+'__'+name:g[name] for key,g in groups.items() for name in ARRAY_NAMES}
    with temporary.open('xb') as stream:
        np.savez_compressed(stream,**arrays);stream.flush();os.fsync(stream.fileno())
    with np.load(temporary,allow_pickle=False) as loaded:
        require(set(loaded.files) == set(arrays) and all(np.array_equal(loaded[key],v) for key,v in arrays.items()),'scientific generation round-trip changed an array')
    temporary.rename(array_path)
    metadata['arrays'] = identity(array_path)
    metadata['groups'] = {key:dict(region=g['region'],field=g['field'],points=len(g['mean']),components=g['mean'].shape[1]) for key,g in groups.items()}
    with path.open('x') as stream:
        json.dump(metadata,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    # The same reader used by the next cold statistics process must accept it.
    reloaded,again = read_generation(path)
    require(reloaded['arrays'] == metadata['arrays'] and all(np.array_equal(groups[key][name],again[key][name]) for key in groups for name in ARRAY_NAMES),
            'persisted scientific generation cannot be resumed exactly')
    return path


def run(manifest,folder,table,prefix,deadline):
    require(manifest['schema'] == 'WJH27ScientificBatchV1','wrong scientific batch schema')
    root = hr.resolve(folder,manifest.get('root','.'));begin,end = map(float,manifest['window'])
    require(math.isfinite(begin) and math.isfinite(end) and 0 <= begin < end,'invalid scientific window')
    require(identity(table)['sha256'] == po.TABLE_SHA,'unqualified scientific reference table');ref = po.Reference(read(table))
    result_path,result_id = history.input_path(manifest['run']['result'],root);packet = read(result_path);control = packet['state']
    require(control['schema'] == 'WJH27NativeCoupledStateV2' and control['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and
            control['manufactured'] is False and control['gauge'] is False and control['basis'] == 'Modified_A', 'unsupported scientific model/pressure level')
    require(control['reference_identity'] == po.TABLE_SHA and control['implicit_solver'] == 'CoupledJFNKV1','wrong scientific recipe')
    layout = {key:control[key] for key in METADATA};layout_hash = digest(layout)
    owners,_ = fields.inventory(control,ref);mesh = dict(r_tolerance=64*np.finfo(float).eps*control['geometry']['top'])
    index_path,index_id = history.input_path(manifest['run']['index'],root);index = read(index_path)
    require(index['schema'] == 'WJH27NativeScientificIndexV1' and index['run_complete'] is True and
            index['coefficients_included'] is True and index['terminal_output'] is True and
            (index['terminal_begin'],index['terminal_end']) == (begin,end),'incomplete or incompatible scientific index')
    require(all(index[key] == control[key] for key in ('source_identity','reference_identity','signature','dt')),'scientific index identity differs')
    session_path,session_id = history.input_path(manifest['run']['session'],root);xml = ET.parse(session_path).getroot()
    parameters = {};solver = {}
    for item in xml.findall('.//PARAMETERS/P'):
        key,value = ''.join(item.itertext()).split('=',1);require(key.strip() not in parameters,'duplicate session parameter');parameters[key.strip()] = value.strip()
    for item in xml.findall('.//SOLVERINFO/I'):
        key = item.get('PROPERTY');require(key not in solver,'duplicate solver info');solver[key] = item.get('VALUE')
    start,finish = index['initial_step'],index['final_step'];dt = control['dt']
    require(type(start) is int and type(finish) is int and finish > start and index['initial_time'] == start*dt and
            index['final_time'] == finish*dt == control['time'] and finish == control['step'],'scientific batch clock differs')
    for key,value in dict(TimeStep=dt,NumSteps=finish-start,IO_InfoSteps=1,WJScientificSteps=1,WJTerminalBegin=begin,WJTerminalEnd=end).items():
        require(float(parameters[key]) == value,'scientific session differs: '+key)
    for key,value in dict(WJNativeMode='Physical',WJTerminalOutput='True',WJOutputCoefficients='True',WJImplicitSolver='CoupledJFNKV1').items():
        require(solver[key] == value,'scientific solver differs: '+key)
    native_prefix = history.replace_suffix(result_path,'.result.json','')
    require(hr.resolve(root,solver['WJOutputPrefix']) == native_prefix and hr.resolve(root,packet['scientific_index']) == index_path,
            'scientific output prefix differs')
    require(index_path == Path(str(native_prefix)+'.science-index.json') and
            hr.identity(hr.resolve(root,solver['WJReferenceTable']))['sha256'] == po.TABLE_SHA,'scientific index/table path differs')
    rows = [hr.native_json(line) for line in Path(str(native_prefix)+'.steps.jsonl').read_text().splitlines()]
    require([row['step'] for row in rows] == list(range(start+1,finish+1)),'scientific accepted log is incomplete')
    by_step = {row['step']:row for row in rows}
    for row in rows:
        clock(deadline);require(row['time'] == row['step']*dt and all(row[key] == control[key] for key in ('source_identity','reference_identity','signature')),'accepted scientific provenance differs')
        require(all(math.isfinite(v) and 0 <= v <= control['nonlinear_tolerance'] <= 1e-8 for v in row['residuals']) and
                all(row[key] > 0 for key in history.GUARDS),'accepted scientific residual/guard failed')
    expected = terminal.expected_schedule(start,finish-start,dt,1,begin,end)
    require([f['step'] for f in index['frames']] == sorted(expected),'scientific every-step/bracket schedule differs')
    groups = None;weight = 0.;covered = begin;last_time = None;last_step = None;last_hashes = None;segments = [];inputs = []
    prior = None
    helper_sources = {key:identity(path) for key,path in dict(accumulator=__file__,observations=terminal.__file__,fields=fields.__file__,
        reference=po.__file__,history=history.__file__,JSON_reader=hr.__file__).items()}
    if manifest.get('resume') is not None:
        previous_path,previous_id = history.input_path(manifest['resume'],root);prior,groups = read_generation(previous_path)
        require(prior['layout_sha256'] == layout_hash and tuple(prior['window']) == (begin,end),'scientific source/owner/window changed on resume')
        require(prior['sampling'] == 'IndependentSignedGauss5EndpointsV1' and prior['math'] == 'WeightedLinearSegmentWelfordV1','incompatible scientific accumulation algorithm')
        require(prior['helper_sources'] == helper_sources,'scientific resume helper identities changed')
        weight,covered,last_time,last_step,last_hashes = (prior[key] for key in ('weight','covered_until','last_time','last_step','last_native_hashes'))
        require(last_step == start and last_time == start*dt,'scientific generation is not the immediate predecessor')
        checkpoint,checkpoint_id = history.input_path(manifest['restart'],root)
        require(hr.resolve(root,solver['WJRestart']) == checkpoint and packet['restart'] is True,'native/statistics restart inputs differ')
        _,binary = history.binary_history(checkpoint,control,dt,deadline,checkpoint_id['sha256'])
        require(binary['step'] == start and all(binary['hashes'][role] == last_hashes[key] for role,key in zip(history.HISTORY_NAMES,HASHES)),
                'scientific generation does not link the complete imported native history')
        require(rows[0]['previous_sha256'] == last_hashes['state_sha256'] and rows[0]['previous_advection_sha256'] == last_hashes['advection_sha256'],
                'native cold step lost scientific predecessor')
        segments = list(prior['segments']);inputs = list(prior['frames'])
    else:
        require(start == 0 and not packet['restart'] and not solver.get('WJRestart'),'fresh scientific accumulation must start with a fresh native epoch')
    for frame in index['frames']:
        clock(deadline);path = hr.resolve(root,frame['path']);require(path == Path(str(native_prefix)+'.science-'+str(frame['step'])+'.json'),'scientific frame owner path differs')
        require(identity(path)['sha256'] == frame['sha256'] and path.stat().st_size == frame['bytes'] and frame['reason_mask'] == expected[frame['step']],
                'scientific frame bytes/reason differ')
        record = read(path);require(digest({key:record[key] for key in METADATA}) == layout_hash,'scientific frame layout/recipe changed')
        require(record['step'] == frame['step'] and record['time'] == frame['time'] == record['step']*dt,'scientific frame clock differs')
        require(transient.state_hash(transient.state(record['coefficients'],record['sizes'])) == record['state_sha256'],'scientific current coefficient hash differs')
        for key,header in [('current_advection','advection_sha256'),('previous_advection','previous_advection_sha256')]:
            require(transient.state_hash(transient.state(record[key],record['sizes'])) == record[header],'scientific explicit coefficient hash differs')
        if record['step'] > start:
            require(frame['origin'] == 'current_process_accepted_step' and all(record[key] == by_step[record['step']][key] for key in HASHES),
                    'scientific frame is not its declared accepted step')
        else:require(frame['origin'] in ('fresh_initial_observation','imported_accepted_history'),'invalid initial scientific provenance')
        hashes = {key:record[key] for key in HASHES}
        if last_step == record['step']:
            require(last_time == record['time'] and last_hashes == hashes,'conflicting duplicate scientific boundary frame');continue
        observation = terminal.observations(record,owners,mesh,ref,deadline)
        if groups is None:
            require(record['time'] <= begin,'scientific history begins after the requested window')
            groups = empty_groups(observation)
        else:
            same_inventory(groups,observation)
            require(record['time'] > last_time,'scientific clock moved backwards')
            lo,hi = max(last_time,begin),min(record['time'],end)
            if hi > lo:
                require(lo == covered,'scientific accumulation has an overlapping/gapped terminal interval')
                w = hi-lo;new_weight = weight+w
                lower = (lo-last_time)/(record['time']-last_time);upper = (hi-last_time)/(record['time']-last_time)
                for key,g in groups.items():
                    a = g['last']+lower*(observation[key]['value']-g['last']);b = g['last']+upper*(observation[key]['value']-g['last'])
                    mean = (a+b)/2;centered = w*(a-b)**2/12
                    if weight == 0:g['mean'] = mean;g['m2'] = centered
                    else:
                        delta = mean-g['mean'];g['mean'] = g['mean']+(w/new_weight)*delta
                        g['m2'] = g['m2']+centered+(weight*w/new_weight)*delta**2
                    require(np.all(g['m2'] >= 0) and np.all(np.isfinite(g['mean'])) and np.all(np.isfinite(g['m2'])),'invalid centered scientific moment')
                segments.append(dict(left_step=last_step,right_step=record['step'],left=last_time,right=record['time'],
                    lo=lo,hi=hi,weight=w,lower_fraction=lower,upper_fraction=upper))
                weight,covered = new_weight,hi
            for key in groups:groups[key]['last'] = observation[key]['value'].copy()
        last_time,last_step,last_hashes = record['time'],record['step'],hashes
        inputs.append(dict(**identity(path),step=last_step,time=last_time,hashes=hashes,origin=frame['origin']))
    require(last_step == finish and last_time == control['time'] and all(last_hashes[key] == control[key] for key in HASHES),'scientific generation does not end at native endpoint')
    require(abs(weight-(covered-begin)) <= 128*np.finfo(float).eps*(end-begin),'scientific window weights do not close')
    metadata = dict(schema=SCHEMA,status='PASS_SCIENTIFIC_ACCUMULATOR_BATCH_ONLY',full_H27_qualification=False,H27_D_pass=False,TASK3_READY=False,
        window=[begin,end],weight=weight,covered_until=covered,window_complete=covered == end,statistics_defined=weight > 0,
        last_time=last_time,last_step=last_step,last_native_hashes=last_hashes,layout_sha256=layout_hash,
        source_identity=control['source_identity'],reference_identity=control['reference_identity'],signature=control['signature'],dt=dt,
        sampling='IndependentSignedGauss5EndpointsV1',math='WeightedLinearSegmentWelfordV1',
        previous_generation=identity(previous_path) if prior is not None else None,
        native_result=result_id,native_index=index_id,session=session_id,segments=segments,frames=inputs,
        variance_policy='M2 accumulates nonnegative centered linear-segment moments; no floor or raw-moment cancellation',
        observations_are_scaled=True,raw_scientific_frames_deleted=False,continuum_supremum_claim=False,
        helper_sources=helper_sources)
    return save_generation(prefix,metadata,groups)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('manifest','table','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True);args=parser.parse_args();began=time.monotonic()
    require(sys.platform.startswith('linux') and sys.byteorder == 'little' and math.isfinite(args.seconds) and 0 < args.seconds <= 120,
            'bounded authoritative little-endian Linux analysis required')
    soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=min([8*1024**3]+[v for v in (soft,hard) if v != resource.RLIM_INFINITY]);resource.setrlimit(resource.RLIMIT_AS,(cap,hard))
    try:
        path=run(read(args.manifest),args.manifest.resolve().parent,args.table.resolve(),args.output.resolve(),began+args.seconds)
        print(json.dumps(dict(status='PASS_SCIENTIFIC_ACCUMULATOR_BATCH_ONLY',generation=identity(path),elapsed_seconds=time.monotonic()-began,
            maximum_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)))
    except Exception as error:
        failure=Path(str(args.output)+'.failure.json')
        if not failure.exists():
            with failure.open('x') as stream:json.dump(dict(schema=SCHEMA,status='FAIL',exception=type(error).__name__,reason=str(error),
                manifest=identity(args.manifest),table=identity(args.table),source=identity(__file__),elapsed_seconds=time.monotonic()-began),stream,indent=2);stream.write('\n')
        raise


if __name__ == '__main__':main()
