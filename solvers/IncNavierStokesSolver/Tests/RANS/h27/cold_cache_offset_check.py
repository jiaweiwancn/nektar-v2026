#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only current-data cache cold restart at an unscheduled refresh time.

Fixed physical700-cell source15 diagnostic: uninterrupted12, prefix4, cold8.
Require control step5 to reuse and cold step5 to rebuild; compare complete
endpoints at unchanged1e-8 scales. No full H27 or production qualification.
"""
import argparse
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
import physical_cache_history_check as history
import physical_oracle as po
import physical_transient_oracle as transient
import history_recovery_check as hr
import implicit_algorithm_check as fields

require, read = po.require, history.read


def clock(deadline): require(time.monotonic() < deadline, 'cold offset checker exhausted explicit cap')


def load(spec, root, control, ref, deadline, begin, end, checkpoint=None):
    path, evidence = history.input_path(spec['result'], root); packet = read(path); record = packet['state']
    session_path, session_id = history.input_path(spec['session'], root)
    xml = ET.parse(session_path).getroot(); parameters = {}; solver = {}
    for p in xml.findall('.//PARAMETERS/P'):
        key, value = ''.join(p.itertext()).split('=',1); require(key.strip() not in parameters,'duplicate parameter');parameters[key.strip()] = value.strip()
    for item in xml.findall('.//SOLVERINFO/I'):
        key = item.get('PROPERTY'); require(key not in solver,'duplicate solver info');solver[key] = item.get('VALUE')
    for key, value in dict(TimeStep=.005, NumSteps=end-begin, IO_InfoSteps=1, WJScientificSteps=3).items():
        require(float(parameters[key]) == value, 'offset session time/output differs: '+key)
    for key, value in dict(SolverType='WJNativeCoupledV2', WJNativeMode='Physical', WJImplicitSolver='CoupledJFNKV1',
            WJStepPreconditioner='CurrentDataDriftV1', WJOutputCoefficients='True', WJNonlinearTrace='True').items():
        require(solver[key] == value, 'offset solver differs: '+key)
    require(solver.get('WJTerminalOutput','False') == 'False', 'unexpected terminal-output policy')
    prefix = history.replace_suffix(path,'.result.json','')
    require(hr.resolve(root,solver['WJOutputPrefix']) == prefix and
            hr.identity(hr.resolve(root,solver['WJReferenceTable']))['sha256'] == po.TABLE_SHA, 'offset prefix/table differs')
    require((hr.resolve(root,solver['WJRestart']) == checkpoint) if checkpoint is not None else not solver.get('WJRestart'), 'offset restart input differs')
    require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True and
            packet['steady_mode'] is False and packet['restart'] is (checkpoint is not None), 'offset operation provenance differs')
    history.configuration(record,control)
    require(record['step'] == end and record['time'] == end*.005, 'offset endpoint differs')
    mesh = history.session_geometry(spec['session'],root,ref); owners = history.inventory(record,control,ref,mesh)
    require(len(packet['complete_step_seconds']) == end-begin and all(math.isfinite(v) and v > 0 for v in packet['complete_step_seconds']),
            'offset accepted timing inventory differs')
    log_path = Path(str(prefix)+'.steps.jsonl'); rows = [hr.native_json(line) for line in log_path.read_text().splitlines()]
    require([r['step'] for r in rows] == list(range(begin+1,end+1)), 'offset accepted log is incomplete')
    previous = None; event = None
    for i,row in enumerate(rows):
        clock(deadline); history.accepted(row,control,False); history.scientific_time(row['step'],row['time'],.005,'offset log')
        _,event = history.cache_events(row,event,i == 0,False,'CurrentDataDriftV1')
        if previous is not None:
            require(row['previous_sha256'] == previous['state_sha256'] and row['previous_advection_sha256'] == previous['advection_sha256'], 'offset history shift differs')
        previous = row
    require(all(rows[-1][key] == record[key] for key in history.HASH_NAMES), 'offset endpoint/log differ')
    states,binary = history.binary_history(Path(str(prefix)+'.final.bin'),record,.005,deadline)
    require(binary['step'] == end,'offset binary clock differs')
    return dict(packet=packet,record=record,rows=rows,states=states,values=history.exported(record,states),owners=owners,
        evidence=dict(result=evidence,session=session_id,log=hr.identity(log_path),binary=binary,whole_native_times={k:packet[k]for k in
            ['construction_seconds','initialisation_seconds','complete_advancing_seconds','complete_step_seconds','final_export_seconds']}))


def run(manifest,folder,table,deadline,out):
    require(manifest['schema'] == 'WJH27ColdCacheOffsetV1','wrong offset manifest')
    root = hr.resolve(folder,manifest.get('root','.')); ref = po.Reference(read(table));require(hr.identity(table)['sha256'] == po.TABLE_SHA,'wrong reference')
    cp,_ = history.input_path(manifest['control']['result'],root); control_record = read(cp)['state']
    require(control_record['source_identity'] == manifest['source_identity'] and
            control_record['coupled_newton']['preconditioner_policy'] == 'CurrentDataDriftV1','offset source/policy differs')
    control = load(manifest['control'],root,control_record,ref,deadline,0,12)
    prefix = load(manifest['prefix'],root,control_record,ref,deadline,0,4)
    checkpoint,checkpoint_id = history.input_path(manifest['cold_input'],root)
    require(checkpoint == Path(prefix['evidence']['binary']['path']) and checkpoint_id['sha256'] == prefix['evidence']['binary']['sha256'],
            'offset cold input is not the pre-launch prefix4 checkpoint')
    cold = load(manifest['cold'],root,control_record,ref,deadline,4,12,checkpoint)
    fourth,fifth = control['rows'][3:5]
    require(all(fourth[key] == prefix['record'][key] for key in history.HASH_NAMES),'offset prefix4 does not match uninterrupted step4')
    before = fifth['coupled_newton']['preconditioner_events'][0]; after = cold['rows'][0]['coupled_newton']['preconditioner_events'][0]
    require(before['reason'] == 'reuse' and before['refreshed'] is False and before['age_before'] > 0,
            'uninterrupted offset step is not a warm inverse reuse')
    require(after['reason'] == 'cold' and after['refreshed'] is True and after['age_before'] == 0,
            'offset cold process did not reconstruct its inverse')
    require(cold['rows'][0]['previous_sha256'] == prefix['record']['state_sha256'] and
            cold['rows'][0]['previous_advection_sha256'] == prefix['record']['advection_sha256'], 'offset cold step lost input history')
    out['runs'] = {name:value['evidence']for name,value in dict(control=control,prefix=prefix,cold=cold).items()}
    out['cache_lineages'] = dict(uninterrupted_first_base=before,cold_first_base=after,different_inverse_preparations=True)
    out['comparison'] = history.compare(control['values'],cold['values'],control_record,control['owners'],ref,deadline)
    out['current_equations'] = history.equation_report(manifest['current_equation_report'],root,
        Path(cold['evidence']['result']['path']),cold['record'],cold['values'],po.TABLE_SHA)
    out.update(status='PASS_COLD_CACHE_OFFSET_FOUNDATION_ONLY',foundation_gate_pass=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['manifest','table','out']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True);args=parser.parse_args();began=time.monotonic()
    if args.out.exists():parser.error('immutable output exists')
    out=dict(schema='WJH27ColdCacheOffsetCheckV1',status='FAIL',foundation_gate_pass=False,full_H27_qualification=False,H27_D_pass=False,
        TASK3_READY=False,workers=1,comparison_limit=1e-8,manifest=hr.identity(args.manifest),table=hr.identity(args.table),
        sources={k:hr.identity(v)for k,v in dict(checker=__file__,history=history.__file__,fields=fields.__file__,reference=po.__file__,
            transient=transient.__file__,JSON_reader=hr.__file__).items()})
    try:
        require(sys.platform.startswith('linux')and sys.byteorder == 'little'and math.isfinite(args.seconds)and 0 < args.seconds <= 120,
                'bounded authoritative little-endian Linux analysis required')
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=min([8*1024**3]+[v for v in [soft,hard]if v != resource.RLIM_INFINITY]);resource.setrlimit(resource.RLIMIT_AS,(cap,hard));out['memory_cap_bytes']=cap
        run(read(args.manifest),args.manifest.resolve().parent,args.table.resolve(),began+args.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds']=time.monotonic()-began;out['maximum_resident_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    with args.out.open('x')as stream:json.dump(out,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({key:out[key]for key in ['status','reason','elapsed_seconds']if key in out}))
    if not out['foundation_gate_pass']:raise SystemExit(2)


if __name__ == '__main__':main()
