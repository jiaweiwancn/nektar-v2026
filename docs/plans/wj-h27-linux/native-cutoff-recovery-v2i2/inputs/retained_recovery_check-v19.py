#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check an advancing cold import and both retained native recoveries."""
import argparse,json,math,os,resource,sys,time
from pathlib import Path
sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import cold_cache_offset_check as cold
import history_recovery_check as hr
import physical_cache_history_check as history
import physical_oracle as po
import terminal_output_check as terminal
require,read = po.require,history.read


def run(manifest,folder,table,deadline,out):
    require(manifest['schema'] == 'WJH27RetainedRecoveryContinuationV1','wrong retained-history manifest')
    root = hr.resolve(folder,manifest['root']);require(hr.identity(table)['sha256'] == po.TABLE_SHA,'reference changed')
    ref = po.Reference(read(table));cp,_ = history.input_path(manifest['control']['result'],root);control_data = read(cp)['state']
    require(control_data['source_identity'] == manifest['source_identity'] and
        control_data['scalar_quadrature_policy'] == 'CubicCutoffFacePanelsV1','wrong source/panel recipe')
    control = cold.load(manifest['control'],root,control_data,ref,deadline,0,12)
    packet = control['packet'];prefix = history.replace_suffix(cp,'.result.json','')
    require(packet['rotated_recovery_count'] == 2 and len(packet['recovery_seconds']) == 4 and
        all(math.isfinite(v) and v > 0 for v in packet['recovery_seconds']),'four writes/two rotations not observed')
    expected = [Path(str(prefix)+'.recovery-'+str(step)+'.bin') for step in (9,12)]
    require([Path(p) for p in packet['retained_recovery']] == expected,'retained recovery names differ')
    removed = [Path(str(prefix)+'.recovery-'+str(step)+'.bin') for step in (3,6)]
    require(all(not p.exists() for p in removed),'native removed recovery reappeared')
    retained = {}
    for step in (9,12):
        path,ident = history.input_path(manifest['inputs'][str(step)],root)
        require(path == expected[(0 if step == 9 else 1)],'prelaunch retained input differs')
        scientific = Path(str(prefix)+'.science-'+str(step)+'.json');record = read(scientific)
        require(all(record[key] == control['rows'][step-1][key] for key in history.HASH_NAMES),'retained science/log history differs')
        states,binary = history.binary_history(path,record,.005,deadline,ident['sha256']);history.exported(record,states)
        require(binary['step'] == step,'retained recovery clock differs')
        retained[str(step)] = dict(input=binary,science=hr.identity(scientific),four_history_binding=True)
    incoming = expected[0];continued = cold.load(manifest['cold'],root,control_data,ref,deadline,9,12,incoming)
    require(continued['rows'][0]['previous_sha256'] == retained['9']['input']['hashes']['current'] and
        continued['rows'][0]['previous_advection_sha256'] == retained['9']['input']['hashes']['advection'],'cold advancement lost input histories')
    event = continued['rows'][0]['coupled_newton']['preconditioner_events'][0]
    require(event['reason'] == 'cold' and event['refreshed'] is True and event['age_before'] == 0,'cold advancing inverse was not rebuilt')
    out['comparison'] = history.compare(control['values'],continued['values'],control_data,control['owners'],ref,deadline)
    out['current_equations'] = history.equation_report(manifest['current_equations'],root,
        Path(continued['evidence']['result']['path']),continued['record'],continued['values'],po.TABLE_SHA)
    zp,zid = history.input_path(manifest['zero']['result'],root);zero = read(zp);z = zero['state'];history.configuration(z,control_data)
    _,session_id,parameters,solver = terminal.session(manifest['zero']['session'],root)
    require(zero['restart'] is True and zero['steady_mode'] is False and z['step'] == 12 and z['time'] == .06 and
        parameters['NumSteps'] == '0' and parameters['IO_CheckSteps'] == '0' and parameters['WJScientificSteps'] == '0' and
        solver.get('WJTerminalOutput','False') == 'False' and hr.resolve(root,solver['WJRestart']) == expected[1],'retained zero recipe differs')
    zprefix = history.replace_suffix(zp,'.result.json','');final = Path(str(zprefix)+'.final.bin')
    states,binary = history.binary_history(final,z,.005,deadline);history.exported(z,states,reload=True)
    require(final.read_bytes() == expected[1].read_bytes() and not Path(str(zprefix)+'.steps.jsonl').read_bytes() and
        zero['complete_step_seconds'] == [] and zero['complete_advancing_seconds'] == 0 and z['rebuilds'] == 0 and
        z['coupled_newton']['refreshes'] == z['coupled_newton']['reuses'] == 0,'retained zero serialization/work differs')
    receipt_path,receipt_id = history.input_path(manifest['negative']['receipt'],root);receipt = read(receipt_path)
    log_path,log_id = history.input_path(manifest['negative']['log'],root)
    old,old_id = history.input_path(manifest['negative']['input'],root)
    require(receipt['status'] == 'FAIL' and receipt['exit'] != 0 and
        'coupled checkpoint scientific identity differs' in log_path.read_text(),'prior-source checkpoint was not explicitly rejected')
    require(old_id['sha256'] != retained['12']['input']['sha256'],'negative input is not a different checkpoint')
    negative_prefix = hr.resolve(root,manifest['negative']['output_prefix'])
    require(not Path(str(negative_prefix)+'.result.json').exists() and not Path(str(negative_prefix)+'.steps.jsonl').exists(),
        'prior-source checkpoint reached advancement/output')
    out.update(status='PASS_RETAINED_RECOVERY_CONTINUATION_ONLY',foundation_gate_pass=True,retained_recoveries=retained,
        observed_writes=[3,6,9,12],observed_removed=[str(p)for p in removed],removed_binary_byte_audit_claimed=False,
        write_seconds=packet['recovery_seconds'],control=control['evidence'],continued=continued['evidence'],
        first_cold_inverse_event=event,uninterrupted_first_event=control['rows'][9]['coupled_newton']['preconditioner_events'][0],
        zero=dict(result=zid,session=session_id,binary=binary,four_binary_histories_exact=True,advanced_steps=0,inverse_factorizations=0),
        incompatible_source_rejection=dict(receipt=receipt_id,log=log_id,input=old_id,pre_advancement=True),
        scope='Short source19 physical cold continuation/retained payloads. Does not qualify monitoring escalation, persisted statistics, selected h/p/q/dt, physical development or complete cost.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','table','out'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True);args=parser.parse_args();began=time.monotonic()
    require(not args.out.exists(),'retained-history output is immutable')
    out=dict(schema='WJH27RetainedRecoveryContinuationCheckV1',status='FAIL',foundation_gate_pass=False,full_H27_qualification=False,TASK3_READY=False,
        manifest=hr.identity(args.manifest),table=hr.identity(args.table),workers=1,
        sources={key:hr.identity(path)for key,path in dict(checker=__file__,cold=cold.__file__,history=history.__file__,reference=po.__file__,sessions=terminal.__file__).items()})
    try:
        require(sys.platform.startswith('linux') and math.isfinite(args.seconds) and 0 < args.seconds <= 60,'bounded Linux assessment required')
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=min([8*1024**3]+[v for v in (soft,hard)if v != resource.RLIM_INFINITY]);resource.setrlimit(resource.RLIMIT_AS,(cap,hard))
        run(read(args.manifest),args.manifest.resolve().parent,args.table.resolve(),began+args.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds']=time.monotonic()-began;out['maximum_resident_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    with args.out.open('x')as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({key:out[key]for key in ('status','reason','elapsed_seconds')if key in out}))
    if not out['foundation_gate_pass']:raise SystemExit(2)


if __name__=='__main__':main()
