#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded transaction and seen-grid endpoint rechecks after WP2 trajectories."""
import argparse,os,subprocess,time
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cao_tamura_budget_wp2 import charge
from run_cls_epawf_picard import require,sha,state_hash,THREADS
from wp2_qualification import physical_gates


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','probe','h1','preflight']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();r=a.root;record=r/'final-checks.json';require(not record.exists(),'final checks already have an owner')
    state=dict(version='WP2FinalChecksV1',status='waiting_for_completed_waves',results={});atomic(record,state)
    try:
        while True:
            require(time.time()-load(r/'origin.json')['started_epoch']<43200,'WP2 wall envelope')
            d=load(r/'completion-waves.json');require(d['status']!='terminal_stop','completion waves stopped: '+str(d.get('failure')))
            if d['status']=='completed_bounded_waves':break
            time.sleep(15)
        jobs=[];env=dict(os.environ,**THREADS)
        for model in ['CLS','standard']:
            command=['/usr/bin/python3',str(Path(__file__).with_name('verify_wp2_long_transaction.py')),
                '--root',str(r),'--solver',str(a.solver),'--h1',str(a.h1),'--preflight',str(a.preflight),
                '--output',str(r/f'transit-transaction-{model}'),'--model',model]
            atomic(r/f'transit-transaction-{model}-command.json',command)
            with (r/f'transit-transaction-{model}.log').open('x') as log:
                job=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            jobs.append((model,job));state['results'][model]=dict(transaction_pid=job.pid)
        state['status']='running_transit_transactions';atomic(record,state)
        while any(job.poll() is None for _,job in jobs):time.sleep(15)
        require(all(job.returncode==0 for _,job in jobs),'transit transaction failed')
        for model,_ in jobs:
            transaction=r/f'transit-transaction-{model}/progress.json'
            require(load(transaction)['status']=='passed_transit_transaction','transaction incomplete')
            state['results'][model]['transaction']=dict(path=str(transaction),sha256=sha(transaction))
            source=r/f'preservation-fine-{model}/step120/case.evolution.json';d=load(source)
            xml=source.with_suffix('').with_suffix('.xml')
            require(xml.name=='case.xml' and xml.exists(),'endpoint XML missing')
            packet=dict(operation='wp2_equilibrium',inner_stress_model='Linear',expected_model='CLSHighRe' if model=='CLS' else 'StandardKEpsilon',
                xml_sha256=sha(xml),state=d['state'],state_sha256=state_hash(d['state']),order=28,boundary_checks=False,
                purpose='terminal repeat on the GL28 grid seen at spatial epoch0; not a fresh sampling claim',
                physical_time=d['time'],parent_report=str(source),parent_report_sha256=sha(source),
                parent_checkpoint=d['checkpoint'],parent_checkpoint_sha256=sha(Path(d['checkpoint'])/'State.bin'),
                fixed_reference_sha256=state_hash(d['fixed_reference']),fixed_defect_sha256=state_hash(d['fixed_defect']),
                native_history_sha256=[state_hash(v) for v in d['native_history']],native_history_times=d['native_history_times'])
            charge(r,'terminal seen GL28 '+model,spatial_cases=1)
            result,receipt=run_probe(a.probe,xml,packet,r/f'terminal-GL28-{model}')
            K=d['parameters'][5];checks=physical_gates(result['spatial'],d['ustar'],K,complete=True)
            require(result['state_sha256']==packet['state_sha256'],'terminal audit changed state')
            require(result['base_sha256']==packet['fixed_reference_sha256'],'terminal audit changed the immutable base')
            require(checks['passed'],'terminal physical/GL28 criteria failed')
            errors=[]
            for witness in result['spatial']['independent_covariance_witnesses']:
                ev=np.linalg.eigvalsh(np.asarray(witness['tensor']));error=float(max(abs(ev-np.asarray(witness['eigenvalues']))));require(error<1e-12,'terminal independent tensor replay');errors.append(error)
            state['results'][model]['terminal_seen_grid_audit']=dict(receipt=receipt,checks=checks,independent_eigenvalue_replay_max=max(errors),
                source_report=str(source),source_report_sha256=sha(source),state_sha256=packet['state_sha256'],base_sha256=result['base_sha256'])
            atomic(record,state)
        state['status']='completed_bounded_final_checks'
    except Exception as e:state.update(status='terminal_stop',failure=str(e));raise
    finally:atomic(record,state)


if __name__=='__main__':main()
