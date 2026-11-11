#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native size/seed, strict pressure-space restart and rollback qualification.

All governing solves in this script are independent verification work. The
production comparison is launched separately from fresh native h0.
"""
import argparse
from datetime import datetime,timezone
import os
from pathlib import Path
import platform
import subprocess
import time
import xml.etree.ElementTree as ET

from analyse_cls_accepted_awf import normalized
from analyse_cls_epawf_height_comparison import differences
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import CampaignLock,atomic,process_identity
from diagnose_cls_accepted_covariance import checkpoint
from run_cls_accepted_awf import case_xml,check_report
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,write


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('solver','baseline-initial-report','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--timeout',type=float,default=1800);a=p.parse_args()
    require(platform.system()=='Linux','authoritative Linux native verification required')
    os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output)
    out=dict(status='started',started_UTC=datetime.now(timezone.utc).isoformat(),solver_sha256=sha(a.solver),
        script_sha256=sha(Path(__file__)),production_governing_solves=0,verification_governing_solves_per_equation=0,
        receipts={},checks={},threads=THREADS,ranks=1)
    atomic(a.output/'summary.json',out)
    def run(name,steps,mode='inner-p2-v1',restart=None,hard_reject=0,expected=True):
        folder=a.output/name
        xml=case_xml('Accepted-AWF-05G',folder,steps,pressure_space=mode,restart=restart,
            companion='regularized',companion_reason='matched inner-pressure qualification; preserve faithful covariance blocker',
            inner_refinement='cutoff-v1',hard_reject=hard_reject,checkpoint_every=1)
        tree=ET.parse(xml).getroot();set_info(tree,'RANSWallEvaluationReuse','ExactStepV1');write(tree,xml)
        start=time.monotonic();env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(a.solver.parent/'lib'),TMPDIR=str(folder))
        with (folder/'solver.log').open('w') as log:
            child=subprocess.Popen(['/usr/bin/time','-f','%e %U %S %M','-o',str(folder/'usage.txt'),str(a.solver),str(xml)],
                cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,pass_fds=(lock.stream.fileno(),))
            atomic(folder/'process.json',dict(process=process_identity(child.pid),command=[str(a.solver),str(xml)]))
            code=child.wait(timeout=a.timeout)
        file=folder/'case.evolution.json';d=load(file) if file.exists() else None
        receipt=dict(exit_code=code,seconds=time.monotonic()-start,xml=str(xml),xml_sha256=sha(xml),
            report=str(file),report_sha256=sha(file) if file.exists() else None,log_sha256=sha(folder/'solver.log'),
            usage=(folder/'usage.txt').read_text().strip())
        out['receipts'][name]=receipt;atomic(folder/'receipt.json',receipt);atomic(a.output/'summary.json',out)
        if expected:
            require(code==0 and d is not None and not d['failure'],'native qualification failed: '+name)
            check_report(d,bool(restart),load(folder/'definition.json'))
            require(d['evaluation_reuse']['mode']=='ExactStepV1' and d['audit_revision']=='NativeElementSideV2','unexpected audit/reuse')
        else:require(code!=0,'incompatible or deliberately rejected request unexpectedly passed')
        return d
    try:
        base=run('p1-initial',0,mode='baseline-v1');old=load(a.baseline_initial_report)
        for key in ('state','fixed_reference','fixed_defect','native_history_times'):
            out['checks']['baseline_archived_'+key]=state_hash(base[key])==state_hash(old[key])
        out['checks']['baseline_archived_history']=[state_hash(v) for v in base['native_history']]==[state_hash(v) for v in old['native_history']]
        require(all(out['checks'].values()),'default baseline initialization changed')
        new=run('p2-initial',0);cp0=Path(new['checkpoint']);cp0sha=sha(cp0/'State.bin');definition=load(cp0.parent/'definition.json')
        sizes=new['accepted_awf_initial']['native_degrees_of_freedom'];prefix=int(sum(sizes[:4]))
        out['checks']['native_sizes']=sizes==[784,784,784,400,187,54,1]
        out['checks']['startup_bound_noop']=new['initial_k_bound_check']['activations']==0
        out['checks']['nonpressure_seed_bitwise']=state_hash(new['state'][:prefix])==state_hash(base['state'][:prefix])
        delta=max(abs(x-y) for x,y in zip(new['state'][:prefix],base['state'][:prefix]));out['nonpressure_seed_coefficient_max_difference']=delta
        require(delta<=1e-12,'pressure control changed its unchanged initial differential fields')
        out['checks']['differential_initial_histories_bitwise']=[state_hash(v) for v in new['native_history']]==[state_hash(v) for v in base['native_history']]
        out['checks']['initial_history_time_labels_bitwise']=state_hash(new['native_history_times'])==state_hash(base['native_history_times'])
        initial_delta=differences(new['accepted_awf_initial'],base['accepted_awf_initial'])
        out['initial_physical_difference']=initial_delta;out['initial_common_normalized_difference']=normalized(initial_delta)
        require(max(out['initial_common_normalized_difference'].values())<1e-10,'initial physical fields/gauged pressure differ above roundoff scale')
        require(checkpoint(cp0)['signature']!=checkpoint(Path(base['checkpoint']))['signature'],'pressure spaces lack strict native identity')
        out['checks']['native_signature_distinct']=True;out['initialization_work']=new['initialization_lift_calls_iterations']
        out['initial_h0_checkpoint']=dict(path=str(cp0),sha256=cp0sha,state_sha256=state_hash(new['state']),
            reference_sha256=state_hash(new['fixed_reference']),history_sha256=[state_hash(v) for v in new['native_history']])
        whole=run('whole2',2,restart=cp0);first=run('split1',1,restart=cp0)
        second=run('restart1',1,restart=Path(first['checkpoint']))
        out['verification_governing_solves_per_equation']=4
        out['checks']['same_dt_restart_exact']=same_history(whole,second)
        out['checks']['restart_checkpoint_payload_exact']=sha(Path(whole['checkpoint'])/'State.bin')==sha(Path(second['checkpoint'])/'State.bin')
        rejected=run('forced-rejection',1,restart=cp0,hard_reject=1,expected=False)
        require(rejected and 'deliberate' in rejected['failure'].lower(),'wrong forced rejection reason')
        out['verification_governing_solves_per_equation']+=1
        out['checks']['rollback_exact']=same_history(new,rejected)
        out['checks']['actual_rejected_candidate_retained']=bool(rejected['rejected_state']) and state_hash(rejected['rejected_state'])!=state_hash(new['state'])
        out['forced_rejection_reason']=rejected['failure']
        for name,mode,cp in [('p1-into-p2','inner-p2-v1',Path(base['checkpoint'])),('p2-into-p1','baseline-v1',cp0)]:
            run(name,0,mode=mode,restart=cp,expected=False)
            text=(a.output/name/'solver.log').read_text();require(any(w in text.lower() for w in ('signature','shape','size','field')),'unexpected incompatible restart error')
            out['checks'][name+'_rejected']=True;out.setdefault('restart_rejection_reasons',{})[name]=text[-1500:]
        out['checks']['original_checkpoint_preserved']=sha(cp0/'State.bin')==cp0sha
        require(all(v for k,v in out['checks'].items() if k not in ('nonpressure_seed_bitwise','differential_initial_histories_bitwise')),'native inner-pressure qualification failed')
        out.update(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat());atomic(a.output/'summary.json',out)
    except BaseException as error:
        out.update(status='failed',failure=str(error));atomic(a.output/'summary.json',out);raise
    print(a.output/'summary.json')


if __name__=='__main__':main()
