#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Replay all retained late-segment transaction packets without a solver."""
import argparse
import os
from pathlib import Path
import time
from wj_immutable import Ledger,read,write_once,artifacts,verify_artifacts,digest as digest_object
from run_cls_epawf_picard import require,state_hash
from compare_cls_epawf_pimple import same_history


def verify(whole,first,second,segment_length=32):
    require(not isinstance(segment_length,bool) and isinstance(segment_length,int) and segment_length>0,
            'transaction segment length must be a positive integer')
    require(len(whole['accepted_states'])==segment_length+1,'transaction accepted-state coverage differs')
    hashes={}
    for key in ('accepted_states','audit_times','state_audits','wj_audit_packets'):
        combined=first[key]+second[key][1:]
        require(digest_object(whole[key])==digest_object(combined),'all-step transaction mismatch: '+key)
        hashes[key]=digest_object(combined)
    joined=first['pimple_attempts']+second['pimple_attempts']
    require(len(whole['pimple_attempts'])==len(joined)==segment_length,'late segment scheduled attempt coverage differs')
    columns=('history_sha256','incoming_sha256','time','alpha','governing_solves_complete','scheduled_complete',
             'momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves',
             'bound_calls','bound_activations','additional_governing_equation_solves','additional_correction_or_restoration_solves')
    for a,b in zip(whole['pimple_attempts'],joined):
        require(digest_object({k:a[k] for k in columns})==digest_object({k:b[k] for k in columns}),'scheduled state/history/time/count mismatch')
        sa,sb=read(a['k_snapshot']),read(b['k_snapshot'])
        require(all(state_hash(sa[k])==state_hash(sb[k]) for k in ('raw','corrected')),'all-step raw/corrected k mismatch')
    require(same_history(whole,second),'final state/history mismatch')
    require(digest_object(whole['physical_samples'])==digest_object(second['physical_samples']),'final physical samples mismatch')
    return dict(status='PASS',accepted_states_including_startup=len(whole['accepted_states']),scheduled_attempts=segment_length,
        exact_all_step_state_time_audits=True,exact_all_step_callback_history_incoming_counts=True,
        exact_all_step_raw_corrected_k=True,exact_final_primary_gauge_glm_history_time_samples=True,sequence_sha256=hashes)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--policy',choices=('EP','Budget'),required=True);a=p.parse_args();root=a.study.resolve()/a.policy
    paths=[root/name/'case.evolution.json' for name in ('longer-coarse-256','transaction-first16','transaction-second16')]
    before=artifacts(paths+[p.parent/'receipt.json' for p in paths])
    require(len(before)==6,'missing native transaction report/receipt')
    for path in paths:verify_artifacts(read(path.parent/'receipt.json')['artifacts'])
    ledger=Ledger(Path(os.environ['WJ4_PACKAGE_ROOT']),a.policy);label=str(a.output.resolve());start=time.time()
    ledger.reserve(label,dict(inputs=before,code=artifacts([Path(__file__)])),0,30)
    try:
        result=verify(*(read(p) for p in paths));verify_artifacts(before)
        write_once(a.output,dict(schema='WJ4AllStepTransactionV1',policy=a.policy,inputs=before,**result,governing_attempts=0))
    finally:ledger.finish(label,0,time.time()-start,str(a.output))
    print(a.output)


if __name__=='__main__':main()
