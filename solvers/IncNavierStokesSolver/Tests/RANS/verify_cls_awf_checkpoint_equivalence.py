#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify every retained per-step native payload, including GLM slots/times."""
import argparse
from pathlib import Path
from continue_cls_epawf_pimple import header,atomic
from run_cls_epawf_picard import load,sha,require


def payloads(folder):
    receipt=load(folder/'receipt.json');last=int(receipt['final']['steps']);first=last-len(receipt['packets'])
    result={}
    for file in folder.glob('checkpoint*/State.bin'):
        n=header(file.parent)['steps']
        if first<n<=last:
            digest=sha(file)
            require(n not in result or result[n]==digest,'conflicting native aliases at step '+str(n));result[n]=digest
    require(sorted(result)==list(range(first+1,last+1)),'missing native per-step checkpoint in '+str(folder))
    return {str(k):v for k,v in sorted(result.items())}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--benchmarks',type=Path,required=True)
    p.add_argument('--reference-results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',nargs='+',default=['BudgetH2','UpstreamOne']);a=p.parse_args()
    out=dict(status='passed',scope='bitwise native payload SHA256; coupled fields, every GLM history slot/time label, reference, defect and physical signature',cases={})
    for name in a.cases:
        root=a.benchmarks/name;data=dict(pairs=[]);reference=payloads(a.reference_results/name/'reference')
        for i in range(1,4):
            left=payloads(root/f'pair-{i}-Off');right=payloads(root/f'pair-{i}-ExactStepV1')
            require(left==right==reference,'per-step payload equality failed')
            data['pairs'].append(dict(pair=i,physical_steps=len(left),payloads=left,all_native_history_and_time_bits_equal=True))
        combined={}
        for i in range(1,5):combined.update(payloads(root/f'blocks8-{i}'))
        long=payloads(root/'long-32');require(combined==long,'8/32 block sizes changed a native history/checkpoint')
        data['block_sizes']=dict(physical_steps=len(long),payloads=long,all_native_history_and_time_bits_equal=True)
        out['cases'][name]=data
    atomic(a.output,out);print(a.output)

if __name__=='__main__':main()
