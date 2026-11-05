#!/usr/bin/env python3
"""Compact cut-exchange diagnosis from a frozen R12 durable progress stream."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--progress',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    records=[json.loads(x) for x in a.progress.read_text().splitlines()]
    cuts={int(x['id']):x for x in records if x['kind']=='r12_original_cut'}
    passes=[x for x in records if x['kind']=='r12_qp_pass']
    oracle=[x for x in records if x['kind']=='r6_cut_audit']
    last=oracle[-1] if oracle else None
    identities=Counter(x['identity'] for x in cuts.values())
    working=Counter(tuple(sorted(cuts[int(i)]['identity'] for i in x['active_ids'])) for x in passes)
    lower=[x['nonnegative_dual_bound'] for x in passes]
    witness=next(x for x in records if x['kind']=='r12_candidate' and x['name']=='L1_witness')
    result=dict(progress=dict(path=str(a.progress.resolve()),bytes=a.progress.stat().st_size,
        sha256=hashlib.sha256(a.progress.read_bytes()).hexdigest()),
        complete_passes=len(passes),admitted_cuts=len(cuts),distinct_original_rows=len(identities),
        exact_repeat_admissions=sum(v-1 for v in identities.values()),
        repeated_working_sets=sum(v-1 for v in working.values()),
        final_stage=last['status'] if last else None,
        final_active_rows=len(last['rows']) if last else None,
        minimum_final_QR=min(last['qr_diagonal'],default=None) if last else None,
        nonnegative_dual_bound_first=lower[0] if lower else None,
        nonnegative_dual_bound_last=lower[-1] if lower else None,
        largest_dual_decrease=max((u-v for u,v in zip(lower,lower[1:])),default=0.),
        feasible_upper_Psi=witness['Psi_Lf'],
        final_upper_lower_gap=witness['Psi_Lf']-lower[-1] if lower else None,
        entering_families=dict(Counter(x['kind'] for x in cuts.values())),
        first_passes=passes[:4],last_passes=passes[-8:],
        outcomes=[x for x in records if x['kind']=='r12_frozen_outcome'])
    a.out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
