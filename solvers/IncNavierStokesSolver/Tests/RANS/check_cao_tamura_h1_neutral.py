#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One frozen neutral restart step: old runtime versus H1 runtime, bitwise."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import platform
from continue_cls_epawf_pimple import CampaignLock,atomic,header
from run_cls_awf_late_ivp import prepare,native,packet
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import require,state_hash,sha
from run_cls_epawf_pimple import identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','baseline','parent','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();require(platform.system()=='Linux','native Linux regression required')
    a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output)
    before=sha(a.parent/'State.bin');q=header(a.parent)
    def one(label,solver):
        xml=prepare(a.output/label,1,q['dt'],restart=a.parent,refinement='cutoff-v1')
        return native(solver.resolve(),xml,1800,lock)
    with ThreadPoolExecutor(max_workers=2) as pool:
        old=pool.submit(one,'baseline',a.baseline);new=pool.submit(one,'H1',a.solver)
        aa,ra=old.result();bb,rb=new.result()
    checks=dict(accepted_state_history_time=same_history(aa,bb),signature=aa['signature']==bb['signature'],
        fixed_reference=state_hash(aa['fixed_reference'])==state_hash(bb['fixed_reference']),
        fixed_defect=state_hash(aa['fixed_defect'])==state_hash(bb['fixed_defect']),
        parent_bytes=sha(a.parent/'State.bin')==before)
    pa,pb=packet(aa),packet(bb)
    for k in ('committed','post_flow','pre_k','raw','corrected','history','fixed_reference'):
        checks[k]=state_hash(pa[k])==state_hash(pb[k])
    checks['audit_values']=aa['state_audits']==bb['state_audits']
    require(all(checks.values()),'H1 changed default neutral closure/stages/residual/history: '+str(checks))
    atomic(a.output/'summary.json',dict(status='passed',checks=checks,parent_sha256=before,
        baseline=identity(a.baseline),H1=identity(a.solver),receipts=[ra,rb],
        physical_steps_per_independent_leg=1,new_governing_solves_per_equation=2,
        interpretation='short frozen neutral regression only; no long neutral trajectory repeated'))
    print(a.output/'summary.json')


if __name__=='__main__':main()
