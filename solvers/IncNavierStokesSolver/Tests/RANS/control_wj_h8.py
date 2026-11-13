#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One-case H8 coordinator. Every consuming child uses the persistent ledger."""
import argparse
from pathlib import Path
import subprocess
import sys
import time
from run_wj_h8 import execute,load,save,sha
from study_wj_h5 import identity
import study_wj_h8 as study


def command(root,label,args,seconds,**kw):
    start=time.monotonic();result=execute(root,label,args,seconds=seconds,**kw)
    # execute includes its child/mapping/output work. Charge parent admission,
    # metadata and serialization overhead too; never reclaim a reservation.
    from run_wj_h7 import update_journal
    elapsed=time.monotonic()-start
    def charge(book):
        v=book['entries'][label];v['seconds']=max(v['seconds'],elapsed)
    update_journal(root/'accounting.json',charge)
    return result


def passive(root,label,mode,seconds=20,**settings):
    args=['python3','-B',study.__file__,'--root',str(root),'--mode',mode]
    for k,v in settings.items():
        args+=['--'+k]if v is True else ['--'+k,str(v)]
    return command(root,label,args,seconds)


def old_ready(root):
    assert load(root/'stage-receipt-v2.json')['status']=='PASS'
    for h in (5.,5.5,6.):
        v=load(root/f'diagnosis-view-H{h:g}.json');assert sha(v['path'])==v['sha256'];s=load(v['path'])
        assert s['status']==dict(coverage='PASS',quality='PASS'),'old-state qualification must pass'
    assert not load(root/'accounting.json').get('coverage_halt')


def views(root,new=False):
    if not new:old_ready(root)
    targets=[('quarter',4),('quarter',12),('quarter',16)]if new else [('coarse',1),('coarse',3),('fine',2),('fine',6)]
    for branch,step in targets:
        if not passive(root,f'request-{branch}-{step}','view-request',branch=branch,step=step):return False
        folder=root/'views'/f'{branch}-{step}';v=load(folder/'context.json')
        if not command(root,f'strict-{branch}-{step}',[str(root/'runtime/RANSWallResolutionProbe'),v['xml']['path'],str(folder/'request.json'),str(folder/'native.json')],120,strict_view=True):return False
        if not passive(root,f'finish-{branch}-{step}','view-finish',seconds=120,branch=branch,step=step):return False
        if load(folder/'result.json')['status']!='PASS':return False
    return True


def admission(root):
    old_ready(root);book=load(root/'accounting.json')
    for name in ('coarse-1','coarse-3','fine-2','fine-6'):assert load(root/'views'/name/'result.json')['status']=='PASS'
    spent=book['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    # H7 4/8-step costs imply an8.9s startup +4.225s/step. Use92s rather
    # than76.5s for16steps, plus current imports/diagnostic/strict costs.
    import_cost=load(root/'operations/old-scan-6.0/receipt.json')['seconds']
    view_cost=max(sum(book['entries'][k]['seconds']for k in (f'request-{b}-{s}',f'strict-{b}-{s}',f'finish-{b}-{s}'))for b,s in (('coarse',1),('coarse',3),('fine',2),('fine',6)))
    costs=dict(cp0=12.,prepare_and_pair=7.,governing16=92.,final_union_scan=max(60.,import_cost*31/14*1.15),strict_views3=3*view_cost,publication_and_closure_algebra=50.,parent_io=15.)
    total=spent+1.25*sum(costs.values())+180
    value=dict(status='PASS'if total<=1200 else 'RESOURCE_LIMIT',spent=spent,remaining_costs=costs,headroom=1.25,untouched_closure=180.,projected_total=total,ceiling=1200,
        sixteen_step_per_invocation_limit=120.,H7_measured_costs=dict(steps4=25.8,steps8=42.7,CP0=9.6,strict_combined=35.2),old_scan_seconds=import_cost,measured_strict_combined_max=view_cost,
        source=identity(__file__),ledger=identity(root/'accounting.json'))
    assert not(root/'admission.json').exists(),'write-once admission; no reset'
    save(root/'admission.json',value);print('H8 admission',value,flush=True);return value['status']=='PASS'


def native(root):
    if not admission(root):return False
    if not passive(root,'prepare-quarter','prepare'):return False
    if not command(root,'quarter-cp0',[str(root/'runtime/IncNavierStokesSolver'),str(root/'case/cp0/case.xml')],60,cp0=True,cwd=root/'case/cp0'):return False
    if not passive(root,'initial-pair','pair'):return False
    ok=command(root,'quarter-block16',[str(root/'runtime/IncNavierStokesSolver'),str(root/'case/run/case.xml')],120,governing=True,cwd=root/'case/run')
    if(root/'case/run/case.evolution.json').exists():
        scan=passive(root,'final-six-mm-union','diagnose',seconds=180,height=6.,version='final-v1',new=True)
        return ok and scan and not load(root/'accounting.json').get('coverage_halt')
    return False

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);p.add_argument('--mode',required=True,choices=('old-views','native','new-views'));a=p.parse_args()
    value=views(a.root)if a.mode=='old-views'else native(a.root)if a.mode=='native'else views(a.root,True)
    raise SystemExit(0 if value else 1)
