#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Sequential H10 build-qualified admission; never retry a native reservation."""
import argparse
from pathlib import Path
import time
from run_wj_h10 import execute,load,save
from run_wj_h7 import update_journal
from wj_h10_analysis import MATRIX
import study_wj_h10 as study


def command(root,label,argv,seconds,**kw):
    t=time.monotonic();ok=execute(root,label,argv,seconds=seconds,**kw)
    def finish(book):book['entries'][label]['seconds']=max(book['entries'][label]['seconds'],time.monotonic()-t)
    update_journal(root/'accounting.json',finish);return ok


def passive(root,label,mode,seconds=45,**args):
    argv=['python3','-B',study.__file__,'--root',str(root),'--mode',mode]
    for k,v in args.items():argv+=['--'+k,str(v)]
    return command(root,label,argv,seconds)


def view(root,branch,phase):
    if not passive(root,f'{branch}-{phase}-view-request','view-request',branch=branch,phase=phase,seconds=10):return False
    p=root/f'science/views/{branch}-{phase}';ctx=load(p/'context.json')
    if not command(root,f'{branch}-{phase}-strict',[str(root/'runtime/RANSWallResolutionProbe'),ctx['xml']['path'],str(p/'request.json'),str(p/'native.json')],90,kind='strict'):return False
    if not passive(root,f'{branch}-{phase}-view-finish','view-finish',branch=branch,phase=phase,seconds=60):return False
    return load(p/'result.json')['status']=='PASS'


def projection(root,tag):
    book=load(root/'accounting.json');entries=book['entries'];spent=sum(v.get('seconds',v['reserved_seconds'])for v in entries.values()if v['phase']=='science')
    # After the retained first CP0, use its actual complete strict/scan costs.
    measured=lambda key:entries[key]['seconds']
    cp=measured('coarse-cp0');scan=measured('coarse-cp0-scan');strict=measured('coarse-cp0-strict')+measured('coarse-cp0-view-finish')+measured('coarse-cp0-view-request')
    left_cp=3-book['cp0_reservations'];left_steps=3-book['governing_attempts'];left_views=6-book['strict_readbacks']
    costs=dict(cp0=left_cp*max(cp,15.),governing=left_steps*max(cp+10,25.),scans=(left_cp+left_steps)*max(scan,10.),strict_views=left_views*max(strict,30.),parent_and_outputs=20.)
    remaining=1.25*sum(costs.values());total=book['setup_seconds']+sum(v.get('seconds',v['reserved_seconds'])for v in entries.values())+remaining+book['closure_reserve']
    value=dict(status='PASS'if spent+remaining<=900 and total<=3600 else 'RESOURCE_LIMIT',science_spent=spent,remaining_costs=costs,headroom=1.25,remaining_with_headroom=remaining,projected_science=spent+remaining,projected_total=total,closure_untouched=300,science_ceiling=900,total_ceiling=3600)
    save(root/f'science/admission-{tag}.json',value);print('H10 remaining admission',tag,value,flush=True);return value['status']=='PASS'


def native(root):
    assert load(root/'science/fixture-v2/native.json')['status']=='PASS'
    assert load(root/'operations/focused-tests-v2/receipt.json')['status']=='PASS'
    for branch,dt in MATRIX:
        if branch!='coarse'and not projection(root,branch+'-cp0'):return False
        folder=root/f'science/cases/{branch}/cp0'
        if not command(root,branch+'-cp0',[str(root/'runtime/IncNavierStokesSolver'),str(folder/'case.xml')],60,kind='cp0',cwd=folder):return False
        if not passive(root,branch+'-cp0-scan','scan',branch=branch,phase='cp0',seconds=60):return False
        if load(root/f'science/analysis/{branch}-cp0/result.json')['status']!='PASS':return False
        if not view(root,branch,'cp0'):return False
        if not projection(root,branch+'-step'):return False
        folder=root/f'science/cases/{branch}/run'
        ok=command(root,branch+'-step',[str(root/'runtime/IncNavierStokesSolver'),str(folder/'case.xml')],90,kind='governing',cwd=folder)
        if not ok:return False
        if not passive(root,branch+'-step-scan','scan',branch=branch,phase='run',seconds=60):return False
        if load(root/f'science/analysis/{branch}-run/result.json')['status']!='PASS':return False
        if not view(root,branch,'run'):return False
    return True

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();raise SystemExit(0 if native(a.root)else 1)
