#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Coarse then fine then quarter, each fully audited before next admission."""
import argparse
from pathlib import Path
from run_wj_h11 import command,load,save,H10
from wj_h11_analysis import MATRIX,VIEWS
import study_wj_h11 as study


def passive(root,label,mode,seconds=90,**kwargs):
    argv=['python3','-B',study.__file__,'--root',str(root),'--mode',mode]
    for k,v in kwargs.items():argv+=['--'+k,str(v)]
    return command(root,label,argv,seconds)


def admission(root,tag):
    b=load(root/'accounting.json');entries=b['entries'];spent=sum(v.get('seconds',v['reserved_seconds'])for v in entries.values())
    left=[v for v in MATRIX if 'block:'+v[0]not in b['native_reservations']]
    costs=dict(governing=sum(max(15+4*n,entries.get('coarse-block',{}).get('seconds',0)*(n+1)/4)for branch,dt,n,target in left),
        scans=sum(max(20+2*n,entries.get('coarse-scan',{}).get('seconds',0)*(n+1)/4)for branch,dt,n,target in left),
        strict=sum(len(VIEWS[branch])*20 for branch,dt,n,target in left),final_union=150.,selected_arithmetic=20.,publication_verification=60.)
    remaining=1.25*sum(costs.values());status='PASS'if spent+remaining<=1200 else 'RESOURCE_LIMIT'
    save(root/f'science/admission-{tag}.json',dict(status=status,spent=spent,costs=costs,headroom=1.25,projected_operations=spent+remaining,projected_charge=60+spent+remaining+240,closure=240))
    print('H11 admission',tag,status,spent+remaining,flush=True);return status=='PASS'


def run(root):
    assert load(root/'operations/focused-tests/receipt.json')['status']=='PASS'
    assert load(root/'operations/imports/receipt.json')['status']=='PASS'
    if not admission(root,'initial'):return False
    for branch,dt,n,target in MATRIX:
        if branch!='coarse'and not admission(root,branch):return False
        case=study.parent(root,branch);folder=root/'science/cases'/branch
        if not command(root,branch+'-block',[str(H10/'runtime/IncNavierStokesSolver'),case['xml']['path']],120,kind='governing',cwd=folder):return False
        if not passive(root,branch+'-scan','scan',branch=branch,seconds=120):return False
        if study.audit(root,branch)['status']!='PASS':return False
        for step in VIEWS[branch]:
            if not passive(root,f'{branch}-{step}-request','view-request',branch=branch,step=step,seconds=10):return False
            out=root/f'science/views/{branch}-{step}'
            if not command(root,f'{branch}-{step}-strict',[str(H10/'runtime/RANSWallResolutionProbe'),case['xml']['path'],str(out/'request.json'),str(out/'native.json')],60,kind='strict'):return False
            if not passive(root,f'{branch}-{step}-finish','view-finish',branch=branch,step=step,seconds=45):return False
            if load(out/'result.json')['status']!='PASS':return False
    return True


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();raise SystemExit(0 if run(a.root)else 1)
