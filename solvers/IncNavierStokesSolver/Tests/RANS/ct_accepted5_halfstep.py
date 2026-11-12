#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One conditional CT5 traction control, preserving the original pair receipt."""
import argparse,copy,os
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted5_run import Campaign,basic_gate,traction_gate,read
from ct_accepted_pilot import artifact_hash
from ct_accepted_open import require
from ct_accepted4_reference import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct4','runtime','diagnosis','contract','campaign','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');prior=read(a.campaign/'results.json')
    require(prior['models']['A95']['status']=='pass','primary_not_admitted')
    require(not (a.root/'initial-quarter').exists(),'single_conditional_halfstep_already_consumed')
    old=prior['models']['A93'];q=old['milestones']['0.008'];prep=old['initial']['initialization'];ih=prep['members']['0.001']['state_sha256']
    require(q['basic']==basic_gate(q['members'],q['temporal'],.008,(.002,.001),ih) and q['basic']['status']=='pass','nontraction_failure_not_eligible')
    require(q['traction_uncertainty']==traction_gate(q['members'],q['case'],.008) and q['traction_uncertainty']['status']=='ambiguous','traction_control_not_triggered')
    m=RequestManifest(a.output,dict(schema='CT5SingleTractionHalfstepV1',model='A93',original_campaign_sha256=sha(a.campaign/'results.json'),
        contract_sha256=sha(a.contract),trigger=q['traction_uncertainty'],runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),
        sources={f:sha(Path(__file__).with_name(f)) for f in ('ct_accepted5_halfstep.py','ct_accepted5_run.py','ct_accepted4_run.py','ct_accepted5_traction.py')},
        scope='one fresh .0005 history; original ambiguous .002/.001 retained; .02 conditional on complete .001/.0005 admission'))
    result=copy.deepcopy(prior);row=result['models']['A93'];row['milestones']['0.008-original-pair']=copy.deepcopy(q)
    camp=Campaign(a);camp.model='A93';camp.c=q['case'];camp.initial_hash=ih;ip=a.root/'initial/A93'
    seal=read(ip/'completed-artifacts.json');require(all(sha(ip/f)==h for f,h in seal['files'].items()),'changed_original_initialization')
    for dt,member in prep['members'].items():
        require(member['state_sha256']==ih,'changed_same_space_initial')
        camp.inits[dt]=(ip,prep)
        r=q['members'][dt];native=read(Path(r['directory'])/'case.evolution.json')
        require(sha(Path(r['directory'])/'case.evolution.json')==r['native_report_sha256'] and artifact_hash(native['checkpoint'])==r['history_sha256'],'changed_original_state_or_history')
        camp.members[(dt,.008)]=r
    row['quarter_initial']=camp.prepare('A93',(.0005,),True)
    fine=camp.pair(.008,(.001,.0005));row['milestones']['0.008']=fine
    write_json(a.output/'progress.json',result)
    if fine['qualification']['status']=='pass':row['milestones']['0.02']=camp.pair(.02,(.001,.0005))
    row['status']='pass' if all(t in row['milestones'] and row['milestones'][t]['qualification']['status']=='pass' for t in ('0.008','0.02')) else 'failed_or_unresolved'
    row['half_step_control']='executed_once_on_A93_for_ambiguous_interface_traction'
    result['single_conditional_halfstep']=dict(model='A93',dt=.0005,original_ambiguous_pair_preserved=True,original_campaign_sha256=sha(a.campaign/'results.json'))
    result['status']='pass' if all(v['status']=='pass' for v in result['models'].values()) else 'failed_or_unresolved'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT5 conditional control',row['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
