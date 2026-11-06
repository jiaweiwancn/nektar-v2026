#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable V1 failures and V2 projection-only coverage, never migration."""
import argparse
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import load,run,write_json,set_info,write,require,state_hash
from run_cls_epawf_pimple import identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--legacy',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    primary=load(a.legacy/'main/steps-000006-000021/case.evolution.json')
    finer=load(a.legacy/'finer-from-h0/case.evolution.json')
    original=load(a.prior/'fixed-horizon/bound-0/case.evolution.json')
    summary=dict(schema=1,probe=identity(a.probe.resolve()),cases={},purpose='projection-only, not a production step or restart migration')
    def save():write_json(a.output/'summary.json',summary)
    def execute(name,state,pre,history,time,alpha,policy,operation='k_bound',witnesses=()):
        folder=a.output/name;folder.mkdir()
        tree=ET.parse(a.prior/'fixed-horizon/bound-0/case.xml').getroot()
        set_info(tree,'RANSWallPimpleKInventory',policy);write(tree,folder/'case.xml')
        request=dict(operation=operation,state=state,pre_k=pre,history=history,time=time,alpha=alpha,
            require_history_h0=False,inner_stress_model='SaturatedR2',projection_only=policy=='ConstitutiveV2',witnesses=list(witnesses))
        write_json(folder/'input.json',request)
        receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],folder,'probe.log',3600)
        write_json(folder/'receipt.json',receipt);require(receipt['exit_code']==0,'native replay failed: '+name)
        return load(folder/'result.json')
    for name,report in [('primary',primary),('finer',finer)]:
        packet=load(Path(report['pimple_attempts'][-1]['k_snapshot']))
        old=execute(name+'-legacy',packet['raw'],packet['pre_k'],packet['history'],packet['time'],packet['alpha'],'LegacyV1')
        require(old['k_bound']['success'] and old['scalar_domain']['admissible'],'legacy limited inventory changed')
        require(state_hash(old['corrected'])==state_hash(packet['corrected']),'legacy corrected candidate changed')
        require('nonnegative transported k' in old['complete_audit_failure'],'legacy missing-ray failure did not reproduce')
        new=execute(name+'-v2',packet['raw'],packet['pre_k'],packet['history'],packet['time'],packet['alpha'],'ConstitutiveV2',witnesses=[old['complete_audit_context']['chart']])
        b=new['k_bound'];require(b['success'] and new['scalar_domain']['admissible'],'V2 projection failed: '+b['failure'])
        require(new['complete_audit_failure']=='','V2 final consumer audit failed: '+new['complete_audit_failure'])
        require(new['non_k_change']==new['essential_k_change']==0 and new['shared_cg_jump']<1e-12,'invariance/continuity changed')
        require(b['stationarity_mass_dual']<1e-10 and b['primal_floor_defect']<1e-12,'V2 QP validation failed')
        # The exact old witness must be in the declared family, not an ad hoc cut.
        witness=old['complete_audit_context']['chart']
        locations=new['scalar_domain']['observations']
        require(any(v['block']==2 and abs(v['coordinate'][0]-witness[0])<1e-12 for v in locations),
            'missing the old normalization-ray family')
        require(new['witness_values'][0]['corrected']>=0,'old physical witness remains negative')
        summary['cases'][name]=dict(legacy_failure=old['complete_audit_failure'],legacy_witness=old['complete_audit_context'],
            raw_sha256=state_hash(packet['raw']),legacy_corrected_sha256=state_hash(packet['corrected']),
            revised_corrected_sha256=state_hash(new['corrected']),bound=b,complete_audit_passed=True,
            pre_k_inventory=new['pre_k_inventory'],retained_physical_witness=new['witness_values'],
            non_k_change=new['non_k_change'],essential_k_change=new['essential_k_change'],shared_cg_jump=new['shared_cg_jump'])
        save()
    s=primary['state'];t=primary['time']
    contact=execute('accepted12-legacy-contact',s,s,s,t,0.,'LegacyV1','contact')
    require(not contact['contacts'][0]['point_source_valid'],'legacy diagnostic-ray failure disappeared')
    new=execute('accepted12-v2',s,s,s,t,0.,'ConstitutiveV2')
    require(new['k_bound']['success'] and not new['complete_audit_failure'],'accepted V1 fixture was not repaired')
    contact_new=execute('accepted12-v2-contact',new['corrected'],new['corrected'],new['corrected'],t,0.,'ConstitutiveV2','contact')
    require(all(v['point_source_valid'] for v in contact_new['contacts']),'V2 ground-minimum contact ray still invalid')
    summary['cases']['legacy_accepted12']=dict(old_contact_failure=contact['contacts'][0],bound=new['k_bound'],
        repaired_contact=contact_new['contacts'][0],production_restart_accepted=False)
    for name,s in [('actual_h0',original['accepted_states'][0]),('neutral_reference',original['fixed_reference'])]:
        r=execute(name,s,s,s,0.,0.,'ConstitutiveV2')
        require(r['k_bound']['success'] and not r['k_bound']['activations'] and state_hash(r['corrected'])==state_hash(s),'inactive V2 changed '+name)
        summary['cases'][name]=dict(sha256=state_hash(s),bitwise_noop=True,bound=r['k_bound'])
    s=original['accepted_states'][0]
    summary['coverage']=execute('consumer-tripwire',s,s,s,0.,0.,'ConstitutiveV2','coverage')
    summary['status']='passed';save();print(a.output/'summary.json')


if __name__=='__main__':main()
