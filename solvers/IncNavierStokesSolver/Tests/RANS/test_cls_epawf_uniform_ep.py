#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused Linux uniform-flat policy/initialization and inherited-prefix checks."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import load,run,write_json,require,state_hash,sha,set_parameter,write
from run_cls_epawf_pimple import prepare,check_report
from continue_cls_epawf_pimple import header
from compare_cls_epawf_pimple import same_history


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--upstream-reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.solver=a.solver.resolve();a.probe=a.probe.resolve();a.output=a.output.resolve()
    a.output.mkdir(parents=True,exist_ok=False);summary=dict(schema=1,status='running',cases={})
    def case(name,steps,policy,restart=None):
        folder=a.output/name
        prepare(folder,.000625,steps,True,restart=restart,checkpoint=True,k_bound='MassProjection',floor_relative=1e-8,
            snapshots=True,inventory='ConstitutiveV2',wall_policy=policy)
        tree=ET.parse(folder/'case.xml').getroot();set_parameter(tree,'IO_CheckSteps',1);write(tree,folder/'case.xml')
        receipt=run([str(a.solver),'case.xml'],folder,'solver.log',14400);write_json(folder/'receipt.json',receipt)
        r=load(folder/'case.evolution.json');summary['cases'][name]=dict(exit_code=receipt['exit_code'],failure=r.get('failure'),steps=r.get('steps'),time=r.get('time'))
        write_json(a.output/'summary.json',summary)
        if r.get('initialization_complete') is not False:check_report(r,bool(restart))
        return r
    upstream=case('upstream6',6,'Upstream');prior=load(a.upstream_reference)
    require(not upstream['failure'] and upstream['steps']==6,'inherited prefix failed');same_history(upstream,prior)
    require(any(v['bound_activations'] for v in upstream['pimple_attempts']),'inherited prefix lacks active correction')
    uniform=case('uniform0',0,'UniformEPFlat')
    require(not uniform['failure'] and uniform['wall_policy']=='UniformEPFlat-v1','uniform native initialization failed')
    require(uniform['initial_k_bound_check']['activations']==0,'control modifies initial k')
    require(uniform['signature']!=upstream['signature'],'wall policy missing from checkpoint identity')
    initial=uniform['state'];old=upstream['accepted_states'][0]
    # All six block sizes come from actual native spaces; full state remains distinct.
    folder=a.output/'initial-contact';folder.mkdir();xml=a.output/'uniform0/case.xml'
    request=dict(operation='contact',state=initial,history=initial,time=0.,alpha=0.,h=[],require_history_h0=True,
        inner_stress_model='SaturatedR2',export_mass=True,isolate=False)
    write_json(folder/'input.json',request)
    rc=run([str(a.probe),str(xml),str(folder/'input.json'),str(folder/'result.json')],folder,'probe.log',3600)
    require(rc['exit_code']==0,'uniform initialization contact readback failed')
    contact=load(folder/'result.json');offsets=list(map(int,contact['offsets']));sizes=[b-a for a,b in zip(offsets,offsets[1:]+[len(initial)])]
    diffs=[max(abs(x-y) for x,y in zip(initial[o:o+n],old[o:o+n])) for o,n in zip(offsets,sizes)]
    require(diffs[0]==diffs[1]==diffs[2]==0,'intended initial velocity/constant k differs')
    require(contact['epsilon_trace_defect']<1e-12,'uniform initialized matching trace inconsistent')
    for operation in ('coverage','source_balance'):
        packet=dict(request,operation=operation,stations=[-3.5,-2.625,-1.,1.6],heights=[0.,.0001,.0025,.01,1.,2.])
        write_json(folder/(operation+'-input.json'),packet)
        receipt=run([str(a.probe),str(xml),str(folder/(operation+'-input.json')),str(folder/(operation+'.json'))],folder,operation+'.log',3600)
        require(receipt['exit_code']==0,operation+' readback failed')
    source=load(folder/'source_balance.json');selected=[v for v in source['points'] if v.get('variant')==-1]
    require(selected and all(v.get('valid') and v['chi']==1 for v in selected),'uniform source/profile fell back')
    require(all(abs(v['budget_identity_defect'])<1e-11 for v in source['ray_budgets']),'normalization budget mismatch')
    full=case('uniform2',2,'UniformEPFlat')
    require(not full['failure'] and full['steps']==2,'short uniform trajectory failed')
    cps={header(p.parent)['steps']:p.parent for p in (a.output/'uniform2').glob('checkpoint*/State.bin')}
    resumed=case('uniform-restart1',1,'UniformEPFlat',cps[1]);same_history(full,resumed)
    checkpoint_sha=sha(cps[1]/'State.bin')
    mismatch=case('uniform-as-upstream-refused',0,'Upstream',cps[1])
    require('checkpoint signature mismatch' in mismatch['failure'],'uniform checkpoint accepted under Upstream')
    require(checkpoint_sha==sha(cps[1]/'State.bin'),'immutable checkpoint changed')
    summary.update(status='passed',upstream_active_prefix_bitwise=True,uniform_initial_sha256=state_hash(initial),
        upstream_initial_sha256=state_hash(old),initial_block_max_differences=diffs,
        native_history_hashes=[state_hash(v) for v in uniform['native_history']],native_history_times=uniform['native_history_times'],
        uniform_initialized_matching_defect=contact['epsilon_trace_defect'],uniform_chi_one_all_queries=True,
        uniform_initial_bound_inactive=True,coverage_tripwire=load(folder/'coverage.json'),
        bitwise_restart_state_time_all_histories=True,incompatible_wall_checkpoint_refused=True)
    write_json(a.output/'summary.json',summary);print(a.output/'summary.json')


if __name__=='__main__':main()
