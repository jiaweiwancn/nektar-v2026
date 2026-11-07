#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native corrector checks; fixture fields are assembled by the native probe."""
import argparse
from pathlib import Path
import platform
from run_cls_epawf_picard import load,write_json,state_hash,require,sha,run
from run_cls_epawf_pimple import prepare
from test_cls_epawf_k_bound import replay


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--solver',type=Path)
    p.add_argument('--startup-report',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();require(platform.system()=='Linux','native Linux checks required')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);a.probe=a.probe.resolve()
    xml=prepare(a.output/'template',.000625,0,True,k_bound='MassProjection')
    if a.startup_report:startup=a.startup_report.resolve()
    else:
        require(a.solver is not None,'solver or startup report required')
        startup_xml=prepare(a.output/'startup',.000625,0,True)
        receipt=run([str(a.solver.resolve()),'case.xml'],startup_xml.parent,'solver.log',900)
        require(receipt['exit_code']==0,'startup initialization failed');startup=startup_xml.with_suffix('.evolution.json')
    base=load(startup)
    # A supplied report may include physical evolution. Its final state is not
    # h0; use the actual native initialized state, independently of reference.
    require(base.get('accepted_states'),'startup report lacks actual initialized state')
    initial=dict(state=base['accepted_states'][0]);results={}
    for fixture in ('constant','reference','smooth','between_nodes'):
        r=replay(a.probe,xml,initial,a.output/fixture,fixture=fixture);b=r['k_bound']
        require(b['success'],fixture+': '+b['failure'])
        if fixture!='between_nodes':require(state_hash(r['raw'])==state_hash(r['corrected']),'no-op failed')
        else:require(r['raw_native_ground_node_minimum']>0 and b['raw_ground']['minimum']<0 and b['corrected_ground']['minimum']>0,'hidden ground undershoot test failed')
        results[fixture]=dict(bound=b,native_node_minimum=r['raw_native_ground_node_minimum'],cg_jump=r['shared_cg_jump'])
    packet=dict(raw=initial['state'],pre_k=initial['state'],history=initial['state'],time=0,alpha=0)
    r=replay(a.probe,xml,packet,a.output/'h0');b=r['k_bound']
    require(b['success'] and state_hash(r['raw'])==state_hash(r['corrected']),'actual h0 no-op failed')
    results['h0']=b
    r=replay(a.probe,xml,initial,a.output/'incompatible',floor=1.,fixture='reference');b=r['k_bound']
    require(not b['success'] and 'prescribed' in b['failure'],'incompatible floor not handled cleanly')
    results['incompatible']=b
    write_json(a.output/'summary.json',dict(status='passed',startup=str(startup),startup_sha256=sha(startup),
        actual_initial_sha256=state_hash(initial['state']),startup_report_final_time=base['time'],cases=results))
    print('native corrector checks passed')


if __name__=='__main__':main()
