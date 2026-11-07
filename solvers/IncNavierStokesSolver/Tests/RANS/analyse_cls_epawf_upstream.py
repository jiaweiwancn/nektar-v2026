#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Physical comparisons for fresh Upstream cases; optional readbacks stay separate."""
import argparse
import datetime
from pathlib import Path
from run_cls_epawf_picard import load,write_json,sha
from run_cls_epawf_height_comparison import UPSTREAM_CASES
from analyse_cls_epawf_height_comparison import differences,evolution_difference,compact_diagnostic,attempts_summary

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--screening',type=Path,required=True)
    p.add_argument('--uniform-screening',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=dict(schema=1,captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases={},paired={},policy_comparisons={},
        norms='dimensional physical L2 at identical physical quadrature; epsilon z>=2; one common weighted pressure gauge',
        positivity_scope='declared Q0 union permitted Q+ and retained/local/mean inventories, plus continuous ground; not volume positivity')
    states={};initial={}
    for c in UPSTREAM_CASES:
        path=a.screening/c/'result.json'
        if not path.exists():continue
        r=load(path);q0=r.get('initialization');out['cases'][c]=case=dict(status=r['status'],failure=r.get('failure'),definition=r['definition'])
        if not q0:continue
        epath=a.screening/c/'trajectory/case.evolution.json';e=load(epath);initial[c]=q0;states[c]={0:q0}
        states[c].update({int(n):load(v['path']) for n,v in r['observations'].items() if v['exit_code']==0})
        if e.get('height_comparison_final'):states[c][e['steps']]=e['height_comparison_final']
        case.update(steps=e['steps'],time=e['time'],native_report=dict(path=str(epath),sha256=sha(epath)),
            actual_h0_sha256=r['actual_h0_sha256'],native_history_hashes=r['initial_native_history_hashes'],
            initialization=compact_diagnostic(q0),initial_error_to_own_reference=differences(q0,q0,True),
            work=attempts_summary(e['pimple_attempts'],e['time']),observations={},qualification=r['qualification_events'])
        if e['failure']:
            case['terminal_attempt']=e['pimple_attempts'][-1]
            case['rejected_state_sha256']=__import__('run_cls_epawf_picard').state_hash(e['rejected_state'])
        for n,q in sorted(states[c].items()):
            if n:case['observations'][str(n)]=dict(time=q['time'],drift_from_actual_h0=differences(q,q0),
                error_to_own_reference=differences(q,q,True),diagnostic=compact_diagnostic(q))
    for left,right in [('U05-graded','U05-one'),('U20-graded','U20-one'),('U05-graded','U20-graded'),('U05-one','U20-one')]:
        if left not in states or right not in states:continue
        out['paired'][left+':'+right]={str(n):dict(time=states[left][n]['time'],fields=differences(states[left][n],states[right][n]),
            evolutionary_difference=evolution_difference(states[left][n],initial[left],states[right][n],initial[right]))
            for n in sorted(states[left].keys() & states[right].keys())}
    if a.uniform_screening:
        for c,u in [('U05-graded','G'),('U05-one','H'),('U20-graded','C'),('U20-one','D')]:
            if c not in states:continue
            r=load(a.uniform_screening/u/'result.json');q0=r['initialization'];uniform={0:q0}
            for n,v in r['observations'].items():uniform[int(n)]=load(v['path'])
            out['policy_comparisons'][c+':'+u]={str(n):dict(time=states[c][n]['time'],fields=differences(states[c][n],uniform[n]),
                evolutionary_difference=evolution_difference(states[c][n],initial[c],uniform[n],q0)) for n in sorted(states[c].keys() & uniform.keys())}
    write_json(a.output,out);print(a.output)

if __name__=='__main__':main()
