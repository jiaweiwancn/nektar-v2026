#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Accepted-step correction budgets on fixed physical windows; no solver work."""
import argparse
import datetime
import json
from pathlib import Path
from analyse_cls_epawf_height_comparison import attempts_summary
from continue_cls_epawf_pimple import load_progress
from run_cls_epawf_picard import load,write_json,require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=dict(captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases={},
        convention='integrated k per unit span; rates divide the same window integral by its physical duration',
        diagnostic_governing_solves=0)
    for c in ('U05-graded','U05-one'):
        root=a.root/('selected-'+c)/'trajectory';progress=load_progress(root/'progress.json')
        prefix=load(Path(progress['prefix_report']));start=int(progress['prefix_steps']);end=int(progress['accepted_steps'])
        attempts={i+1:q for i,q in enumerate(prefix['pimple_attempts']) if i<start}
        text=(root/'steps.jsonl').read_text()
        if text and not text.endswith('\n'):text=text[:text.rfind('\n')+1]
        for row in map(json.loads,text.splitlines()):
            n=int(row['global_step'])
            if n<=end:
                require(n not in attempts,'duplicate accepted correction budget');attempts[n]=row['attempt']
        require(sorted(attempts)==list(range(1,end+1)),'noncontiguous accepted ledger')
        dt=progress['dt']
        def budget(first,last):
            work=attempts_summary([attempts[n] for n in range(first+1,last+1)],last*dt)
            rows=work.pop('increments');norms={k:[q[k] for q in rows if q.get(k) is not None]
                for k in ('raw_step_mass_norm','correction_mass_norm','correction_to_step_ratio',
                    'k_added_storage_residual','k_added_frozen_residual','k_corrected_frozen_linear_residual')}
            work['increment_diagnostics']={k:dict(minimum=min(v),maximum=max(v),sum=sum(v)) for k,v in norms.items() if v}
            duration=(last-first)*dt;work.update(first_step=first,last_step=last,start=first*dt,end=last*dt,duration=duration,
                signed_integral_rate=work['correction_integral']/duration)
            return work
        complete=budget(0,end)
        require(abs(complete['correction_integral']-progress['cumulative_bound_integral'])<1e-13,
            'native prefix/restart correction budgets were duplicated or lost')
        windows=[budget(n-32,n) for n in range(32,end+1,32)]
        if end>=160:windows.append(budget(80,160))
        out['cases'][c]=dict(status=progress['status'],accepted_steps=end,accepted_time=progress['accepted_time'],
            cumulative=complete,physical_windows=windows,
            milestone_cumulative={str(n):budget(0,n) for n in (32,64,80,160) if n<=end})
    write_json(a.output,out);print(a.output)


if __name__=='__main__':main()
