#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Accepted-step correction budgets on fixed physical windows; no solver work."""
import argparse
import datetime
import json
from pathlib import Path
from analyse_cls_epawf_height_comparison import attempts_summary
from continue_cls_epawf_pimple import load_progress
from run_cls_epawf_picard import load,write_json,require,sha


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
        # Sparse immutable raw/corrected packets quantify cancellation, rather
        # than inferring it merely from a correction/raw norm ratio near one.
        import numpy as np
        mass_file=a.root/'discrete-budget'/c/'mass-and-bands.json';mass=load(mass_file)
        M=np.array(mass['free_k_mass']);free=np.array(mass['free_k_rows'],dtype=int);samples={}
        norm=lambda v:float(np.sqrt(max(0.,v@M@v)))
        for n in sorted({32,64,80,160,32*(end//32)}):
            if n>end:continue
            attempt=attempts[n];packet=Path(attempt['k_snapshot'])
            require(sha(packet)==attempt['k_snapshot_sha256'],'retained increment packet changed')
            q=load(packet);require(np.array_equal(free,np.array(q['transport_rows'],dtype=int)),'consistent-mass free rows differ')
            pre=np.array(q['pre_k']);raw=np.array(q['raw']);corrected=np.array(q['corrected'])
            r=(raw-pre)[free];b=(corrected-raw)[free];net=(corrected-pre)[free]
            rn,bn,nn=map(norm,(r,b,net));defect=norm(net-r-b)
            require(defect<1e-13,'raw+bound increment identity failed')
            require(abs(rn-q['k_bound']['mass_norm_raw_step'])<1e-13 and
                abs(bn-q['k_bound']['mass_norm_correction'])<1e-13,'mass replay differs from native correction accounting')
            samples[str(n)]=dict(time=q['time'],packet=str(packet),packet_sha256=sha(packet),
                mass_file_sha256=sha(mass_file),raw_mass_norm=rn,bound_mass_norm=bn,accepted_mass_norm=nn,
                mass_angle_cosine=float(r@M@b/(rn*bn)) if rn*bn else None,
                accepted_over_raw=nn/rn if rn else None,increment_identity_mass_norm=defect)
        out['cases'][c]=dict(status=progress['status'],accepted_steps=end,accepted_time=progress['accepted_time'],
            cumulative=complete,physical_windows=windows,
            retained_increment_samples=samples,
            milestone_cumulative={str(n):budget(0,n) for n in (32,64,80,160) if n<=end})
    write_json(a.output,out);print(a.output)


if __name__=='__main__':main()
