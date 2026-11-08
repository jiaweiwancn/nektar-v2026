#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Plot retained native cubic fields at a declared common physical time."""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path

from analyse_cls_accepted_awf import dense_column, analytical_column
from continue_cls_epawf_pimple import atomic
from run_cls_accepted_awf import CASES
from run_cls_epawf_picard import THREADS, load, require, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--suite',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',nargs='+',choices=tuple(CASES),default=['Accepted-AWF-05G','Accepted-AWF-05O'])
    p.add_argument('--step',type=int,default=80);p.add_argument('--station',type=float,default=1.6)
    a=p.parse_args();os.environ.update(THREADS)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    require(a.step>0 and len(set(a.cases))==len(a.cases),'invalid matched observations')
    require(not a.output.exists(),'use a new plot path rather than overwrite evidence')
    fig,axes=plt.subplots(2,4,figsize=(15,8),layout='constrained')
    sources=[];times=[];heights=[]
    for case in a.cases:
        source=a.suite/case/'result.json';result=load(source);initial=result['initialization']
        entry=result['snapshots'][str(a.step)];snapshot=Path(entry['path'])
        require(entry['exit_code']==0 and sha(snapshot)==entry['sha256'],'retained observation identity differs')
        selected=load(snapshot)
        require(selected.get('valid') is not False and selected.get('common_reporting_queries_valid') is not False,
            'optional diagnostic invalid; never repair a retained field for plotting')
        times.append(selected['time']);heights.append(selected['H'])
        sources.append(dict(case=case,result=str(source.resolve()),result_sha256=sha(source),
            snapshot=str(snapshot),snapshot_sha256=sha(snapshot),actual_h0_sha256=result['actual_h0_sha256']))
        for state,when in [(initial,'native h0'),(selected,'T='+format(selected['time'],'.6g'))]:
            x,values=dense_column(state,a.station)
            for axis,field in zip(axes[0],('k','u','du_dz','dk_dz')):
                axis.plot(values[field],values['z'],label=case+' '+when)
            for axis,field in zip(axes[1],('k','u','du_dz','dk_dz')):
                axis.plot(values[field],[z/state['H'] for z in values['z']],label=case+' '+when)
        x,values=dense_column(initial,a.station);reference=analytical_column(initial,values['z'])
        for line,scale in zip(axes,[1.,initial['H']]):
            for axis,field in zip(line,('k','u','du_dz','dk_dz')):
                axis.plot(reference[field],[z/scale for z in values['z']],':',color='black')
    require(max(times)-min(times)<1e-12 and max(heights)-min(heights)<1e-12,
        'matched plot requires one common physical time and wall-model height')
    for line,ylabel in zip(axes,('physical z','z/H')):
        for axis,field in zip(line,('transported k','velocity u','du/dz','dk/dz')):
            axis.set(xlabel=field,ylabel=ylabel)
    for axis in axes[0]:axis.set_ylim(0,.1)
    for axis in axes[1]:axis.set_ylim(0,1)
    axes[0,2].set_xscale('symlog');axes[1,2].set_xscale('symlog')
    axes[0,0].legend(fontsize=7)
    fig.suptitle(f'Matched native cubic fields at T={times[0]:.6g}, x={x:.6g}; dotted: compatible analytical reference')
    a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=160);plt.close(fig)
    atomic(a.output.with_suffix('.json'),dict(captured_UTC=datetime.now(timezone.utc).isoformat(),
        command_purpose='read-only matched-time profile comparison; no native solves or history changes',
        step=a.step,time=times[0],station=x,H=heights[0],sources=sources,
        plot_sha256=sha(a.output),script_sha256=sha(Path(__file__)),threads=THREADS))
    print(a.output)


if __name__=='__main__':main()
