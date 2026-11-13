#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H7 results; large state/row packets stay in the isolated Linux root."""
import argparse
import csv
import gc
from pathlib import Path
import numpy as np
from run_wj_h7 import load,save,sha,MATRIX
from study_wj_h5 import identity
from study_wj_h6 import configuration
from study_wj_h7 import reports,accounting,paths
from wj_h7_analysis import Rows


def write_csv(p,rows):
    if not rows:return
    with p.open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def publish(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=root/'publication';out.mkdir();gates={};timeline=[];errors=[];pairs=[];interfaces=[];witnesses=[];norms={};heavy=[];details={};steps=attempts=0
    for h in (5.5,6.,5.):
        view=root/f'scan-view-H{h:g}.json'
        if not view.exists():gates[str(h)]=dict(status='NOT_RUN');continue
        v=load(view);assert sha(v['path'])==v['sha256'];summary=load(v['path']);folder=Path(v['path']).parent
        inv=load(summary['inventory']['path']);rs=reports(root,h);r0=rs[0][1];c=configuration(h);cache=Rows(r0,[p['physical']for p in inv['points']],c)
        rec=load(summary['records']['path']);pair=load(summary['pairs']['path']);details[h]=rec
        native={};endpoints={}
        for branch in ('coarse','fine'):
            rp=root/'cases'/f'H{h:g}mm'/branch/'case.evolution.json';ep=root/'observations'/f'H{h:g}mm'/branch/'result.json'
            if rp.exists():
                r=load(rp);a=accounting(r);native[branch]=dict(a,report=identity(rp),time=r['time'],initialization_work=r['initialization_lift_calls_iterations']);attempts+=a['attempted'];steps+=a['accepted']
                heavy.append(identity(rp));heavy.append(identity(rp.with_name('case.xml')))
                for f in rp.parent.glob('CP*/State.bin'):heavy.append(identity(f))
            if ep.exists():
                e=load(ep);norms[f'{h:g}-{branch}']=dict(status=e['status'],role=e.get('role'),time=e.get('time'),alpha=e.get('alpha'),norms=e.get('spatial_and_completed_norms'),
                    independent_weak=e.get('independent_weak'),quadrature=e.get('quadrature'),producer=e.get('producer'),checkpoint=e.get('checkpoint'))
                endpoints[branch]=e['status'];heavy.append(identity(ep.parent/'native.json'))
                save(out/f'endpoint-H{h:g}-{branch}.json',{k:v for k,v in e.items()if k not in ('scalar_domain',)})
        accepted=[r for r in rec if r['role']=='accepted'and r['branch']!='H6'];expanded=all(r['status']in ('PASS','QUALITY_LOST')for r in accepted)
        quality=all(r['status']=='PASS'for r in accepted)
        reached=all(b in native and native[b]['time']==.04 and native[b]['accepted']==n and not native[b]['failure']for b,n in (('coarse',4),('fine',8)))
        expected={(b,s)for b,n in (('coarse',4),('fine',8))for s in range(n+1)};actual={(r['branch'],r['step'])for r in accepted}
        complete=expected==actual and len([r for r in pair if r['status']=='PASS'])==5
        weak=all(endpoints.get(b)=='PASS'for b in ('coarse','fine'))
        status='OBSTRUCTED'if not expanded or any(a['failure']for a in native.values())else 'NOT_RUN'if not(reached and complete and weak)else 'QUALITY_LOST'if not quality else 'PASS'
        gates[str(h)]=dict(status=status,reached_target=reached,complete_accepted_and_common_times=complete,native=native,expanded_admissibility='PASS'if expanded else 'OBSTRUCTED',
            quality='PASS'if quality else 'QUALITY_LOST',endpoint_weak=endpoints,inventory=summary['inventory'],source_view=v)
        for r in rec:
            timeline.append(dict(height_mm=h,branch=r['branch'],role=r['role'],step=r['step'],time=r['time'],alpha=r['alpha'],status=r['status'],min_in_plane=r['min_in_plane'],min_transverse=r['min_transverse'],min_covariance=r['min_covariance'],
                max_divergence=r['raw_divergence_max'],state_sha256=r['state_sha256'],history_sha256=r['history_sha256'],inventory_sha256=r['inventory_sha256']))
            errors.append(dict(height_mm=h,branch=r['branch'],role=r['role'],step=r['step'],time=r['time'],**r['errors']))
            witnesses.append(dict(height_mm=h,branch=r['branch'],role=r['role'],step=r['step'],time=r['time'],witness=r['witness'],independent_compiled_replay=r['expanded_compiled_replay'],native_replays=r['native_witness_replays']))
            for kind in ('reference_error','initial_drift'):
                for scope,fields in r[kind].items():
                    for field,m in fields.items():
                        scale=1 if field in ('U','G')else c.K if field in ('k','R')else c.A/c.H if field=='epsilon'else c.u_tau**2
                        pairs.append(dict(height_mm=h,branch=r['branch'],step=r['step'],time=r['time'],kind=kind,scope=scope,field=field,L2=m['L2'],sampled_max=m['sampled_max'],fixed_scale=scale,L2_over_fixed_scale=m['L2']/scale,max_over_fixed_scale=m['sampled_max']/scale))
            with np.load(r['arrays']['path'])as data:
                # The original v1 array lacks only the newly published inner
                # turbulent-production contraction; derive it without a solve.
                d=dict(data);d['turbulent_production']=d['nu_eff']*np.sum((d['G']+d['G'].transpose(0,2,1))*d['G'],axis=(1,2))
                inter=cache.interface_history(d);interfaces.append(dict(height_mm=h,branch=r['branch'],role=r['role'],step=r['step'],time=r['time'],**inter))
            heavy.append(r['arrays'])
        for p in pair:
            if p['status']!='PASS':continue
            for scope,fields in p['difference'].items():
                for field,m in fields.items():
                    scale=1 if field in ('U','G')else c.K if field in ('k','R')else c.A/c.H if field=='epsilon'else c.u_tau**2
                    pairs.append(dict(height_mm=h,branch='coarse-fine',step=p['coarse_step'],time=p['time'],kind='timestep_difference',scope=scope,field=field,L2=m['L2'],sampled_max=m['sampled_max'],fixed_scale=scale,L2_over_fixed_scale=m['L2']/scale,max_over_fixed_scale=m['sampled_max']/scale))
        save(out/f'pairs-H{h:g}.json',pair);del cache;gc.collect()
    overall='PASS'if len(gates)==3 and all(g['status']=='PASS'for g in gates.values())else 'QUALITY_LOST'if all(g['status']in ('PASS','QUALITY_LOST')for g in gates.values())else 'OBSTRUCTED'
    save(out/'component-gates.json',dict(status=overall,cases=gates,accepted_steps=steps,attempts=attempts,fresh_cp0=load(root/'accounting.json')['cp0_reservations'],
        scope='finite inventory/T=.04 idealized neutral transport; no order, long-time, developed-flow or experimental claim'))
    write_csv(out/'timeline.csv',timeline);write_csv(out/'reference-max-errors.csv',errors);write_csv(out/'norms-and-differences.csv',pairs)
    save(out/'principal-witnesses.json',witnesses);save(out/'interface-budget-history.json',interfaces);save(out/'residual-weak-checks.json',norms)
    for name in ('imports.json','runtime-equivalence.json','private-runtime.json','runtime-differences-v1.json','setup.json','allocation-denial-v1.json','parent-thread-fix.json'):
        if (root/name).exists():save(out/name,load(root/name))
    save(out/'accounting-snapshot.json',load(root/'accounting.json'))
    commands={}
    for f in sorted((root/'operations').glob('*/receipt.json')):
        commands[f.parent.name]=load(f);log=f.with_name('command.log');(out/'logs').mkdir(exist_ok=True)
        (out/'logs'/f'{f.parent.name}.log').write_bytes(log.read_bytes())
    save(out/'commands.json',commands)
    fig,ax=plt.subplots(1,3,figsize=(11,3.4),constrained_layout=True)
    for a,h in zip(ax,(5.,5.5,6.)):
        for b,mark in (('coarse','o'),('fine','s')):
            rr=[v for v in timeline if v['height_mm']==h and v['branch']==b and v['role']=='accepted']
            a.plot([v['time']for v in rr],[v['min_in_plane']for v in rr],mark+'-',ms=3,label=b)
        a.axhline(.004137527846908,color='.4',ls=':',label='exact neutral');a.axhline(.004137527846908/2,color='r',ls='--',label='quality threshold')
        a.set(title=f'Hw={h:g}mm',xlabel='time [Hh/Uinf]',ylabel='min m/nu_eff');a.grid(alpha=.2)
    ax[0].legend(fontsize=7);fig.savefig(out/'margin-history.png',dpi=160);fig.savefig(out/'margin-history.pdf');plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(11,3.4),constrained_layout=True)
    for a,field in zip(ax,('U','k','epsilon')):
        for h in (5.,5.5,6.):
            for kind,ls in (('reference_error','-'),('initial_drift','--'),('timestep_difference',':')):
                pp=[p for p in pairs if p['height_mm']==h and p['field']==field and p['scope']=='full'and p['kind']==kind and p['branch']==('coarse-fine'if kind=='timestep_difference'else 'fine')]
                a.plot([p['time']for p in pp],[p['L2_over_fixed_scale']for p in pp],ls,label=f'{h:g}mm {kind}')
        a.set(title=field,xlabel='time [Hh/Uinf]',ylabel='physical L2 / fixed field scale');a.grid(alpha=.2)
    ax[0].legend(fontsize=6);fig.savefig(out/'error-drift-timestep.png',dpi=170);fig.savefig(out/'error-drift-timestep.pdf');plt.close(fig)
    # Exact physical profiles are overlaid with both H-sided/face owners.
    fig,ax=plt.subplots(1,3,figsize=(11,3.5),constrained_layout=True)
    for a,h in zip(ax,(5.,5.5,6.)):
        if h not in details:continue
        for branch in ('coarse','fine'):
            rr=[r for r in details[h]if r['branch']==branch and r['role']=='accepted'];r=rr[-1]
            view=load(root/f'scan-view-H{h:g}.json');summary=load(view['path']);inv=load(summary['inventory']['path']);base=reports(root,h)[0][1];c=configuration(h);cache=Rows(base,[p['physical']for p in inv['points']],c)
            with np.load(r['arrays']['path'])as d:
                ii=[i for i,(p,o)in enumerate(cache.rows[:cache.point_count])if p[0]==8 and o['region']=='outer'and p[2]<1]
                ii=sorted(ii,key=lambda i:cache.rows[i][0][2]);z=[cache.rows[i][0][2]*40 for i in ii]
                a.plot(z,d['G'][ii,0,2],'.-',ms=2,label=branch)
                if branch=='fine':a.plot(z,cache.exact['G'][ii,0,2],'k:',label='exact')
            del cache;gc.collect()
        a.set(title=f'Hw={h:g}mm, x=8, terminal',xlabel='physical z [mm]',ylabel='Gxz [Uinf/Hh]');a.grid(alpha=.2)
    ax[0].legend();fig.savefig(out/'terminal-gradients.png',dpi=160);plt.close(fig)
    save(out/'manifest.json',dict(root=str(root),publisher=identity(__file__),heavy=heavy,producer_scope='governing source unchanged c344dc2; native observer e62b5f4; Python versions in receipts',
        note='native full scalar/state/history rows remain heavy; reference/CP0 and actual completed callback roles distinct'))
    save(out/'seal.json',dict(files=[dict(name=str(p.relative_to(out)),**identity(p))for p in sorted(out.rglob('*'))if p.is_file()]))
    print('H7 compact publication',overall,'attempts',attempts,'accepted',steps,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);publish(p.parse_args().root)
