#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H8 measured results; full owned fields/checkpoints remain outside Git."""
import argparse
import csv
import math
from pathlib import Path
import numpy as np
from run_wj_h8 import load,save,sha,coverage_quality
from study_wj_h5 import identity
from study_wj_h6 import configuration
from study_wj_h7 import accounting
from wj_h7_analysis import Rows,EXACT_MARGIN
from wj_h8_analysis import FACE,ratio
from wj_h4_points import aligned


def csvfile(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('x',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def comparison(root,out,folder):
    summary=load(folder/'summary.json');records=load(summary['records']['path']);inv=load(summary['inventory']['path']);c=configuration(6.)
    r=load(root/'imports-v2/H6mm/initial.evolution.json');cache=Rows(r,[p['physical']for p in inv['points']],c)
    index={(v['branch'],v['step']):v for v in records if v['role']=='accepted'};metrics=[];local=[];envelopes=[];pairs=[]
    point=np.arange(len(cache.rows))<cache.point_count;outer=point & cache.outer
    fixed=[i for i,(p,o)in enumerate(cache.rows[:cache.point_count])if abs(p[0]-8)<1e-14 and abs(p[2]-FACE)<1e-14 and o['region']=='outer']
    assert len(fixed)==2,'both normal owners at the fixed face required'
    variables={'m':('margin',None,1.),'transverse':('transverse',None,1.),'Gxx':('G',(0,0),1.),'Gxz':('G',(0,2),1.),'Gzx':('G',(2,0),1.),'Gzz':('G',(2,2),1.),
        'U':('U',(0,),1.),'W':('U',(2,),1.),'k':('k',None,c.K),'epsilon':('epsilon',None,c.A/c.H),'pi':('pi',None,c.u_tau**2),'p_mechanical':('mechanical_pressure',None,c.u_tau**2),
        'divergence':('divergence',None,1.),'raw_production':('production',None,c.A/c.H),'Rxx':('R',(0,0),c.K),'Rxz':('R',(0,2),c.K),'Rzz':('R',(2,2),c.K)}
    for s in range(5):
        keys=[('coarse',s),('fine',2*s),('quarter',4*s)]
        if not all(k in index for k in keys):pairs.append(dict(time=s*.01,status='NOT_RUN',missing=[list(k)for k in keys if k not in index]));continue
        rr=[index[k]for k in keys]
        assert len({aligned(v['step'],dict(coarse=.01,fine=.005,quarter=.0025)[v['branch']],v['time'])for v in rr})==1
        with np.load(rr[0]['arrays']['path'])as a,np.load(rr[1]['arrays']['path'])as b,np.load(rr[2]['arrays']['path'])as d:
            D1=cache.norms(a,b);D2=cache.norms(b,d)
            for scope,fields in D1.items():
                mask=np.ones(len(cache.rows),bool)if scope=='full'else cache.outer if scope=='outer'else ~cache.outer
                for field,norms in fields.items():
                    scale=1. if field in ('U','G')else c.K if field in ('k','R')else c.A/c.H if field=='epsilon'else c.u_tau**2
                    active=mask & (cache.outer if field=='epsilon'else True);measure=float(sum(cache.weights[active]))
                    for norm in norms:
                        floor=4e-11*scale*(math.sqrt(measure)if norm=='L2'else 1.)
                        metrics.append(dict(time=s*.01,scope=scope,field=field,norm=norm,fixed_scale=scale,measure=measure,**ratio(D1[scope][field][norm],D2[scope][field][norm],floor)))
            for name,(field,comp,scale)in variables.items():
                aa=a[field]if comp is None else a[field][(slice(None),)+comp];bb=b[field]if comp is None else b[field][(slice(None),)+comp];dd=d[field]if comp is None else d[field][(slice(None),)+comp]
                x=bb-aa;y=dd-bb;ii=np.flatnonzero(outer);i1=ii[np.argmax(abs(x[ii]))];i2=ii[np.argmax(abs(y[ii]))]
                envelopes.append(dict(time=s*.01,variable=name,D1_max_owner=cache.rows[i1],D2_max_owner=cache.rows[i2],**ratio(float(abs(x[i1])),float(abs(y[i2])),4e-11*scale)))
                for i in fixed:
                    p,o=cache.rows[i];local.append(dict(time=s*.01,geometry=o['geometry'],normal_side=o['normal_side'],x=p[0],physical_z=p[2],variable=name,
                        coarse=float(aa[i]),fine=float(bb[i]),quarter=float(dd[i]),signed_D1=float(x[i]),signed_D2=float(y[i]),**ratio(float(abs(x[i])),float(abs(y[i])),4e-11*scale)))
            for i in fixed:
                # Direction is pi-periodic; a change is a diagnostic, not a state correction.
                angles=[float(v['angle'][i])for v in (a,b,d)];changes=[abs((angles[j+1]-angles[j]+math.pi/2)%math.pi-math.pi/2)for j in (0,1)]
                local.append(dict(time=s*.01,geometry=cache.rows[i][1]['geometry'],normal_side=cache.rows[i][1]['normal_side'],variable='minimizing_angle',coarse=angles[0],fine=angles[1],quarter=angles[2],D1=changes[0],D2=changes[1],status='DIRECTION_DIAGNOSTIC'))
        pairs.append(dict(time=s*.01,status='PASS',recorded_times=[v['time']for v in rr],states=[v['state_sha256']for v in rr],producers=[v['producer']for v in rr]))
    csvfile(out/'three-level-norms.csv',metrics);csvfile(out/'fixed-face-three-level.csv',local);save(out/'same-owner-envelopes.json',envelopes);save(out/'common-time-bindings.json',pairs)
    return dict(common_times=pairs,primary_face=[v for v in local if v['variable']=='m'],margin_envelopes=[v for v in envelopes if v['variable']=='m'])


def publish(root,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out.mkdir();timeline=[];error=[];interfaces=[];projection=[];witnesses=[];heavy=[];gates={};localfixed=[];views={};summaries={}
    for h in (5.,5.5,6.):
        binding=root/f'diagnosis-view-H{h:g}.json'
        if not binding.exists():gates[str(h)]={'status':'NOT_RUN'};continue
        v=load(binding);assert sha(v['path'])==v['sha256'];summary=load(v['path']);summaries[h]=summary;records=load(summary['records']['path']);folder=Path(v['path']).parent
        gates[str(h)]=dict(summary['status'],inventory=summary['inventory'],states=len(records),diagnostic_mass_solves=summary['diagnostic_mass_solves'],diagnostic_factorizations=summary['factorizations'])
        for r in records:
            timeline.append(dict(height_mm=h,branch=r['branch'],step=r['step'],time=r['time'],alpha=r['alpha'],role=r['role'],expanded_status=r['status'],min_in_plane=r['min_in_plane'],min_transverse=r['min_transverse'],min_covariance=r['min_covariance'],
                max_divergence=r['raw_divergence_max'],state_sha256=r['state_sha256'],callback_history_sha256=r['callback_history_sha256'],compiled_minimum=r['compiled_minimum'],inventory_sha256=r['inventory_sha256']))
            witnesses.append(dict(height_mm=h,branch=r['branch'],step=r['step'],time=r['time'],role=r['role'],witness=r['witness'],compiled_minimum=r['compiled_minimum'],producer=r['producer'],checkpoint=r['checkpoint']))
            error.append(dict(height_mm=h,branch=r['branch'],step=r['step'],time=r['time'],role=r['role'],**r['errors']))
            for region,p in r['pressure_projection']['regions'].items():
                projection.append(dict(height_mm=h,branch=r['branch'],step=r['step'],time=r['time'],region=region,**{k:p[k]for k in ('divergence_L2','continuous_Q1_projection_L2','continuous_Q1_complement_L2','broken_Q1_projection_L2','broken_Q1_complement_L2','maximum_pressure_test_moment','integrated_net_flux','identity_roundoff','broken_identity_roundoff','gauge_contribution')}))
            interfaces.append(dict(height_mm=h,branch=r['branch'],step=r['step'],time=r['time'],role=r['role'],**r['interface_budget']))
            heavy.append(r['arrays'])
        heavy.extend([summary['records'],summary['local'],summary['inventory']])
        local=load(summary['local']['path'])
        for r in local:
            if abs(r['physical'][0]-8)<1e-14 and abs(r['physical'][2]-FACE)<1e-14 and r['owner']['region']=='outer':
                localfixed.append(dict(height_mm=h,branch=r['branch'],step=r['step'],time=r['time'],role=r['role'],geometry=r['owner']['geometry'],normal_side=r['owner']['normal_side'],x=r['physical'][0],z=r['physical'][2],
                    Gxx=r['G'][0][0],Gxz=r['G'][0][2],Gzx=r['G'][2][0],Gzz=r['G'][2][2],U=r['U'][0],W=r['U'][2],k=r['k'],epsilon=r['epsilon'],pi=r['pi'],mechanical_pressure=r['mechanical_pressure'],
                    margin=r['margin'],transverse=r['transverse'],angle=r['angle'],divergence=r['divergence'],production=r['production'],root=r['root'],F_N=r['F_N'],scaled_root_residual=r['scaled_root_residual']))
    comparison_result=None
    if 6. in summaries:
        folder=Path(summaries[6.]['records']['path']).parent;comparison_result=comparison(root,out,folder)
    for v in (root/'views').glob('*/result.json')if(root/'views').exists()else []:
        value=load(v);views[v.parent.name]=value;heavy.extend([identity(v),value['native']])
    native={}
    rp=root/'case/run/case.evolution.json'
    if rp.exists():
        r=load(rp);native=dict(accounting(r),report=identity(rp),time=r['time'],initialization_work=r['initialization_lift_calls_iterations']);heavy.extend([identity(rp),identity(rp.with_name('case.xml'))]);heavy.extend(identity(f)for f in rp.parent.glob('CP*/State.bin'))
    expected_views={'coarse-1','coarse-3','fine-2','fine-6','quarter-4','quarter-12','quarter-16'}
    weak='PASS'if expected_views==set(views)and all(v['status']=='PASS'for v in views.values())else 'UNRESOLVED'if views and any(v['status']!='PASS'for v in views.values())else 'NOT_RUN'
    book=load(root/'accounting.json');overall='OBSTRUCTED'if book.get('coverage_halt')or native.get('failure')else 'NOT_RUN'if native.get('accepted')!=16 else 'UNRESOLVED'if weak!='PASS'else 'QUALITY_LOST'if any(g.get('quality')!='PASS'for g in gates.values())else 'PASS'
    save(out/'component-gates.json',dict(status=overall,old_reuse='H7 closed36/36 unchanged',new_native=native,heights=gates,strict_weak=weak,strict_views=sorted(views),required_strict_views=sorted(expected_views),
        coverage_halt=book.get('coverage_halt'),coverage_earliest=book.get('coverage_earliest'),target=.04,new_cp0=book['cp0_reservations'],strict_view_reservations=book.get('strict_view_reservations',0),
        predictor_pressure_causality='UNRESOLVED:no genuine separated predictor/correction vectors; local action does not take pressure',contraction=comparison_result))
    csvfile(out/'timeline.csv',timeline);csvfile(out/'reference-errors.csv',error);csvfile(out/'pressure-complement.csv',projection);csvfile(out/'fixed-face-history.csv',localfixed)
    save(out/'principal-witnesses.json',witnesses);save(out/'interface-budget-history.json',interfaces);save(out/'strict-weak-stage-views.json',views)
    for name in ('setup.json','imports-v2.json','runtime-equivalence-v2.json','stage-receipt-v2.json','admission.json','case/intent.json','case/initial-pair.json'):
        if(root/name).exists():save(out/Path(name).name,load(root/name));heavy.append(identity(root/name))
    # All operation receipts, including failed imports, survive. Native maps are
    # compared against the verified private closure rather than host libc.
    commands={};maps={}
    expected={v['sha256']for v in load(root/'runtime-equivalence-v2.json')['qualified_ELFs']}
    for f in sorted((root/'operations').glob('*/receipt.json')):
        commands[f.parent.name]=load(f);heavy.extend([identity(f),identity(f.with_name('command.log')),identity(f.with_name('loaded-maps.json'))])
        loaded=load(f.with_name('loaded-maps.json'))
        if commands[f.parent.name].get('cp0')or commands[f.parent.name].get('units')or commands[f.parent.name].get('strict_view'):
            bypath={v['path']:v for v in loaded};process=load(f.with_name('process-maps.json'));native=[]
            for pid,v in process.items():
                if any(Path(p).name in ('IncNavierStokesSolver.elf','RANSWallResolutionProbe.elf')for p in v['files']):
                    native.extend(bypath[p]for p in v['files']if p in bypath)
            native=list({v['path']:v for v in native}.values());unexpected=[v for v in native if v['sha256']not in expected]
            maps[f.parent.name]=dict(status='PASS'if native and not unexpected else 'UNRESOLVED',native_loaded=native,unexpected=unexpected,process_maps=identity(f.with_name('process-maps.json')),scope='actual native PID; host time/shell harness excluded')
    save(out/'commands.json',commands);save(out/'actual-native-runtime.json',maps)
    save(out/'accounting-snapshot.json',book)
    fig,axs=plt.subplots(1,3,figsize=(11,3.4),constrained_layout=True)
    for h,ax in zip((5.,5.5,6.),axs):
        for branch in ('coarse','fine','quarter'):
            rr=[r for r in timeline if r['height_mm']==h and r['branch']==branch and r['role']=='accepted']
            if rr:ax.plot([v['time']for v in rr],[v['min_in_plane']for v in rr],'.-',label=branch)
        ax.axhline(EXACT_MARGIN/2,color='k',ls='--',lw=.8);ax.set_title(f'Hw={h:g} mm');ax.set_xlabel('time Hh/Uinf');ax.grid(alpha=.2)
    axs[0].set_ylabel('minimum in-plane principal / nu_eff');axs[-1].legend();fig.savefig(out/'margins.png',dpi=160);plt.close(fig)
    if comparison_result:
        rr=[v for v in comparison_result['primary_face']if v['geometry']==17 and v['time']>0];fig,axs=plt.subplots(1,2,figsize=(9,3.4),constrained_layout=True)
        axs[0].plot([v['time']for v in rr],[v['D1']for v in rr],'o-',label='dt .01/.005');axs[0].plot([v['time']for v in rr],[v['D2']for v in rr],'s-',label='dt .005/.0025');axs[0].set_ylabel('same-owner principal difference');axs[0].legend();axs[1].plot([v['time']for v in rr],[v['rho']for v in rr],'o-');axs[1].axhline(1,color='k',ls='--');axs[1].set_ylabel('D2/D1')
        for ax in axs:ax.set_xlabel('time Hh/Uinf');ax.grid(alpha=.2)
        fig.suptitle('x=8, physical z=10.243007 mm, geometry17 upper trace');fig.savefig(out/'local-timestep.png',dpi=160);plt.close(fig)
    # Runtime heavy records are explicitly remote verification, not downloaded fields.
    save(out/'manifest.json',dict(heavy=list({v['path']:v for v in heavy}.values()),scope='hashes verified on Linux; compact transport separately verified on Mac',source=identity(__file__)))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()and p.name!='seal.json'],status=overall))
    print('H8 publication',overall,'accepted',native.get('accepted',0),'strict',len(views),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();publish(a.root,a.output)
