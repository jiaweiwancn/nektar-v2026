#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Physical common-rule comparisons for accepted-AWF native case suites.

Initial representation errors, evolution from actual native h0, cross-height
reference changes and one/graded differences remain separate quantities.
"""
import argparse
import datetime
import gzip
import json
import math
import os
from pathlib import Path

from run_cls_accepted_awf import CASES,KN,X,FAMILY,SPEC_SHA256
from run_cls_epawf_picard import load,write_json,require,sha,THREADS
from analyse_cls_epawf_height_comparison import differences,evolution_difference,norm,attempts_summary


def report(path):
    if path.suffix=='.gz':
        with gzip.open(path,'rt') as f:return json.load(f,parse_int=float)
    return load(path)


def valid_common_query(q):
    return bool(q and q.get('common_samples')) and q.get('valid') is not False and \
        q.get('common_reporting_queries_valid') is not False


def comparison_states(result):
    """Invalid optional queries cannot enter a cross-suite field comparison."""
    initial=result.get('initialization')
    if not valid_common_query(initial):return {}
    states={0:initial}
    if valid_common_query(result.get('final')):states[int(result['steps'])]=result['final']
    for n,entry in result.get('snapshots',{}).items():
        if entry['exit_code']==0:
            q=load(Path(entry['path']))
            if valid_common_query(q):states[int(n)]=q
    return states


def compact(q):
    vectors=('same_state_accepted_source_load','same_state_calibrated_EP_source_load')
    d={k:v for k,v in q.items() if k not in ('common_samples','weak_k_component_loads',*vectors)}
    if 'weak_k_component_loads' in q:
        d['weak_load_vector_receipt']={name:dict(coefficients=len(v),maximum_abs=max(map(abs,v),default=0.),
            nonzero_modes=sum(x!=0 for x in v),interpretation='coefficient load diagnostic; not a physical norm')
            for name,v in q['weak_k_component_loads'].items()}
    d['same_state_source_load_receipts']={name:dict(coefficients=len(q[name]),maximum_abs=max(map(abs,q[name]),default=0.),
        nonzero_modes=sum(x!=0 for x in q[name]),interpretation='native all-mode signed source loads; full vectors retained in the native artifact')
        for name in vectors if name in q}
    return d


def reference_errors(q):
    r=differences(q,q,True);s=q['common_samples'];r['v']=norm(s['weights'],s['v'])
    r['velocity']=math.hypot(r['u'],r['v']);return r


def normalized(change,K0=.51050801305988613):
    scales=dict(u=.4*math.sqrt(24),v=.4*math.sqrt(24),k=K0*math.sqrt(24),
        epsilon=(.4**3/.4405539051866587/2)*math.sqrt(12),pressure=.16*math.sqrt(24))
    return {k:change[k]/v for k,v in scales.items() if k in change}


def band_key(row):return row['name']+':'+row.get('x_region','all')


def estimated_integrals(by_step,through):
    """Trapezoidal sampled continuum rates, separate from exact step repairs."""
    estimates={};old=by_step[0]
    for n,current in sorted(by_step.items()):
        if not n or n>through:continue
        dt=current['time']-old['time'];before={band_key(x):x for x in old['fixed_bands']['regions']}
        for y in current['fixed_bands']['regions']:
            name=band_key(y);x=before[name]
            budget=estimates.setdefault(name,dict(start=0.,end=0.,production=0.,destruction=0.,valid=True))
            for key in ('production','destruction'):
                if x[key] is None or y[key] is None:budget['valid']=False
                else:budget[key]+=.5*dt*(x[key]+y[key])
            budget['end']=current['time']
        old=current
    return estimates


def work_summary(attempts,time):
    out=attempts_summary(attempts,time)
    chosen={a['time']:a for a in attempts if not a['failure'] and a['scheduled_complete'] and a['time']<=time+1e-13}
    for row in out['increments']:
        b=chosen[row['time']].get('k_bound',{})
        row.update(accepted_step_mass_norm=b.get('mass_norm_accepted_step'),
            raw_bound_inner_product=b.get('mass_inner_product_raw_bound'),
            increment_norm_identity_defect=b.get('mass_norm_increment_identity_defect'))
    return out


def reference_offsets(a,b):
    x,y=a['common_samples'],b['common_samples'];require(x['weights']==y['weights'],'common physical rules differ')
    return {field:norm(x['epsilon_weights'] if field=='epsilon' else x['weights'],
        [v-w for v,w in zip(x['reference_'+field],y['reference_'+field])]) for field in ('u','k','epsilon')}


def dense_column(q,target=1.6):
    """Recover native cubic values/derivatives from common physical GL samples."""
    import numpy as np
    gx,_=np.polynomial.legendre.leggauss(8);gz,_=np.polynomial.legendre.leggauss(16)
    xs=[a+(v+1)*(b-a)/2 for a,b in zip(X,X[1:]) for v in gx];column=min(range(len(xs)),key=lambda i:abs(xs[i]-target))
    knots=q.get('common_reporting_knots',KN)
    size=16*(len(knots)-1);s=q['common_samples'];require(len(s['u'])==len(xs)*size,'native common column layout changed')
    u=np.array(s['u']).reshape((len(xs),size))[column];k=np.array(s['k']).reshape((len(xs),size))[column]
    values=dict(z=[],u=[],du_dz=[],k=[],dk_dz=[])
    for panel,(a,b) in enumerate(zip(knots,knots[1:])):
        if b>q['H']+1e-12:continue
        sample=a+(gz+1)*(b-a)/2;z=np.linspace(a,b,40)
        up=np.polynomial.Legendre.fit(sample,u[16*panel:16*(panel+1)],3,domain=[a,b])
        kp=np.polynomial.Legendre.fit(sample,k[16*panel:16*(panel+1)],3,domain=[a,b])
        for name,v in dict(z=z,u=up(z),du_dz=up.deriv()(z),k=kp(z),dk_dz=kp.deriv()(z)).items():values[name].extend(v)
    return xs[column],values


def analytical_column(q,z):
    """Independent rough-reference integral, for plot lines only.

    Numerical validation of this construction belongs to the kernel oracle;
    this routine never changes a native state or supplies a transport gradient.
    """
    from scipy.integrate import quad
    ref=q['reference'];K=ref['K0'];us=.4;Q=ref.get('q_star',K/us**2);kap=ref.get('kappa',.4405539051866587)
    H=q['H'];nu=1.5e-5;h=q['h'];r=h*math.sqrt(K)/(70*nu)
    b=10.7*(1-r**max(.5-.4*r**.7,1-.79*r**(-.28)))
    zv=nu*b/math.sqrt(K);a=kap*math.sqrt(K/Q);d=nu-a*zv
    require(zv<=0 and d>=nu,'plot reference routine is restricted to the tested rough branch')
    t0=K/Q;A=t0**1.5/kap;c=ref['cutoff']
    def flux(s):
        Rc=s/c if s<=c else 1+math.log(s/c)
        return A*(Rc-math.log1p(a*s/d)+nu*(1/d-1/(d+a*s)))
    def kval(s):
        split=[v for v in (s,c,H) if s<=v<=H];split=sorted(set(split))
        return K-sum(quad(lambda r:flux(r)/(d+a*r),x,y,epsabs=2e-12,epsrel=2e-12)[0] for x,y in zip(split,split[1:]))
    return dict(u=[t0/a*math.log1p(a*s/d) for s in z],du_dz=[t0/(d+a*s) for s in z],
        k=[kval(s) for s in z],dk_dz=[flux(s)/(d+a*s) for s in z])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--suite',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--compare-suite',type=Path)
    a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=True);states={};native={}
    out=dict(schema=1,captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),family=FAMILY,
        specification_sha256=SPEC_SHA256,suite=str(a.suite.resolve()),cases={},one_graded={},cross_height={},
        normalization='dimensional physical L2 and common U=.4,K=K0,H=2,V=24,outerV=12; epsilon z>=2; one pressure gauge',
        budget_units='repair budgets are integrals of k; sampled production/destruction are rates integrated over physical time',
        caveat='finite-inventory/ground admissibility and survival do not establish global positivity, temporal accuracy or steady convergence')
    for case in CASES:
        path=a.suite/case/'result.json'
        if not path.exists():continue
        result=load(path);initial=result.get('initialization');row={k:v for k,v in result.items() if k not in ('initialization','final','snapshots')}
        out['cases'][case]=row
        if not initial or not result.get('actual_h0_sha256'):continue
        if not valid_common_query(initial):
            row['invalid_optional_initial_query']=initial
            continue
        states[case]={0:initial};row['initial']=compact(initial);row['actual_h0_error_to_own_reference']=reference_errors(initial)
        nf=a.suite/case/'trajectory/case.evolution.json'
        if nf.exists():native[case]=load(nf)
        for step,entry in result.get('snapshots',{}).items():
            file=Path(entry['path'])
            if entry['exit_code']==0 and file.exists():
                d=load(file)
                if valid_common_query(d):states[case][int(step)]=d
                else:row.setdefault('invalid_optional_queries',[]).append(dict(entry,payload=d))
            else:row.setdefault('invalid_optional_queries',[]).append(entry)
        if result.get('final'):
            if valid_common_query(result['final']):
                states[case][int(result['steps'])]=result['final']
            else:row.setdefault('invalid_optional_queries',[]).append(dict(time=result.get('time'),payload=result['final']))
        row['observations']={};attempts=native.get(case,{}).get('pimple_attempts',[])
        s0=initial['common_samples'];bulk=sum(w*u for w,u in zip(s0['weights'],s0['u']))/sum(s0['weights']);flow=6/bulk
        row['initial_flow_through_time']=flow;row['wall_time_scale']=initial['H']/.4
        for step,q in sorted(states[case].items()):
            if not step:continue
            change=differences(q,initial);budget=work_summary(attempts,q['time'])
            row['observations'][str(step)]=dict(time=q['time'],diagnostic=compact(q),drift_from_actual_h0=change,
                common_normalized_drift=normalized(change),error_to_own_reference=reference_errors(q),
                exact_discrete_work_and_correction=budget,sampled_physical_time_integrals=estimated_integrals(states[case],step),
                elapsed_wall_times=q['time']/row['wall_time_scale'],elapsed_flow_through_times=q['time']/flow)
    for left,right in [('Accepted-AWF-2G','Accepted-AWF-2O'),('Accepted-AWF-05G','Accepted-AWF-05O')]:
        if left not in states or right not in states:continue
        out['one_graded'][left+' vs '+right]={str(n):dict(time=states[left][n]['time'],
            represented_difference=differences(states[left][n],states[right][n]),
            evolution_difference=evolution_difference(states[left][n],states[left][0],states[right][n],states[right][0]))
            for n in sorted(states[left].keys()&states[right].keys())}
    for left,right in [('Accepted-AWF-2G','Accepted-AWF-05G'),('Accepted-AWF-2O','Accepted-AWF-05O')]:
        if left not in states or right not in states:continue
        out['cross_height'][left+' vs '+right]=dict(analytical_reference_offset=reference_offsets(states[left][0],states[right][0]),
            matched_times={str(n):dict(time=states[left][n]['time'],represented_difference=differences(states[left][n],states[right][n]),
                evolution_difference=evolution_difference(states[left][n],states[left][0],states[right][n],states[right][0]))
                for n in sorted(states[left].keys()&states[right].keys())})
    if a.compare_suite:
        out['cross_suite']={}
        for case,by_step in states.items():
            f=a.compare_suite/case/'result.json'
            if not f.exists():continue
            result=load(f)
            other=comparison_states(result)
            if not other:continue
            out['cross_suite'][case]=dict(other=str(a.compare_suite.resolve()),
                same_actual_h0=result['actual_h0_sha256']==out['cases'][case]['actual_h0_sha256'],
                common_times={str(q['time']):differences(q,next(v for v in other.values() if abs(v['time']-q['time'])<1e-12))
                    for q in by_step.values() if any(abs(v['time']-q['time'])<1e-12 for v in other.values())})
    write_json(a.output/'native.json',out)
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for pair in [('Accepted-AWF-2G','Accepted-AWF-2O'),('Accepted-AWF-05G','Accepted-AWF-05O')]:
            included=[name for name in pair if name in states]
            if not included:continue
            for station in (-3.5,1.6):
                fig,axes=plt.subplots(2,4,figsize=(15,8),layout='constrained')
                for name in included:
                    companion=out['cases'][name]['definition']['companion']
                    label=name+(' ['+companion+']' if companion!='faithful' else '')
                    for n in sorted({0,max(states[name])}):
                        q=states[name][n];x,v=dense_column(q,station)
                        for axis,field in zip(axes[0],('k','u','du_dz','dk_dz')):axis.plot(v[field],v['z'],label=f'{label} t={q["time"]:.4g}')
                        for axis,field in zip(axes[1],('k','u','du_dz','dk_dz')):axis.plot(v[field],[z/q['H'] for z in v['z']],label=f'{label} t={q["time"]:.4g}')
                    q=states[name][0];x,v=dense_column(q,station);ref=analytical_column(q,v['z'])
                    for line,scale in zip(axes,[1.,q['H']]):
                        for axis,field in zip(line,('k','u','du_dz','dk_dz')):axis.plot(ref[field],[z/scale for z in v['z']],':',color='black')
                for line,label in zip(axes,('physical z','z/H')):
                    for axis,field in zip(line,('transported k','velocity u','du/dz','dk/dz')):axis.set(xlabel=field,ylabel=label)
                for axis in axes[0]:axis.set_ylim(0,.1)
                axes[0,2].set_xscale('symlog');axes[1,2].set_xscale('symlog');axes[0,0].legend(fontsize=7)
                fig.suptitle(f'Actual native cubic fields, x={x:.6g}; dotted: accepted analytical reference')
                fig.savefig(a.output/(pair[0]+'-profiles-x'+str(station)+'.png'),dpi=160);plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
        for name,by_step in states.items():
            companion=out['cases'][name]['definition']['companion'];label=name+(' ['+companion+']' if companion!='faithful' else '')
            points=sorted(by_step.values(),key=lambda q:q['time']);times=[q['time'] for q in points]
            for axis,band in zip(axes[:2],('z_le_0.01','z_le_0.1')):
                axis.plot(times,[next(v['k_integral'] for v in q['fixed_bands']['regions'] if v['name']==band and v.get('x_region','all')=='all') for q in points],'o-',label=label)
                axis.set(xlabel='physical time',ylabel='integral k, '+band)
            axes[2].plot(times,[q['mean_wall_stress'] for q in points],'o-',label=label)
        axes[2].set(xlabel='physical time',ylabel='mean wall stress');axes[0].legend(fontsize=7)
        fig.savefig(a.output/'fixed-bands-wall-stress.png',dpi=160);plt.close(fig)
    except ImportError:out['plots']='Matplotlib/SciPy unavailable; numerical evidence remains in native.json';write_json(a.output/'native.json',out)
    print(a.output/'native.json')


if __name__=='__main__':main()
