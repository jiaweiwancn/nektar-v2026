# SPDX-License-Identifier: MIT
"""Common physical composite observations and separate engineering pilot gates."""
import math
import numpy as np
from numpy.testing import assert_allclose
from wj2_cases import normal_knots
from wj3_measure import FIELDS,K0,arrays,wall,arithmetic
from wj_qualification import component,aggregate
from test_wj2_native import independent,CSTAR,KAPPA

AMP=dict(u=.5,v=.5,velocity=.5,k=K0,outer_epsilon=.125/KAPPA,reduced_pressure=.25,mechanical_pressure=.25)
BUDGET=dict(L2=1e-4,sampled_Linf=1e-3,traction=1e-3,k_flux=1e-4,pressure_offset=1e-4,divergence=1e-3)

def rule(xmeshes=((0.,1.),(0.,.5,1.)),order=6,xorder=4):
    _,zv,_,_,zc=wall(K0)
    xs=sorted(set([x for knots in xmeshes for x in knots]+[.25,.75]))
    zs=sorted(set(normal_knots(2)+[z for z in (zc,zv,.1,.5,1) if 0<z<2]))
    gx,wx=np.polynomial.legendre.leggauss(xorder);gz,wz=np.polynomial.legendre.leggauss(order)
    rows=[]
    for xl,xh in zip(xs,xs[1:]):
        for zl,zh in zip(zs,zs[1:]):
            for x,w in zip(xl+(gx+1)*(xh-xl)/2,wx*(xh-xl)/2):
                for z,v in zip(zl+(gz+1)*(zh-zl)/2,wz*(zh-zl)/2):rows.append([float(x),0.,float(z),float(w*v),0,0])
    # Complete face intersections plus common dense horizontal maximum samples.
    # Explicit x-side is distinct from the radial side at H.
    samplex=sorted(set(xs+np.linspace(xs[0],xs[-1],33).tolist()))
    for x in samplex:
        for z in zs:
            for side in ((-1,1) if z==1 else (0,)):
                for xside in ((-1,1) if x in xs[1:-1] else (0,)):
                    rows.append([float(x),0.,float(z),0.,side,xside])
    return rows

def owned_rule(physical,level,xknots):
    zknots=normal_knots(level);nx=len(xknots)-1;result=[]
    for x,y,z,w,side,xside in physical:
        def choose(knots,value,selection):
            candidates=[i for i,(a,b) in enumerate(zip(knots,knots[1:])) if a-1e-14<=value<=b+1e-14]
            if not candidates:raise ValueError('observation outside physical domain')
            # Negative selection approaches the face from smaller coordinates.
            return min(candidates) if selection<0 else max(candidates)
        ix=choose(xknots,x,xside);iz=choose(zknots,z,side)
        result.append([x,y,z,w,side,iz*nx+ix,xside])
    return result

def same_rule(a,b):
    key=lambda r:(r['position'],r['weight'],r['side'],r.get('requested_x_side',0))
    if [key(r) for r in a['observations']]!=[key(r) for r in b['observations']]:raise ValueError('unequal common physical observation rule')

def physical_norms(d,w,outer,rows,mask=None,gauge=None):
    d=np.array(d,copy=True);w=np.array(w);mask=np.ones(len(w),bool) if mask is None else np.array(mask)
    # Gauge is determined over the FULL physical domain before any extraction.
    shift=float(w@d[:,4]/sum(w)) if gauge is None else gauge
    d[:,4]-=shift;d[:,5]-=shift
    result={}
    for j,name in enumerate(FIELDS+('velocity',)):
        active=mask & outer if name=='outer_epsilon' else mask
        measure=float(sum(w[active]));v=np.hypot(d[:,0],d[:,1]) if name=='velocity' else abs(d[:,j])
        ids=np.flatnonzero(active);which=ids[np.argmax(v[active])]
        result[name]=dict(L2=float(np.sqrt(w[active]@(v[active]**2))),sampled_Linf=float(v[which]),measure=measure,
            normalized_L2=float(np.sqrt(w[active]@(v[active]**2))/(AMP[name]*math.sqrt(measure))) if measure else None,
            normalized_sampled_Linf=float(v[which]/AMP[name]),maximum_provenance={k:rows[which][k] for k in
                ('position','geometry','element','side','requested_x_side','face_sides_xz','provenance')})
    return dict(fields=result,overall_pressure_gauge_shift=shift)

def comparison(a,b,policy):
    same_rule(a,b);av,_,_,w,outer=arrays(a,policy);bv,_,_,_,_=arrays(b,policy)
    return regional_norms(av-bv,w,outer,a['observations'])

def regional_norms(d,w,outer,rows):
    allnorm=physical_norms(d,w,outer,rows);g=allnorm['overall_pressure_gauge_shift']
    x=np.array([r['position'][0] for r in rows]);z=np.array([r['position'][2] for r in rows]);
    central=(x>=.25)&(x<=.75)
    # Boundary is a sampled trace norm, not a zero-measure volume integral.
    boundary=(x==min(x))|(x==max(x))|(z==2)
    bw=np.array([1. if b else 0. for b in boundary])
    return dict(full=allnorm,central=physical_norms(d,w,outer,rows,central,g),
        constrained_boundary_sample=physical_norms(d,bw,outer,rows,boundary,g),
        boundary_scope='sampled RMS and maximum on both lateral boundaries and top; own sample count, not a volume L2')

def evaluate(m,policy):
    v,initial,ref,w,outer=arrays(m,policy)
    metrics={key:regional_norms(d,w,outer,m['observations']) for key,d in
        (('continuum_error',v-ref),('initial_drift',v-initial),('initial_representation',initial-ref))}
    flux=m['interfaces'];jump=-K0/(81/20)
    metrics['interface']=dict(traction=max(np.linalg.norm(np.array(f['inner_outward_traction'])+f['outer_outward_traction']) for f in flux)/.25,
        k_flux=max(abs(sum(f['opposing_k_diffusive_flux'])) for f in flux)/(K0*.5),
        pressure_offset=max(abs(f['inner_outer_mechanical_pressure'][0]-f['inner_outer_mechanical_pressure'][1]-jump) for f in flux)/.25,
        k_trace=max(abs(np.diff(f['inner_outer_k'])[0]) for f in flux),epsilon_trace=max(abs(np.diff(f['epsilon_trace_target'])[0]) for f in flux),
        sampled_K=[min(f['sampled_K'] for f in flux),max(f['sampled_K'] for f in flux)],
        cutoff=[min(f['cutoff'] for f in flux),max(f['cutoff'] for f in flux)])
    metrics['strong_divergence']=m['strong_divergence_sampled_Linf']*2
    # Horizontal mean/spread at each common height/side, include face traces.
    spread=[]
    for z,side in sorted(set((r['position'][2],r['side']) for r in m['observations'])):
        ids=[i for i,r in enumerate(m['observations']) if r['position'][2]==z and r['side']==side]
        a=v[ids];spread.append(dict(z=z,side=side,sampled_horizontal_mean=np.mean(a,axis=0).tolist(),
            sampled_horizontal_spread=np.ptp(a,axis=0).tolist(),vertical_velocity_max=float(max(abs(a[:,1])))))
    metrics['horizontal_samples']=spread
    tests={}
    for kind in ('continuum_error','initial_drift'):
        for field,row in metrics[kind]['full']['fields'].items():
            for quantity in ('L2','sampled_Linf'):
                value=row['normalized_'+quantity];ok=math.isfinite(value) and value<=BUDGET[quantity]
                tests[kind+'/'+field+'/'+quantity]=component('PASS' if ok else 'OBSTRUCTED','engineering accuracy budget',value=value,limit=BUDGET[quantity])
    for name in ('traction','k_flux','pressure_offset'):
        value=metrics['interface'][name];tests[name]=component('PASS' if math.isfinite(value) and value<=BUDGET[name] else 'OBSTRUCTED',value=value,limit=BUDGET[name])
    value=metrics['strong_divergence'];tests['divergence']=component('PASS' if math.isfinite(value) and value<=BUDGET['divergence'] else 'OBSTRUCTED',value=value,limit=BUDGET['divergence'])
    metrics['accuracy_gate']=aggregate(tests)
    return metrics

def sensitivity_gate(diff,uncertainty,volume=2.):
    tests={}
    for field,row in diff['full']['fields'].items():
        measure=row['measure'];floor=uncertainty[field]/(AMP[field]*math.sqrt(measure))
        # Offline physical L2 arithmetic estimate is not a sampled-maximum bound.
        # Linf uncertainty separately conservatively uses physical floor * 100
        # on this finite observation inventory; it is labelled an estimate.
        for norm,uf in (('L2',floor),('sampled_Linf',100*floor)):
            v=row['normalized_'+norm];limit=.1*BUDGET[norm]
            tests[field+'/'+norm]=component('PASS' if v+uf<=limit else 'UNRESOLVED','difference plus finite diagnostic uncertainty <=10% budget',
                value=v,uncertainty=uf,limit=limit,resolved=v>uf)
    return aggregate(tests)

def independent_weak(checks,native):
    lo,H,hi=native['support_z'];xl,xh=native['x_bounds'];L=xh-xl
    expected=np.zeros((2,2,4,5))
    for q in checks['independent_quadrature']:
        x,_,z=q['coordinate']
        if z<lo or z>hi:continue
        inner=q['inner'];side=0 if inner else 1;k=q['k'];G=np.array(q['G']).reshape(3,3);u=np.array(q['velocity']);pi=q['reduced_pressure'];w=q['weight']
        gk=np.array(q['grad_k']);ge=np.array(q['grad_epsilon']);e=0 if inner else q['epsilon']
        if inner:
            a,zv,t,A,c=wall(q['sampled_K']);nut=a*max(z-zv,0);dk=1.5e-5+nut;de=0
            stress=dk*(G+G.T)-pi*np.eye(3)
            P=.5*nut*np.sum((G+G.T)**2)-2*k*np.trace(G)/3
            sk=P-(t*t*nut/dk**2 if checks['policy']=='EP' else A/max(z,c));se=0
        else:
            c=independent(k,e,G);dk=CSTAR*k*k/e;de=dk/1.3;stress=-c['R']+(2*k/3-pi)*np.eye(3)
            sk=c['scalar'][5];se=c['scalar'][8]
        h=q['accepted_history_values'];alpha=q['accepted_stage_alpha']
        for mode in (0,1):
            t=(x-xl)/L;b=t*(1-t);db=(1-2*t)/L;xx=b if mode==0 else b*(2*t-1);dx=db if mode==0 else db*(2*t-1)+2*b/L
            r=(z-lo)/(H-lo) if inner else (hi-z)/(hi-H);dr=1/(H-lo) if inner else -1/(hi-H)
            phi=xx*r;dphi=np.array([dx*r,0,xx*dr]);a=expected[mode,side]
            for f,value in enumerate([u[0],u[2],k]+([] if inner else [e])):
                a[f,0]+=w*phi*alpha*(value-(h[f] if h else value))
                a[f,1]-=w*float(dphi@u)*value
                a[f,2]+=w*float(dphi@(stress[(0,2)[f]] if f<2 else dk*gk if f==2 else de*ge))
                a[f,3]-=w*phi*(sk if f==2 else se if f==3 else 0)
                # Both support elements are above h: roughness force is zero.
                if z<=.1:raise ValueError('weak support entered roughness force region')
    errors={}
    for mode,row in enumerate(native['modes']):
        actual=np.array(row['regional_loads']);delta=float(np.max(abs(actual-expected[mode])))
        assert_allclose(actual,expected[mode],rtol=3e-9,atol=3e-10,err_msg='independent regional weak loads')
        assert_allclose(row['signed_sum'],expected[mode].sum(axis=(0,2)),rtol=3e-8,atol=5e-10)
        errors[str(mode)]=dict(max_absolute_error=delta,regional_loads=actual.tolist(),signed_sum=row['signed_sum'],prescribed_exterior_reaction=row['prescribed_exterior_reaction'])
    return dict(status='PASS',support_z=[lo,H,hi],modes=errors,
        scope='two shared tensor-cubic interface test modes; storage/convection/stress/diffusion/source from actual represented state and accepted callback; exterior reactions vanish with tests')
