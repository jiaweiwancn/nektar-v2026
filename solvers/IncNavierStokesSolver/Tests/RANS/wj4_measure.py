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
    for xl,xh in zip(xs,xs[1:]):
        for x in xl+(gx+1)*(xh-xl)/2:
            for side in (-1,1):rows.append([float(x),0.,1.,0.,side,0])
    for zl,zh in zip(zs,zs[1:]):
        for z in zl+(gz+1)*(zh-zl)/2:
            for x in (xs[0],xs[-1]):rows.append([float(x),0.,float(z),0.,0,0])
    for xl,xh in zip(xs,xs[1:]):
        for x in xl+(gx+1)*(xh-xl)/2:rows.append([float(x),0.,2.,0.,0,0])
    # Complete face intersections plus common dense horizontal maximum samples.
    # Explicit x-side is distinct from the radial side at H.
    samplex=sorted(set(xs+np.linspace(xs[0],xs[-1],33).tolist()))
    for x in samplex:
        for z in zs:
            for side in ((-1,1) if z==1 else (0,)):
                for xside in ((-1,1) if x in xs[1:-1] else (0,)):
                    rows.append([float(x),0.,float(z),0.,side,xside])
    return rows

def boundary_rule(order=6,xorder=4):
    _,zv,_,_,zc=wall(K0)
    xs=[0,.25,.5,.75,1];zs=sorted(set(normal_knots(2)+[z for z in (zc,zv,.1,.5,1) if 0<z<2]))
    gx,wx=np.polynomial.legendre.leggauss(xorder);gz,wz=np.polynomial.legendre.leggauss(order);weights={}
    for lo,hi in zip(zs,zs[1:]):
        for z,w in zip(lo+(gz+1)*(hi-lo)/2,wz*(hi-lo)/2):
            for x in (0.,1.):weights[(x,float(z))]=float(w)
    for lo,hi in zip(xs,xs[1:]):
        for x,w in zip(lo+(gx+1)*(hi-lo)/2,wx*(hi-lo)/2):weights[(float(x),2.)]=float(w)
    return weights

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

def regional_norms(d,w,outer,rows,xorder=4):
    allnorm=physical_norms(d,w,outer,rows);g=allnorm['overall_pressure_gauge_shift']
    x=np.array([r['position'][0] for r in rows]);z=np.array([r['position'][2] for r in rows]);
    central=(x>=.25)&(x<=.75)
    boundary=(x==min(x))|(x==max(x))|(z==2)
    weights=boundary_rule(xorder=xorder)
    bw=np.array([weights.get((r['position'][0],r['position'][2]),0.) for r in rows])
    return dict(full=allnorm,central=physical_norms(d,w,outer,rows,central,g),
        constrained_boundary=physical_norms(d,bw,outer,rows,boundary,g),
        boundary_scope='physical trace L2 on both lateral boundaries and top; own surface measure (5; epsilon3); pressure sampled there is an algebraic multiplier')

def evaluate(m,policy,xorder=4):
    v,initial,ref,w,outer=arrays(m,policy)
    metrics={key:regional_norms(d,w,outer,m['observations'],xorder) for key,d in
        (('continuum_error',v-ref),('initial_drift',v-initial),('initial_representation',initial-ref))}
    flux=m['interfaces'];jump=-K0/(81/20)
    metrics['interface']=dict(traction=max(np.linalg.norm(np.array(f['inner_outward_traction'])+f['outer_outward_traction']) for f in flux)/.25,
        k_flux=max(abs(sum(f['opposing_k_diffusive_flux'])) for f in flux)/(K0*.5),
        pressure_offset=max(abs(f['inner_outer_mechanical_pressure'][0]-f['inner_outer_mechanical_pressure'][1]-jump) for f in flux)/.25,
        k_trace=max(abs(np.diff(f['inner_outer_k'])[0]) for f in flux),epsilon_trace=max(abs(np.diff(f['epsilon_trace_target'])[0]) for f in flux),
        sampled_K=[min(f['sampled_K'] for f in flux),max(f['sampled_K'] for f in flux)],
        cutoff=[min(f['cutoff'] for f in flux),max(f['cutoff'] for f in flux)])
    metrics['strong_divergence']=m['strong_divergence_sampled_Linf']*2
    rows=m['observations'];div=np.array([abs(r['raw_G'][0]+r['raw_G'][8])*2 for r in rows])
    masks=dict(full=np.ones(len(rows),bool),central=np.array([.25<=r['position'][0]<=.75 for r in rows]),
        constrained_boundary=np.array([r['position'][0] in (0.,1.) or r['position'][2]==2 for r in rows]))
    metrics['divergence_by_scope']={}
    for scope,mask in masks.items():
        ids=np.flatnonzero(mask);i=int(ids[np.argmax(div[mask])]);r=rows[i]
        metrics['divergence_by_scope'][scope]=dict(normalized_maximum=float(div[i]),
            witness={k:r[k] for k in ('position','geometry','element','side','requested_x_side','face_sides_xz','raw_G','value')})
    # Complement interface Gauss points with both H sides at every common
    # horizontal face/corner and maximum-observation point.
    bykey={}
    for r in rows:
        if r['position'][2]==1:bykey.setdefault((r['position'][0],r['requested_x_side']),{})[r['side']]=r
    face=[]
    for key,pair in bykey.items():
        if -1 not in pair or 1 not in pair:raise ValueError('missing H observation side')
        i,o=pair[-1],pair[1];K=i['value'][2];a,zv,t,A,c=wall(K)
        fi=(1.5e-5+a*max(1-zv,0))*i['grad_k'][2];fo=-CSTAR*o['value'][2]**2/o['value'][3]*o['grad_k'][2]
        traction=np.array(i['stress'])[:,2]-np.array(o['stress'])[:,2]
        face.append(dict(x=key[0],requested_x_side=key[1],inner_geometry=i['geometry'],outer_geometry=o['geometry'],
            traction=float(np.linalg.norm(traction))/.25,traction_jump=traction.tolist(),k_flux_jump=fi+fo,
            k_flux=abs(fi+fo)/(K0*.5),pressure_offset=abs(i['value'][5]-o['value'][5]-jump)/.25,
            sampled_K=K,cutoff=c,epsilon_trace=abs(o['value'][3]-A)))
    for field in ('traction','k_flux','pressure_offset','epsilon_trace'):
        metrics['interface'][field]=max(metrics['interface'][field],max(r[field] for r in face))
    metrics['face_interfaces']=face
    gx,wx=np.polynomial.legendre.leggauss(xorder);cuts=[0,.25,.5,.75,1]
    weights={float(x):float(w) for lo,hi in zip(cuts,cuts[1:]) for x,w in zip(lo+(gx+1)*(hi-lo)/2,wx*(hi-lo)/2)}
    moments=[]
    for mode in (0,1):
        independent_value=np.zeros(4);native_value=np.zeros(4)
        for r in face:
            if r['x'] in weights and r['requested_x_side']==0:
                x=r['x'];phi=x*(1-x)*(2*x-1)**mode
                independent_value+=weights[x]*phi*np.array(r['traction_jump']+[r['k_flux_jump']])
        for r in m['interfaces']:
            x=r['position'][0];phi=x*(1-x)*(2*x-1)**mode
            native_value+=r['weight']*phi*np.array((np.array(r['inner_outward_traction'])+r['outer_outward_traction']).tolist()+[sum(r['opposing_k_diffusive_flux'])])
        assert_allclose(native_value,independent_value,rtol=2e-7,atol=5e-10,err_msg='independent interface trace moments')
        moments.append(dict(mode=mode,native_Gauss8=native_value.tolist(),independent_composite=independent_value.tolist(),max_absolute_difference=float(max(abs(native_value-independent_value)))))
    metrics['interface_trace_moments']=dict(order_x=xorder,panels=cuts,moments=moments,
        scope='pointwise constitutive traction/diffusive flux moments; separate from signed regional volume loads and algebraic value trace')
    # Horizontal mean/spread at each common height/side, include face traces.
    spread=[]
    for z,side in sorted(set((r['position'][2],r['side']) for r in m['observations'])):
        ids=[i for i,r in enumerate(m['observations']) if r['position'][2]==z and r['side']==side]
        a=v[ids];spread.append(dict(z=z,side=side,sampled_horizontal_mean=np.average(a,axis=0,weights=w[ids] if sum(w[ids])>0 else None).tolist(),
            mean_rule='composite physical x Gauss weights' if sum(w[ids])>0 else 'zero-weight trace inventory arithmetic mean',
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
    for scope in ('full','central','constrained_boundary'):
        for field,row in diff[scope]['fields'].items():
            measure=row['measure'];entry=uncertainty[field]
            if not isinstance(entry,dict):entry=dict(L2=entry,sampled_Linf=entry)
            sampled=entry['sampled_Linf']/AMP[field]
            # A trace L2 is bounded using the finite point envelope and its own
            # surface measure, never by reinterpreting a volume L2 uncertainty.
            floor=sampled if scope=='constrained_boundary' else entry['L2']/(AMP[field]*math.sqrt(measure))
            for norm,uf in (('L2',floor),('sampled_Linf',sampled)):
                v=row['normalized_'+norm];limit=.1*BUDGET[norm]
                tests[scope+'/'+field+'/'+norm]=component('PASS' if v+uf<=limit else 'UNRESOLVED','difference plus finite diagnostic uncertainty <=10% budget',
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
                test=phi if f!=3 else phi*(1-r);gradient=dphi if f!=3 else np.array([dx*r*(1-r),0,xx*dr*(1-2*r)])
                a[f,0]+=w*test*alpha*(value-(h[f] if h else value))
                a[f,1]-=w*float(gradient@u)*value
                a[f,2]+=w*float(gradient@(stress[(0,2)[f]] if f<2 else dk*gk if f==2 else de*ge))
                a[f,3]-=w*test*(sk if f==2 else se if f==3 else 0)
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
