#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Four-owner saved-phase algebra. No volume rows, mass system or PDE entrypoint."""
import hashlib
import math
import numpy as np
from validate_cls_epawf_m3 import basis
from wj_h4_points import mapping,owner_point,field,aligned,vector_hash
from wj_h7_analysis import context
from wj_h11_analysis import selected,callback_relation
from wj_h9_local import phase_bindings,same,actions,DIRECTION,delta_B
from wj_h3_analysis import algebra,principal
from wj_h5_neutral import CD

EPS=np.finfo(float).eps
BRANCHES=(('coarse',.01),('fine',.005),('quarter',.0025))
WINDOWS=((0.,.01),(.03,.04))
ATTEMPTS=tuple((w,b,dt,s)for w,(lo,hi)in enumerate(WINDOWS)for b,dt in BRANCHES
               for s in range(round(lo/dt)+1,round(hi/dt)+1))
OWNER_KEYS=((17,'upper',.25607517888672315),(18,'lower',.25607517888672315),
            (18,'upper',.31946961248954786),(19,'lower',.31946961248954786))
PHASES=('P','C','F','K','A')
OBSERVABLES=('U','G','devG','k','epsilon','tau','B','R','pi','mechanical_pressure',
             'd','nu_eff','a','g','mhat','transverse','Praw','N','F_N')

def scale(name):
    return .001 if name in ('d','nu_eff')else .008809015145515192 if name=='R'else 1.

def comparison_scale(name):
    if name in ('k','R'):return .008809015145515192
    if name=='epsilon':return .051**3/(.381917863546733*.15)
    if name in ('pi','mechanical_pressure'):return .051**2
    return scale(name)

def identity(residual,terms,amplitude=1.):
    """Componentwise roundoff bound, independent of physical comparison floors."""
    residual=np.asarray(residual);bound=128*EPS*(amplitude+sum(abs(np.asarray(x))for x in terms))
    if np.any(~np.isfinite(residual))or np.any(abs(residual)>bound):
        raise ValueError('algebraic identity exceeds componentwise roundoff bound')
    return dict(maximum_error=float(np.max(abs(residual))),maximum_bound=float(np.max(bound)),
                residual=residual.tolist(),bound=bound.tolist(),status='PASS')

def checked_context(r,step):
    c=context(r,step)
    if c['step']!=step:raise ValueError('context has a different integer step')
    aligned(step,r['dt'],c['time'])
    return c

def role_context(binding,phase):
    if phase not in PHASES:raise ValueError('unknown saved phase')
    return dict(attempt_time=binding['time'],time=binding['predecessor_time']if phase=='P'else binding['time'],
        state_absolute_step=binding['predecessor_step']if phase=='P'else binding['step']if phase=='A'else None,
        accepted_role=phase in ('P','A'),callback_history_sha256=binding['stage_history_sha256'],alpha=binding['alpha'])

def endpoint_margin(record):
    """H11 fixed-owner records export margin, not H4 scan's principal dict."""
    value=float(record['margin'])
    if not math.isfinite(value):raise ValueError('nonfinite published margin')
    return value

def bind(r,packet,step,producer,expected_producer):
    if producer!=expected_producer:raise ValueError('wrong original attempt producer')
    c,attempt,local=selected(r,step);p=checked_context(r,step-1)
    checked_context(r,step)
    states,v=phase_bindings(r,packet,step)
    if not v['P_C_bitwise']:raise ValueError('unexpected ExternalLift change')
    if vector_hash(packet['history'])!=attempt['history_sha256']or attempt['incoming_sha256']!=v['states']['C']:
        raise ValueError('attempt/snapshot/history binding differs')
    if attempt['scalar_order']!='k_then_epsilon':raise ValueError('wrong native phase order')
    if attempt['bound_activations']!=0:raise ValueError('material bound event')
    v.update(local_attempt=local,producer=producer,predecessor_step=p['step'],predecessor_time=p['time'],
             callback_relation=callback_relation(c['history'],p['state']),
             fixed_reference_sha256=vector_hash(r['fixed_reference']),fixed_defect_sha256=vector_hash(r['fixed_defect']))
    return states,v

def validate_owner(r,record,index):
    p=record['physical'];o=record['owner'];geo,side,z=OWNER_KEYS[index]
    if p!=[8.,0.,z]or o['geometry']!=geo or o['normal_side']!=side or o['x_side']!='upper' or o['region']!='outer':
        raise ValueError('wrong exact physical owner')
    el=next(v for v in r['fields'][0]['mesh']if v['geometry']==geo)
    fresh=owner_point(r,el,p)
    for k in fresh:
        if fresh[k]!=o[k]:raise ValueError('owner map/snap provenance differs: '+k)

class SparseRows:
    """Only four tensor-element rows; cache native signs and immutable map once."""
    def __init__(self,r,records):
        if len(records)!=4:raise ValueError('exactly four owners required')
        self.r=r;self.records=records;self.rows={}
        for i,rec in enumerate(records):
            validate_owner(r,rec,i);o=rec['owner']
            for block in range(6):
                el=next((e for e in r['fields'][block]['mesh']if e['geometry']==o['geometry']),None)
                if el is None:continue
                lo,hi=el['lower'],el['upper'];x,_,rho=o['chart']
                px,dx=basis(el['degree'],[2*(x-lo[0])/(hi[0]-lo[0])-1])
                pr,dr=basis(el['degree'],[2*(rho-lo[2])/(hi[2]-lo[2])-1])
                J=mapping(r)[3]
                rows=np.array([np.outer(pr[0],px[0]).ravel(),np.outer(pr[0],dx[0]).ravel()*2/(hi[0]-lo[0]),
                               np.outer(dr[0],px[0]).ravel()*2/(hi[2]-lo[2])/J])*np.asarray(el['sign'])
                self.rows[i,block]=(np.asarray(el['ids'],int)+int(r['offsets'][block]),rows)
    def block(self,state,block,i):
        ids,rows=self.rows[i,block];v=rows@np.asarray(state)[ids]
        return float(v[0]),np.array([v[1],0.,v[2]])
    def linear(self,state):
        out=[]
        for i,rec in enumerate(self.records):
            U=np.zeros(3);G=np.zeros((3,3))
            for block,axis in ((0,0),(1,2)):U[axis],G[axis]=self.block(state,block,i)
            k,gk=self.block(state,2,i);e,ge=self.block(state,3,i);pi,gpi=self.block(state,5,i)
            out.append(dict(U=U,G=G,trace_G=float(np.trace(G)),devG=G-np.trace(G)*np.eye(3)/3,k=k,epsilon=e,pi=pi,
                            mechanical_pressure=pi-2*k/3,grad_k=gk,grad_epsilon=ge,grad_pi=gpi))
        return out
    def reference_check(self,state):
        errors=[]
        for i,rec in enumerate(self.records):
            for b in (0,1,2,3,5):
                a,g=self.block(state,b,i);v,h=field(self.r,state,b,rec['owner'])
                errors.extend([abs(a-v),*abs(g-h)])
        return float(max(errors))

def local(linear):
    """All closure calculations use actual phases, never scalar/gradient hybrids."""
    v=dict(linear);k,e,G=v['k'],v['epsilon'],v['G']
    if not np.isfinite(k+e)or min(k,e)<=0:
        return dict(v,constitutive_status='NOT_EVALUABLE',reason='nonpositive/nonfinite scalar')
    try:
        raw=actions(np.array([k]),np.array([e]),G[None,:,:]);a={key:np.asarray(x)[0]for key,x in raw.items()}
        if any(np.any(~np.isfinite(x))for x in a.values())or a['nu_eff']<=0 or a['F_N']<=0:raise ValueError('undefined closure/tangent')
        v.update(R=a['R'],nu_eff=a['nu_eff'],N=a['N'],F_N=a['F_N'],root_residual=a['root_residual'],
                 tau=k/e,B=a['B_t'],nu_d=CD*k*k/e,Praw=a['raw_production'],implicit_dN=a['implicit_dN'],
                 implicit_derivative_error=a['implicit_derivative_error'],d=a['dimensional_fixed_action'],a=a['normalized_fixed_action'],
                 mhat=a['in_plane'],transverse=a['transverse'],minimizing_angle=a['minimizing_angle'])
        v['g']=v['a']-v['mhat'];v['covariance']=np.linalg.eigvalsh(v['R'])
        angle=v['minimizing_angle'];xi=np.array([math.cos(angle),0,math.sin(angle)]);h=np.outer([-xi[2],0,xi[0]],xi)
        v['direction_projector']=np.outer(h.ravel(),h.ravel());v['minimizing_xi_projector']=np.outer(xi,xi)
        O=(v['B']-v['B'].T)/2;c=3*v['N']/(5*(v['N']**2-2*np.trace(O@O)))
        v['c_B']=c;v['r_B']=v['R']/k
        v['scaling_nu_identity']=identity(v['nu_eff']-k*k/e*c,[v['nu_eff'],k*k/e*c],.001)
        # Homogeneity is checked by differentiation of the very same r(B)
        # formula at its actual B. This is a variable change, not a saved hybrid.
        devH=DIRECTION-np.trace(DIRECTION)*np.eye(3)/3
        RB=algebra(1.,1.,v['B'].astype(complex)+1e-27j*devH)[0]
        DB=RB.imag/1e-27
        dscaled=-k*k/e*np.sum(DIRECTION*DB)
        v['scaling_stress_identity']=identity(v['R']-k*RB.real,[v['R'],k*RB.real],.008809015145515192)
        v['scaling_tangent_identity']=identity(v['d']-dscaled,[v['d'],dscaled],.001)
        v['denominator_Q']=5/6*(v['N']**2-2*np.trace(O@O))*(2*v['N']**2-np.trace(O@O))
        v['physical_checks']=dict(k_positive=k>0,epsilon_positive=e>0,physical_root=v['N']>=1.8,
            root_regular=v['F_N']>0,denominator_positive=v['denominator_Q']>0,nu_eff_positive=v['nu_eff']>0,
            nu_d_positive=v['nu_d']>0,covariance_min=float(min(v['covariance'])),raw_production=float(v['Praw']),
            in_plane=float(v['mhat']),transverse=float(v['transverse']))
        v['constitutive_status']='DEFINED'
        v['physical_status']='WARNING'if min(v['covariance'])<0 or min(v['mhat'],v['transverse'])<0 else 'SAMPLED_POSITIVE'
        v['gap_status']='UNRESOLVED'if abs(v['g'])<=2e-11 else 'NONNEGATIVE'if v['g']>0 else 'INCONSISTENT'
        if v['g']< -2e-11:raise ValueError('negative optimizer gap')
    except (ValueError,ArithmeticError,np.linalg.LinAlgError)as err:
        if 'algebraic identity'in str(err)or 'optimizer gap'in str(err):raise
        return dict(v,constitutive_status='NOT_EVALUABLE',reason=str(err))
    return v

def decompositions(before,after):
    if any(v['constitutive_status']!='DEFINED'for v in (before,after)):
        return dict(status='NOT_EVALUABLE',reason='undefined actual phase')
    bt=delta_B(before,after)
    mb=identity(bt['residual'],[bt['total'],bt['gradient_term'],bt['timescale_term']])
    dd=after['d']-before['d'];dn=1/after['nu_eff']-1/before['nu_eff']
    dimensional=(1/after['nu_eff']+1/before['nu_eff'])/2*dd
    normalization=(after['d']+before['d'])/2*dn
    da=after['a']-before['a'];dg=after['g']-before['g'];dm=after['mhat']-before['mhat']
    return dict(status='PASS',delta_B=bt,B_identity=mb,delta_d=dd,delta_inverse_nu=dn,
                dimensional_term=dimensional,normalization_term=normalization,delta_a=da,delta_g=dg,delta_mhat=dm,
                normalized_identity=identity(da-dimensional-normalization,[da,dimensional,normalization]),
                optimizer_identity=identity(dm-da+dg,[dm,da,dg]))

def increments(v,name):
    if any(name not in v[p]for p in PHASES):return None
    return {tag:np.asarray(v[b][name])-v[a][name]for tag,a,b in
            (('lift','P','C'),('flow','C','F'),('k','F','K'),('epsilon','K','A'))}

def window(rows,name):
    pieces=[increments(v,name)for v in rows]
    if any(x is None for x in pieces):return dict(status='NOT_EVALUABLE')
    total={p:sum(v[p]for v in pieces)for p in ('lift','flow','k','epsilon')}
    start=np.asarray(rows[0]['P'][name]);end=np.asarray(rows[-1]['A'][name]);terms=list(total.values())
    return dict(status='PASS',start=start,end=end,cumulative=total,
                identity=identity(end-start-sum(terms),[start,end,*terms],scale(name)))

def cancellation(start,contributions,end,amplitude):
    S=abs(np.asarray(start))+sum(abs(np.asarray(v))for v in contributions)
    floor=4e-11*amplitude;ratio=np.full(np.shape(S),np.nan)
    np.divide(abs(end),S,out=ratio,where=S>floor)
    return dict(S=S,arithmetic_floor=floor,cancellation=np.where(S>floor,1-ratio,np.nan),resolved=S>floor)

def paired_window(a,b,name):
    if a['status']!='PASS'or b['status']!='PASS':return dict(status='NOT_EVALUABLE')
    start=a['start']-b['start'];end=a['end']-b['end'];terms={p:a['cumulative'][p]-b['cumulative'][p]for p in a['cumulative']}
    return dict(status='PASS',incoming=start,endpoint=end,contributions=terms,
                identity=identity(end-start-sum(terms.values()),[end,start,*terms.values()],scale(name)),
                cancellation=cancellation(start,list(terms.values()),end,comparison_scale(name)))

def kernel_input(v):
    return hashlib.sha256(np.asarray([v['k'],v['epsilon'],*np.asarray(v['G']).ravel()],dtype='<f8').tobytes()).hexdigest()

def scatter(rows,values,size):
    ids=np.asarray(rows,int);v=np.asarray(values)
    if len(ids)!=len(v)or len(set(ids))!=len(ids)or min(ids)<0 or max(ids)>=size:raise ValueError('invalid native velocity map')
    out=np.zeros(size);out[ids]=v;return out

def predictor_eligible(branch,step):
    return (branch,step)in (('coarse',1),('fine',2),('quarter',4))

def sparse_flow(cache,states,packet,branch,step):
    if not predictor_eligible(branch,step):raise ValueError('no saved predictor coverage for this attempt')
    if packet['binding']['states']!={p:vector_hash(v)for p,v in states.items()}:
        raise ValueError('arithmetic phases are from a different attempt')
    if packet['status']!='PASS' or packet['checks']['velocity_update']['status']!='PASS':raise ValueError('unqualified saved arithmetic')
    pred=scatter(packet['velocity_rows'],packet['predictor'],len(states['C']))
    correction=scatter(packet['velocity_rows'],packet['pressure_correction_velocity'],len(states['C']))
    bg=scatter(packet['velocity_rows'],packet['BGp_delta_pi'],len(states['C']))
    # No governing matvec: these are already sealed signed vectors.
    sign=identity(correction+bg,[correction,bg]);u=cache.linear(states['F']-states['C'])
    vp=cache.linear(pred);vc=cache.linear(correction);out=[]
    base=cache.linear(states['C'])
    for j,(total,p,c)in enumerate(zip(u,vp,vc)):
        tau=base[j]['k']/base[j]['epsilon'];dev=lambda G:G-np.trace(G)*np.eye(3)/3
        out.append(dict(predictor_U=p['U'],pressure_U=c['U'],total_U=total['U'],predictor_G=p['G'],pressure_G=c['G'],total_G=total['G'],
                        U_identity=identity(total['U']-p['U']-c['U'],[total['U'],p['U'],c['U']]),
                        G_identity=identity(total['G']-p['G']-c['G'],[total['G'],p['G'],c['G']]),
                        predictor_delta_B=tau*dev(p['G']),pressure_delta_B=tau*dev(c['G']),
                        total_delta_B=tau*dev(total['G']),scope='actual increments only; no constitutive state at predictor or correction'))
    return dict(status='PASS',sealed_sign_identity=sign,owners=out)
