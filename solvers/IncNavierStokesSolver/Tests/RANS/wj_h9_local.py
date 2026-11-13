#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only saved phases, signed moments and exact local change identities."""
import math
import numpy as np
from scipy.linalg import cho_solve
from wj_h4_points import vector_hash,aligned
from wj_h7_analysis import Rows,principal_batch,context
from wj_h8_analysis import PressureProjection,local_rows
from wj_h3_analysis import algebra

MATRIX=[('coarse',.01,1),('fine',.005,1),('fine',.005,2)]+[('quarter',.0025,s)for s in range(1,5)]
XI=np.array([1.,0.,-7.])/math.sqrt(50);POL=np.array([7.,0.,1.])/math.sqrt(50);DIRECTION=np.outer(POL,XI)

def same(a,b):return np.asarray(a,dtype='<f8').tobytes()==np.asarray(b,dtype='<f8').tobytes()

def phase_bindings(report,packet,step):
    ctx=context(report,step);prev=context(report,step-1);offset=report['offsets']
    assert packet['time']==ctx['time']and packet['alpha']==ctx['alpha']
    assert vector_hash(packet['history'])==vector_hash(ctx['history'])
    assert same(packet['fixed_reference'],report['fixed_reference'])
    assert same(packet['post_flow'],packet['pre_k'])and same(packet['raw'],packet['corrected'])
    assert packet['k_bound']['activations']==0 and packet['k_bound']['mass_norm_correction']==0
    states=dict(P=np.array(prev['state']),C=np.array(packet['committed']),F=np.array(packet['post_flow']),K=np.array(packet['raw']),A=np.array(ctx['state']))
    masks={}
    for before,after,allowed in (('C','F',(0,1,4,5)),('F','K',(2,)),('K','A',(3,))):
        for b in range(6):
            sl=slice(offset[b],offset[b+1]);changed=not same(states[before][sl],states[after][sl])
            assert b in allowed or not changed,f'{before}->{after} changes field {b}'
        masks[before+'->'+after]=[b for b in range(6)if not same(states[before][offset[b]:offset[b+1]],states[after][offset[b]:offset[b+1]])]
    for b in (0,1,4,5):assert same(states['F'][offset[b]:offset[b+1]],states['A'][offset[b]:offset[b+1]])
    return states,dict(step=step,time=ctx['time'],alpha=ctx['alpha'],stage_history_sha256=vector_hash(ctx['history']),states={k:vector_hash(v)for k,v in states.items()},
        P_C_bitwise=same(states['P'],states['C']),lift_max=float(max(abs(states['C']-states['P']))),field_masks=masks,
        post_flow_pre_k_bitwise=True,raw_corrected_bitwise=True,flow_A_U_pi_bitwise=True,accepted_role='A only; P is predecessor; C/F/K diagnostic aliases',
        initial_differential_storage_max=float(max(abs(states['C'][:offset[4]]-np.array(packet['history'])[:offset[4]]))))


def external_ids(report):
    result=set();H=report['parameters'][1];top=max(e['upper'][2]for e in report['fields'][0]['mesh'])
    for b in range(4):
        for el in report['fields'][b]['mesh']:
            m=el['degree']+1;ids=np.array(el['ids']).reshape(m,m);off=report['offsets'][b]
            if el['lower'][0]==0:result.update((off+ids[:,0]).tolist())
            if b in (0,1)and el['lower'][2]==0:result.update((off+ids[0,:]).tolist())
            if b==1 and el['upper'][2]==top:result.update((off+ids[1,:]).tolist())
        if b==3:
            for el in report['fields'][b]['mesh']:
                if el['lower'][2]==H:
                    ids=np.array(el['ids']).reshape(el['degree']+1,-1);result.difference_update((report['offsets'][b]+ids[0,:]).tolist())
    return result


def delta_B(a,b):
    def value(d):
        G=np.asarray(d['G']);return np.asarray(d['k'])/np.asarray(d['epsilon']),G-np.trace(G,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)/3
    ta,ga=value(a);tb,gb=value(b)
    gradient=(ta+tb)[...,None,None]/2*(gb-ga);timescale=(tb-ta)[...,None,None]*(ga+gb)/2;total=tb[...,None,None]*gb-ta[...,None,None]*ga
    return dict(gradient_term=gradient,timescale_term=timescale,total=total,residual=total-gradient-timescale)


def telescope(values):
    a=np.asarray(values);pieces=np.diff(a,axis=0);return pieces,a[-1]-a[0]-np.sum(pieces,axis=0)


def safe_ratio(n,d,floor=1e-15):
    if not all(np.isfinite([n,d,floor]))or min(n,d,floor)<0:raise ValueError('finite nonnegative norms required')
    return None if d<=floor else float(n/d)


class MomentAudit:
    """Full regional pressure tests. Only diagnostic physical mass factors/solves."""
    def __init__(self,rows):
        self.rows=rows;self.projection=PressureProjection(rows);self.solves=0
    def moments(self,state,solve=True):
        _,gx=self.rows.field(state,0);_,gz=self.rows.field(state,1);div=gx[:,0]+gz[:,2];out={}
        for region,(idx,w,B,factor,elements)in self.projection.regions.items():
            rhs=B.T@(w*div[idx]);v=dict(moments=rhs,integral=float(sum(rhs)))
            if solve:
                coeff=cho_solve(factor,rhs);self.solves+=1;v['projection_L2']=math.sqrt(max(0,float(rhs@coeff)))
            out[region]=v
        return out
    def identity(self,C,F):
        a=self.moments(C);b=self.moments(F);d=self.moments(F-C,False)
        return {r:dict(C=a[r]['moments'].tolist(),F=b[r]['moments'].tolist(),D_delta_U=d[r]['moments'].tolist(),
            linear_identity_max=float(max(abs(b[r]['moments']-a[r]['moments']-d[r]['moments']))),repair_max=float(max(abs(d[r]['moments']+a[r]['moments']))),
            post_flow_max=float(max(abs(b[r]['moments']))),C_projection_L2=a[r]['projection_L2'],F_projection_L2=b[r]['projection_L2'])for r in a}


def actions(k,e,G):
    k=np.asarray(k);e=np.asarray(e);G=np.asarray(G);R,nt,N,m,t,angles=principal_batch(k,e,G)
    action=-algebra(k,e,G.astype(complex)+1e-27j*DIRECTION)[0].imag/1e-27
    dim=np.sum(DIRECTION*action,axis=(-2,-1));tau=k/e;devG=G-np.trace(G,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)/3;B=tau[...,None,None]*devG;S=(B+B.swapaxes(-1,-2))/2;O=(B-B.swapaxes(-1,-2))/2
    s=np.trace(S@S,axis1=-2,axis2=-1);w=np.trace(O@O,axis1=-2,axis2=-1);FN=3*N*N-3.6*N-2*w-2.7*s
    db=tau[...,None,None]*(DIRECTION-np.trace(DIRECTION)*np.eye(3)/3);ds=2*np.trace(S@((db+db.swapaxes(-1,-2))/2),axis1=-2,axis2=-1);dw=2*np.trace(O@((db-db.swapaxes(-1,-2))/2),axis1=-2,axis2=-1)
    dN=(2.7*N*ds+2*(N-1.8)*dw)/FN;csN=algebra(k,e,G.astype(complex)+1e-27j*DIRECTION)[2].imag/1e-27
    return dict(R=R,nu_eff=nt,N=N,F_N=FN,root_residual=(N-1.8)*(N*N-2*w)-2.7*s*N,implicit_dN=dN,implicit_derivative_error=dN-csN,
        B_t=B,dimensional_fixed_action=dim,normalized_fixed_action=dim/nt,in_plane=m,transverse=t,minimizing_angle=angles,raw_production=-np.sum(R*G,axis=(-2,-1)))


def matrix_check(packet,path):
    rows=np.asarray(packet['transport_rows'],int);n=len(rows);M=np.fromfile(path,dtype='<f8');assert M.size==n*n;M=M.reshape(n,n)
    step=np.array(packet['raw_increment']);res=np.array(packet['transport_residual']);saved=np.array(packet['raw_linear_residual']);d=np.array(packet['raw'])[rows]-np.array(packet['pre_k'])[rows]
    defect=M@step+res;from_state=M@d+res;bound=64*np.finfo(float).eps*(abs(M)@abs(step)+abs(res))+2e-20
    assert np.all(abs(defect-saved)<=bound);assert max(abs(d-step))<2e-17
    assert same(packet['raw'],packet['corrected'])
    return dict(shape=[n,n],matvecs=3,factorizations=0,solves=0,recorded_defect_max=float(max(abs(saved))),matvec_recorded_error=float(max(abs(defect-saved))),
        state_increment_error=float(max(abs(d-step))),state_matvec_residual_max=float(max(abs(from_state))),raw_corrected_bitwise=True)
