# SPDX-License-Identifier: MIT
"""Independent flat ground/interface traces from represented native polynomials.

One-sided element polynomials are evaluated separately; stresses are recomputed
from their actual gradients/scalars. No stress interpolation or side averaging.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss,legvander
from cao_tamura_h2g_constitutive import parts
from wp2_all_state import wall,US
from run_cls_epawf_picard import require


def traces(rows,state,report,order=36):
    d=rows.tuple(state,report);z=np.asarray(report['wp2_common_samples']['physical_z_partition'])
    nz=len(z)-1;raw=d['values'].reshape(6,8,nz,8,14).transpose(0,2,1,3,4)
    old,w=leggauss(8);v=legvander(old,3);projection=(2*np.arange(4)+1)[:,None]/2*v.T*w[None,:]
    x,wx=leggauss(order);ix=legvander(x,3)@projection
    xknots=np.array([-.2,-.1,-.05,0,.05,.1,.6]);X=(xknots[:-1,None]+(x[None,:]+1)*np.diff(xknots)[:,None]/2).reshape(-1)
    weight=(np.diff(xknots)[:,None]/2*wx).reshape(-1)
    nu,H,hs=report['parameters'][:3];q,kap=report['parameters'][6:8]
    h=np.flatnonzero(z==H);require(len(h)==1 and h[0]>0,'missing exact physical interface')
    kh=(d['endpoint_K'].reshape(6,8)@ix.T).reshape(-1)
    gk=d['gk'].reshape(6,8,nz,8,3).transpose(0,2,1,3,4)
    output={}
    for label,iz,side,height,inner in [('ground',0,-1.,0.,True),('H_minus',h[0]-1,1.,H,True),('H_plus',h[0],-1.,H,False)]:
        endpoint=(legvander([side],3)@projection)[0]
        value=np.einsum('ai,xijf,j->xaf',ix,raw[:,iz],endpoint,optimize=True).reshape(-1,14)
        gradk=np.einsum('ai,xijf,j->xaf',ix,gk[:,iz],endpoint,optimize=True).reshape(-1,3)
        G=np.zeros((len(value),3,3));G[:,0,0]=value[:,10];G[:,0,2]=value[:,11];G[:,2,0]=value[:,12];G[:,2,2]=value[:,13]
        physical=np.column_stack((X,np.zeros(len(X)),np.full(len(X),height)))
        auxiliary=wall(kh,np.full(len(X),height),nu,H,hs,q,kap)
        epsilon=auxiliary['E'] if inner else value[:,3]
        mask=np.full(len(X),inner)
        law='calibrated_standard' if report['model_specification']['name']=='StandardKEpsilon' else 'regularized_CLS'
        closure=parts(G,value[:,2],epsilon,mask,auxiliary['nt'],physical,US,kap,law)
        T=-closure['tensor']+(2*value[:,2]/3-value[:,4])[:,None,None]*np.eye(3)
        if inner:T+=nu*(G+G.transpose(0,2,1))
        # Only primitive coefficients/gradients determine the result. Replace
        # every derived sample column; no interpolated stress is exported.
        value[:,3]=epsilon;value[:,5]=value[:,4]-2*value[:,2]/3
        value[:,6]=closure['nu_t'];value[:,7]=T[:,0,2];value[:,8]=T[:,2,2];value[:,9]=closure['production']
        output[label]=dict(values=value,G=G,gk=gradk,R=closure['tensor'],T=T,epsilon=epsilon,native_z_panel=iz,
            k_flux_z=(nu+closure['nu_t'] if inner else closure['nu_t'])*gradk[:,2],physical=physical)
    jump=output['H_minus']['T'][:,:,2]-output['H_plus']['T'][:,:,2]
    flux=output['H_minus']['k_flux_z']-output['H_plus']['k_flux_z']
    K0=report['parameters'][5];rms=lambda a:float(np.sqrt(np.dot(weight,a*a)/sum(weight)))
    summary=dict(order=order,wall_mean_shear=float(np.dot(weight,output['ground']['T'][:,0,2])/sum(weight)),
        wall_shear_RMS_error_over_ustar2=rms(output['ground']['T'][:,0,2]-US**2)/US**2,
        interface_normal_traction_jump_RMS=rms(jump[:,2]),interface_full_traction_jump_RMS=float(np.sqrt(np.dot(weight,np.sum(jump*jump,axis=1))/sum(weight))),
        interface_k_flux_jump_RMS=rms(flux),minimum_actual_outer_trace_epsilon=float(min(output['H_plus']['epsilon'])),
        actual_dependent_trace_epsilon_relative_max=float(max(abs(output['H_plus']['epsilon']-output['H_minus']['epsilon'])/output['H_minus']['epsilon'])),
        covariance_minima={n:float(np.linalg.eigvalsh(o['R'])[:,0].min()) for n,o in output.items()},
        scope='represented physical one-sided affine polynomials; immutable analytical reference is separate; sigmaK=1 in both selected models; finite surface sampling complements native sided hard audit')
    summary['covariance_witnesses']={}
    geometry=rows.a['geometry'][:,0].reshape(6,8,nz,8)[:,0,:,0]
    elements=rows.a['elements'][:,0].reshape(6,8,nz,8)[:,0,:,0]
    for label,o in output.items():
        eigen,vector=np.linalg.eigh(o['R']);j=int(np.argmin(eigen[:,0]));panel=j//order;iz=o['native_z_panel']
        summary['covariance_witnesses'][label]=dict(physical=o['physical'][j].tolist(),geometry=int(geometry[panel,iz]),element=int(elements[panel,iz]),
            radial_side=-1 if label=='H_minus' else 1 if label=='H_plus' else 0,tangential_side=0,
            k=float(o['values'][j,2]),epsilon=float(o['epsilon'][j]),epsilon_semantics='actual retained outer' if label=='H_plus' else 'dependent E_C',
            gradient=o['G'][j].tolist(),tensor=o['R'][j].tolist(),eigenvalues=eigen[j].tolist(),critical_eigenvector=vector[j,:,0].tolist())
    summary['covariance_margin_passed']=min(summary['covariance_minima'].values())>=1e-4*K0
    summary['positive_scalars']=all(np.all(o['values'][:,2]>0) and np.all(o['epsilon']>0) and np.all(np.isfinite(o['R'])) for o in output.values())
    summary['passed']=summary['covariance_margin_passed'] and summary['positive_scalars']
    return dict(summary=summary,surfaces=output,weights=weight)


def compare(left,right,K0):
    require(np.array_equal(left['weights'],right['weights']),'surface integration differs')
    w=left['weights'];out={}
    for label in ['ground','H_minus','H_plus']:
        a,b=left['surfaces'][label],right['surfaces'][label]
        require(np.array_equal(a['physical'],b['physical']),'surface physical points differ')
        out[label+'_traction_RMS_over_ustar2']=float(np.sqrt(np.sum(w[:,None]*(a['T'][:,:,2]-b['T'][:,:,2])**2)/sum(w))/US**2)
        out[label+'_k_flux_RMS_over_ustar_K0']=float(np.sqrt(np.dot(w,(a['k_flux_z']-b['k_flux_z'])**2)/sum(w))/(US*K0))
    return dict(values=out,limit=.0025,passed=max(out.values())<=.0025)
