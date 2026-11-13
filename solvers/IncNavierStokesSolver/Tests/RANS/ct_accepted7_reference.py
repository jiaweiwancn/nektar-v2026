#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Displaced experimental fit and an EXPLICIT proposed wall-policy oracle.

This module never calls the production cutoff kernel or changes native data.
The proposed shifted sink is not the original Accepted appendix policy.
"""
import argparse,csv,itertools,json,math,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.stats import qmc
from ct_accepted4_reference import observations,resistance
from ct_accepted_feasibility import Column,parameters,require
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha

MODELS=('A93','A95','E94E95')
HEIGHTS=(5.5,6.,6.5,6.99)


def fit(model,z,v,normal,d=.0035):
    p=parameters(model);us=.051/normal
    C=float(np.mean(v[:4]/normal-us/p['kappa']*np.log((z[:4]-d)/.04)))
    z0=.04*math.exp(-C*p['kappa']/us)
    predicted=us/p['kappa']*np.log((z-d)/z0)
    return dict(C=C,z0_m=z0,ustar=us,K=p['q']*us*us,model_parameters=p,
        predicted_U_Hh=us/p['kappa']*math.log((.04-d)/z0),
        predicted=predicted.tolist(),residuals=(predicted-v/normal).tolist())


def wall_root(model,H,f,nu,disp):
    target=f['ustar']/parameters(model)['kappa']*math.log((H-disp)/(f['z0_m']/.04))
    fun=lambda hs:f['ustar']**2*resistance(model,f['K'],H,hs,nu)-target
    end=np.nextafter(H,0.);left,right=fun(0.),fun(end)
    if left*right>0:return dict(status='unattainable_original_roughness_interval',target=target,bracket_residuals=[left,right])
    hs=brentq(fun,0.,end,xtol=2e-14)
    return dict(status='root',h_s=hs,h_s_mm=hs*40,target=target,root_defect=fun(hs),bracket_residuals=[left,right])


def height_limit(model,f,nu,disp):
    def residual(H):
        return f['ustar']**2*resistance(model,f['K'],H,H,nu)-f['ustar']/parameters(model)['kappa']*math.log((H-disp)/(f['z0_m']/.04))
    return brentq(residual,disp+1e-5,.5,xtol=2e-14)*40


class DisplacedColumn(Column):
    """Candidate only; J=1+gamma*r, fixed physical wall at r=0."""
    def __init__(self,model,K,H,hr,nu,disp,order=48,gamma=0.):
        self.disp=disp
        # Initialize viscosity/resistance without invoking the base class's
        # original curved cutoff solve with this different sink policy.
        super().__init__(model,K,H,hr,nu,curvature=0.,order=order)
        self.curvature=-gamma
        self.cutoff=None
        self.PJ=self.integrate(lambda r:self.J(r)*self.production(r))
        require(0<=disp<H and min(1.,1+gamma*H)>0,'invalid shifted reference geometry')
        self.cutoff=None
        self.minimum_budget=self.A*(H+gamma*H*H/2)/(H-disp)
        if self.PJ>=self.minimum_budget:
            if disp==0 and gamma==0:self.cutoff=H*math.exp(1-self.PJ/self.A)
            else:
                self.cutoff=brentq(lambda c:float(self.sink_integral(H,c))-self.PJ,
                                  np.nextafter(disp,H),H,xtol=2e-14*H)
    def sink_integral(self,r,c=None):
        c=self.cutoff if c is None else c;d=self.disp;r=np.asarray(r);g=-self.curvature
        inside=np.minimum(r,c)
        plateau=(inside+g*inside*inside/2)/(c-d)
        tail=(c+g*c*c/2)/(c-d)+(1+g*d)*(np.log(np.maximum(r,c)-d)-math.log(c-d))+g*(np.maximum(r,c)-c)
        return self.A*np.where(r<=c,plateau,tail)
    def k_extrema(self):
        c=self.cutoff;points=set(self.breaks());ac=self.a*(c-self.disp);disc=ac*ac-4*ac*self.nu
        if disc>=0:
            for D in ((ac+math.sqrt(disc))/2,(ac-math.sqrt(disc))/2):
                r=(D-self.d)/self.a
                if max(0.,self.zv)<r<c:points.add(r)
        b=self.d+self.a*self.disp
        if b+self.nu!=0:
            r=(b*self.nu/(b+self.nu)-self.d)/self.a
            if max(c,self.zv)<r<self.H:points.add(r)
        points=sorted(points);roots=[]
        for lo,hi in zip(points[:-1],points[1:]):
            if float(self.Q(lo))*float(self.Q(hi))<0:
                roots.append(brentq(lambda s:float(self.Q(s)),lo,hi,xtol=2e-14*self.H))
        places=[0.,self.H,*roots];values=[self.k(x)/self.K for x in places];i=int(np.argmin(values))
        return dict(minimum_over_K=values[i],location=places[i],stationary_Q_roots=roots,candidates=places,values_over_K=values)
    def result(self):
        r=dict(policy='SamePhysicalWallDisplacedAcceptedCandidateV1',model=self.model,H=self.H,h_s=self.hr,d=self.disp,
            gamma=-self.curvature,J0=self.J0,J2=self.J2,P0=self.PJ,A=self.A,minimum_required_budget=self.minimum_budget,
            cutoff_margin=self.PJ/self.minimum_budget-1.,cutoff=self.cutoff,original_roughness_interval=0<=self.hr<self.H,
            matching_velocity=self.t0*self.J0,metric_minimum=min(1.,1-self.curvature*self.H))
        if self.cutoff is None:return dict(r,status='cutoff_budget_failure')
        ex=self.k_extrema();cov=[]
        if self.curvature==0:
            fun=lambda x:2*self.k(x)/3-self.t0*float(self.nt(x)/self.D(x))
            cuts=sorted(set([*self.breaks(),*ex['candidates']]))
            cov=[(fun(x)/self.K,x) for x in cuts]
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                v=minimize_scalar(fun,bounds=(lo,hi),method='bounded',options={'xatol':1e-11*self.H})
                cov.append((float(v.fun)/self.K,float(v.x)))
        budget=self.integrate(lambda z:self.J(z)*self.A/(np.maximum(z,self.cutoff)-self.disp))-self.PJ
        r.update(reference_k=ex,epsilon_H=self.A/(self.H-self.disp),epsilon_plateau=self.A/(self.cutoff-self.disp),
            budget_defect=budget,Q_ground=float(self.Q(0.)),Q_H=float(self.Q(self.H)),
            ground_k_gradient=float(self.Q(0.)/self.D(0.)),interface_k_gradient=float(self.Q(self.H)/(self.J(self.H)*self.D(self.H))),
            counterfactual_flat_mechanical_work_defect=self.t0**2*self.J0-self.nu*self.t0**2*self.J2-self.integrate(self.production))
        if cov:
            r['flat_covariance_minimum_over_K'],r['flat_covariance_location']=min(cov)
            G=np.zeros((3,3));G[0,2]=self.um/(self.p['kappa']*(self.H-self.disp))
            law=oracle(self.model,self.K,self.A/(self.H-self.disp),G)
            r.update(outer_covariance_minimum=float(np.linalg.eigvalsh(law['R'])[0]),Q=law['Q'].tolist(),Pi_outer=-float(law['Q'][2,2]),
                outer_shear=-float(law['R'][0,2]),outer_normal_traction=-float(law['R'][2,2])+2*self.K/3+float(law['Q'][2,2]))
        else:r['curved_covariance']='undefined without a declared actual curved velocity/strain field'
        r['status']='negative_reference_k' if ex['minimum_over_K']<=0 else 'covariance_failure' if cov and min(cov)[0]<0 else 'admissible_frozen_column'
        return r


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();m=RequestManifest(a.output,dict(schema='CT7DisplacedReferenceStudyV1',contract=sha(a.contract),dataset={n:sha(a.data/n) for n in ('snapshot.json','profiles.csv','extraction.json','source_manifest.json')},source=sha(Path(__file__))))
    rows,z,v,dz,dv=observations(a.data);normal=float(np.interp(.04,z,v));nu=1/10500;disp=.0035/.04
    raw=list(csv.DictReader((a.data/'profiles.csv').open()))
    shear=[dict(point_id=r['point_id'],physical_height_m=.25*float(r['axis_y']),stress_over_utau2=float(r['value'])/.051**2)
           for r in raw if r['series_id']=='f3b_R' and r['quality']=='visible' and .009<=.25*float(r['axis_y'])<=.031]
    # Membership is frozen before the first call to wall_root, not selected by
    # whether a candidate returns positive k or an admissible h_s.
    require([int(r['point_id']) for r in rows[:4]]==[0,1,2,3],'changed frozen fit membership')
    result=dict(schema='CT7DisplacedReferenceStudyV1',fits={},cases={},extraction={},shear_support=shear,
        normalization=dict(U_Hh_over_Uinf=normal,ustar=.051/normal,nu=nu,policy='fixed original normalization, u_tau and nu in every extraction scenario; no fitted-U(Hh) renormalization'),
        source_roles_unchanged=True,primary_fit='equal weights f3a_R IDs0--3; CT7 proposed interval, author exact fit membership unavailable')
    scenarios=np.vstack([np.array(list(itertools.product((-1.,1.),repeat=12))),2*qmc.Sobol(12,scramble=True,seed=4).random_base2(8)-1,np.zeros((1,12))])
    for model in MODELS:
        f=fit(model,z,v,normal);old=fit(model,z,v,normal,0.);f['limiting_H_mm']=height_limit(model,f,nu,disp)
        f['archived_d0_control']=old
        f['displacement_plane_limiting_physical_H_mm']=height_limit(model,f,nu,0.)+3.5
        f['point_observations']=[dict(id=int(r['point_id']),z_m=float(zz),used=i<4,observed_U_over_UH=float(vv/normal)) for i,(r,zz,vv) in enumerate(zip(rows,z,v,strict=True))]
        result['fits'][model]=f
        f['duplicate_mean_checks']=[dict(series=r['series_id'],point_id=r['point_id'],physical_height_m=.0035+float(r['axis_x'])/1000,
            observed_U_over_utau=float(r['value']),predicted_U_over_utau=math.log((float(r['axis_x'])/1000)/f['z0_m'])/parameters(model)['kappa'],
            role=r['role'],coordinate='Fig3e x axis is physical y-d in mm; y axis U/u_tau, not y/z0') for r in raw if r['series_id']=='f3e_R' and r['quality']=='visible']
        f['full_height_outer_reference']=[]
        for z_m in (.007,.04,.25,.8):
            zeta=(z_m-.0035)/.04;p=parameters(model);A=f['ustar']**3/p['kappa'];G=np.zeros((3,3));G[0,2]=f['ustar']/(p['kappa']*zeta)
            law=oracle(model,f['K'],A/zeta,G)
            f['full_height_outer_reference'].append(dict(physical_z_m=z_m,zeta_over_Hh=zeta,
                U=f['ustar']/p['kappa']*math.log((z_m-.0035)/f['z0_m']),K=f['K'],epsilon=A/zeta,
                nut=law['nut'],Pi=-float(law['Q'][2,2]),R=law['R'].tolist(),
                k_source=law['production']-A/zeta,epsilon_balance=law['eps_source']+f['ustar']**4/(p['sigma_epsilon']*zeta**2),
                prescribed_top_diffusive_epsilon=-f['ustar']**4/(p['sigma_epsilon']*zeta),
                no_cap=True,scope='outer idealization only; complete inner compatibility is a separate gate'))
        for mm in HEIGHTS:
            H=mm/40;p=parameters(model)
            root=wall_root(model,H,f,nu,disp);row=dict(model=model,H_mm=mm,roughness_match=root,
                H_plus=H*f['ustar']/nu,zeta_H_plus=(H-disp)*f['ustar']/nu,cube_clearance_mm=mm-5,
                outer_nu_over_nut_H=nu/(p['kappa']*f['ustar']*(H-disp)),
                curvature_height_ratio=math.pi**2/12.5*H,historical_thin_point1=math.pi**2/12.5*H<=.1,
                below_nominal_Fig12_R_height=mm<7.,below_flat_stress_extraction_lower=mm<9.121)
            row['alternative_displacement_plane_velocity_match']=wall_root(model,H-disp,f,nu,0.)
            if root['status']=='root':
                coarse=DisplacedColumn(model,f['K'],mm/40,root['h_s'],nu,disp).result()
                fine=DisplacedColumn(model,f['K'],mm/40,root['h_s'],nu,disp,96).result()
                row.update(candidate=fine,quadrature_uncertainty={k:abs(fine[k]-coarse[k]) for k in ('J0','P0')})
                row['velocity_match_defect']=fine['matching_velocity']-root['target']
            result['cases'][model+'-'+str(mm)]=row
        vals=[];brackets=set();u40=[];z0=[];failed={mm:0 for mm in HEIGHTS}
        for s in scenarios:
            zz=z.copy();vv=v.copy();zz[:6]+=s[:6]*dz[:6];vv[:6]+=s[6:]*dv[:6]
            sf=fit(model,zz,vv,normal);vals.append(height_limit(model,sf,nu,disp));z0.append(sf['z0_m'])
            j=int(np.searchsorted(zz,.04));brackets.add((j-1,j));u40.append(float(np.interp(.04,zz,vv)/normal))
            for mm in HEIGHTS:failed[mm]+=wall_root(model,mm/40,sf,nu,disp)['status']!='root'
        result['extraction'][model]=dict(scenarios=len(scenarios),z0_m_range=[min(z0),max(z0)],limiting_H_mm_range=[min(vals),max(vals)],
            failed_velocity_roots={str(k):int(v) for k,v in failed.items()},observed_U40_over_original_scale_range=[min(u40),max(u40)],normalization_brackets=sorted(brackets),
            scope='joint deterministic extraction evaluations; fixed membership and units; correlated fit/40mm check includes shared point3; not statistical or a certified box extremum')
    regression=[]
    for model in MODELS:
        K=parameters(model)['q']*.16
        for gamma in (-.25,0.,.25):
            old=Column(model,K,.5,.05,1.5e-5,curvature=-gamma,order=96)
            new=DisplacedColumn(model,K,.5,.05,1.5e-5,0.,96,gamma)
            defect=abs(old.cutoff-new.cutoff)
            require(defect<2e-11,'d0 metric cutoff regression')
            shifted=DisplacedColumn(model,K,.5,.05,1.5e-5,.05,96,gamma).result()
            regression.append(dict(model=model,gamma=gamma,d0_cutoff_defect=defect,shifted_budget=shifted))
    result['candidate_metric_controls']=regression
    result['decision']=dict(execution='pass',physical_family='not_admitted_pending_complete_fit_wall_height_and_domain_capacity',native_displaced_policy='not_implemented_or_admitted; independent proposed specification extension only',
        smooth='CT4/5 smooth preparation limitations unchanged',modeled_displacement_plane='alternative changes physical wall/domain/H interpretation; not silently used')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT7 displaced fit and candidate column study complete',flush=True)
    return 0

if __name__=='__main__':raise SystemExit(main())
