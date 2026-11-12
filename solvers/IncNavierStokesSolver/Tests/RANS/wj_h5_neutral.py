#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Distinct exact neutral WJ/Budget reference and mean-only calibration.

Lengths use H_h=40mm unless an explicit *_mm field says otherwise. The fixed
experimental viscosity uses reported z0=.2mm, never the fitted intercept.
"""
import csv
from dataclasses import dataclass,field
import hashlib
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from wj_h2_profiles import wall

NE=81/20;QE=81/(2*math.sqrt(143));CD=572/6561
KAPPA=math.sqrt(1.30*(1.82-1.44)*math.sqrt(CD))
SCALE_MM=40.;DEFAULT_US=.051;DEFAULT_NU=.051*.005/3.97
DATA=Path(__file__).with_name('cao-tamura-fig3')/'rough-fig3.csv'
DIGITIZATION=DATA.with_name('digitization.json')


def positive(**values):
    if any(isinstance(v,bool)or not math.isfinite(v)or v<=0 for v in values.values()):
        raise ValueError('positive finite parameters required: '+','.join(values))


def measurements(d0_mm=3.5,delta_mm=250.):
    positive(d0_mm=d0_mm,delta_mm=delta_mm)
    import json
    resolution=json.loads(DIGITIZATION.read_text())['resolution'];result=[]
    for row in csv.DictReader(DATA.open()):
        panel=row['panel'];raw=float(row['raw_height'])
        z=raw if panel=='a' else raw*delta_mm if panel=='b' else d0_mm+raw*delta_mm
        result.append(dict(panel=panel,point=int(row['point']),z_mm=z,r_mm=z-d0_mm,value=float(row['value']),
            height_origin=row['height_origin'],reading_height_mm=resolution[panel]['height_mm_reading_resolution'],
            reading_value=resolution[panel]['value_reading_resolution'],page=7,figure='3('+panel+')',
            kind='digitized hollow rough-case symbol'))
    return result


def mean_fit(*,u_tau=DEFAULT_US,d0_mm=3.5,delta_mm=250.,window_mm=25.,omit_lowest=False):
    positive(u_tau=u_tau,window_mm=window_mm)
    points=[r for r in measurements(d0_mm,delta_mm)if r['panel']=='a'and r['z_mm']<=window_mm]
    if omit_lowest:points=points[1:]
    if not points:raise ValueError('empty outer mean-fit window')
    if any(r['r_mm']<=r['reading_height_mm']for r in points):raise ValueError('reading window crosses effective origin')
    logz=sum(math.log(r['r_mm']/SCALE_MM)-KAPPA*r['value']/u_tau for r in points)/len(points)
    z0=math.exp(logz)
    def bound(height_sign,value_sign):
        return SCALE_MM*math.exp(sum(math.log((r['r_mm']+height_sign*r['reading_height_mm'])/SCALE_MM)-
            KAPPA*(r['value']+value_sign*r['reading_value'])/u_tau for r in points)/len(points))
    residuals=[]
    for r in points:
        predicted=u_tau/KAPPA*math.log((r['r_mm']/SCALE_MM)/z0)
        residuals.append(dict(r,predicted=predicted,residual=predicted-r['value']))
    return dict(schema='Fig3FixedSlopeMeanV1',u_tau=u_tau,d0_mm=d0_mm,delta_mm=delta_mm,window_mm=window_mm,
        omit_lowest=omit_lowest,z0_e=z0,z0_e_mm=SCALE_MM*z0,points=residuals,
        rms=math.sqrt(sum(r['residual']**2 for r in residuals)/len(residuals)),
        reading_bounds_mm=dict(height_only=[bound(-1,0),bound(1,0)],value_only=[bound(0,1),bound(0,-1)],
            combined=[bound(-1,1),bound(1,-1)]),
        reading_scope='fixed center-selected membership; bounded reading effects, not statistical confidence; d0/delta/u_tau uncertainty unquantified',
        data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest())


def correlation_b(x):
    positive(x=x)
    m=max(.5-.4*x**.7,1-.79*x**(-.28))
    return -10.7*math.expm1(m*math.log(x))


def required_roughness(H,z0_e,*,u_tau=DEFAULT_US,nu=DEFAULT_NU):
    """Required inverse, also returned when h>=H to expose infeasibility.

    This is not an admissible model constructor. The reference class below
    enforces h<H and all cutoff/reference conditions independently.
    """
    positive(H=H,z0_e=z0_e,u_tau=u_tau,nu=nu)
    if H<=z0_e:raise ValueError('intercept must lie below matching height')
    a=KAPPA*u_tau;d=a*H*z0_e/(H-z0_e)
    if d<=nu:raise ValueError('negative_intercept_branch_required')
    alpha=KAPPA/math.sqrt(QE);b=(1-d/nu)/alpha
    hi=2.
    while correlation_b(hi)>b:
        hi*=2
        if hi>1e12:raise ValueError('roughness_inverse_absent')
    x=brentq(lambda t:correlation_b(t)-b,1.,hi,xtol=5e-15,rtol=1e-14)
    h=x*70*nu/math.sqrt(QE*u_tau*u_tau)
    return dict(h=h,d_nu=d,zv=(nu-d)/a,b=b,x=x,a=a)


@dataclass(frozen=True)
class NeutralReference:
    H:float
    z0_e:float
    u_tau:float=DEFAULT_US
    nu:float=DEFAULT_NU
    d0:float=.0875
    top:float=20.
    fit_identity:str='primary fixed-slope Fig.3 mean fit'
    w:dict=field(init=False,repr=False,compare=False)
    h:float=field(init=False)
    d_nu:float=field(init=False)

    def __post_init__(self):
        positive(H=self.H,z0_e=self.z0_e,u_tau=self.u_tau,nu=self.nu,top=self.top)
        if not math.isfinite(self.d0)or self.d0<0 or self.d0+self.H>=self.top:raise ValueError('invalid physical geometry')
        v=required_roughness(self.H,self.z0_e,u_tau=self.u_tau,nu=self.nu)
        if v['h']>=self.H:raise ValueError('h_s_not_below_matching_height')
        object.__setattr__(self,'h',v['h']);object.__setattr__(self,'d_nu',v['d_nu'])
        w=wall(self.K,self.H,self.h,nu=self.nu)
        if w['zv']>=0:raise ValueError('negative_intercept_branch_required')
        if w['Lambda']<1 or not 0<w['rc']<=self.H:raise ValueError('budget_cutoff_inadmissible')
        object.__setattr__(self,'w',w)

    @property
    def K(self):return QE*self.u_tau**2
    @property
    def A(self):return self.u_tau**3/KAPPA
    @property
    def a(self):return KAPPA*self.u_tau
    @property
    def Rc(self):
        return self.K*np.array([[2/3+1/NE,0,-1/QE],[0,2/3,0],[-1/QE,0,2/3-1/NE]])

    def D(self,r):return self.nu+self.a*max(r-self.w['zv'],0)
    def nt(self,r):return self.a*max(r-self.w['zv'],0)
    def U_inner(self,r):return self.u_tau**2*float(self.w['J'](r))
    def U_outer(self,r):return self.u_tau/KAPPA*math.log(r/self.z0_e)
    def epsilon(self,r):return self.A/r
    def k_inner(self,r):return self.w['kref'](float(r))
    def production(self,r):return self.u_tau**4*self.nt(r)/self.D(r)**2
    def destruction(self,r):return self.A/max(r,self.w['rc'])
    def covariance_inner(self,r):
        R=2*self.k_inner(r)/3*np.eye(3);R[0,2]=R[2,0]=-self.nt(r)*self.u_tau**2/self.D(r);return R

    def summary(self):
        w=self.w;rc=w['rc'];P0=self.production(0);plateau=self.A/rc
        peak=(2*self.nu-self.d_nu)/self.a
        # Positive tail defect and a unimodal plateau production with P(0)>E
        # extend the sufficient d>=2nu proposition without assuming it applies.
        conditions=dict(negative_intercept=w['zv']<0,cutoff_gate=w['Lambda']>=1,
            production_unimodal=self.d_nu>self.nu,initial_negative_defect=P0>plateau,
            peak_before_cutoff=peak<rc,positive_tail=(self.d_nu>0 and self.nu>0),
            zero_total_budget=abs(w['Q'](self.H))<=2e-12*self.A)
        return dict(family='WJNeutralBudgetReferenceV1',H_w=self.H,H_w_mm=SCALE_MM*self.H,
            physical_z_H_mm=SCALE_MM*(self.d0+self.H),h_s=self.h,h_s_mm=SCALE_MM*self.h,z0_e=self.z0_e,
            z0_e_mm=SCALE_MM*self.z0_e,u_tau=self.u_tau,nu=self.nu,d0=self.d0,top=self.top,K=self.K,A=self.A,
            Ne=NE,qe=QE,Cd=CD,kappa=KAPPA,c_mu_wall=CD,c_length_wall=KAPPA*QE**1.5,alpha=KAPPA/math.sqrt(QE),
            a=self.a,d_nu=self.d_nu,d_over_nu=self.d_nu/self.nu,zv=w['zv'],Lambda=w['Lambda'],rc=rc,rc_mm=SCALE_MM*rc,
            k_ground_over_K=self.k_inner(0)/self.K,reference_min_k=self.K,
            reference_covariance_min_over_utau2=float(np.linalg.eigvalsh(self.covariance_inner(self.H)).min()/self.u_tau**2),
            Q_endpoints=[w['Q'](0),w['Q'](self.H)],positivity_extension=conditions,
            positivity_scope='unimodal P, P(0)>plateau, positive tail and zero budget imply one crossing, Q<0, k>K; stationary reference only',
            H_plus=self.H*self.u_tau/self.nu,outer_molecular_turbulent_ratio=self.nu/(self.a*self.H),
            inner_molecular_fraction=self.nu/self.D(self.H),h_over_H=self.h/self.H,
            physical_plane_over_cube_height=SCALE_MM*(self.d0+self.H)/5,
            gradient_ratio=self.u_tau**2/self.D(self.H)/(self.u_tau/KAPPA/self.H),
            outer_R_over_utau2=(self.Rc/self.u_tau**2).tolist(),
            full_top_reference_traction=[self.u_tau**2,0,self.K/NE],top_normal_role='impermeability reaction, not prescribed',
            finite_depth_mismatch={str(z):self.U_outer(z/SCALE_MM-self.d0)for z in (250,800)},
            fit_identity=self.fit_identity,missing_inputs='no measured-rms/spanwise k completion; no Blackadar/finite-depth tail')


def feasibility(fit):
    z0=fit['z0_e'];us=fit['u_tau'];d0=fit['d0_mm'];rows=measurements(d0,fit['delta_mm'])
    minima={p:min((r for r in rows if r['panel']==p),key=lambda r:r['z_mm'])for p in 'abcd'}
    ceiling=min(minima[p]['r_mm']-minima[p]['reading_height_mm']for p in 'bcd')
    def threshold(z):
        return SCALE_MM*brentq(lambda H:required_roughness(H,z,u_tau=us)['h']-H,
            z*1.01,1.,xtol=1e-14,rtol=1e-14)
    primary=threshold(z0);printed=threshold(.2/SCALE_MM)
    lower3=SCALE_MM*DEFAULT_NU/(KAPPA*us*.03);lower2=SCALE_MM*DEFAULT_NU/(KAPPA*us*.02)
    return dict(minima=minima,fig3_conservative_ceiling_mm=ceiling,fig12=dict(page=15,figure='12',height_mm=7.,height_source_page=13,
        r_mm=7-d0,requirement='H_w<7mm-d0 if all rms samples must be outside AWF',digitized_values='not extracted in H5'),
        primary_h_equals_H_mm=primary,printed_h_equals_H_mm=printed,
        primary_fig3_strict_interval_mm=[primary,ceiling],printed_fig3_strict_interval_mm=[printed,ceiling],
        three_percent_lower_mm=lower3,two_percent_lower_mm=lower2,
        primary_three_percent_interval_mm=[max(primary,lower3),ceiling],primary_two_percent_interval='EMPTY'if lower2>=ceiling else [max(primary,lower2),ceiling],
        all_7mm_outside=dict(primary='EMPTY',printed='EMPTY',required_upper_mm=7-d0,
            primary_required_h_at_3p5mm=SCALE_MM*required_roughness((7-d0)/SCALE_MM,z0,u_tau=us)['h']),
        qualification='Fig.3-only reference candidates; no whole-experiment matching-height selection')


def observation_regions(config):
    rows=measurements(SCALE_MM*config.d0)
    rows.append(dict(panel='Fig12',point='all downstream rms heights',z_mm=7.,r_mm=7-SCALE_MM*config.d0,value=None,
        height_origin='physical floor; 2mm above5mm cubes',reading_height_mm=None,reading_value=None,page=15,figure='12',kind='reported height; rms values not digitized'))
    return [dict(r,H_w_mm=SCALE_MM*config.H,region='outer WJ'if r['r_mm']>SCALE_MM*config.H else 'inner AWF',
        prediction_role='full3D WJ stress'if r['r_mm']>SCALE_MM*config.H else 'local k_w/full-strain AWF stress; not outer WJ extrapolation')for r in rows]
