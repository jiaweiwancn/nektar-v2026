#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT-Accepted-1: independent flat/metric column oracle, never a hill solver."""
import argparse,csv,json,math,os,subprocess
from functools import lru_cache
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from numpy.polynomial.legendre import leggauss,Legendre
from scipy.optimize import brentq,minimize_scalar
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import write_json


def require(condition,reason):
    if not condition:raise ValueError(reason)


def parameters(model):
    # Defining simple-shear roots, independent of the production descriptor.
    q=(57+math.sqrt(6249))/40 if model=='A93' else (3+math.sqrt(61 if model=='A95' else 41))/(2*math.sqrt(2))
    delta=.48 if model!='E94E95' else 1.9-.43*q
    sigma=1.3 if model!='E94E95' else 1.2
    kap=math.sqrt(sigma*delta/q)
    return dict(q=q,kappa=kap,sigma_epsilon=sigma,c_mu_wall=q**-2,
                c_l_wall=kap*q**1.5,alpha=kap/math.sqrt(q))


@lru_cache(None)
def gauss(order):return leggauss(order)


def incident(data,surface,height,logarithmic=False):
    rows=[r for r in csv.DictReader((data/'profiles.csv').open())
          if r['series_id']=='f3a_'+surface and float(r['weight'])>0]
    rows.sort(key=lambda r:float(r['y_over_H']))
    z=np.array([.04*float(r['y_over_H']) for r in rows]);v=np.array([float(r['value']) for r in rows])
    bound=np.array([float(r['value_bound']) for r in rows]);d=.0035 if surface=='R' else 0.
    require(height>d,'observation_height_not_above_assumed_displacement')
    x=np.log(z-d) if logarithmic else z;t=math.log(height-d) if logarithmic else height
    index=max(0,min(len(x)-2,int(np.searchsorted(x,t))-1));a=(t-x[index])/(x[index+1]-x[index])
    value=(1-a)*v[index]+a*v[index+1];envelope=abs(1-a)*bound[index]+abs(a)*bound[index+1]
    outside=height<z[0] or height>z[-1]
    method_difference=0.
    if outside:
        ids=list(range(3)) if height<z[0] else list(range(len(z)-3,len(z)))
        third=float(np.polyval(np.polyfit(x[ids],v[ids],1),t));method_difference=abs(third-value)
    return dict(value=float(value),extraction_envelope=float(envelope),
        two_vs_three_point_method_difference=method_difference,extrapolated=outside,
        visible_height_range_m=[float(z[0]),float(z[-1])],height_m=height,
        point_ids=[rows[index]['point_id'],rows[index+1]['point_id']],
        method='linear in log(y-d)' if logarithmic else 'linear in physical height',
        uncertainty='extraction envelope and extrapolation-method difference separate; displacement/facility/measurement uncertainty remains unquantified')


class Column:
    """All lengths nondimensionalized by Hh, U(Hh)=1; J=1-curvature*r."""
    def __init__(self,model,K,H,hr,nu,curvature=0.,order=48):
        self.model=model;self.p=parameters(model);self.K=K;self.H=H;self.hr=hr;self.nu=nu;self.curvature=curvature;self.order=order
        require(K>0 and H>0 and nu>0 and hr>=0,'invalid_column_scalar')
        self.um=math.sqrt(K/self.p['q']);self.t0=K/self.p['q'];self.a=self.p['kappa']*self.um;self.A=self.um**3/self.p['kappa']
        x=hr*math.sqrt(K)/(70*nu)
        b=10.7 if x==0 else 10.7*(1-x**max(.5-.4*x**.7,1-.79*x**(-.28)))
        self.zv=nu*b/math.sqrt(K);self.d=nu-self.a*self.zv
        self.cutoff=None
        self.J0=self.integrate(lambda r:1/self.D(r))
        self.J2=self.integrate(lambda r:1/self.D(r)**2)
        self.PJ=self.integrate(lambda r:self.J(r)*self.production(r))
        self.minimum_budget=self.A*(1-curvature*H/2)
        self.Lambda=self.PJ/self.A
        if self.PJ>=self.minimum_budget and min(1.,1-curvature*H)>0:
            if curvature==0:self.cutoff=H*math.exp(1-self.Lambda)
            else:self.cutoff=math.exp(brentq(lambda logc:self.sink_integral(H,math.exp(logc))-self.PJ,
                math.log(H)-max(64.,2*self.Lambda),math.log(H),xtol=5e-14))

    def J(self,r):return 1-self.curvature*np.asarray(r)
    def nt(self,r):return self.a*np.maximum(np.asarray(r)-self.zv,0.)
    def D(self,r):return self.nu+self.nt(r)
    def production(self,r):return self.t0**2*self.nt(r)/self.D(r)**2
    def breaks(self,lower=0.,upper=None):
        upper=self.H if upper is None else upper
        points={lower,upper}
        for value in (self.zv,self.hr,self.cutoff):
            if value is not None and lower<value<upper:points.add(value)
        # Resolve the viscosity length at its actual positive-viscosity start.
        start=max(0.,self.zv);length=self.D(start)/self.a
        for j in range(80):
            value=start+length*(1.5**j-1)
            if value>=upper:break
            if value>lower:points.add(float(value))
        return sorted(points)
    def integrate(self,function,lower=0.,upper=None):
        x,w=gauss(self.order);cuts=self.breaks(lower,upper);total=0.
        for a,b in zip(cuts[:-1],cuts[1:]):
            r=(a+b)/2+(b-a)*x/2;total+=(b-a)/2*float(np.dot(w,function(r)))
        return total
    def production_integral(self,r):
        # Independent closed primitive, checked against split numerical quadrature.
        lo=max(0.,self.zv);r=np.asarray(r);active=np.maximum(r,lo)
        dl=self.D(lo);D=self.D(active);change=D-dl;log=np.log1p(change/dl);inverse=1/D-1/dl
        i0=self.t0**2/self.a*(log+self.nu*inverse)
        i1=self.t0**2/self.a**2*(change-(self.d+self.nu)*log-self.nu*self.d*inverse)
        return np.where(r<=lo,0.,i0-self.curvature*i1)
    def sink_integral(self,r,c=None):
        c=self.cutoff if c is None else c;r=np.asarray(r)
        f0=np.where(r<=c,r/c,1+np.log(np.maximum(r,c)/c))
        f1=np.where(r<=c,r*r/(2*c),r-c/2)
        return self.A*(f0-self.curvature*f1)
    def Q(self,r):return self.sink_integral(r)-self.production_integral(r)
    def k(self,r):return self.K-self.integrate(lambda s:self.Q(s)/(self.J(s)*self.D(s)),float(r),self.H)
    def k_extrema(self):
        # Q'=J(E-p0). On each of these intervals Q is monotone, so every
        # interior Q zero is bracketed; a uniform sample is not used.
        c=self.cutoff;points=set(self.breaks());ac=self.a*c;disc=ac*ac-4*ac*self.nu
        if disc>=0:
            for D in ((ac+math.sqrt(disc))/2,(ac-math.sqrt(disc))/2):
                r=(D-self.d)/self.a
                if max(0.,self.zv)<r<c:points.add(r)
        if self.d+self.nu!=0:
            r=-self.d**2/(self.a*(self.d+self.nu))
            if max(c,self.zv)<r<self.H:points.add(r)
        points=sorted(points);roots=[]
        for a,b in zip(points[:-1],points[1:]):
            qa,qb=float(self.Q(a)),float(self.Q(b))
            if qa*qb<0:roots.append(brentq(lambda r:float(self.Q(r)),a,b,xtol=2e-14*self.H))
        candidates=[0.,self.H,*roots];values=[self.k(r) for r in candidates]
        i=int(np.argmin(values));return dict(minimum_over_K=values[i]/self.K,location=candidates[i],
            stationary_Q_roots=roots,candidates=candidates,values_over_K=[k/self.K for k in values])
    def result(self):
        geom=self.hr<self.H and min(1.,1-self.curvature*self.H)>=.2
        result=dict(model=self.model,coefficients=self.p,K=self.K,H=self.H,hr=self.hr,nu=self.nu,
            counterfactual=dict(u=self.um,traction=self.t0,a=self.a,epsilon_amplitude=self.A),
            J0=self.J0,J2=self.J2,finite_nu_production=self.PJ,Lambda=self.Lambda,
            cutoff_existence_margin=(self.PJ-self.minimum_budget)/self.minimum_budget,
            cutoff=self.cutoff,zv=self.zv,d_over_nu=self.d/self.nu,
            unshifted_compatible_z0=self.H*math.exp(-self.a*self.J0),
            matching_velocity=self.t0*self.J0,total_reference_traction=self.t0,
            geometry_valid=geom,curvature=self.curvature,minimum_stretch=min(1.,1-self.curvature*self.H))
        if self.cutoff is None:result.update(classification='cutoff_failure' if geom else 'geometric_exclusion');return result
        extrema=self.k_extrema();minimum=extrema['minimum_over_K']
        moments=[]
        for degree in range(7):
            phi=Legendre.basis(degree);der=phi.deriv();source=self.integrate(lambda r:(self.A/np.maximum(r,self.cutoff)-self.production(r))*self.J(r)*phi(2*r/self.H-1))
            diffusion=self.integrate(lambda r:self.Q(r)*2/self.H*der(2*r/self.H-1))
            moments.append(dict(degree=degree,source=source,diffusion=diffusion,defect=source+diffusion))
        result.update(cutoff_over_H=self.cutoff/self.H,dissipation_plateau=self.A/self.cutoff,
            reference_k=extrema,Q_end=float(self.Q(self.H)),weak_moments=moments,
            budget_defect=self.integrate(lambda r:self.J(r)*self.A/np.maximum(r,self.cutoff))-self.PJ,
            production_primitive_defect=abs(float(self.production_integral(self.H))-self.PJ))
        if self.curvature==0:
            locations=sorted(set([*extrema['candidates'],0.,self.H,self.cutoff,
                *([self.zv] if 0<self.zv<self.H else [])]))
            covariance=lambda r:2*self.k(r)/3-self.t0*float(self.nt(r)/self.D(r))
            minima=[(covariance(r),r) for r in locations]
            for a,b in zip(locations[:-1],locations[1:]):
                opt=minimize_scalar(covariance,bounds=(a,b),method='bounded',options={'xatol':1e-11*self.H})
                minima.append((float(opt.fun),float(opt.x)))
            margin,point=min(minima)
            result.update(flat_covariance_margin_over_K=margin/self.K,flat_covariance_minimum_location=point,
                flat_covariance_sufficient_lower_bound_over_K=2*minimum/3-self.t0/self.K*float(self.nt(self.H)/self.D(self.H)),
                mechanical_work_defect=self.t0*result['matching_velocity']-self.nu*self.t0**2*self.J2-self.PJ)
        result['classification']='geometric_exclusion' if not geom else 'negative_k' if minimum<0 else 'covariance_failure' if result.get('flat_covariance_margin_over_K',1)<0 else 'admissible_reference'
        return result


def fit_roughness(model,nu,us,H,target,kind):
    q=parameters(model)['q'];K=q*us*us
    def residual(hr):
        c=Column(model,K,H,hr,nu)
        return c.t0*c.J0-target if kind=='velocity' else c.H*math.exp(-c.a*c.J0)-target
    a,b=0.,.0144/.04;fa,fb=residual(a),residual(b)
    if fa*fb>0:return dict(status='unattainable_roughness_target',bracket=[a,b],residuals=[fa,fb],target=target)
    hr=brentq(residual,a,b,xtol=5e-14)
    return dict(status='fit',hr=hr,target=target,residual=residual(hr),condition='16mm, K=q_m*utau^2; held fixed thereafter',kind=kind)


def native_compare(probe,c,row):
    if probe is None:return dict(status='not_attempted')
    points=sorted(set([0.,c.H]+([] if c.cutoff is None else [c.cutoff,*row['reference_k']['candidates']])))
    values=[c.K,c.H,c.hr,c.H,c.H,1.,96,c.t0*c.J0,0.,0.,0.,len(points),*points]
    probe.stdin.write('wall '+c.model+' Natural Accepted '+str(c.nu)+' '+' '.join(map(str,values))+'\n');probe.stdin.flush()
    native=json.loads(probe.stdout.readline())
    if c.hr>=c.H:
        require(not native['ok'] and '0<=h<H' in native['error'],'native_geometric_rejection_reason')
        return dict(status='pass',expected='geometric_rejection',native=native)
    if c.cutoff is None:
        require(not native['ok'] and 'cutoff_budget_inadmissible' in native['error'],'native_cutoff_rejection_reason')
        return dict(status='pass',expected='cutoff_rejection',native=native)
    require(native['ok'],'native_forward_column_failure')
    defects=dict(J0=abs(native['J'][0]-c.J0)/c.J0,J2=abs(native['J'][2]-c.J2)/c.J2,
        cutoff=abs(native['cutoff']-c.cutoff)/c.H,U=abs(native['UH']-c.t0*c.J0),
        production=abs(native['production_budget']-c.PJ)/max(c.PJ,1e-30))
    for point in native['samples']:
        z=point[0];defects['k']=max(defects.get('k',0.),abs(point[6]-c.k(z))/c.K)
        defects['Q']=max(defects.get('Q',0.),abs(point[7]-float(c.Q(z)))/c.A)
    require(max(defects.values())<=2e-9,'native_flat_column_disagrees_with_independent_oracle')
    return dict(status='pass',expected='forward mathematical profile; field admissibility classified separately',defects=defects,native=native)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--contract',type=Path,required=True);p.add_argument('--data',type=Path,required=True);p.add_argument('--native-probe',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    contract=json.loads(a.contract.read_text());require(contract['schema']=='CT-Accepted-1-v1','contract_version')
    for name,digest in contract['frozen_dataset']['hashes'].items():require(sha(a.data/name)==digest,'frozen_dataset_changed:'+name)
    intent=dict(schema='CTAcceptedIndependentColumnsV1',contract_sha256=sha(a.contract),source_sha256=sha(Path(__file__)),
        dataset_snapshot_sha256=sha(a.data/'snapshot.json'),native_probe=runtime_identity(a.native_probe) if a.native_probe else None)
    seal=RequestManifest(a.output,intent);process=subprocess.Popen([str(a.native_probe)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True) if a.native_probe else None
    incident_data={};fits={};flat=[];curved=[]
    try:
        for surface in ('S','R'):
            uh=incident(a.data,surface,.04);uinf=1/uh['value'];nu=1/contract['benchmark']['Re_H'][surface];us=contract['benchmark']['utau_over_Uinf'][surface]*uinf
            anchors={str(H):incident(a.data,surface,H,True) if H>(.0035 if surface=='R' else 0) else dict(classification='inadequate_observational_support',reason='below assumed observational displacement; physical wall unchanged') for H in contract['feasibility']['heights_m']}
            incident_data[surface]=dict(U_H_over_Uinf=uh,Uinf_over_U_H=uinf,nu=nu,ustar=us,anchors=anchors)
            for model in contract['physics']['models']:
                q=parameters(model)['q'];K=q*us*us;scenarios={'smooth':0.} if surface=='S' else {'common-1mm':.001/.04,'cube-proxy-5mm':.005/.04}
                if surface=='R':
                    target=anchors['0.016']['value']*uinf
                    for name,kind,value in [('model-velocity-fit','velocity',target),('reported-z0-unshifted-fit','z0',.0002/.04)]:
                        fit=fit_roughness(model,nu,us,.016/.04,value,kind);fits[model+':'+name]=fit
                        if fit['status']=='fit':scenarios[name]=fit['hr']
                        else:flat.append(dict(model=model,surface=surface,scenario=name,classification=fit['status'],fit=fit))
                for name,hr in scenarios.items():
                    for Hm in contract['feasibility']['heights_m']:
                        for factor in contract['feasibility']['endpoint_K_factors']:
                            c=Column(model,factor*K,Hm/.04,hr,nu);row=c.result();tight=Column(model,factor*K,Hm/.04,hr,nu,order=96).result()
                            uncertainty={key:abs(row[key]-tight[key]) for key in ('J0','J2','finite_nu_production','matching_velocity')}
                            if c.cutoff is not None:uncertainty['minimum_k_over_K']=abs(row['reference_k']['minimum_over_K']-tight['reference_k']['minimum_over_K'])
                            row.update(surface=surface,scenario=name,H_m=Hm,K_factor=factor,quadrature_uncertainty=uncertainty)
                            row['native']=native_compare(process,c,row);flat.append(row)
        for model in contract['physics']['models']:
            for surface in ('S','R'):
                ref=incident_data[surface];fit=fits.get(model+':model-velocity-fit')
                if surface=='R' and fit['status']!='fit':continue
                hr=0. if surface=='S' else fit['hr'];K=parameters(model)['q']*ref['ustar']**2
                for Hm in contract['feasibility']['heights_m']:
                    H=Hm/.04;stations=sorted(set(np.linspace(-2.5,2.5,33).tolist()+[-3.,-2.5-1e-8,-2.5+1e-8,-1.25,0.,1.25,2.5-1e-8,2.5+1e-8,3.]))
                    for x in stations:
                        slope=-math.pi/5*math.sin(math.pi*x/2.5) if abs(x)<2.5 else 0.
                        second=-math.pi**2/12.5*math.cos(math.pi*x/2.5) if abs(x)<2.5 else 0.
                        curvature=second/(1+slope*slope)**1.5
                        c=Column(model,K,H,hr,ref['nu'],curvature);row=c.result()
                        thin=abs(curvature)*H<=.1;reach=H<12.5/math.pi**2
                        row.update(surface=surface,H_m=Hm,x_over_Hh=x,thin_layer_eligible=thin,global_normal_reach_bound=12.5/math.pi**2,within_sufficient_global_reach=reach,
                            geometric_classification='admissible_geometry' if thin and reach and row['geometry_valid'] else 'geometric_exclusion',covariance='not_certified: frozen normal column has no declared complete terrain velocity/strain')
                        if curvature==0 and c.cutoff is not None:
                            f=Column(model,K,H,hr,ref['nu']).result();require(abs(row['cutoff']-f['cutoff'])<1e-13,'zero_curvature_cutoff_limit')
                        curved.append(row)
        # Deliberate molecular-only cutoff rejection, not a failed physical case.
        controls=[]
        for model in contract['physics']['models']:
            ref=incident_data['S'];K=parameters(model)['q']*ref['ustar']**2;H=.5*ref['nu']/math.sqrt(K)
            c=Column(model,K,H,0.,ref['nu']);row=c.result();require(c.cutoff is None,'molecular_only_control_not_rejected');row['native']=native_compare(process,c,row);controls.append(row)
        smooth=[r for r in flat if r.get('surface')=='S' and r.get('H_m')==.008 and r.get('K_factor')==1]
        require(len(smooth)==3 and all(r['Lambda']>1 and r['reference_k']['values_over_K'][0]<0 and r['classification']=='negative_k' for r in smooth),'source_scaled_smooth_counterexample_not_reproduced')
        result=dict(status='pass',execution_status='completed',verification_status='pass',incident=incident_data,fits=fits,flat=flat,curved=curved,negative_controls=controls,
            smooth_8mm_counterexample=[dict(model=r['model'],Lambda=r['Lambda'],ground_k_over_K=r['reference_k']['values_over_K'][0]) for r in smooth],
            decision='feasibility classifications are physical/reference outcomes, not numerical verification failures or experimental validation',global_reach_justification='For the compact C1,1 graph, h(y)-h(x)-hprime(x)*(y-x)<=C*(y-x)^2/2, C=pi^2/12.5 in Hh units. At r<1/C, every normal footpoint is the unique nearest ground point, so distinct normal rays cannot meet within the layer. This is a sufficient global bound, not an exact reach estimate.')
    except Exception as error:result=dict(status='fail',failure=repr(error),incident=incident_data,fits=fits,flat=flat,curved=curved)
    finally:
        if process is not None:process.stdin.close();process.wait()
    write_json(a.output/'results.json',result);seal.seal(a.output,['results.json']);print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
