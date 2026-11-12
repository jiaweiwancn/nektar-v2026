#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT4 independent neutral fit, joint extraction scenarios and reference oracle."""
import argparse,csv,itertools,json,math,os,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from scipy.optimize import brentq
from scipy.stats import qmc
from ct_accepted_feasibility import Column,parameters,native_compare,write_json,require
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha


def observations(data,surface='R'):
    rows=[r for r in csv.DictReader((data/'profiles.csv').open()) if r['series_id']=='f3a_'+surface and r['quality']=='visible']
    rows.sort(key=lambda r:int(r['point_id']))
    panel=next(p for p in json.loads((data/'extraction.json').read_text())['panels'] if p['id']=='f3a')
    (x0,v0),(x1,v1)=panel['axis_x']['anchors'];(y0,z0),(y1,z1)=panel['axis_y']['anchors']
    z=np.array([((float(r['py'])-y0)*(z1-z0)/(y1-y0)+z0)/1000 for r in rows])
    v=np.array([(float(r['px'])-x0)*(v1-v0)/(x1-x0)+v0 for r in rows])
    dz=np.full(len(z),.5*abs((z1-z0)/(y1-y0))/1000);dv=np.full(len(v),.5*abs((v1-v0)/(x1-x0)))
    return rows,z,v,dz,dv


def fit(model,z,v,ratio=.051):
    p=parameters(model);normal=float(np.interp(.04,z,v));us=ratio/normal
    B=float(np.mean(v[:4]-ratio/p['kappa']*np.log(z[:4]/.04))/normal)
    return dict(U_H_over_Uinf=normal,Uinf=1/normal,ustar=us,B=B,K0=p['q']*us*us,
        normalization_bracket=[int(i) for i in (max(0,np.searchsorted(z,.04)-1),min(len(z)-1,np.searchsorted(z,.04)))])


def resistance(model,K,H,hr,nu):
    p=parameters(model);us=math.sqrt(K/p['q']);a=p['kappa']*us;x=hr*math.sqrt(K)/(70*nu)
    b=10.7 if hr==0 else 10.7*(1-x**max(.5-.4*x**.7,1-.79*x**(-.28)))
    zv=nu*b/math.sqrt(K);start=max(0.,zv)
    if start>=H:return H/nu
    return start/nu+math.log1p(a*(H-start)/(nu+a*(start-zv)))/a


def roughness(model,H,f,nu,branches=False):
    kap=parameters(model)['kappa'];target=f['B']+f['ustar']/kap*math.log(H)
    fun=lambda hr:f['ustar']**2*resistance(model,f['K0'],H,hr,nu)-target
    end=np.nextafter(H,0.);a,b=fun(0.),fun(end)
    if a*b>0:return dict(status='unattainable_roughness_target',residual_bracket=[a,b],target=target)
    hr=brentq(fun,0.,end,xtol=2e-14)
    out=dict(status='fit',h_s=hr,target=target,residual=fun(hr),residual_bracket=[a,b])
    if branches:
        points=np.linspace(0.,end,129);values=[fun(x) for x in points]
        count=sum(a*b<0 for a,b in zip(values[:-1],values[1:]))
        require(count==1 and 0<=hr<H,'nonunique_or_invalid_roughness_branch')
        out['branch_crossings_on_128_panels']=count
    return out


def reference(c,z,outer=None):
    outer=z>c.H if outer is None else outer
    if outer:
        U=c.t0*c.J0+c.um/c.p['kappa']*math.log(z/c.H);du=c.um/(c.p['kappa']*z)
        G=np.zeros((3,3));G[0,2]=du;law=oracle(c.model,c.K,c.A/z,G)
        return np.array([U,c.K,c.A/z,-law['Q'][2,2],du,0.,-c.A/z**2,0.,c.p['kappa']*c.um*z])
    return np.array([c.t0*c.integrate(lambda s:1/c.D(s),upper=z) if z>0 else 0.,c.k(z),c.A/max(z,c.cutoff),0.,
        c.t0/c.D(z),float(c.Q(z)/c.D(z)),0. if z<c.cutoff else -c.A/z**2,float(c.Q(z)),float(c.D(z))])


def normal_distance(x,z):
    def h(t):return math.cos(math.pi*t/5)**2 if abs(t)<2.5 else 0.
    def hp(t):return -math.pi/5*math.sin(math.pi*t/2.5) if abs(t)<2.5 else 0.
    if z<h(x):return dict(status='inside_terrain',distance=None)
    f=lambda t:t-hp(t)*(z-h(t))-x
    ts=np.linspace(-2.5+1e-12,2.5-1e-12,129);roots=[]
    for a,b in zip(ts[:-1],ts[1:]):
        if f(a)*f(b)<0:roots.append(brentq(f,a,b,xtol=1e-13))
    candidates=[-2.5,2.5,*roots]
    if abs(x)>=2.5:candidates.append(x)
    t=min(candidates,key=lambda t:(t-x)**2+(h(t)-z)**2)
    return dict(status='closest_normal' if t in roots or t==x else 'foot_limit',footpoint_x=t,
        distance=math.hypot(t-x,z-h(t)),vertical=z-h(x),normal_root_count=len(roots))


def heights(data):
    records=list(csv.DictReader((data/'profiles.csv').open()));flat={};hill=[]
    for series in sorted({r['series_id'] for r in records if r['series_id'].startswith(('f3b','f3c','f3d'))}):
        rows=[r for r in records if r['series_id']==series and r['quality']=='visible'];surface=series[-1];d=.0035 if surface=='R' and not series.startswith('f3b') else 0.
        z=[d+.25*float(r['axis_y']) for r in rows];lower=[d+.25*(float(r['axis_y'])-float(r['axis_y_bound'])) for r in rows]
        flat[series]=dict(minimum_height_m=min(z),minimum_extraction_lower_height_m=min(lower),
            assumed_displacement_m=d,coordinate='paper physical height; delta=.25m nominal')
    for r in records:
        if r['quality']!='visible' or not r['series_id'].startswith(('f4','f5','f8','f10','f12')):continue
        if not r['x_over_H']:continue
        z=.002/.04 if r['series_id'].startswith('f12') and '_SS' in r['series_id'] else .007/.04 if r['series_id'].startswith('f12') else float(r['y_over_H']) if r['y_over_H'] else None
        if z is None:continue
        v=normal_distance(float(r['x_over_H']),z)
        hill.append(dict(series=r['series_id'],point_id=r['point_id'],x_over_Hh=float(r['x_over_H']),global_z_over_Hh=z,**v))
    return dict(flat=flat,hill=hill,Fig12=dict(S_height_m=.002,R_height_m=.007,uncertainty='fixed nominal probe heights; experimental position uncertainty not reported'),
        Fig11='physical height unavailable: delta75 unspecified; no physical height invented',
        uncertainty='frozen deterministic extraction envelopes, not experimental/probe-position uncertainty; nominal delta and assumed displacement separate')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--native-probe',type=Path);a=p.parse_args()
    m=RequestManifest(a.output,dict(schema='CT4IndependentReferenceV1',contract_sha256=sha(a.contract),dataset_sha256=sha(a.data/'snapshot.json'),source_sha256=sha(Path(__file__))))
    rows,z,v,dz,dv=observations(a.data);nu=1/10500;central={};fits={};cases={}
    proc=subprocess.Popen([str(a.native_probe)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True) if a.native_probe else None
    try:
        for model in ('A95','A93','E94E95'):
            f=fit(model,z,v);fits[model]=f;kap=parameters(model)['kappa'];fitU=f['B']+f['ustar']/kap*np.log(z/.04)
            fits[model]['mean_residuals']=[dict(point_id=int(r['point_id']),used=int(r['point_id'])<4,z_m=float(zz),observed=float(vv*f['Uinf']),predicted=float(uu),residual=float(uu-vv*f['Uinf'])) for r,zz,vv,uu in zip(rows,z,v,fitU,strict=True)]
            fits[model].update(unshifted_z0_m=.04*math.exp(-f['B']*kap/f['ustar']),predicted_U_Hh=f['B'],normalization_preserved=True)
            for mm in (6.,5.5,6.5):
                H=mm/40;fitr=roughness(model,H,f,nu,True);require(fitr['status']=='fit','central_roughness_unattainable')
                c=Column(model,f['K0'],H,fitr['h_s'],nu);q=c.result();tight=Column(model,c.K,H,c.hr,nu,order=96).result()
                q.update(mean_fit=fitr,utau=f['ustar'],B=f['B'],Uinf=f['Uinf'],H_mm=mm,
                    H_plus=H*c.um/nu,H_minus_dobs_plus=(H-.0035/.04)*c.um/nu,cube_clearance_mm=mm-5,
                    nu_over_nut_H=nu/float(c.nt(H)),finite_Re_indicator=math.sqrt(nu*c.A/H)/c.K,
                    curvature_ratio=math.pi**2/12.5*H,thin_point1=math.pi**2/12.5*H<=.1,thin_point25=math.pi**2/12.5*H<=.25,
                    global_reach=H<12.5/math.pi**2,quadrature_uncertainty={k:abs(q[k]-tight[k]) for k in ('J0','J2','matching_velocity','budget_defect')})
                q['native']=native_compare(proc,c,q) if proc else dict(status='not_attempted')
                G=np.zeros((3,3));G[0,2]=c.um/(kap*H);law=oracle(model,c.K,c.A/H,G)
                q['outer_equilibrium']=dict(R=law['R'].tolist(),Q=law['Q'].tolist(),covariance_minimum=float(np.linalg.eigvalsh(law['R'])[0]),
                    RMS_over_utau=(np.sqrt(np.diag(law['R']))/c.um).tolist(),Pi=-float(law['Q'][2,2]),outlet_traction=[float(law['Q'][2,2]-law['Q'][0,0]),c.um*c.um],
                    Qxz=float(law['Q'][0,2]),Rxz=float(law['R'][0,2]),production=float(law['production']))
                key=model+'-'+str(mm);central[key]=q
                q['reference_samples']=[dict(z=r,outer=side,values=reference(c,r,side).tolist()) for r,side in [(0.,False),(H/32,False),(c.cutoff,False),(c.hr,False),(H,False),(H,True),(1.,True),(2.,True)]]
        # The SAME source coordinate/velocity realization is used in the fit
        # and normalization, including point3 and a changed 40mm bracket.
        vertices=np.array(list(itertools.product((-1.,1.),repeat=12)))
        interior=2*qmc.Sobol(12,scramble=True,seed=4).random_base2(8)-1
        scenarios=np.vstack([vertices,interior,np.zeros((1,12))]);uncertainty={}
        for model in fits:
            out={mm:[] for mm in (5.5,6.,6.5)};brackets=set()
            for s in scenarios:
                zz=z.copy();vv=v.copy();zz[:6]+=s[:6]*dz[:6];vv[:6]+=s[6:]*dv[:6]
                f=fit(model,zz,vv);brackets.add(tuple(f['normalization_bracket']))
                for mm in out:
                    q=roughness(model,mm/40,f,nu)
                    out[mm].append(q['h_s']*40 if q['status']=='fit' else None)
            uncertainty[model]=dict(h_s_mm={str(mm):dict(minimum=min(x for x in vals if x is not None),maximum=max(x for x in vals if x is not None),failed_scenarios=sum(x is None for x in vals)) for mm,vals in out.items()},
                normalization_brackets=sorted(brackets),scenarios=len(scenarios),fit_point_ids=[0,1,2,3],
                method='4096 joint vertices + 256 deterministic Sobol interior + central; fixed fit memberships; correlated shared fit/normalization values; evaluated extraction envelope, not a confidence interval or a proof of extrema over the continuum box')
        smooth=[];sr,sz,sv,_,_=observations(a.data,'S')
        for model in fits:
            f=fit(model,sz,sv,.033)
            for mm in (.5,1.,1.5,2.,8.):
                c=Column(model,f['K0'],mm/40,0.,1/12000);q=c.result();q.update(H_mm=mm,H_plus=c.H*c.um/c.nu,
                    mean_matching_defect=c.t0*c.J0-(f['B']+c.um/c.p['kappa']*math.log(c.H)),scope='diagnostic unshifted four-point smooth log intercept only; no repair or smooth native admission')
                smooth.append(q)
        result=dict(schema='CT4IndependentReferenceV1',diagnostic_execution='pass',fits=fits,cases=central,extraction_sensitivity=uncertainty,
            observation_heights=heights(a.data),smooth=smooth,scales='Hh=.04m,U_H from unchanged40mm interpolation; nu*=1/Re_H; no measured dimensional U or nu claimed',
            decisions=dict(mean_fit='computed; unused data discrepancies retained',mathematical='per-case cutoff/k/covariance',native_numerical='separate',future_terrain='not_admitted: thin0.1 and physical roughness-sublayer adequacy unresolved'))
        write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT4 mathematical cases',len(central),'extraction scenarios',len(scenarios))
    finally:
        if proc:proc.stdin.close();proc.wait()
    return 0
if __name__=='__main__':raise SystemExit(main())
