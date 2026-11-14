#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT8 independent three-point fit and displaced Accepted oracle.

No production kernel supplies expected values. Native comparison is a separate
sampling adapter. All lengths in a column are normalized by Hh=.04m.
"""
import argparse,itertools,json,math,os,subprocess
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from scipy.optimize import brentq
from ct_accepted_feasibility import parameters,write_json,require,Column
from ct_accepted4_reference import observations,resistance
from ct_accepted7_reference import DisplacedColumn,wall_root
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha,runtime_identity

POLICY='SamePhysicalWallDisplacedAcceptedV1'
PROVIDER='DisplacedNeutralCompositeV1'
FIT='CT8PhysicalBelow25mmEqualIds012V1'


def fit(z,v,normal):
    ids=np.where(z<=.025)[0].tolist()
    require(ids==[0,1,2],'changed CT8 physical-height membership')
    p=parameters('A95');us=.051/normal;d=.0035
    C=float(np.mean(v[ids]/normal-us/p['kappa']*np.log((z[ids]-d)/.04)))
    z0=.04*math.exp(-C*p['kappa']/us)
    pred=us/p['kappa']*np.log((z-d)/z0)
    return dict(model='A95',C=C,z0_m=z0,ustar=us,K=p['q']*us*us,model_parameters=p,
        ids=ids,weights=[1,1,1],predicted=pred.tolist(),residuals=(pred-v/normal).tolist(),
        predicted_U_Hh=us/p['kappa']*math.log((.04-d)/z0),normalization=normal)


@lru_cache(None)
def column(K,H,hr,nu,d,order=96,gamma=0.):
    require(np.isfinite([K,H,hr,nu,d,gamma]).all() and K>0 and nu>0 and 0<=d<H and 0<=hr<H,
        'invalid displaced physical scalar, displacement, height or roughness')
    require(min(1.,1+gamma*H)>=.1,'invalid metric')
    return DisplacedColumn('A95',K,H,hr,nu,d,order,gamma)


def reference(col,z,outer=None):
    outer=z>col.H if outer is None else outer;d=col.disp;kap=col.p['kappa']
    if outer:
        v=col.t0*col.J0+col.um/kap*math.log((z-d)/(col.H-d));du=col.um/(kap*(z-d))
        G=np.zeros((3,3));G[0,2]=du;law=oracle('A95',col.K,col.A/(z-d),G)
        return np.array([v,col.K,col.A/(z-d),-law['Q'][2,2],du,0.,-col.A/(z-d)**2,0.,kap*col.um*(z-d)])
    require(col.cutoff is not None,'cutoff budget inadmissible')
    return np.array([col.t0*col.integrate(lambda r:1/col.D(r),upper=z) if z else 0.,col.k(z),col.A/(max(z,col.cutoff)-d),0.,
        col.t0/col.D(z),float(col.Q(z)/(col.J(z)*col.D(z))),0. if z<col.cutoff else -col.A/(z-d)**2,float(col.Q(z)),float(col.D(z))])


def make_case(f,root):
    H=.1875;d=.0875;nu=1/10500;col=column(f['K'],H,root['h_s'],nu,d);r=col.result()
    r.update(H=H,H_mm=7.5,nu=nu,hr=root['h_s'],zv=col.zv,
        counterfactual=dict(epsilon_amplitude=col.A,traction=col.t0,u=col.um),
        outer_equilibrium=dict(outlet_traction=[float(reference_stress(col)[2,2]-reference_stress(col)[0,0]),col.t0]))
    # One frozen representation recipe; no outcome-driven split search.
    inner=[0.,H/64,H/32,H/16,H/8,H/4,H/2,H,root['h_s'],col.cutoff]
    if 0<col.zv<H:inner.append(col.zv)
    outer=[H];z=H
    while z*1.5<20.:z*=1.5;outer.append(z)
    return dict(model='A95',policy='Natural',provider=PROVIDER,wall_policy=POLICY,fit_policy=FIT,
        H=H,H_mm=7.5,Z=20.,displacement=d,nu=nu,ustar=f['ustar'],B=f['C'],K0=f['K'],z0_m=f['z0_m'],
        roughness_input=root['h_s'],reference_column=r,x=[0.,.125,.25],knots=sorted(set(inner+outer+[20.])),degree=4,pressure_degree=2,
        coordinate='physical floor r=z=0; zeta=z-d only in outer log/epsilon; fixed Hh40mm normalization',
        integration=dict(normal=32,tangential=8,mass=12,volume=12,face=64))


def reference_stress(col):
    G=np.zeros((3,3));G[0,2]=col.um/(col.p['kappa']*(col.H-col.disp))
    return oracle('A95',col.K,col.A/(col.H-col.disp),G)['Q']


def native_checks(rt,out,f,root):
    out.mkdir();H=.1875;d=.0875;nu=1/10500;K=f['K'];hr=root['h_s'];rows={}
    def sample(label,request,expected=None,reason=None):
        inp=out/(label+'-input.json');dest=out/(label+'.json');write_json(inp,request)
        code=subprocess.run([rt/'RANSCurvedAcceptedKernelProbe',inp,dest],capture_output=True,text=True)
        q=json.loads(dest.read_text());row=dict(exit_code=code.returncode,packet_sha256=sha(dest),request_sha256=sha(inp))
        if reason:
            require(code.returncode!=0 and reason in q['reason'],'wrong rejection:'+label+str(q));row.update(status='pass',expected_rejection=q)
        else:
            require(code.returncode==0,'kernel forward failed:'+str(q));col=expected;defects={}
            defects['cutoff']=abs(q['cutoff']-col.cutoff)/col.H
            defects['budget']=abs(q['budget']-col.PJ)/col.A
            defects['minimum_budget']=abs(q['minimum_budget']-col.minimum_budget)/col.A
            for actual in q['samples']:
                z=actual[0];truth=[z,col.J(z),col.production(z),col.A/(max(z,col.cutoff)-col.disp),col.k(z),col.Q(z)/(col.J(z)*col.D(z)),col.Q(z),col.D(z)]
                scales=[1.,1.,col.A/col.H,col.A/col.H,col.K,col.K/col.H,col.A,col.D(col.H)]
                for i,(a,b,s) in enumerate(zip(actual,truth,scales,strict=True)):defects[str(i)]=max(defects.get(str(i),0.),abs(a-b)/s)
            require(max(defects.values())<=2e-9,'independent displaced kernel mismatch:'+str(defects));row.update(status='pass',defects=defects)
        rows[label]=row;return q
    for gamma in (-math.pi**2/12.5,0.,math.pi**2/12.5):
        col=column(K,H,hr,nu,d,96,gamma);require(col.cutoff is not None,'central metric cutoff not admissible')
        pts=sorted(set([0.,d,col.cutoff,H,*col.k_extrema()['candidates']]))
        req=dict(model='A95',K=K,H=H,hr=hr,nu=nu,displacement=d,displacement_policy=POLICY,gamma=gamma,r=pts,order=96)
        sample('central-gamma-'+str(gamma),req,col)
    for gamma in (-.25,0.,.25):
        col=column(K,H,hr,nu,0.,96,gamma);req=dict(model='A95',K=K,H=H,hr=hr,nu=nu,displacement=0.,displacement_policy=POLICY,gamma=gamma,r=[0.,H/2,H],order=96)
        new=sample('d0-'+str(gamma),req,col);req['displacement_policy']='UnshiftedV1'
        old=sample('original-'+str(gamma),req,col)
        require(new==old,'d0 profile arithmetic changed')
    base=dict(model='A95',K=K,H=H,hr=hr,nu=nu,displacement=d,displacement_policy=POLICY,r=[0.,H],order=96)
    for label,change,why in [('negative-d',dict(displacement=-.01),'0<=d<H'),('height',dict(displacement=H),'0<=d<H'),
        ('roughness',dict(hr=H),'0<=h<H'),('metric',dict(gamma=-2/H),'metric'),('budget',dict(hr=0.,displacement=H*.999),'cutoff_budget_inadmissible'),
        ('foreign-model',dict(model='A93'),'Natural A95'),('unshifted-d',dict(displacement_policy='UnshiftedV1'),'UnshiftedV1')]:sample(label,dict(base,**change),reason=why)
    return dict(status='pass',checks=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--runtime',type=Path);a=p.parse_args()
    m=RequestManifest(a.output,dict(schema='CT8IndependentReferenceV1',contract=sha(a.contract),data={n:sha(a.data/n) for n in ('profiles.csv','snapshot.json','extraction.json','source_manifest.json')},source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'RANSCurvedAcceptedKernelProbe') if a.runtime else None))
    rows,z,v,dz,dv=observations(a.data);normal=float(np.interp(.04,z,v));f=fit(z,v,normal);nu=1/10500;d=.0875;H=.1875
    root=wall_root('A95',H,f,nu,d);require(root['status']=='root','central root failed:'+str(root));col=column(f['K'],H,root['h_s'],nu,d)
    central=col.result();coarse=column(f['K'],H,root['h_s'],nu,d,48).result()
    require(central['status']=='admissible_frozen_column' and abs(central['Q_H'])/col.A<2e-10,'central reference inadmissible')
    ex=[]
    for signs in [np.zeros(6),*itertools.product((-1.,1.),repeat=6)]:
        zz=z.copy();vv=v.copy();zz[:3]+=np.array(signs[:3])*dz[:3];vv[:3]+=np.array(signs[3:])*dv[:3]
        sf=fit(zz,vv,normal);sr=wall_root('A95',H,sf,nu,d);item=dict(signs=list(signs),z0_m=sf['z0_m'],fit=sf,roughness=sr)
        if sr['status']=='root':item['column']=column(sf['K'],H,sr['h_s'],nu,d).result()
        ex.append(item)
    c=make_case(f,root);result=dict(schema='CT8IndependentReferenceV1',central_fit=f,central_root=root,central_column=central,
        quadrature_difference={k:abs(central[k]-coarse[k]) for k in ('J0','P0','cutoff')},case=c,extraction_scenarios=ex,
        extraction_summary=dict(cases=65,failed_roots=sum(q['roughness']['status']!='root' for q in ex),z0_mm_range=[min(q['z0_m'] for q in ex)*1000,max(q['z0_m'] for q in ex)*1000],scope='joint deterministic evaluated extraction corners, no probability; fixed normalization/displacement/friction/viscosity and membership'),
        central_status='admissible',native='not_attempted',physical_scope='central fit/reference, not experimental validation or native preservation')
    if a.runtime:result['native']=native_checks(a.runtime,a.output/'kernel',f,root)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT8 independent central reference and65 extraction scenarios complete');return 0
if __name__=='__main__':raise SystemExit(main())
