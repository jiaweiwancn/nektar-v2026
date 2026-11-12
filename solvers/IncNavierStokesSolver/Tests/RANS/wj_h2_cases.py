# SPDX-License-Identifier: MIT
"""Production session/Functions cases. Physical heights, no probe-built state."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline
from numpy.polynomial import Chebyshev, Polynomial
import xml.etree.ElementTree as ET
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_info,set_parameter,write
from wj_h2_profiles import Completion,wall,QE,CD,KAPPA,NU,US,D0,TOP


def num(x):return format(float(x),'.17g')
def poly(c,x):
    s=num(c[0])
    for v in c[1:]:s=f'({s}*({x})+({num(v)}))'
    return s


def piecewise(p,lower,upper):
    """PPoly local polynomials with finite inactive branches (no if division)."""
    x=p.x;c=p.c;terms=[f'(y<{num(x[0])})*({num(lower)})']
    for j in range(len(x)-1):
        terms.append(f'((y>={num(x[j])})*(y<{num(x[j+1])}))*({poly(c[:,j],f"y-({num(x[j])})")})')
    terms.append(f'(y>={num(x[-1])})*({num(upper)})')
    return '('+'+'.join(terms)+')'


def tail_expression(c,panel,target):
    rows=[v for v in c.rows if v['panel']==panel]
    x=np.array([float(v['y_over_H']) for v in rows]);v=np.array([float(v['value']) for v in rows])
    if panel!='a':v*=US
    p=PchipInterpolator(x,v);t=CubicHermiteSpline([x[-1],12.5],[v[-1],target],[p.derivative()(x[-1]),0])
    from scipy.interpolate import PPoly
    return piecewise(PPoly(np.column_stack([p.c,t.c]),np.r_[p.x,12.5]),v[0],target)


def approximate(fn,knots,degree=12):
    """Declared independent expression approximation; native Q3 projection separate."""
    expressions=[];maxerr=0.
    for a,b in zip(knots[:-1],knots[1:]):
        cp=Chebyshev.interpolate(lambda t:np.vectorize(fn,otypes=[float])((a+b)/2+(b-a)/2*t),degree)
        coefficients=cp.convert(kind=Polynomial).coef
        t=f'(2*(y-({num(D0+a)}))/({num(b-a)})-1)'
        expression=poly(coefficients[::-1],t)
        mask=f'((y>={num(D0+a)})*(y<{num(D0+b)}))'
        expressions.append(f'{mask}*({expression})')
        test=np.linspace(-1,1,101);maxerr=max(maxerr,float(np.max(abs(cp(test)-np.array([fn((a+b)/2+(b-a)/2*s) for s in test])))))
    expressions.append(f'(y>={num(D0+knots[-1])})*({num(fn(knots[-1]))})')
    return '('+'+'.join(expressions)+')',maxerr


def continuous_outer(c,z):
    """Physical parallel seed gradient/covariance; spanwise stress retained."""
    from wj_h2_profiles import constitutive
    G=np.zeros((3,3));G[0,2]=c.mean(z,1)
    R,_,_=constitutive(float(c.k(z)),float(c.epsilon(z)),G)
    return G,R


def reduced_normal_offset(K,Rnn):
    """T_w,nn=-pi_w, T_o,nn=2K/3-pi_o-Rnn for M_w,nn=0."""
    return float(Rnn-2*K/3)


def expressions(c,kind,pressure_seed='SignCorrectedV1',outlet='ZeroTractionV1'):
    H=c.H;r=f'(y-({num(D0)}))';K=c.K;w=c.w
    if kind=='analytic':
        K=QE*US**2;w=wall(K,H,c.h);A=US**3/KAPPA
        nt=f'({num(w["a"])})*(({r})-({num(w["zv"])}))'
        # Include the physical molecular-only segment for positive intercept.
        d=NU+w['a']*max(-w['zv'],0);start=max(0,w['zv'])
        J=f'((({r})<{num(start)})*({r})/{num(NU)}+(({r})>={num(start)})*({num(start/NU)}+log(1+({num(w["a"])})*(({r})-({num(start)}))/({num(NU+w["a"]*max(start-w["zv"],0))}))/{num(w["a"])}))'
        # Negative log argument in inactive branches is avoided by separate
        # approximation of the finite-nu inner velocity on physical panels.
        cuts=sorted(set([0.,H]+[v for v in (w['rc'],c.h,w['zv'],.001,.005,.02,.05,.1,.25) if 0<v<H]))
        ui,ue=approximate(lambda z:US**2*float(w['J'](z)),cuts)
        ki,ke=approximate(w['kref'],cuts)
        inner=dict(u=ui,v='0',k=ki,p=num(-K/(81/20)))
        outer=dict(u=f'({num(US**2*w["J0"])})+({num(US/KAPPA)})*log(({r})/{num(H)})',v='0',k=num(K),epsilon=f'{num(A)}/({r})',p='0')
        from wj_h2_profiles import constitutive
        G=np.zeros((3,3));G[0,2]=US/KAPPA/H;R,_,_=constitutive(K,A/H,G)
        outlet_inner=dict(u=num(K/(81/20)),v=num(US**2))
        outlet_outer=dict(u=num(2*K/3-R[0,0]),v=num(US**2))
        top=dict(u=num(US**2));topflux=dict(k='0',epsilon=num(CD*K*K/1.3/(TOP-D0)))
        return inner,outer,outlet_inner,outlet_outer,top,topflux,dict(expression_errors=dict(inner_U=ue,inner_k=ke),K=K,wall=w)
    U=tail_expression(c,'a',1.);su=tail_expression(c,'c',.003);sw=tail_expression(c,'d',.003)
    k=f'(.5*((1+({num(c.span)}))*({su})^2+(2-({num(c.span)}))*({sw})^2))'
    L=c.length_fraction*6.25
    mixing=f'({num(KAPPA)}*({r})*({num(1+KAPPA*H/L)})/(1+({num(KAPPA/L)})*({r})))'
    epsilon=f'({num(CD**.75)})*({k})^1.5/({mixing})'
    cuts=sorted(set([0.,H]+[v for v in (w['zv'],w['rc'],c.h,.001,.005,.02,.05,.1,.25) if 0<v<H]))
    ui,ue=approximate(lambda z:float(c.U(D0+z)),cuts)
    # Near-wall measured stresses are extrapolated constant, so modeled ground
    # k flux is initially compatible; no k(0)=0 and no neutral-k substitution.
    ki,ke=approximate(lambda z:float(c.k(D0+z)),cuts)
    _,R=continuous_outer(c,D0+H)
    C=reduced_normal_offset(K,R[2,2])
    inner=dict(u=ui,v='0',k=ki,p=num(C))
    outer=dict(u=U,v='0',k=k,epsilon=epsilon,p='0')
    errors=dict(inner_U=ue,inner_k=ke)
    if pressure_seed not in ('SignCorrectedV1','VerticalBalanceV1'):raise ValueError('physical pressure seed')
    if outlet not in ('ZeroTractionV1','FrozenSeedTractionV1'):raise ValueError('physical outlet diagnostic')
    if outlet=='FrozenSeedTractionV1' and pressure_seed!='VerticalBalanceV1':raise ValueError('frozen traction requires the balanced seed')
    # All continuous profile breakpoints, independent of the native mesh.
    cuts=sorted(set([H,12.5-D0,TOP-D0]+[float(v['y_over_H'])-D0 for v in c.rows if D0+H<float(v['y_over_H'])<TOP]))
    # Only the new nonlinear pressure/traction expressions need extra panels.
    # Keep the historical velocity/k/epsilon completion expressions unchanged.
    cuts=sorted(set(cuts+[(a+b)/2 for a,b in zip(cuts[:-1],cuts[1:])]))
    if pressure_seed=='VerticalBalanceV1':
        outer['p'],errors['outer_pi']=approximate(lambda r:2*c.k(D0+r)/3-continuous_outer(c,D0+r)[1][2,2]+C,cuts,18)
    oi,oo=dict(u='0',v='0'),dict(u='0',v='0')
    if outlet=='FrozenSeedTractionV1':
        oi=dict(u=num(-C),v=num(c.UH/w['J0']))
        oo['u'],errors['outlet_Txx']=approximate(lambda r:continuous_outer(c,D0+r)[1][2,2]-continuous_outer(c,D0+r)[1][0,0]-C,cuts,18)
        oo['v'],errors['outlet_Tzx']=approximate(lambda r:-continuous_outer(c,D0+r)[1][2,0],cuts,18)
    return inner,outer,oi,oo,None,None,dict(expression_errors=errors,K=K,wall=w,
        seed=dict(version=pressure_seed,offset=C,outlet=outlet,scope='initialization only; no imposed interface pressure row or restoring source'))


def case_xml(folder,kind='physical',steps=4,dt=.02,refined=False,reuse='Off',span=.5,length=.09,restart=None,reject=0,
             pressure_seed='SignCorrectedV1',outlet='ZeroTractionV1',xknots=None,tangential_order=8,extra_physical_z=None):
    if kind not in ('physical','analytic') or not math.isfinite(dt) or dt<=0 or not 0<=steps<=32:raise ValueError('bounded H2 case')
    c=Completion(span,length);H=c.H
    extras=[] if extra_physical_z is None else list(extra_physical_z)
    if len(extras)>2 or any(isinstance(v,bool) or not math.isfinite(v) or not D0+H<v<TOP for v in extras) or any(a>=b for a,b in zip(extras[:-1],extras[1:])):
        raise ValueError('at most two strictly increasing finite outer PHYSICAL normal knots')
    xknots=[0.,8.] if xknots is None else list(xknots)
    if len(xknots)<2 or xknots[0]!=0 or xknots[-1]!=8 or not all(math.isfinite(v) for v in xknots) or any(a>=b for a,b in zip(xknots[:-1],xknots[1:])):raise ValueError('strict physical x knots in [0,8]')
    if isinstance(tangential_order,bool) or tangential_order not in (8,12,16):raise ValueError('declared H3 tangential orders8/12/16')
    inner,outer,oi,oo,top,tf,metadata=expressions(c,kind,pressure_seed,outlet)
    w=metadata.pop('wall');K=metadata['K']
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    root=ET.parse(Path(__file__).with_name('EquilibriumWallFlat.xml')).getroot();co=root.find('CONDITIONS')
    for child in list(co):
        if child.tag!='VARIABLES':co.remove(child)
    ET.SubElement(co,'SOLVERINFO');ET.SubElement(co,'PARAMETERS')
    scheme=ET.SubElement(co,'TIMEINTEGRATIONSCHEME');ET.SubElement(scheme,'METHOD').text='IMEX';ET.SubElement(scheme,'ORDER').text='1'
    # Fit the exact mean-data knots; r is physical normal distance in flat inner.
    inner_knots=sorted(set([0.,.001,.005,.02,.05,.1,.25,H]+[v for v in (w['zv'],w['rc'],c.h) if 0<v<H]))
    # Avoid a sub-1e-4 metre cell at the virtual viscosity origin: quadrature
    # already splits it; the first element resolves molecular diffusion.
    inner_knots=[v for v in inner_knots if v==0 or v>=.001]
    physical_outer=sorted(set([H,12.5-D0,TOP-D0]+[float(v)-D0 for v in c.Udata[0] if v>D0+H]))
    if refined:
        inner_knots=sorted(set(inner_knots+[(a+b)/2 for a,b in zip(inner_knots[:-1],inner_knots[1:]) if b<=.25]))
        physical_outer=sorted(set(physical_outer+[(a+b)/2 for a,b in zip(physical_outer[:-1],physical_outer[1:]) if a<6.3]))
    if extras:
        if any(abs(z-D0-r)<1e-12 for z in extras for r in physical_outer):raise ValueError('extra physical knot duplicates an existing face')
        physical_outer=sorted(physical_outer+[z-D0 for z in extras])
    rho=lambda r:H+(r-H)*(TOP-H)/(TOP-D0-H)
    z=inner_knots+[rho(v) for v in physical_outer[1:]]
    tensor_mesh(root,z,3,2,8.,1.,points=[8,8],xknots=xknots)
    for node in root.findall('./EXPANSIONS/E'):node.set('FIELDS',node.get('FIELDS').replace('eps_tilde','epsilon'))
    for node in co.findall('./VARIABLES/V'):
        if node.text=='eps_tilde':node.text='epsilon'
    info=dict(EQTYPE='UnsteadyNavierStokes',SolverType='RANSEquilibriumWallEvolution',Projection='Galerkin',GlobalSysSoln='DirectFull',
        TurbulenceModel='WJHighRe',RANSWJGradientPolicy='DevStrainRawWorkV1',RANSWallFormulation='WJCalibratedBudgetAWF',
        RANSWallNativeGeometry='BudgetTerrainV1',RANSWallMapContinuation='CompositeTunnelV1',RANSWallMap='Plane',
        RANSWallTemporalCoupling='PimpleSingleSweep',RANSWallNonlinearPolicy='Legacy',RANSWallVariant='Original',
        RANSWallInnerStress='Linear',RANSWallDestruction='Original',RANSWallTransition='Fixed',RANSWallKSampling='InterfaceTrace',
        RANSWallPhysicalRangeGuard='WJNativeV1',RANSWallCondense='False',RANSWallReference='Functions',RANSWallFunctionCoordinates='Physical',
        RANSWallEvaluationReuse=reuse,RANSWallPimpleVerify='Off',RANSWallComparison='Off',RANSWallBoundaryMode='OpenTunnelV1',
        RANSWallOpenTop='ReferenceTractionV1' if kind=='analytic' else 'SlipV1',RANSWallPimpleStopOnBound='On',
        RANSWallPimpleKBound='MassProjection',RANSWallPimpleKInventory='WJNativeV1',RANSWallPimpleKBoundSolver='RetainedPool',
        RANSWallPimpleSnapshots='On',RANSWallCheckpoint=str(folder/'CP'))
    if restart:info['RANSWallRestart']=str(Path(restart).resolve())
    for key,value in info.items():set_info(root,key,value)
    for key,value in dict(Time=0.,TimeStep=dt,NumSteps=steps,IO_CheckSteps=1,RANSWallDiagnostics=1,Kinvis=NU,RANSWallHeight=H,
        RANSWallRoughness=c.h,RANSWallSampleHeight=H,RANSWallTailStart=H,RANSWallChi=0.,RANSWallFrictionVelocity=US,
        RANSWallMapDisplacement=D0,RANSWallMapTunnelTop=TOP,RANSWallMapAmplitude=0.,RANSWallMapLengthX=2.5,
        RANSWallEta=0.,RANSWallTangentialOrder=tangential_order,RANSWallNormalOrder=16,RANSWallPimpleKBoundMaxSeconds=60.,
        RANSWallPimpleKFloorRelative=1e-8,RANSWallPimpleHardAuditRejectStep=reject).items():set_parameter(root,key,value)
    for family in ('InitialConditions','BoundaryValues'):
        for side,values in (('Inner',inner),('Outer',outer)):
            f=ET.SubElement(co,'FUNCTION',{'NAME':family+side})
            for variable,value in values.items():
                if family=='BoundaryValues' and variable=='p':continue
                ET.SubElement(f,'E',{'VAR':variable,'VALUE':value})
    for name,values in (('OutletTractionInner',oi),('OutletTractionOuter',oo),('TopTractionOuter',top),('TopDiffusiveFluxOuter',tf)):
        if values is None:continue
        f=ET.SubElement(co,'FUNCTION',{'NAME':name})
        for variable,value in values.items():ET.SubElement(f,'E',{'VAR':variable,'VALUE':value})
    path=folder/'case.xml';write(root,path)
    metadata.update(kind=kind,dt=dt,steps=steps,reuse=reuse,refined=refined,chart_z=z,x=xknots,tangential_order=tangential_order,
        physical_outer_r=physical_outer,completion=c.summary(),K=K,
        reference_cutoff=w['rc'],reference_k_min=min(w['kref'](r) for r in np.linspace(0,H,41)))
    if extras:metadata['extra_physical_z']=extras
    (folder/'case-design.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    return path

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--kind',choices=('physical','analytic'),default='physical')
    p.add_argument('--steps',type=int,default=4);p.add_argument('--dt',type=float,default=.02);p.add_argument('--refined',action='store_true')
    p.add_argument('--reuse',choices=('Off','ExactStepV1'),default='Off');p.add_argument('--restart',type=Path)
    p.add_argument('--pressure-seed',choices=('SignCorrectedV1','VerticalBalanceV1'),default='SignCorrectedV1')
    p.add_argument('--outlet',choices=('ZeroTractionV1','FrozenSeedTractionV1'),default='ZeroTractionV1')
    p.add_argument('--x-knots',type=lambda s:[float(v) for v in s.split(',')]);p.add_argument('--tangential-order',type=int,default=8)
    p.add_argument('--extra-physical-z',type=lambda s:[float(v) for v in s.split(',')])
    a=p.parse_args()
    print(case_xml(a.output,a.kind,a.steps,a.dt,a.refined,a.reuse,restart=a.restart,pressure_seed=a.pressure_seed,
        outlet=a.outlet,xknots=a.x_knots,tangential_order=a.tangential_order,extra_physical_z=a.extra_physical_z))
