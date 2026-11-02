"""Physical upstream masks and strong scalar MMS for M6, run on Linux.

The source is an independently prescribed differential expression, not an
assembled native load. U=pi=0 and k=k(t). E=B(ell,t)/r has a spatially varying
EP/original matching trace. Broad state thresholds deliberately select their
unit plateau for this MMS; separate probes exercise every switching branch.
"""
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from cls_epawf_m4_cases import tensor_mesh, _evolution
from cls_epawf_m3_cases import set_info, set_parameter, write


def arc_offset(amplitude=.02, length=1.):
    # Independent quadrature; upstream ell=q1-offset for the cosine hill.
    import numpy as np
    z,w=np.polynomial.legendre.leggauss(96)
    x=(z+1)*length/2
    return float(length/2*np.dot(w,np.sqrt(1+(amplitude*math.pi/(2*length)*np.sin(math.pi*x/length))**2))-length)


def transition_session(template,destination,family='flat',degree=3,dimension=2,
                       steps=0,order=2,dt=.005,time=.3,condensed=False,
                       restart=None,checkpoint=None,ownership_shift=0,
                       quadrature=12,normal_order=16,refine=1,probe=False,
                       active=True,stationary=False,model='CLSHighRe'):
    if family not in ('flat','feet'):raise ValueError('M6 requires an exactly flat upstream fetch')
    root=ET.parse(template).getroot()
    offset=arc_offset() if family=='feet' else 0.
    x=[-4+offset,-3+offset,-2+offset,-1.,1.,2.]
    if refine>1:x=sorted(set(x+[a+(b-a)*i/refine for a,b in zip(x,x[1:]) for i in range(1,refine)]))
    z=[0.,2.,4.]
    if refine>1:z=[0.,1.,2.,3.,4.]
    tensor_mesh(root,z,degree,dimension,1.,1.,xknots=x)
    _evolution(root,steps,order,dt,time,condensed,restart,checkpoint)
    set_info(root,'RANSWallMap','CosineHill' if family=='feet' else 'Flat')
    set_info(root,'RANSWallFunctionCoordinates','Chart')
    set_info(root,'RANSWallTransition','Upstream')
    set_info(root,'RANSWallReference','Functions')
    set_info(root,'TurbulenceModel',model)
    parameters={'RANSWallHeight':2.,'RANSWallEta':0,'RANSWallChi':int(active),
        'RANSWallTangentialOrder':quadrature,'RANSWallNormalOrder':normal_order,
        'RANSWallOwnershipShift':ownership_shift,'RANSWallTransitionStart':-3.,
        'RANSWallTransitionEnd':-2.25,'RANSWallBaselineStart':-2.,
        'RANSWallSensorDistanceA':2.5,'RANSWallSensorDistanceB':3.5,
        'RANSWallSensorOn':1. if probe else 100.,'RANSWallSensorOff':30. if probe else 200.,
        'RANSWallExpressionDimension':dimension,
        'RANSWallMapAmplitude':.02,'RANSWallMapLengthX':1.}
    if model=='StandardKEpsilon':
        parameters['Qstar']=1/math.sqrt(.09)
        parameters['Kappa']=math.sqrt(.48*1.3/parameters['Qstar'])
    for name,value in parameters.items():set_parameter(root,name,repr(value))
    conditions=root.find('CONDITIONS')
    for node in list(conditions.findall('FUNCTION')):conditions.remove(node)
    if not probe:
        for name,values in expressions(dimension,offset,active,stationary,model).items():
            node=ET.SubElement(conditions,'FUNCTION',{'NAME':name})
            for var,value in values.items():ET.SubElement(node,'E',{'VAR':var,'VALUE':value})
    return write(root,destination)


def expressions(dimension,offset,active=True,stationary=False,model='CLSHighRe'):
    r='y' if dimension==2 else 'z'
    k='K0' if stationary else '(K0*(1+0.02*sin(t)))'
    kd='0' if stationary else '(0.02*K0*cos(t))'
    # The installed interpreter omits loading a constant operand of a binary
    # function when only its other operand is variable. Explicit zero*t keeps
    # both operands dynamic, without changing the prescribed mathematics.
    # The expression probe independently verifies the actual evaluated values.
    q=f'min((1+0*t),max((0*t),((x-({offset:.17g}))+3)/0.75))'
    chi=f'(1-10*({q})^3+15*({q})^4-6*({q})^5)' if active else '0'
    first=f'(-30*({q})^2*(1-({q}))^2/0.75)' if active else '0'
    second=f'(-60*({q})*(1-({q}))*(1-2*({q}))/(0.75^2))' if active else '0'
    original=f'(({k})^1.5/2.55)'
    corrected=f'((({k})/Qstar)^1.5/Kappa)'
    delta=f'(({corrected})-({original}))'
    B=f'(({original})+({delta})*({chi}))'
    E=f'(({B})/({r}))'
    Ed='0' if stationary else f'(1.5*({kd})/({k})*({E}))'
    ze=f'(5.1*Kinvis/sqrt({k}))'
    roughx=f'(RANSWallRoughness*sqrt({k})/(70*Kinvis))'
    roughm=f'max(0.5-0.4*({roughx})^0.7,1-0.79*({roughx})^(-0.28))'
    zv=f'(10.7*Kinvis*(1-({roughx})^({roughm}))/sqrt({k}))'
    nt=f'(0.2295*sqrt({k})*max((0*t),({r})-({zv})))'
    p0=f'((({k})/Qstar)^2*({nt})/(Kinvis+({nt}))^2)'
    ea=f'(({original})/max(({r}),(({ze})+0*t)))'
    logBxx=f'(({delta})*({second})/({B})-(({delta})*({first})/({B}))^2)'
    if offset:
        fx='(-0.02*PI/2*sin(PI*x))';fxx='(-0.02*PI^2/2*cos(PI*x))'
        curvature=f'(-({fxx})/(1+({fx})^2)^1.5*(1-sign(fabs(x)-1))/2)'
        trBr=f'(({curvature})/(1+({r})*({curvature})))'
    else:trBr='0'
    cmu='(0.3*(1-exp(-0.36)))' if model=='CLSHighRe' else '0.09'
    diffusion=f'(({cmu})*({k})^2/1.3*(1/({r})^2-({trBr})/({r})+({logBxx})))'
    ratio=f'(({k})^1.5/(2.5*({E})*({r})))'
    yap=f'(0.83*({E})^2/({k})*max(({ratio})-1,(0*t))*({ratio})^2)' if model=='CLSHighRe' else '0'
    source_e=f'(({Ed})-({diffusion})+1.92*({E})^2/({k})-({yap}))'
    source_i=f'(({kd})+(1-({chi}))*({ea})+({chi})*({p0}))'
    source_o=f'(({kd})+({E}))'
    names=['u','v'] if dimension==2 else ['u','v','w']
    inner={**{name:'0' for name in names},'k':k,'p':'0'}
    outer={**inner,'eps_tilde':E}
    functions={}
    for side,values in [('Inner',inner),('Outer',outer)]:
        for prefix in ['InitialConditions','ExactSolution']:
            functions[prefix+side]=values
        functions['BoundaryValues'+side]={key:value for key,value in values.items() if key!='p'}
        functions['BodyForce'+side]={name:'0' for name in names}
    functions['ScalarSourceInner']={'k':source_i}
    functions['ScalarSourceOuter']={'k':source_o,'eps_tilde':source_e}
    return functions


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('destination',type=Path)
    p.add_argument('--family',choices=['flat','feet'],default='feet');p.add_argument('--degree',type=int,default=3)
    p.add_argument('--steps',type=int,default=2);p.add_argument('--dt',type=float,default=.005)
    a=p.parse_args()
    transition_session(Path(__file__).with_name('EquilibriumWallEvolution.xml'),a.destination,
                       family=a.family,degree=a.degree,steps=a.steps,dt=a.dt)
