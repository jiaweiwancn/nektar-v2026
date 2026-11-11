# SPDX-License-Identifier: MIT
"""Fresh, bounded WJ serial-flat tensor-cubic pilots; no CLS model choices."""
from pathlib import Path
import xml.etree.ElementTree as ET
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_info,set_parameter,write

KN=[0.,.0025,.005,.01,.02,.04,.08,.1,.2,.4,.7,1.,1.15,1.35,1.6,2.]


def normal_knots(level=0):
    if level not in (0,1,2):raise ValueError('WJ normal refinement level must be 0, 1 or 2')
    level=int(level)
    knots=KN
    for _ in range(level):knots=sorted(set(knots+[.5*(a+b) for a,b in zip(knots[:-1],knots[1:])]))
    return knots


def validated_knots(knots,lo,hi,required=()):
    import math
    values=list(map(float,knots))
    if len(values)<2 or not all(math.isfinite(v) for v in values) or values[0]!=lo or values[-1]!=hi:
        raise ValueError('invalid WJ knot endpoints/nonfinite values')
    if any(a>=b for a,b in zip(values,values[1:])) or any(v not in values for v in required):
        raise ValueError('WJ knots must be strictly increasing and retain physical interfaces')
    return values


def local_normal_knots(mesh):
    if mesh=='L1':return normal_knots(1)
    if mesh not in ('M1','M2'):raise ValueError('unknown WJ local mesh')
    knots=normal_knots(1)+[.085,.095,.125,.175]
    if mesh=='M2':knots += [.0825,.0875,.0925,.0975,.1125,.1375,.1625,.1875]
    return sorted(set(knots))


def case_xml(folder,policy='EP',steps=0,restart=None,refined=False,normal=16,dt=1e-6,hard_reject=0,level=None,xknots=None,zknots=None):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    root=ET.parse(Path(__file__).with_name('EquilibriumWallFlat.xml')).getroot()
    conditions=root.find('CONDITIONS')
    for node in list(conditions):
        if node.tag!='VARIABLES':conditions.remove(node)
    ET.SubElement(conditions,'SOLVERINFO');ET.SubElement(conditions,'PARAMETERS')
    scheme=ET.SubElement(conditions,'TIMEINTEGRATIONSCHEME')
    ET.SubElement(scheme,'METHOD').text='IMEX';ET.SubElement(scheme,'ORDER').text='1'
    if level is None:level=int(refined)
    if refined and level!=1:raise ValueError('refined=True is the level-1 compatibility spelling')
    knots=normal_knots(level) if zknots is None else validated_knots(zknots,0.,2.,(.1,1.))
    xknots=[0.,1.] if xknots is None else validated_knots(xknots,float(xknots[0]),float(xknots[-1]))
    tensor_mesh(root,knots,3,2,1.,1.,points=[8,8],xknots=xknots)
    for e in root.findall('./EXPANSIONS/E'):e.set('FIELDS',e.get('FIELDS').replace('eps_tilde','epsilon'))
    for v in conditions.findall('./VARIABLES/V'):
        if v.text=='eps_tilde':v.text='epsilon'
    for key,value in dict(EQTYPE='UnsteadyNavierStokes',SolverType='RANSEquilibriumWallEvolution',
        Projection='Galerkin',GlobalSysSoln='DirectFull',TurbulenceModel='WJHighRe',
        RANSWJGradientPolicy='DevStrainRawWorkV1',RANSWallFormulation='WJCalibrated'+('EPAWF' if policy=='EP' else 'BudgetAWF'),
        RANSWallTransition='UniformEPFlat' if policy=='EP' else 'Fixed',RANSWallTemporalCoupling='PimpleSingleSweep',
        RANSWallNonlinearPolicy='Legacy',RANSWallPhysicalRangeGuard='WJNativeV1',RANSWallVariant='Original',
        RANSWallInnerStress='Linear',RANSWallDestruction='Original',RANSWallKSampling='InterfaceTrace',
        RANSWallReference='Neutral',RANSWallMap='Flat',RANSWallCondense='False',RANSWallComparison='Off',
        RANSWallEvaluationReuse='Off',RANSWallPimpleKBound='MassProjection',RANSWallPimpleKInventory='WJNativeV1',
        RANSWallPimpleKBoundSolver='RetainedPool',RANSWallPimpleSnapshots='On',RANSWallPimpleVerify='Off',
        RANSWallCheckpoint=str(folder.resolve()/'checkpoint'),RANSWallRestart=str(Path(restart).resolve()) if restart else 'None').items():
        set_info(root,key,value)
    for key,value in dict(NumSteps=steps,TimeStep=dt,Time=0.,Kinvis=1.5e-5,RANSWallHeight=1.,RANSWallRoughness=.1,
        RANSWallSampleHeight=1.,RANSWallTailStart=.5 if policy=='EP' else 1.,RANSWallChi=1. if policy=='EP' else 0.,
        RANSWallFrictionVelocity=.5,RANSWallEta=0.,RANSWallTangentialOrder=8,RANSWallNormalOrder=normal,
        RANSWallDiagnostics=1,IO_CheckSteps=1,IO_InfoSteps=1,RANSWallPimpleHardAuditRejectStep=hard_reject,
        RANSWallPimpleKBoundMaxAdd=2000,RANSWallPimpleKBoundMaxTotal=2000,RANSWallPimpleKBoundMaxSeparations=2000,
        RANSWallPimpleKBoundMaxSeconds=120.,RANSWallPimpleKFloorRelative=1e-8).items():set_parameter(root,key,value)
    write(root,folder/'case.xml')
    return folder/'case.xml'
