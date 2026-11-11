# SPDX-License-Identifier: MIT
"""Fresh budget-only H1 verification cases; old flat checkpoints are incompatible."""
from pathlib import Path
import math
import functools
import xml.etree.ElementTree as ET
from scipy.optimize import brentq
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_info, set_parameter, write
from validate_wj1 import QE, KAPPA

NU=.051*.005/3.97
US=.051
HW=.375
X=[-5.,-2.5,0.,2.5,5.]
Z=[0.,HW/2,HW,2.,6.25,20.]

@functools.lru_cache(maxsize=16)
def provisional_roughness(H=HW):
    """Conditional flat inversion from a WJ-slope log anchor, not measured U_m."""
    K=QE*US**2
    ell=H/(H/.005-1)
    target=NU/(KAPPA*US)-ell
    def zv(h):
        x=h*math.sqrt(K)/(70*NU)
        b=10.7 if h==0 else 10.7*(1-x**max(.5-.4*x**.7,1-.79*x**(-.28)))
        return NU*b/math.sqrt(K)
    return brentq(lambda h:zv(h)-target,0.,H*(1-1e-12),xtol=1e-15)

def case_xml(folder,kind='hill',steps=3,dt=1e-6):
    if kind not in ('hill','zero','plane'):raise ValueError('unknown H1 geometry')
    if steps not in (0,1,2,3) or not math.isfinite(dt) or dt<=0:raise ValueError('bounded H1 step design')
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    root=ET.parse(Path(__file__).with_name('EquilibriumWallFlat.xml')).getroot()
    conditions=root.find('CONDITIONS')
    for child in list(conditions):
        if child.tag!='VARIABLES':conditions.remove(child)
    ET.SubElement(conditions,'SOLVERINFO');ET.SubElement(conditions,'PARAMETERS')
    scheme=ET.SubElement(conditions,'TIMEINTEGRATIONSCHEME')
    ET.SubElement(scheme,'METHOD').text='IMEX';ET.SubElement(scheme,'ORDER').text='1'
    tensor_mesh(root,Z,3,2,10.,1.,points=[8,8],xknots=X)
    for node in root.findall('./EXPANSIONS/E'):node.set('FIELDS',node.get('FIELDS').replace('eps_tilde','epsilon'))
    for node in conditions.findall('./VARIABLES/V'):
        if node.text=='eps_tilde':node.text='epsilon'
    for key,value in dict(EQTYPE='UnsteadyNavierStokes',SolverType='RANSEquilibriumWallEvolution',
        Projection='Galerkin',GlobalSysSoln='DirectFull',TurbulenceModel='WJHighRe',
        RANSWJGradientPolicy='DevStrainRawWorkV1',RANSWallFormulation='WJCalibratedBudgetAWF',
        RANSWallNativeGeometry='BudgetTerrainV1',RANSWallMapContinuation='CompositeTunnelV1',
        RANSWallMap='Plane' if kind=='plane' else 'CosineHill',RANSWallTemporalCoupling='PimpleSingleSweep',
        RANSWallNonlinearPolicy='Legacy',RANSWallVariant='Original',RANSWallInnerStress='Linear',
        RANSWallDestruction='Original',RANSWallTransition='Fixed',
        RANSWallKSampling='InterfaceTrace',RANSWallPhysicalRangeGuard='WJNativeV1',
        RANSWallCondense='False',RANSWallReference='Functions',RANSWallFunctionCoordinates='Chart',
        RANSWallEvaluationReuse='Off',RANSWallPimpleVerify='Off',RANSWallComparison='Off',
        RANSWallPimpleKBound='MassProjection',RANSWallPimpleKInventory='WJNativeV1',
        RANSWallPimpleKBoundSolver='RetainedPool',RANSWallPimpleSnapshots='On').items():set_info(root,key,value)
    # No production inlet/outlet is selected by this controlled native probe.
    for key,value in dict(Time=0.,TimeStep=dt,NumSteps=steps,Kinvis=NU,RANSWallHeight=HW,
        RANSWallRoughness=provisional_roughness(),RANSWallSampleHeight=HW,RANSWallTailStart=HW,RANSWallChi=0.,
        RANSWallFrictionVelocity=US,RANSWallMapDisplacement=.0875,RANSWallMapTunnelTop=20.,
        RANSWallMapAmplitude=0. if kind=='zero' else 1.,RANSWallMapLengthX=2.5,
        RANSWallMapSlopeX=.2 if kind=='plane' else 0.,RANSWallEta=0.,RANSWallTangentialOrder=16,
        RANSWallNormalOrder=24,RANSWallPimpleKBoundMaxSeconds=60.,RANSWallPimpleKFloorRelative=1e-8).items():set_parameter(root,key,value)
    path=folder/'case.xml';write(root,path);return path
