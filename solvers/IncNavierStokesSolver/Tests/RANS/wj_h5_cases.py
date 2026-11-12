#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Zero-step physical-Functions sessions for WJNeutralBudgetReferenceV1."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from wj_h5_neutral import NeutralReference,mean_fit,measurements,SCALE_MM
from wj_h2_cases import analytic_expressions,inner_normal_knots
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_info,set_parameter,write


def mesh_knots(c):
    inner=inner_normal_knots(c.H,c.w,c.h)
    # Data heights only define observation/mesh panels, never an old U/k tail.
    outer=sorted(set([c.H,12.5-c.d0,c.top-c.d0]+[
        r['z_mm']/SCALE_MM-c.d0 for r in measurements(SCALE_MM*c.d0)if r['panel']=='a'and c.d0+c.H<r['z_mm']/SCALE_MM<c.top]))
    rho=lambda r:c.H+(r-c.H)*(c.top-c.H)/(c.top-c.d0-c.H)
    chart=inner+[rho(r)for r in outer[1:]]
    if any(a>=b for a,b in zip(chart[:-1],chart[1:]))or chart.count(c.H)!=1:
        raise ValueError('invalid low-H interface/mesh ordering')
    return inner,outer,chart


def case_xml(folder,c):
    if not isinstance(c,NeutralReference):raise TypeError('distinct neutral configuration required; no Completion adapter')
    inner,outer,oi,oo,top,tf,meta=analytic_expressions(c.H,c.h,u_tau=c.u_tau,nu=c.nu,d0=c.d0,top=c.top,z0_e=c.z0_e)
    meta.pop('wall');ik,physical_r,z=mesh_knots(c);folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    root=ET.parse(Path(__file__).with_name('EquilibriumWallFlat.xml')).getroot();co=root.find('CONDITIONS')
    for child in list(co):
        if child.tag!='VARIABLES':co.remove(child)
    ET.SubElement(co,'SOLVERINFO');ET.SubElement(co,'PARAMETERS')
    scheme=ET.SubElement(co,'TIMEINTEGRATIONSCHEME');ET.SubElement(scheme,'METHOD').text='IMEX';ET.SubElement(scheme,'ORDER').text='1'
    tensor_mesh(root,z,3,2,8.,1.,points=[8,8],xknots=[0.,8.])
    for e in root.findall('./EXPANSIONS/E'):e.set('FIELDS',e.get('FIELDS').replace('eps_tilde','epsilon'))
    for e in co.findall('./VARIABLES/V'):
        if e.text=='eps_tilde':e.text='epsilon'
    info=dict(EQTYPE='UnsteadyNavierStokes',SolverType='RANSEquilibriumWallEvolution',Projection='Galerkin',GlobalSysSoln='DirectFull',
        TurbulenceModel='WJHighRe',RANSWJGradientPolicy='DevStrainRawWorkV1',RANSWallFormulation='WJCalibratedBudgetAWF',
        RANSWallNativeGeometry='BudgetTerrainV1',RANSWallMapContinuation='CompositeTunnelV1',RANSWallMap='Plane',
        RANSWallTemporalCoupling='PimpleSingleSweep',RANSWallNonlinearPolicy='Legacy',RANSWallVariant='Original',
        RANSWallInnerStress='Linear',RANSWallDestruction='Original',RANSWallTransition='Fixed',RANSWallKSampling='InterfaceTrace',
        RANSWallPhysicalRangeGuard='WJNativeV1',RANSWallCondense='False',RANSWallReference='Functions',RANSWallFunctionCoordinates='Physical',
        RANSWallEvaluationReuse='ExactStepV1',RANSWallPimpleVerify='Off',RANSWallComparison='Off',RANSWallBoundaryMode='OpenTunnelV1',
        RANSWallOpenTop='ReferenceTractionV1',RANSWallPimpleStopOnBound='On',RANSWallPimpleKBound='MassProjection',
        RANSWallPimpleKInventory='WJNativeV1',RANSWallPimpleKBoundSolver='RetainedPool',RANSWallPimpleSnapshots='On',RANSWallCheckpoint=str(folder/'CP'))
    for key,value in info.items():set_info(root,key,value)
    params=dict(Time=0.,TimeStep=.01,NumSteps=0,IO_CheckSteps=1,RANSWallDiagnostics=1,Kinvis=c.nu,RANSWallHeight=c.H,
        RANSWallRoughness=c.h,RANSWallSampleHeight=c.H,RANSWallTailStart=c.H,RANSWallChi=0.,RANSWallFrictionVelocity=c.u_tau,
        RANSWallMapDisplacement=c.d0,RANSWallMapTunnelTop=c.top,RANSWallMapAmplitude=0.,RANSWallMapLengthX=2.5,
        RANSWallEta=0.,RANSWallTangentialOrder=12,RANSWallNormalOrder=16,RANSWallPimpleKBoundMaxSeconds=60.,RANSWallPimpleKFloorRelative=1e-8)
    for key,value in params.items():set_parameter(root,key,value)
    for family in ('InitialConditions','BoundaryValues'):
        for side,values in (('Inner',inner),('Outer',outer)):
            f=ET.SubElement(co,'FUNCTION',{'NAME':family+side})
            for variable,value in values.items():
                if family=='BoundaryValues'and variable=='p':continue
                ET.SubElement(f,'E',{'VAR':variable,'VALUE':value})
    for name,values in (('OutletTractionInner',oi),('OutletTractionOuter',oo),('TopTractionOuter',top),('TopDiffusiveFluxOuter',tf)):
        f=ET.SubElement(co,'FUNCTION',{'NAME':name})
        for variable,value in values.items():ET.SubElement(f,'E',{'VAR':variable,'VALUE':value})
    path=folder/'case.xml';write(root,path)
    meta.update(reference=c.summary(),inner_knots=ik,physical_outer_r=physical_r,chart_knots=z,
        physical_knots=[c.d0+r for r in ik]+[c.d0+r for r in physical_r[1:]],
        xknots=[0.,8.],normal_cells=len(z)-1,steps=0,scope='IDEALIZED EQUILIBRIUM CONTROL; not finite-depth experimental inlet')
    (folder/'case-design.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n');return path


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--height-mm',type=float,choices=(5.,5.5,6.),required=True)
    a=p.parse_args();print(case_xml(a.output,NeutralReference(a.height_mm/SCALE_MM,mean_fit()['z0_e'])))
