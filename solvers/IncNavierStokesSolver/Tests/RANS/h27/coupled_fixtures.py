#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Create opt-in standard-main H27 coupled fixtures on Linux, create-only."""
import argparse,hashlib,json,sys,collections
from pathlib import Path
import xml.etree.ElementTree as E
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'h26'))
from fixtures import mesh

def generate(path,table,prefix,mode='Manufactured',nx=2,nr=2,dt=.02,steps=1,
             quadrature_x=8,quadrature_r=12,restart='',recovery=0,science=0,info=1,
             linear=1e-12,nonlinear=1e-10,retry=-1,coefficients=True,template=None,
             backend='Historical',mms_body_amplitude=0,mms_variant='Historical',
             solve='Advance',steady_seconds=0,steady_pc_mass=0,steady_seed=0,
             anderson_depth=0,nonlinear_trace=False,maximum_iterations=60,initial_fields='',initial_source_identity=''):
    path=Path(path);assert not path.exists() and mode in ('Manufactured','Physical')
    if template:
        import shutil
        shutil.copy2(template,path)
    else:
        assert mode=='Manufactured'
        xs=[i/nx for i in range(nx+1)]
        rs=[.15*i/nr for i in range(nr)]+[.15+.35*i/nr for i in range(nr+1)]
        mesh(path,xs,rs,quadrature=12)
    root=E.parse(path).getroot()
    if not template:
        for vertex in root.find('GEOMETRY/VERTEX'):
            x,z,d=map(float,vertex.text.split());vertex.text=f'{x:.17g} {z+.0875:.17g} {d:.17g}'
    expansions=root.find('EXPANSIONS');original=dict(expansions[0].attrib)
    for child in list(expansions):expansions.remove(child)
    for names,modes in [('u,v,k,epsilon',4),('p',3)]:
        attrs=dict(original,FIELDS=names,NUMMODES=f'{modes},{modes}')
        E.SubElement(expansions,'E',attrs)
    conditions=root.find('CONDITIONS');si=conditions.find('SOLVERINFO')
    settings=dict(SolverType='WJNativeCoupledV2',Driver='Standard',
                  WJNativeMode=mode,WJReferenceTable=str(table),WJOutputPrefix=str(prefix),
                  WJOutputCoefficients='True'if coefficients else'False',Projection='Continuous',
                  WJConstitutiveBackend=backend,WJMMSVariant=mms_variant,WJNativeSolve=solve,
                  WJNonlinearTrace='True'if nonlinear_trace else'False')
    if restart:settings['WJRestart']=str(restart)
    if initial_fields:settings['WJInitialFields']=str(initial_fields)
    if initial_source_identity:settings['WJInitialSourceIdentity']=str(initial_source_identity)
    for key in ('WJRestart','WJInitialFields','WJInitialSourceIdentity'):
        if key not in settings:
            for old in list(si.findall(f"I[@PROPERTY='{key}']")):si.remove(old)
    for key,value in settings.items():
        item=si.find(f"I[@PROPERTY='{key}']")
        if item is None:item=E.SubElement(si,'I',PROPERTY=key)
        item.set('VALUE',value)
    parameters=conditions.find('PARAMETERS')
    values=dict(Time=0,TimeStep=dt,NumSteps=steps,Kinvis=6.423173803526447e-5,
        IO_CheckSteps=recovery,IO_InfoSteps=info,WJScientificSteps=science,
        WJLinearTolerance=linear,WJCoupledTolerance=nonlinear,WJMaximumIterations=maximum_iterations,
        WJQuadratureX=quadrature_x,WJQuadratureR=quadrature_r,WJPhysicalDisplacement=.0875,
        WJInjectRetryStep=retry,NumQuadPointsError=0,WJMMSBodyAmplitude=mms_body_amplitude,
        WJSteadySeconds=steady_seconds,WJSteadyPreconditionerMass=steady_pc_mass,
        WJSteadySeed=steady_seed,WJAndersonDepth=anderson_depth)
    for child in list(parameters):
        if child.text.split('=')[0].strip() in values:parameters.remove(child)
    for key,value in values.items():E.SubElement(parameters,'P').text=f'{key} = {value}'
    variables=conditions.find('VARIABLES');variables.clear()
    for i,name in enumerate(('u','v','k','epsilon','p')):E.SubElement(variables,'V',ID=str(i)).text=name
    # Identify physical exterior sides from native edge/vertex IDs. A retained
    # template's boundary-region numbering is not a physical boundary law.
    geometry=root.find('GEOMETRY');vertices={int(v.attrib['ID']):tuple(map(float,v.text.split()))for v in geometry.find('VERTEX')}
    edges={int(e.attrib['ID']):tuple(map(int,e.text.split()))for e in geometry.find('EDGE')}
    owners=collections.Counter(e for element in geometry.find('ELEMENT')for e in map(int,element.text.split()))
    left=min(v[0]for v in vertices.values());right=max(v[0]for v in vertices.values())
    bottom=min(v[1]for v in vertices.values());top=max(v[1]for v in vertices.values());side_edges=[[]for _ in range(4)]
    for edge,count in sorted(owners.items()):
        if count!=1:continue
        endpoints=[vertices[v]for v in edges[edge]]
        tests=[all(abs(v[0]-left)<1e-11 for v in endpoints),all(abs(v[0]-right)<1e-11 for v in endpoints),
               all(abs(v[1]-bottom)<1e-11 for v in endpoints),all(abs(v[1]-top)<1e-11 for v in endpoints)]
        assert sum(tests)==1,'unsupported exterior edge geometry';side_edges[tests.index(True)].append(edge)
    composite=geometry.find('COMPOSITE');first=max(int(c.attrib['ID'])for c in composite)+1
    regions=conditions.find('BOUNDARYREGIONS');regions.clear()
    boundary=conditions.find('BOUNDARYCONDITIONS');boundary.clear()
    for bid,ids in enumerate(side_edges):
        assert ids
        E.SubElement(composite,'C',ID=str(first+bid)).text='E['+','.join(map(str,ids))+']'
        E.SubElement(regions,'B',ID=str(bid)).text=f'C[{first+bid}]'
        region=E.SubElement(boundary,'REGION',REF=str(bid))
        for name in ('u','v'):
            essential=mode=='Manufactured'or bid in (0,2)or(bid==3 and name=='v')
            E.SubElement(region,'D'if essential else'N',VAR=name,VALUE='0')
        for name in ('k','epsilon'):
            E.SubElement(region,'D'if bid==0 else'N',VAR=name,VALUE='0')
        E.SubElement(region,'N',VAR='p',VALUE='0')
    initial=conditions.find("FUNCTION[@NAME='InitialConditions']")
    if initial is None:initial=E.SubElement(conditions,'FUNCTION',NAME='InitialConditions')
    initial.clear();initial.set('NAME','InitialConditions')
    for name in ('u','v','k','epsilon','p'):E.SubElement(initial,'E',VAR=name,VALUE='0')
    E.indent(root);path.write_text(E.tostring(root,encoding='unicode')+'\n')
    Path(str(path)+'.opt').write_text('<NEKTAR><COLLECTIONS DEFAULT="NoCollection"/></NEKTAR>\n')
    return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),mode=mode,
        table=dict(path=str(table),sha256=hashlib.sha256(Path(table).read_bytes()).hexdigest()),
        nx=nx,nr=nr,dt=dt,steps=steps,template=str(template)if template else None)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--table',type=Path,required=True)
    p.add_argument('--steps',type=int,default=1);p.add_argument('--nx',type=int,default=2);p.add_argument('--nr',type=int,default=2);p.add_argument('--dt',type=float,default=.02)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    result=generate(a.out/'mms-smoke.xml',a.table,a.out/'mms-smoke',nx=a.nx,nr=a.nr,dt=a.dt,steps=a.steps)
    (a.out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(a.out)
