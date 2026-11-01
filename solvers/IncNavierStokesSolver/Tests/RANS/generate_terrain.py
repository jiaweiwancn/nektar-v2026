#!/usr/bin/env python3
"""Reproducible terrain-wall verification cases; not calibrated hill data."""
import argparse
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
from generate_m6 import generate_abl, original_j0, structured_hex
from generate_m4 import info, parameter
from generate_m3 import element, mesh


def plane(path, model='standard_fixed', angle=0., degree=4, layers=8, nx=4,
          dt=.002, steps=2, legacy=False, augmentation=2., diffusion='HighRe', dimension=2,
          boundary='periodic', height=.008):
    """Rigid rotation of a complete, matched-boundary neutral wall problem."""
    if legacy and angle != 0:
        raise ValueError('The legacy flat adapter cannot be rotated')
    H, P, N, top, length, nu, roughness, star = height, height/2, 2*height, .08, .4, 1.5e-5, .002, .3
    generate_abl(path, model=model, dimension=dimension, degree=degree, layers=layers, nx=nx, width=.032, ny=1,
                 height=H, top=top, length=length, roughness=roughness, dt=dt, steps=steps,
                 boundary='matched', diffusion=diffusion, monitor=False, backend='matrixfree')
    root = ET.parse(path).getroot()
    filters = root.find('FILTERS')
    if filters is not None:
        root.remove(filters)
    standard, log_wall = model.startswith('standard'), model == 'standard'
    q = 10/3 if standard else 3.190675081624288
    kappa = math.sqrt(.48*1.3/q) if standard else .440553905186659
    z0 = .0002
    matching = star/kappa*math.log(H/z0) if log_wall else star**2*original_j0(q*star**2,H,roughness,nu)
    c, s = math.cos(angle), math.sin(angle)
    vertical = 'y' if dimension == 2 else 'z'
    eta = f'({c:.17g}*{vertical}-({s:.17g})*x)'
    velocity = f'(MatchingVelocity+Ustar/Kappa*log({eta}/H))'
    exact = [f'{c:.17g}*{velocity}']+(['0'] if dimension == 3 else [])+[f'{s:.17g}*{velocity}',
             'Qstar*Ustar^2',f'Ustar^3/(Kappa*{eta})','0']
    names = ['u','v']+(['w'] if dimension == 3 else [])+['k','eps' if standard else 'eps_tilde','p']
    for vertex in root.findall('GEOMETRY/VERTEX/V'):
        x, y, z = map(float,vertex.text.split())
        height = y if dimension == 2 else z
        a, b = c*x-s*height, s*x+c*height
        vertex.text = f'{a:.17g} {b:.17g} 0' if dimension == 2 else f'{a:.17g} {y:.17g} {b:.17g}'
    for group in root.findall('CONDITIONS/FUNCTION'):
        if group.attrib['NAME'] in ('InitialConditions','ExactSolution'):
            for entry in group:
                name = entry.attrib.get('VAR')
                if name in names:
                    entry.attrib['VALUE'] = exact[names.index(name)]
                elif 'VALUE' in entry.attrib:
                    entry.attrib['VALUE'] = re.sub(r'\b'+vertical+r'\b',eta,entry.attrib['VALUE'])
    wall_tag = 'RANSLogLaw' if log_wall else ('AWFOriginal' if model.endswith('_original') else 'AWFEquilibrium')
    # Rebuild all four boundary regions: the other three prescribe matched data,
    # so the plane test needs no axis-aligned artificial outlet or roof policy.
    boundaries = root.find('CONDITIONS/BOUNDARYCONDITIONS'); boundaries.clear()
    for b in range(4 if dimension == 2 else 6):
        region = element(boundaries,'REGION',REF=b)
        for i, name in enumerate(names):
            if dimension == 3 and b in (4,5):
                element(region,'P',VAR=name,VALUE='[5]' if b == 4 else '[4]')
            elif boundary == 'periodic' and b in (1,3):
                element(region,'P',VAR=name,VALUE='[3]' if b == 1 else '[1]')
            elif name == 'p':
                element(region,'N',VAR=name,VALUE='0',USERDEFINEDTYPE='RANSPressure')
            elif b != 0:
                element(region,'D',VAR=name,VALUE=exact[i])
            elif i < dimension:
                if legacy and i == dimension-1:
                    element(region,'D',VAR=name,VALUE='0')
                elif legacy and log_wall:
                    element(region,'N',VAR=name,VALUE='0',USERDEFINEDTYPE=wall_tag)
                else:
                    element(region,'R',VAR=name,VALUE='0',PRIMCOEFF='0',USERDEFINEDTYPE=wall_tag)
            elif name == 'k':
                element(region,'N' if log_wall else 'R',VAR=name,VALUE='0',
                        **({'USERDEFINEDTYPE':wall_tag,'PRIMCOEFF':'0'} if not log_wall else {}))
            else:
                element(region,'D',VAR=name,VALUE=exact[i])
    for name, value in {'Ustar':star,'MatchingVelocity':matching,'RANSZ0':z0,
            'Uref':matching+star/kappa*math.log(top/H),
            'RANSStripHeight':H,'RANSStripPointP':P,'RANSStripPointN':N,
            'RANSNormalAugmentation':augmentation,'RANSNormalRelaxation':.5,
            'RANSPicardRelaxation':.1,'RANSPressureRelaxation':.1,
            'RANSNonlinearTolerance':2e-9,'RANSPressureTolerance':5e-8,
            'RANSMaxIterations':4000,'RANSAWFPointDiagnostics':1,'RANSStageAudit':1,
            'RANSReferenceMonitor':0,'IterativeSolverTolerance':1e-14,'IO_CheckSteps':1}.items():
        if name != 'RANSZ0' or log_wall:
            parameter(root,name,value)
    if not legacy:
        info(root,'RANSWallGeometry','Plane')
        parameter(root,'RANSTerrainSlope',math.tan(angle))
    info(root,'RANSAWFFailurePolicy','RequirePreservation')
    ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
    return dict(kind='plane',model=model,angle=angle,height=H,roughness=roughness,dimension=dimension,
                degree=degree,layers=layers,nx=nx,dt=dt,steps=steps,legacy=legacy,boundary=boundary)


def terrain_mesh(root, xs, etas, amplitude, half_length=.1, top=.3, dimension=2, span=.032, normal_n=.016, span_layers=1):
    """Normal columns through N, then a linear blend to a flat outer roof.

    Edge/face geometry interpolation is independent of solution degree. The
    terrain feet are always column boundaries, since the paper's graph is C1.
    """
    N = normal_n
    root.remove(root.find('GEOMETRY'))
    if dimension == 2:
        mesh(root,xs,etas)
    else:
        if span_layers < 1:
            raise ValueError('Use positive span layers')
        structured_hex(root,xs,[span*j/span_layers for j in range(span_layers+1)],etas)
    geometry = root.find('GEOMETRY')
    seeds = {int(v.attrib['ID']): tuple(map(float,v.text.split())) for v in geometry.find('VERTEX')}
    def mapping(seed):
        x = seed[0]
        eta = seed[1] if dimension == 2 else seed[2]
        lateral = 0 if dimension == 2 else seed[1]
        h = first = 0.
        if abs(x) < half_length:
            a = math.pi*x/half_length
            h = amplitude*.5*(1+math.cos(a))
            first = -amplitude*math.pi/(2*half_length)*math.sin(a)
        norm = math.hypot(1,first)
        nx, nz = -first/norm, 1/norm
        if eta <= N:
            px, pz = x+eta*nx, h+eta*nz
        else:
            r = (eta-N)/(top-N)
            px = (1-r)*(x+N*nx)+r*x
            pz = (1-r)*(h+N*nz)+r*top
        return (px,pz,0) if dimension == 2 else (px,lateral,pz)
    for vertex in geometry.find('VERTEX'):
        vertex.text = ' '.join(f'{v:.17g}' for v in mapping(seeds[int(vertex.attrib['ID'])]))
    curves = element(geometry,'CURVED')
    edges = {int(e.attrib['ID']):tuple(map(int,e.text.split())) for e in geometry.find('EDGE')}
    points = 12
    edge_count = 0
    for eid, pair in edges.items():
        a, b = (seeds[v] for v in pair)
        if dimension == 2 and a[0] == b[0]:
            continue  # Normal/span edges are straight within each mesh layer.
        coordinates = []
        for i in range(points):
            r = i/(points-1)
            coordinates += mapping(tuple((1-r)*v+r*w for v,w in zip(a,b)))
        element(curves,'E',' '.join(f'{v:.17g}' for v in coordinates),ID=edge_count,
                EDGEID=eid,NUMPOINTS=points,TYPE='PolyEvenlySpaced')
        edge_count += 1
    if dimension == 3:
        for face in geometry.find('FACE'):
            eids = list(map(int,face.text.split()))
            e0, e1, e2 = (edges[e] for e in eids[:3])
            v0, v1 = e0 if e0[1] in e1 else e0[::-1]
            v2 = next(v for v in e1 if v != v1)
            v3 = next(v for v in e2 if v != v2)
            a,b,c,d = (seeds[v] for v in (v0,v1,v2,v3))
            coordinates = []
            for j in range(points):
                s = j/(points-1)
                for i in range(points):
                    r = i/(points-1)
                    seed = tuple((1-r)*(1-s)*a[k]+r*(1-s)*b[k]+r*s*c[k]+(1-r)*s*d[k] for k in range(3))
                    coordinates += mapping(seed)
            element(curves,'F',' '.join(f'{v:.17g}' for v in coordinates),ID=int(face.attrib['ID']),
                    FACEID=int(face.attrib['ID']),NUMPOINTS=points*points,TYPE='PolyEvenlySpaced')


def hill(path, model='standard_fixed', degree=4, layers=6, refinement=1, dt=.0002,
         steps=2, amplitude=.04, width_ratio=2., empty=False, dimension=2, properties=False,
         height=.008, diffusion='HighRe', periodic=False):
    """Paper-shaped short numerical pilot and its identical empty-floor control."""
    if refinement < 1 or layers < 3 or degree < 3:
        raise ValueError('Use refinement>=1, layers>=3 and degree>=3')
    plane(path,model=model,degree=degree,layers=layers,nx=2,dt=dt,steps=steps,dimension=dimension,
          boundary='matched',height=height,diffusion=diffusion)
    root = ET.parse(path).getroot()
    top, H, N, half_length = .3,height,2*height,.1
    hill_height = 0. if empty else amplitude
    # Fixed downstream ramp endpoint, including identical mesh coordinates for
    # every width. Subdivide all segments for an independent h-refinement study.
    xs = [-.8,-.68,-.56,-.52,-.48,-.44,-.40,-.36,-.32,-.28,-.24,-.20,-.16,-.12,-.1]
    xs += [-.075,-.05,-.025,0.,.025,.05,.075,.1,.14,.2,.3,.45,.6,.8]
    if periodic:
        # A full cosine period has matching one-sided curvature at the seam.
        # It separates smooth-wall consistency from the compact C1 feet.
        xs = [x for x in xs if -.1 <= x <= .1]
    xs = [a+(b-a)*i/refinement for a,b in zip(xs,xs[1:]) for i in range(refinement)]+[xs[-1]]
    etas = [H+(N-H)*i/refinement for i in range(refinement)]
    outer = [N*(top/N)**(i/(layers-1)) for i in range(layers)]
    # Refine the complete resolved domain. Keeping the outer cells fixed
    # leaves the pressure/initial-projection error unchanged and is not an
    # h-convergence test of the coupled strip/resolved problem.
    etas += [a+(b-a)*i/refinement for a,b in zip(outer,outer[1:]) for i in range(refinement)]+[outer[-1]]
    terrain_mesh(root,xs,etas,hill_height,half_length,top,dimension,normal_n=N)
    standard, log_wall = model.startswith('standard'), model == 'standard'
    names = ['u','v']+(['w'] if dimension == 3 else [])+['k','eps' if standard else 'eps_tilde','p']
    vertical = 'y' if dimension == 2 else 'z'
    # This unary form also avoids a nested binary-function expression in the
    # native vector interpreter. Native initial-field receipts verify it.
    a = '((abs(x/0.1+1)-abs(x/0.1-1))/2)'
    surface = f'({hill_height:.17g}/2*(1+cos(PI*{a})))'
    slope = f'(-{hill_height:.17g}*PI/0.2*sin(PI*{a}))'
    # A streamfunction F(eta) gives a divergence-free start, tangent to the
    # physical terrain and with zero normal velocity and unchanged speed at
    # the flat outer roof. This is a pilot initial condition, not measured data.
    gap = f'({top:.17g}-{surface})'
    eta = f'({vertical}-{surface}*(({top:.17g}-{vertical})/{gap})^2)'
    speed = f'(MatchingVelocity+Ustar/Kappa*log({eta}/H))'
    u = f'{speed}*(1+2*{surface}*({top:.17g}-{vertical})/{gap}^2)'
    v = f'{speed}*{slope}*({top:.17g}-{vertical})^2*({top:.17g}+{surface})/{gap}^3'
    exact = [u]+(['0'] if dimension == 3 else [])+[v,'Qstar*Ustar^2',f'Ustar^3/(Kappa*{eta})','0']
    conditions = root.find('CONDITIONS')
    for group in list(conditions.findall('FUNCTION')):
        if group.attrib.get('NAME') in ('InitialConditions','ExactSolution'):
            conditions.remove(group)
    for label in ('InitialConditions',):
        group = element(conditions,'FUNCTION',NAME=label)
        for name,value in zip(names,exact):
            element(group,'E',VAR=name,VALUE=value)
    for region in root.findall('CONDITIONS/BOUNDARYCONDITIONS/REGION'):
        if periodic and region.attrib['REF'] in ('1','3'):
            ref = region.attrib['REF']
            region.clear();region.attrib['REF'] = ref
            for name in names:
                element(region,'P',VAR=name,VALUE='[3]' if ref == '1' else '[1]')
            continue
        if dimension == 3 and region.attrib['REF'] in ('4','5'):
            continue
        if region.attrib['REF'] == '1':
            region.clear(); region.attrib['REF'] = '1'
            for i,name in enumerate(names):
                if name == 'p':
                    element(region,'N',VAR=name,VALUE='0',USERDEFINEDTYPE='RANSPressure')
                elif i < dimension:
                    traction = ('0' if standard else '-(.92-2/3)*Qstar*Ustar^2') if i == 0 else ('Ustar^2' if i == dimension-1 else '0')
                    element(region,'N',VAR=name,VALUE=traction,USERDEFINEDTYPE='RANSOutflow')
                else:
                    element(region,'N',VAR=name,VALUE='0')
        elif region.attrib['REF'] != '0':
            for entry in region:
                if entry.tag == 'D':
                    entry.attrib['VALUE'] = exact[names.index(entry.attrib['VAR'])]
        else:
            for entry in region:
                if entry.attrib['VAR'] == names[-2]:
                    entry.attrib['VALUE'] = 'Ustar^3/(Kappa*H)'
    info(root,'RANSWallGeometry','Plane' if empty else 'CosineHill')
    for key,value in {'RANSTerrainHeight':hill_height,'RANSTerrainHalfLength':half_length,
            'RANSTerrainSlope':0,'L':.2 if periodic else 1.6,'Z':top,'RANSStageAudit':0,'RANSWallGeometryTolerance':1e-7,
            'RANSPicardRelaxation':.1,'RANSPressureRelaxation':.1,
            'RANSNormalAugmentation':2,'RANSNormalRelaxation':.5,
            'RANSNonlinearTolerance':1e-8,'RANSPressureTolerance':1e-8,
            'IO_CheckSteps':max(1,steps//4)}.items():
        parameter(root,key,value)
    q = 10/3 if standard else 3.190675081624288
    kappa = math.sqrt(.48*1.3/q) if standard else .440553905186659
    matching = .3/kappa*math.log(H/.0002) if log_wall else .3**2*original_j0(q*.3**2,H,.002,1.5e-5)
    parameter(root,'Uref',matching+.3/kappa*math.log(top/H))
    if model.endswith(('_fixed','_switched')):
        info(root,'RANSAWFSpatialWeight','Ramp')
        parameter(root,'RANSAWFTransitionStart',-.2-width_ratio*.04)
        parameter(root,'RANSAWFTransitionEnd',-.2)
    if properties:
        group = element(conditions,'FUNCTION',NAME='RANSWallProperties')
        element(group,'E',VAR='z0' if log_wall else 'roughness',
                VALUE='0.0002+0.0001*(1+tanh(x/0.1))/2' if log_wall else '0.002+0.001*(1+tanh(x/0.1))/2')
    ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
    return dict(kind='hill' if not empty else 'paired-floor',model=model,degree=degree,layers=layers,
                refinement=refinement,dt=dt,steps=steps,amplitude=hill_height,width_ratio=width_ratio,
                dimension=dimension,height=H,length=.2 if periodic else 1.6,periodic=periodic)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('session',type=Path)
    p.add_argument('--kind',choices=('plane','hill'),default='plane')
    p.add_argument('--model',default='standard_fixed',choices=(
        'standard','standard_original','standard_fixed','standard_switched',
        'cls_original','cls_fixed','cls_switched'))
    p.add_argument('--angle',type=float,default=0.)
    p.add_argument('--degree',type=int,default=4)
    p.add_argument('--layers',type=int,default=8)
    p.add_argument('--nx',type=int,default=4)
    p.add_argument('--dt',type=float)
    p.add_argument('--steps',type=int,default=2)
    p.add_argument('--legacy',action='store_true')
    p.add_argument('--dimension',type=int,choices=(2,3),default=2)
    p.add_argument('--boundary',choices=('periodic','matched'),default='periodic')
    p.add_argument('--height',type=float,default=.008)
    p.add_argument('--diffusion',choices=('HighRe','Molecular'),default='HighRe')
    p.add_argument('--amplitude',type=float,default=.04)
    p.add_argument('--refinement',type=int,default=1)
    p.add_argument('--width-ratio',type=float,default=2.)
    p.add_argument('--empty',action='store_true')
    p.add_argument('--properties',action='store_true')
    args = p.parse_args()
    args.session.parent.mkdir(parents=True,exist_ok=True)
    values = vars(args).copy(); session = values.pop('session'); kind = values.pop('kind')
    if values['dt'] is None:
        values.pop('dt')
    for name in (('angle','nx','legacy','boundary') if kind == 'hill' else
                 ('amplitude','refinement','width_ratio','empty','properties')):
        values.pop(name)
    print((hill if kind == 'hill' else plane)(session,**values))


if __name__ == '__main__':
    main()
