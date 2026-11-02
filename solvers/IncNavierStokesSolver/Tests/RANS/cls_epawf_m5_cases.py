"""Canonical-chart native-CG sessions for the approved M5 geometry ladder.

All curved runs use eta=chi=0.  Geometry maps are analytic; .chart.fld stores
chart coefficients, and solver VTU output stores actual Cartesian positions.
The scalar MMS has U=0, spatially constant k(t), pi=c*z, and outer E=A(t)/r.
Its diffusion divergence is computed analytically in normal coordinates. A
separate circular MMS includes layer swirl, its centrifugal acceleration and
full rotating-frame Cartesian strain. Sources are strong prescribed functions,
never assembled discrete loads. Solver execution belongs on Linux.
"""
from functools import lru_cache
from pathlib import Path
import math
import xml.etree.ElementTree as ET
from cls_epawf_m4_cases import tensor_mesh, _evolution, neutral_session
from cls_epawf_m3_cases import set_info, set_parameter, write

FAMILIES = {
    'plane': dict(kind='Plane', SlopeX=.25),
    'circle': dict(kind='Cylinder', Radius=12.),
    'smooth': dict(kind='Gaussian', Amplitude=.025, LengthX=2., LengthY=2.),
    'variable': dict(kind='Gaussian', Amplitude=.025, LengthX=2., LengthY=2., HeightSlopeX=.04, HeightAmplitude=.025, HeightWaveX=1.3),
    'surface3d': dict(kind='Gaussian', Amplitude=.025, LengthX=2., LengthY=1.7, HeightSlopeX=.025, HeightSlopeY=-.015, HeightAmplitude=.02, HeightWaveX=1.3, HeightWaveY=1.1),
    'feet': dict(kind='CosineHill', Amplitude=.02, LengthX=1.),
}


def configure_map(root, family, dimension):
    spec = FAMILIES[family]
    set_info(root, 'RANSWallMap', spec['kind'])
    set_info(root, 'RANSWallFunctionCoordinates', 'Chart')
    for key, value in spec.items():
        if key != 'kind': set_parameter(root, 'RANSWallMap'+key, repr(value))
    if dimension == 2 and spec.get('HeightSlopeY', 0):
        raise ValueError('3D surface requires dimension=3')


@lru_cache(None)
def expressions(family, dimension, stationary=False, swirl=False):
    import sympy as sp
    from sympy.printing.str import StrPrinter
    class Printer(StrPrinter):
        def _print_Pow(self, e):
            a,b=e.as_base_exp(); return '('+self._print(a)+')^('+self._print(b)+')'
        def _print_Abs(self,e): return 'fabs('+self._print(e.args[0])+')'
        def _print_Max(self,e):
            args=list(map(self._print,e.args))
            value=args[0]
            for a in args[1:]: value='max('+value+','+a+')'
            return value
    printer=Printer()
    x,y,z,t=sp.symbols('x y z t', real=True)
    rho=y if dimension==2 else z
    q2=sp.Integer(0) if dimension==2 else y
    H0=sp.Integer(2); spec=FAMILIES[family]
    H=H0+sp.Float(spec.get('HeightSlopeX',0))*x+sp.Float(spec.get('HeightSlopeY',0))*q2
    H+=sp.Float(spec.get('HeightAmplitude',0))*sp.sin(sp.Float(spec.get('HeightWaveX',1))*x)*sp.cos(sp.Float(spec.get('HeightWaveY',1))*q2)
    r=rho*H/H0
    if spec['kind']=='Cylinder':
        R=sp.Float(spec['Radius']); angle=x/R
        n=sp.Matrix([sp.sin(angle),0,sp.cos(angle)])
        tangent=sp.Matrix([sp.cos(angle),0,-sp.sin(angle)])
        physical_z=(R+r)*sp.cos(angle)-R
        trBr=1/(R+r)
    else:
        if spec['kind']=='Plane': f=sp.Float(spec.get('SlopeX',0))*x+sp.Float(spec.get('SlopeY',0))*q2
        elif spec['kind']=='CosineHill':
            # Differentiate the smooth interior jet first. The curvature has
            # one-sided limits at explicit native faces; it is zero outside.
            L=sp.Float(spec['LengthX'])
            f=sp.Float(spec['Amplitude'])/2*(1+sp.cos(sp.pi*x/L))
        else:
            f=sp.Float(spec['Amplitude'])*sp.exp(-(x/sp.Float(spec['LengthX']))**2-(q2/sp.Float(spec.get('LengthY',1)))**2)
        fx=sp.diff(f,x); fy=sp.diff(f,y) if dimension==3 else 0
        fxx=sp.diff(fx,x); fxy=sp.diff(fx,y) if dimension==3 else 0
        fyy=sp.diff(fy,y) if dimension==3 else 0
        g=1+fx**2+fy**2
        n=sp.Matrix([-fx,-fy,1])/sp.sqrt(g)
        physical_z=f+r/sp.sqrt(g)
        # B has the opposite sign to the upward graph second fundamental form.
        trB=-((1+fy**2)*fxx-2*fx*fy*fxy+(1+fx**2)*fyy)/g**sp.Rational(3,2)
        detB=(fxx*fyy-fxy**2)/g**2
        trBr=(trB+2*r*detB)/(1+r*trB+r*r*detB)
        if spec['kind']=='CosineHill':
            trBr*=(1-sp.sign(sp.Abs(x)-L))/2
    K0=sp.Symbol('K0', positive=True)
    k=K0 if stationary else K0*(1+sp.sin(t)/50)
    kd=sp.diff(k,t); A=k**sp.Rational(3,2)/sp.Rational(255,100)
    E=A/r
    cmu=sp.Rational(3,10)*(1-sp.exp(-sp.Rational(9,25)))
    diffusion=cmu*k*k/sp.Rational(13,10)*(1/r**2-trBr/r)
    ratio=sp.Rational(255,250)
    yap=sp.Rational(83,100)*E*E/k*(ratio-1)*ratio**2
    source_e=sp.diff(E,t)-diffusion+sp.Rational(48,25)*E*E/k-yap
    nu=sp.Symbol('Kinvis', positive=True); rough=sp.Symbol('RANSWallRoughness', positive=True)
    X=rough*sp.sqrt(k)/(70*nu)
    m=sp.Max(sp.Rational(1,2)-sp.Rational(2,5)*X**sp.Rational(7,10),1-sp.Rational(79,100)*X**sp.Rational(-7,25))
    zv=nu*sp.Rational(107,10)*(1-X**m)/sp.sqrt(k)
    a=sp.Rational(2295,10000)*sp.sqrt(k)
    ze=sp.Rational(51,10)*nu/sp.sqrt(k)
    ea=A/sp.Max(r,ze)
    pressure_scale=sp.Integer(0) if family=='feet' else sp.Rational(1,500)
    pi=pressure_scale*physical_z
    ui=sp.zeros(3,1); bi=sp.Matrix([0,0,pressure_scale]); production=0
    if swirl:
        if family!='circle' or dimension!=2: raise ValueError('swirl MMS is an extruded constant-H circular case')
        # Differentiate at fixed physical normal distance before substituting r.
        rr=sp.Symbol('rr', positive=True); ss=rr/H0
        speed=sp.Rational(3,10)*(1+sp.sin(t)/10)*ss**2*(1-ss)**2
        shear=sp.diff(speed,rr)-speed/(R+rr)
        nt=a*(rr-zv); D=nu+nt
        ft=sp.diff(speed,t)-sp.diff(D*shear,rr)-2*D*shear/(R+rr)
        fn=-speed**2/(R+rr)
        ui=(speed*tangent).subs(rr,r)
        bi+=(ft*tangent+fn*n).subs(rr,r)
        production=(nt*shear**2).subs(rr,r)
    names=['u','v'] if dimension==2 else ['u','v','w']
    axes=[0,2] if dimension==2 else [0,1,2]
    inner={name:ui[axis] for name,axis in zip(names,axes)}
    outer={name:sp.Integer(0) for name in names}
    inner.update(k=k,p=pi);outer.update(k=k,eps_tilde=E,p=pi)
    bodies_i={name:bi[axis] for name,axis in zip(names,axes)}
    bodies_o={name:pressure_scale if axis==2 else 0 for name,axis in zip(names,axes)}
    functions={}
    for side,values in [('Inner',inner),('Outer',outer)]:
        for prefix in ['InitialConditions','ExactSolution']:
            functions[prefix+side]=values
        functions['BoundaryValues'+side]={key:val for key,val in values.items() if key!='p'}
    functions.update(BodyForceInner=bodies_i,BodyForceOuter=bodies_o,
                     ScalarSourceInner={'k':kd-production+ea},
                     ScalarSourceOuter={'k':kd+E,'eps_tilde':source_e})
    return {name:{var:printer.doprint(sp.sympify(value)) for var,value in data.items()} for name,data in functions.items()}


def mapped_session(template, destination, family='circle', dimension=2, degree=3,
                   steps=0, order=2, dt=.025, time=.3, condensed=False,
                   restart=None, checkpoint=None, ownership_shift=0, split=1,
                   stationary=False, swirl=False, quadrature=12, radial_split=1):
    if family=='surface3d' and dimension!=3: raise ValueError('surface3d needs full 3D')
    root=ET.parse(template).getroot()
    xknots=[i/split for i in range(split+1)]
    if family=='feet':xknots=[-1.5,-1.,0.,1.,1.5]
    knots=[0.,1.,2.,2.1,4.]
    if radial_split>1:
        knots=sorted(set(knots+[a+(b-a)*j/radial_split for a,b in zip(knots[:-1],knots[1:]) for j in range(1,radial_split)]))
    tensor_mesh(root,knots,degree,dimension,1.,1.,xknots=xknots)
    _evolution(root,steps,order,dt,time,condensed,restart,checkpoint)
    set_info(root,'RANSWallReference','Functions');configure_map(root,family,dimension)
    for name,value in [('RANSWallHeight',2.),('RANSWallEta',0),('RANSWallChi',0),
                       ('RANSWallTangentialOrder',quadrature),('RANSWallNormalOrder',16),
                       ('RANSWallOwnershipShift',ownership_shift)]:set_parameter(root,name,value)
    conditions=root.find('CONDITIONS')
    for node in list(conditions.findall('FUNCTION')):conditions.remove(node)
    for name,values in expressions(family,dimension,stationary,swirl).items():
        node=ET.SubElement(conditions,'FUNCTION',{'NAME':name})
        for var,value in values.items():ET.SubElement(node,'E',{'VAR':var,'VALUE':value})
    return write(root,destination)


def plane_reference(template,destination,**kwargs):
    destination=neutral_session(template,destination,**kwargs)
    root=ET.parse(destination).getroot();configure_map(root,'plane',kwargs.get('dimension',2))
    return write(root,destination)
