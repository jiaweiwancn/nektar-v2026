#!/usr/bin/env python3
"""Independent high-Re CLS manufactured and original-AWF transport fixtures."""
import argparse
import base64
from functools import lru_cache
import math
from pathlib import Path
import re
import struct
import xml.etree.ElementTree as ET
import zlib
from generate_m3 import element, generate, mesh

MODEL_FIELDS = ("rans_nut", "rans_P", "rans_Yap", "rans_eps_total",
                "rans_R_xx", "rans_R_yy", "rans_R_zz", "rans_R_xy", "rans_R_xz", "rans_R_yz",
                "rans_quadratic_norm", "rans_cubic_norm")


def parameter(root, name, value):
    parameters = root.find("CONDITIONS/PARAMETERS")
    for p in parameters:
        if p.text.split("=")[0].strip() == name:
            p.text = f"{name} = {value}"; return
    element(parameters, "P", f"{name} = {value}")


def info(root, name, value):
    infos = root.find("CONDITIONS/SOLVERINFO")
    for i in infos:
        if i.attrib["PROPERTY"] == name:
            i.attrib["VALUE"] = value; return
    element(infos, "I", PROPERTY=name, VALUE=value)


@lru_cache(maxsize=None)
def manufactured(dimension=2, omit_cubic=False):
    import sympy as s
    x, y, z, t, nu = s.symbols("x y z t ResolvedNu", real=True)
    coords = [x, y] if dimension == 2 else [x, y, z]
    height = coords[-1]
    # A linear velocity activates strain, rotation AND cubic stress without a
    # max-invariant branch crossing. Both k/eps and C(q) vary in space/time.
    a = (1+s.sin(t)/5)/5
    matrix = s.Matrix([[s.Rational(1,5), 1], [-s.Rational(2,5), -s.Rational(1,5)]]) if dimension == 2 else s.Matrix([
        [s.Rational(1,5), 1, s.Rational(3,10)],
        [-s.Rational(2,5), -s.Rational(1,10), s.Rational(1,5)],
        [s.Rational(1,10), -s.Rational(3,10), -s.Rational(1,10)]])
    centered = [x-s.Rational(1,2), y-s.Rational(3,2)] if dimension == 2 else [x-s.Rational(1,2), y-s.Rational(1,2), z-s.Rational(3,2)]
    velocity = list(a*matrix*s.Matrix(centered))
    k = 1+x/10+(height-1)/20+s.sin(t)/50
    eps = s.Rational(3,20)+x/50+(height-1)/100+s.cos(t)/200
    if dimension == 3:
        k += y/25; eps += y/200
    pi = (x-s.Rational(1,2))*(s.Rational(1,5)+s.cos(t)/10)+(height-s.Rational(3,2))*s.sin(t)/20
    gradient = s.zeros(3)
    axes = [0, 2] if dimension == 2 else [0, 1, 2]
    for i in range(dimension):
        for j in range(dimension): gradient[axes[i], axes[j]] = a*matrix[i,j]
    strain, rotation = gradient+gradient.T, gradient-gradient.T
    ss, oo = s.trace(strain*strain.T), s.trace(rotation*rotation.T)
    base_s, base_o = matrix+matrix.T, matrix-matrix.T
    q_constant = s.sqrt(max(s.trace(base_s*base_s.T)/2, s.trace(base_o*base_o.T)/2))
    q = a*k/eps*q_constant
    cmu = s.Rational(3,10)/(1+s.Rational(35,100)*q**s.Rational(3,2))*(1-s.exp(-s.Rational(36,100)*s.exp(s.Rational(75,100)*q)))
    nut = cmu*k*k/eps
    # Direct dimensional matrix law, independent of production C++ evaluation.
    quadratic = k*nut/eps*(-(strain*strain-ss*s.eye(3)/3)/10+
                         (rotation*strain-strain*rotation)/10+
                         s.Rational(26,100)*(rotation*rotation.T-oo*s.eye(3)/3))
    cubic = nut*k*k/eps**2*(-10*cmu**2*(strain*strain*rotation-rotation*strain*strain)+
                          (-5*cmu**2*ss+5*cmu**2*oo)*strain)
    stress = 2*k*s.eye(3)/3-nut*strain+quadratic+cubic
    production = -sum(stress[i,j]*gradient[i,j] for i in range(3) for j in range(3))
    ratio = k**s.Rational(3,2)/(s.Rational(5,2)*eps*height)
    # r>1 throughout these fixtures: no nondifferentiable positive-part branch.
    yap = eps**2/k*s.Rational(83,100)*(ratio-1)*ratio**2
    adv = lambda f: sum(velocity[i]*s.diff(f,c) for i,c in enumerate(coords))
    diff = lambda f,d: sum(s.diff(d*s.diff(f,c),c) for c in coords)
    nonlinear = quadratic if omit_cubic else quadratic+cubic
    forces = [s.diff(velocity[i],t)+adv(velocity[i])+s.diff(pi,coords[i])-
              sum(s.diff((nu+nut)*strain[axes[i],axes[j]]-nonlinear[axes[i],axes[j]],coords[j]) for j in range(dimension)) for i in range(dimension)]
    forces += [s.diff(k,t)+adv(k)-diff(k,nu+nut)-production+eps,
               s.diff(eps,t)+adv(eps)-diff(eps,nu+nut/s.Rational(13,10))-
               s.Rational(144,100)*eps/k*production+s.Rational(192,100)*eps**2/k-yap]
    output = [nut,production,yap,eps]+[stress[i,j] for i,j in ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))]
    output += [s.sqrt(sum(m[i,j]**2 for i in range(3) for j in range(3))) for m in (quadratic,cubic)]
    expression = lambda f: s.sstr(f).replace("**", "^").replace("Abs(", "abs(")
    return [expression(f) for f in velocity+[k,eps,pi]], [expression(f) for f in forces], [expression(f) for f in output]


def hex_mesh(root, heights, length=1., width=1.):
    """Expand the repository's eight-hex fixture; preserve native face ordering."""
    original = ET.fromstring(re.sub(r"<\?xml[^>]*\?>", "", (Path(__file__).parents[1]/"TGV3d_2_2.xml").read_text()))
    source = original.find("GEOMETRY")
    def records(node, fmt):
        raw = zlib.decompress(base64.b64decode(node.text))
        return list(struct.iter_unpack("<"+fmt, raw))
    vertices = {int(v[0]):tuple(round(c*2)/2 for c in v[1:]) for v in records(source.find("VERTEX"), "qddd")}
    edges = {int(v[0]):v[1:] for v in records(source.find("EDGE"), "qqq")}
    faces = {int(v[0]):v[1:] for v in records(source.find("FACE/Q"), "qqqqq")}
    cells = {int(v[0]):v[1:] for v in records(source.find("ELEMENT/H"), "qqqqqqq")}
    geometry = element(root,"GEOMETRY",DIM=3,SPACE=3)
    vertex = element(geometry,"VERTEX")
    for i,(x,y,z) in vertices.items():
        element(vertex,"V",f"{x*length:.17g} {y*width:.17g} {heights[round(2*z)]:.17g}",ID=i)
    for label,tag,items in (("EDGE","E",edges),("FACE","Q",faces),("ELEMENT","H",cells)):
        group=element(geometry,label)
        for i,values in items.items():element(group,tag," ".join(map(str,values)),ID=i)
    sides=[[] for _ in range(6)]
    for f,es in faces.items():
        points=[vertices[v] for e in es for v in edges[e]]
        for b,(d,value) in enumerate(((2,0),(0,1),(2,1),(0,0),(1,0),(1,1))):
            if all(p[d]==value for p in points):sides[b].append(f)
    composites=element(geometry,"COMPOSITE")
    element(composites,"C","H[0-7]",ID=0)
    for b,fs in enumerate(sides,1):element(composites,"C","F["+",".join(map(str,fs))+"]",ID=b)
    element(geometry,"DOMAIN","C[0]")


def awf_reference(k=.5105080130598861, height=2., roughness=.1, nu=1.5e-5):
    """Independent quadrature used to anchor the continuous CLS log profile."""
    from scipy.integrate import quad
    root=math.sqrt(k); c_length=2.55; alpha=.09*c_length
    x=roughness*root/(70*nu)
    intercept=10.7 if roughness==0 else 10.7*(1-x**max(.5-.4*x**.7,1-.79*x**(-.28)))
    zv=nu*intercept/root; ze=2*c_length*nu/root
    nut=lambda z: alpha*root*max(z-zv,0)
    d=lambda z:nu+nut(z)
    points=sorted(set([0.,height]+[p for p in (zv,roughness,ze) if 0<p<height]))
    integral=lambda f:sum(quad(f,a,b,epsabs=1e-12,epsrel=1e-12)[0] for a,b in zip(points,points[1:]))
    return {"j0":integral(lambda z:1/d(z)),"j1":integral(lambda z:max(z-roughness,0)/d(z)),
            "production_weight":integral(lambda z:nut(z)/d(z)**2),
            "production_weight1":integral(lambda z:nut(z)*max(z-roughness,0)/d(z)**2),
            "production_weight2":integral(lambda z:nut(z)*max(z-roughness,0)**2/d(z)**2),
            "epsilon_mean":k**1.5/(c_length*ze) if height<ze else k**1.5/(c_length*height)*(1+math.log(height/ze)),
            "epsilon_p":k**1.5/(c_length*max(1.,ze)),"epsilon_h":k**1.5/(c_length*max(height,ze))}


def generate_cls(path, kind="mms", dimension=2, dt=.025, steps=4, order=2, modes=5,
                 points=12, diffusion="HighRe", omit_cubic=False, roughness=.1,
                 perturbation=0., span_angle=0., physical_nu=1.5e-5):
    generate(path,kind="decay",dt=dt,steps=steps,order=order,modes=modes,points=points,diffusion=diffusion)
    tree=ET.parse(path);root=tree.getroot();c=root.find("CONDITIONS")
    info(root,"TurbulenceModel","CLSHighRe")
    parameter(root,"RANSNonlinearTolerance","1e-10");parameter(root,"RANSMaxIterations",2000)
    names=["u","v"]+(["w"] if dimension==3 else [])+["k","eps_tilde","p"]
    variables=c.find("VARIABLES");variables.clear()
    for i,name in enumerate(names):element(variables,"V",name,ID=i)
    root.remove(root.find("GEOMETRY"))
    if dimension==3:hex_mesh(root,[1.,1.5,2.] if kind=="mms" else [2.,4.,30.],length=1 if kind=="mms" else 300,width=1 if kind=="mms" else 100)
    else:mesh(root,[0,.5,1] if kind=="mms" else [0,150,300],[1,1.5,2] if kind=="mms" else [2,4,8,16,30])
    for expansion,fields,n in zip(root.find("EXPANSIONS"),(",".join(names[:-1]),"p"),(modes,modes-2)):
        expansion.attrib.update(NUMMODES=",".join([str(n)]*dimension),NUMPOINTS=",".join([str(points)]*dimension),
            BASISTYPE=",".join(["Modified_A"]*dimension),POINTSTYPE=",".join(["GaussLobattoLegendre"]*dimension),FIELDS=fields)
    regions=c.find("BOUNDARYREGIONS");regions.clear()
    for b in range(4 if dimension==2 else 6):element(regions,"B",f"C[{b+1}]",ID=b)
    boundaries=c.find("BOUNDARYCONDITIONS");boundaries.clear()
    extras=[];forces=None
    if kind=="mms":
        exact,forces,extras=manufactured(dimension,omit_cubic)
        parameter(root,"Kinvis",.025);parameter(root,"ResolvedNu",0 if diffusion=="HighRe" else .025)
    else:
        q,kappa,star=3.190675081624288,.440553905186659,.4
        k=q*star*star;ref=awf_reference(k,roughness=roughness,nu=physical_nu)
        parameter(root,"Kinvis",physical_nu)
        parameter(root,"ResolvedNu",physical_nu if diffusion=="Molecular" else 0)
        parameter(root,"Ustar",star);parameter(root,"Kappa",kappa);parameter(root,"Qstar",q)
        parameter(root,"MatchingVelocity",star*star*ref["j0"])
        height="y" if dimension==2 else "z"
        u=f"MatchingVelocity+Ustar/Kappa*log({height}/2)"
        if perturbation:u+=f"+({perturbation})*sin(2*PI*x/300)"
        horizontal = [u,"0"] if dimension==2 else [f"({math.cos(span_angle):.17g})*({u})",f"({math.sin(span_angle):.17g})*({u})","0"]
        exact=horizontal+["Qstar*Ustar^2",f"Ustar^3/(Kappa*{height})","0"]
        if kind=="awf":
            info(root,"RANSWallModel","AWFOriginal");parameter(root,"RANSRoughnessHeight",roughness)
            parameter(root,"RANSStripHeight",2);parameter(root,"RANSStripPointP",1);parameter(root,"RANSStripPointN",4)
            parameter(root,"RANSPicardRelaxation",.1)
    for b in range(4 if dimension==2 else 6):
        region=element(boundaries,"REGION",REF=b)
        for i,name in enumerate(names):
            if kind!="mms" and b in (1,3,4,5):
                partner={1:3,3:1,4:5,5:4}[b];element(region,"P",VAR=name,VALUE=f"[{partner}]");continue
            if name=="p":element(region,"N",VAR=name,USERDEFINEDTYPE="RANSPressure",VALUE="0");continue
            if kind=="mms":element(region,"D",VAR=name,USERDEFINEDTYPE="TimeDependent",VALUE=exact[i]);continue
            if b==0 and kind=="awf":
                if i<dimension-1 or name=="k":element(region,"R",VAR=name,USERDEFINEDTYPE="AWFOriginal",VALUE="0",PRIMCOEFF="0")
                else:element(region,"D",VAR=name,VALUE="0" if i==dimension-1 else exact[i])
            elif i==dimension-1:element(region,"D",VAR=name,VALUE="0")
            elif name=="eps_tilde":
                if b==0:element(region,"D",VAR=name,VALUE=exact[i])
                else:element(region,"N",VAR=name,VALUE="-Ustar^4/(1.3*30)")
            elif name=="k":element(region,"N",VAR=name,VALUE="0")
            else:
                direction = 1. if dimension==2 else (math.cos(span_angle) if i==0 else math.sin(span_angle))
                element(region,"N",VAR=name,VALUE=f"({direction:.17g})*("+("Ustar^2" if b==2 else "-Ustar^2")+")")
    for f in c.findall("FUNCTION"):c.remove(f)
    for label in ("InitialConditions","ExactSolution"):
        function=element(c,"FUNCTION",NAME=label)
        for name,value in zip(names,exact):element(function,"E",VAR=name,VALUE=re.sub(r"\bt\b","0",value) if label=="InitialConditions" else value)
        if label=="ExactSolution":
            for name,value in zip(MODEL_FIELDS,extras):element(function,"E",VAR=name,VALUE=value)
    if forces:
        function=element(c,"FUNCTION",NAME="MMSForce")
        for name,value in zip(names,forces):element(function,"E",VAR=name,VALUE=value)
        force=element(element(root,"FORCING"),"FORCE",TYPE="Body");element(force,"BODYFORCE","MMSForce")
    root.find("FILTERS/FILTER/PARAM[@NAME='Points']").text=".25 1.25 0\n.5 1.5 0" if dimension==2 and kind=="mms" else (".25 .25 1.25" if kind=="mms" else ("75 2.1 0\n150 4 0\n225 10 0" if dimension==2 else "75 25 2.1\n150 50 4\n225 75 10"))
    ET.indent(root,space="  ");tree.write(path,encoding="utf-8",xml_declaration=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("path",type=Path)
    parser.add_argument("--kind",choices=("mms","bulk_abl","awf"),default="mms")
    parser.add_argument("--dimension",type=int,choices=(2,3),default=2)
    parser.add_argument("--dt",type=float,default=.025);parser.add_argument("--steps",type=int,default=4)
    parser.add_argument("--modes",type=int,default=5);parser.add_argument("--points",type=int,default=12)
    args=vars(parser.parse_args());path=args.pop("path");generate_cls(path,**args)
