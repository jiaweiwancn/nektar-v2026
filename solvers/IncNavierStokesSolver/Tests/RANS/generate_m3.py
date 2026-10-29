#!/usr/bin/env python3
"""Generate small, auditable standard-RANS transport fixtures on Linux.

SymPy derives independent manufactured forcing. Generated meshes and fields
belong in a dedicated run directory, never in the source checkout.
"""
import argparse
import copy
from functools import lru_cache
import math
import re
from pathlib import Path
import xml.etree.ElementTree as ET


def element(parent, tag, text=None, **attributes):
    child = ET.SubElement(parent, tag, {k: str(v) for k, v in attributes.items()})
    child.text = str(text) if text is not None else None
    return child


def mesh(root, xs, ys):
    g = element(root, "GEOMETRY", DIM=2, SPACE=2)
    vertices, edges, elements = (element(g, t) for t in ("VERTEX", "EDGE", "ELEMENT"))
    nx, ny = len(xs)-1, len(ys)-1
    vertex = lambda i, j: j*(nx+1)+i
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            element(vertices, "V", f"{x:.17g} {y:.17g} 0", ID=vertex(i, j))
    edge_ids = {}
    def edge(a, b):
        key = tuple(sorted((a, b)))
        if key not in edge_ids:
            edge_ids[key] = len(edge_ids)
            element(edges, "E", f"{a} {b}", ID=edge_ids[key])
        return edge_ids[key]
    for j in range(ny):
        for i in range(nx):
            a,b,c,d = (vertex(i,j),vertex(i+1,j),vertex(i+1,j+1),vertex(i,j+1))
            ids = [edge(a,b),edge(b,c),edge(c,d),edge(d,a)]
            element(elements, "Q", " ".join(map(str,ids)), ID=j*nx+i)
    sides = [[edge(vertex(i,0),vertex(i+1,0)) for i in range(nx)],
             [edge(vertex(nx,j),vertex(nx,j+1)) for j in range(ny)],
             [edge(vertex(i,ny),vertex(i+1,ny)) for i in range(nx)],
             [edge(vertex(0,j),vertex(0,j+1)) for j in range(ny)]]
    composites = element(g,"COMPOSITE")
    element(composites,"C",f"Q[0-{nx*ny-1}]",ID=0)
    for i, side in enumerate(sides,1):
        element(composites,"C","E["+",".join(map(str,side))+"]",ID=i)
    element(g,"DOMAIN","C[0]")


@lru_cache(maxsize=None)
def manufactured(missing_transpose=False):
    import sympy as s
    x,y,t,nu = s.symbols("x y t ResolvedNu",real=True)
    Y=y-1
    psi=(1+s.sin(t)/5)*x**2*(1-x)**2*Y**2*(1-Y)**2
    u,v=s.diff(psi,y),-s.diff(psi,x)
    k=1+x/5+Y/10+s.sin(t)/20
    eps=s.Rational(7,10)+x/10+Y/10+s.cos(t)*s.Rational(3,100)
    p=(x-s.Rational(1,2))*(s.Rational(1,5)+s.cos(t)/10)+(Y-s.Rational(1,2))*s.sin(t)/10
    velocity=[u,v]; coords=[x,y]
    nut=s.Rational(9,100)*k*k/eps
    strain=[[s.diff(velocity[i],coords[j])+s.diff(velocity[j],coords[i]) for j in range(2)] for i in range(2)]
    production=nut*sum(strain[i][j]*s.diff(velocity[i],coords[j]) for i in range(2) for j in range(2))
    adv=lambda f: u*s.diff(f,x)+v*s.diff(f,y)
    diff=lambda f,d: sum(s.diff(d*s.diff(f,c),c) for c in coords)
    forces=[s.diff(velocity[i],t)+adv(velocity[i])+s.diff(p,coords[i])-
            sum(s.diff((nu+nut)*(s.diff(velocity[i],coords[j]) if missing_transpose else strain[i][j]),coords[j]) for j in range(2)) for i in range(2)]
    forces += [s.diff(k,t)+adv(k)-diff(k,nu+nut)-production+eps,
               s.diff(eps,t)+adv(eps)-diff(eps,nu+nut/s.Rational(13,10))-
               s.Rational(144,100)*eps/k*production+s.Rational(192,100)*eps*eps/k]
    def expression(f):
        return s.sstr(f).replace("**","^")
    return [expression(f) for f in [u,v,k,eps,p]], [expression(f) for f in forces]


@lru_cache(maxsize=None)
def manufactured_3d():
    import sympy as s
    x,y,z,t,nu=s.symbols("x y z t ResolvedNu",real=True)
    coords=[x,y,z]
    psi=s.prod(c**2*(1-c)**2 for c in coords)
    a=10*(1+s.sin(t)/5);b=5*(1+s.cos(t)/10)
    velocity=[a*s.diff(psi,y),-a*s.diff(psi,x)+b*s.diff(psi,z),-b*s.diff(psi,y)]
    k=1+x/10+y/20+z*s.Rational(3,100)+s.sin(t)/50
    eps=s.Rational(7,10)+x/25+y*s.Rational(3,100)+z/50+s.cos(t)/100
    p=(x-s.Rational(1,2))*(s.Rational(1,5)+s.cos(t)/10)+(y-s.Rational(1,2))*s.sin(t)/20+(z-s.Rational(1,2))*s.cos(2*t)/50
    nut=s.Rational(9,100)*k*k/eps
    strain=[[s.diff(velocity[i],coords[j])+s.diff(velocity[j],coords[i]) for j in range(3)] for i in range(3)]
    production=nut*sum(strain[i][j]*s.diff(velocity[i],coords[j]) for i in range(3) for j in range(3))
    adv=lambda f:sum(velocity[i]*s.diff(f,c) for i,c in enumerate(coords))
    diff=lambda f,d:sum(s.diff(d*s.diff(f,c),c) for c in coords)
    forces=[s.diff(velocity[i],t)+adv(velocity[i])+s.diff(p,coords[i])-
            sum(s.diff((nu+nut)*strain[i][j],coords[j]) for j in range(3)) for i in range(3)]
    forces += [s.diff(k,t)+adv(k)-diff(k,nu+nut)-production+eps,
               s.diff(eps,t)+adv(eps)-diff(eps,nu+nut/s.Rational(13,10))-
               s.Rational(144,100)*eps/k*production+s.Rational(192,100)*eps*eps/k]
    return [s.sstr(f).replace("**","^") for f in velocity+[k,eps,p]], [s.sstr(f).replace("**","^") for f in forces]


def generate(path, kind="mms", dt=.05, steps=4, order=2, modes=7, points=None,
             diffusion="Molecular", tolerance=1e-10, top_flux=True,
             missing_transpose=False):
    root=ET.Element("NEKTAR")
    exp=element(root,"EXPANSIONS")
    n=points or modes+3
    element(exp,"E",COMPOSITE="C[0]",NUMMODES=f"{modes},{modes}",
            BASISTYPE="Modified_A,Modified_A",POINTSTYPE="GaussLobattoLegendre,GaussLobattoLegendre",
            NUMPOINTS=f"{n},{n}",FIELDS="u,v,k,eps")
    element(exp,"E",COMPOSITE="C[0]",NUMMODES=f"{modes-2},{modes-2}",
            BASISTYPE="Modified_A,Modified_A",POINTSTYPE="GaussLobattoLegendre,GaussLobattoLegendre",
            NUMPOINTS=f"{n},{n}",FIELDS="p")
    conditions=element(root,"CONDITIONS")
    info=element(conditions,"SOLVERINFO")
    for name,value in {"SolverType":"RANSVelocityCorrectionScheme","EQTYPE":"UnsteadyNavierStokes",
                       "Projection":"Galerkin","TurbulenceModel":"StandardKEpsilon",
                       "RANSResolvedDiffusion":diffusion,"GlobalSysSoln":"DirectStaticCond"}.items():
        element(info,"I",PROPERTY=name,VALUE=value)
    time=element(conditions,"TIMEINTEGRATIONSCHEME")
    element(time,"METHOD","IMEX");element(time,"ORDER",order)
    parameters=element(conditions,"PARAMETERS")
    kappa=math.sqrt(.48*1.3*.3); star=kappa*10/math.log(10/.02)
    physical_nu=.025 if kind=="mms" else 1.5e-5
    for name,value in {"TimeStep":dt,"NumSteps":steps,"IO_InfoSteps":max(steps,1),
                       "Kinvis":physical_nu,"ResolvedNu":physical_nu if diffusion=="Molecular" else 0,
                       "IterativeSolverTolerance":1e-13,
                       "RANSNonlinearTolerance":tolerance,
                       "RANSMaxIterations":500,"Ustar":star,"Kappa":kappa,"Z0":.02}.items():
        element(parameters,"P",f"{name} = {value:.17g}")
    variables=element(conditions,"VARIABLES")
    names=["u","v","k","eps","p"]
    for i,name in enumerate(names):element(variables,"V",name,ID=i)
    regions=element(conditions,"BOUNDARYREGIONS")
    for i in range(4):element(regions,"B",f"C[{i+1}]",ID=i)
    bcs=element(conditions,"BOUNDARYCONDITIONS")
    forces=None
    if kind=="mms":
        exact,forces=manufactured(missing_transpose);xs=[0,.5,1];ys=[1,1.5,2]
    elif kind in ("decay","stiff"):
        initial_rate=3 if kind=="decay" else 50
        rate=f"(1+0.92*{initial_rate}*t)"
        exact=["0","0",f"{rate}^(-1/0.92)",f"{initial_rate}*{rate}^(-1.92/0.92)","0"]
        xs=[0,.5,1];ys=[1,1.5,2]
    elif kind in ("abl","abl_shift"):
        shifted=kind=="abl_shift"
        if shifted:
            star=kappa*10/math.log((10+.02)/.02)
            for p in parameters:
                if p.text.startswith("Ustar ="):p.text=f"Ustar = {star:.17g}"
            element(parameters,"P","RANSHeightShift = Z0")
        height="(y+Z0)" if shifted else "y"
        exact=[f"Ustar/Kappa*log({height}/Z0)","0","Ustar^2/0.3",f"Ustar^3/(Kappa*{height})","0"]
        xs=[0,150,300];ys=[2,4,8,16,30]
    else:raise ValueError(kind)
    for i in range(4):
        region=element(bcs,"REGION",REF=i)
        for j,name in enumerate(names):
            if kind.startswith("abl") and i in (1,3):
                element(region,"P",VAR=name,VALUE=f"[{4-i}]");continue
            if name=="p":
                element(region,"N",VAR=name,USERDEFINEDTYPE="RANSPressure",VALUE="0");continue
            if kind in ("decay","stiff"):
                element(region,"D" if j<2 else "N",VAR=name,VALUE="0");continue
            if kind.startswith("abl"):
                if i==0:
                    if name=="u":element(region,"N",VAR=name,USERDEFINEDTYPE="RANSLogLaw",VALUE="0")
                    elif name=="k":element(region,"N",VAR=name,VALUE="0")
                    else:element(region,"D",VAR=name,VALUE=exact[j])
                else:
                    val={"u":"Ustar^2","v":"0","k":"0","eps":f"-Ustar^4/(1.3*{height})" if top_flux else "0"}[name]
                    element(region,"D" if name=="v" else "N",VAR=name,VALUE=val)
            else:
                element(region,"D",VAR=name,USERDEFINEDTYPE="TimeDependent",VALUE=exact[j])
    for label in ("InitialConditions","ExactSolution"):
        function=element(conditions,"FUNCTION",NAME=label)
        for name,value in zip(names,exact):
            element(function,"E",VAR=name,VALUE=re.sub(r"\bt\b","0",value) if label=="InitialConditions" else value)
    if forces:
        function=element(conditions,"FUNCTION",NAME="MMSForce")
        for name,value in zip(names,forces):element(function,"E",VAR=name,VALUE=value)
        force=element(element(root,"FORCING"),"FORCE",TYPE="Body")
        element(force,"BODYFORCE","MMSForce")
    mesh(root,xs,ys)
    history=element(element(root,"FILTERS"),"FILTER",TYPE="HistoryPoints")
    element(history,"PARAM",1,NAME="OutputFrequency")
    locations=[(75,2.1),(150,4),(225,10),(150,20),(75,29)] if kind.startswith("abl") else [(0.25,1.25),(.5,1.5),(.75,1.75)]
    element(history,"PARAM","\n"+"\n".join(f"{x} {y} 0" for x,y in locations)+"\n",NAME="Points")
    ET.indent(root,space="  ")
    path.parent.mkdir(parents=True,exist_ok=True)
    ET.ElementTree(root).write(path,encoding="utf-8",xml_declaration=True)


def generate_3d(path, kind="decay", **options):
    """Eight-hex periodic decay or a full three-component manufactured case."""
    generate(path, kind="decay", **options)
    tree=ET.parse(path);root=tree.getroot()
    original_text=(Path(__file__).resolve().parents[1]/"TGV3d_2_2.xml").read_text()
    # This historical fixture places a comment before its XML declaration.
    original=ET.fromstring(re.sub(r"<\?xml[^>]*\?>","",original_text))
    root.remove(root.find("GEOMETRY"));root.append(copy.deepcopy(original.find("GEOMETRY")))
    exp=root.find("EXPANSIONS");exp.clear()
    modes=options.get("modes",5);points=options.get("points") or modes+3
    names=["u","v","w","k","eps","p"]
    for fields,n in [("u,v,w,k,eps",modes),("p",modes-2)]:
        element(exp,"E",COMPOSITE=original.find("EXPANSIONS/E").attrib["COMPOSITE"],NUMMODES=",".join([str(n)]*3),
                BASISTYPE=",".join(["Modified_A"]*3),POINTSTYPE=",".join(["GaussLobattoLegendre"]*3),
                NUMPOINTS=",".join([str(points)]*3),FIELDS=fields)
    c=root.find("CONDITIONS");variables=c.find("VARIABLES");variables.clear()
    for i,name in enumerate(names):element(variables,"V",name,ID=i)
    regions=c.find("BOUNDARYREGIONS");regions.clear()
    for b in original.findall("CONDITIONS/BOUNDARYREGIONS/B"):regions.append(copy.deepcopy(b))
    boundaries=c.find("BOUNDARYCONDITIONS");boundaries.clear()
    for r in original.findall("CONDITIONS/BOUNDARYCONDITIONS/REGION"):
        region=element(boundaries,"REGION",REF=r.attrib["REF"])
        if kind=="mms":
            for name in names:
                element(region,"N" if name=="p" else "D",VAR=name,
                        USERDEFINEDTYPE="RANSPressure" if name=="p" else "TimeDependent",
                        VALUE="0")
        else:
            value=r.find("P").attrib["VALUE"]
            for name in names:element(region,"P",VAR=name,VALUE=value)
    if kind=="mms":
        vertex=root.find("GEOMETRY/VERTEX")
        for axis in "XYZ":vertex.attrib[axis+"SCALE"]="1"
        for p in c.findall("PARAMETERS/P"):
            if p.text.startswith("Kinvis ="):p.text="Kinvis = 0.025"
            if p.text.startswith("ResolvedNu ="):p.text="ResolvedNu = "+("0.025" if options.get("diffusion","Molecular")=="Molecular" else "0")
        exact,forces=manufactured_3d()
        for f in c.findall("FUNCTION"):c.remove(f)
        for label in ("InitialConditions","ExactSolution"):
            f=element(c,"FUNCTION",NAME=label)
            for name,value in zip(names,exact):element(f,"E",VAR=name,VALUE=re.sub(r"\bt\b","0",value) if label=="InitialConditions" else value)
        for region in boundaries:
            for boundary in region:
                if boundary.attrib["VAR"]!="p":boundary.attrib["VALUE"]=exact[names.index(boundary.attrib["VAR"])]
        f=element(c,"FUNCTION",NAME="MMSForce")
        for name,value in zip(names,forces):element(f,"E",VAR=name,VALUE=value)
        force=element(element(root,"FORCING"),"FORCE",TYPE="Body");element(force,"BODYFORCE","MMSForce")
    else:
        for f in c.findall("FUNCTION"):element(f,"E",VAR="w",VALUE="0")
    history=root.find("FILTERS/FILTER/PARAM[@NAME='Points']");history.text="0.1 0.2 0.3\n"
    ET.indent(root,space="  ");tree.write(path,encoding="utf-8",xml_declaration=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path",type=Path)
    parser.add_argument("--kind",choices=["mms","decay","stiff","abl","abl_shift"],default="mms")
    parser.add_argument("--dt",type=float,default=.05)
    parser.add_argument("--steps",type=int,default=4)
    parser.add_argument("--order",type=int,default=2)
    parser.add_argument("--modes",type=int,default=7)
    parser.add_argument("--points",type=int)
    parser.add_argument("--diffusion",choices=["HighRe","Molecular"],default="Molecular")
    args=vars(parser.parse_args())
    path=args.pop("path");generate(path,**args)
