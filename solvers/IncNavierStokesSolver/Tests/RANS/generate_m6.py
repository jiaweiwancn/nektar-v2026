#!/usr/bin/env python3
"""Neutral-ABL preservation meshes, profiles and boundary-control sessions."""
import argparse
import base64
import math
from pathlib import Path
import re
import struct
import xml.etree.ElementTree as ET
import zlib
from generate_m3 import element, generate, mesh
from generate_m4 import info, parameter


def original_j0(k, height, roughness, nu):
    """Closed-form integral of 1/(nu+nu_t) for the original AWF anchor.

    This generator needs only the standard library. The independent SciPy
    quadrature in generate_m4 remains the validation oracle.
    """
    root=math.sqrt(k);slope=.09*2.55*root
    rh=roughness*root/(70*nu)
    intercept=10.7 if roughness==0 else 10.7*(1-rh**max(.5-.4*rh**.7,1-.79*rh**(-.28)))
    zv=nu*intercept/root
    if height<=zv:return height/nu
    if zv>0:return zv/nu+math.log1p(slope*(height-zv)/nu)/slope
    return math.log1p(height/(nu/slope-zv))/slope


def structured_hex(root, xs, ys, zs):
    """Tile the tracked native hex's ordered face template, sharing topology."""
    template=ET.fromstring(re.sub(r"<\?xml[^>]*\?>", "", (Path(__file__).parents[1]/"TGV3d_2_2.xml").read_text())).find("GEOMETRY")
    unpack=lambda node,fmt:list(struct.iter_unpack("<"+fmt,zlib.decompress(base64.b64decode(node.text))))
    vertices={int(v[0]):tuple(round(c*2)/2 for c in v[1:]) for v in unpack(template.find("VERTEX"),"qddd")}
    edges={int(v[0]):v[1:] for v in unpack(template.find("EDGE"),"qqq")}
    faces={int(v[0]):v[1:] for v in unpack(template.find("FACE/Q"),"qqqqq")}
    first=unpack(template.find("ELEMENT/H"),"qqqqqqq")[0][1:]
    ids={v for f in first for e in faces[f] for v in edges[e]}
    low=[min(vertices[v][d] for v in ids) for d in range(3)]
    high=[max(vertices[v][d] for v in ids) for d in range(3)]
    local={v:tuple(round((vertices[v][d]-low[d])/(high[d]-low[d])) for d in range(3)) for v in ids}
    nx,ny,nz=len(xs)-1,len(ys)-1,len(zs)-1
    vertex=lambda i,j,k:(k*(ny+1)+j)*(nx+1)+i
    geometry=element(root,"GEOMETRY",DIM=3,SPACE=3);vg=element(geometry,"VERTEX")
    coordinates={}
    for k,z in enumerate(zs):
        for j,y in enumerate(ys):
            for i,x in enumerate(xs):
                n=vertex(i,j,k);coordinates[n]=(x,y,z);element(vg,"V",f"{x:.17g} {y:.17g} {z:.17g}",ID=n)
    edge_data=[];face_data=[];cell_data=[];edge_map={};face_map={}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                vmap={v:vertex(i+a,j+b,k+c) for v,(a,b,c) in local.items()}
                mapped=[]
                for f in first:
                    ev=[];fv=set()
                    for e in faces[f]:
                        pair=tuple(vmap[v] for v in edges[e]);key=tuple(sorted(pair));fv.update(pair)
                        if key not in edge_map:edge_map[key]=len(edge_data);edge_data.append(pair)
                        ev.append(edge_map[key])
                    key=tuple(sorted(fv))
                    if key not in face_map:face_map[key]=len(face_data);face_data.append(tuple(ev))
                    mapped.append(face_map[key])
                cell_data.append(mapped)
    for label,tag,data in (("EDGE","E",edge_data),("FACE","Q",face_data),("ELEMENT","H",cell_data)):
        node=element(geometry,label)
        for i,values in enumerate(data):element(node,tag," ".join(map(str,values)),ID=i)
    sides=[[] for _ in range(6)]
    for key,f in face_map.items():
        for b,(axis,value) in enumerate(((2,zs[0]),(0,xs[-1]),(2,zs[-1]),(0,xs[0]),(1,ys[0]),(1,ys[-1]))):
            if all(coordinates[v][axis]==value for v in key):sides[b].append((tuple(sum(coordinates[v][d] for v in key)/4 for d in range(3) if d!=axis),f))
    composites=element(geometry,"COMPOSITE");element(composites,"C",f"H[0-{len(cell_data)-1}]",ID=0)
    for b,entries in enumerate(sides,1):element(composites,"C","F["+",".join(str(f) for _,f in sorted(entries))+"]",ID=b)
    element(geometry,"DOMAIN","C[0]")


def generate_abl(path, model="cls_switched", dimension=2, degree=4, layers=12,
                 dt=1., steps=4, flow_times=None, length=300., width=20., nx=2, ny=1,
                 height=2., top=30., roughness=.1, boundary="periodic", diffusion="HighRe",
                 backend="static", perturbation=0., tolerance=1e-9, checkpoint=0,
                 monitor=False, wrong_top=False, wrong_kappa=False):
    if model not in ("standard","cls_original","cls_fixed","cls_switched"):
        raise ValueError("Unknown neutral-ABL model")
    if dimension not in (2,3) or degree<3 or layers<1 or nx<1 or ny<1:
        raise ValueError("Use full 2D/3D, degree >=3 and positive mesh counts")
    if not (dt>0 and length>0 and width>0 and height>0 and 2*height<top):
        raise ValueError("Require positive dimensions/timestep and 0 < H < N=2H < Z")
    if backend not in ("static","full","matrixfree") or boundary not in ("periodic","matched","outflow"):
        raise ValueError("Unknown linear backend or streamwise boundary")
    standard=model=="standard";q=10/3 if standard else 3.190675081624288
    kappa=math.sqrt(.48*1.3/q) if standard else .440553905186659
    star=kappa*10/math.log(10/.02) if standard else .4
    k=q*star**2;nu=1.5e-5
    matching=star/kappa*math.log(height/.02) if standard else star**2*original_j0(k,height,roughness,nu)
    uref=matching+star/kappa*math.log(10/height)
    if flow_times is not None:steps=math.ceil(flow_times*length/uref/dt)
    generate(path,kind="decay",dt=dt,steps=steps,order=2,modes=degree+1,points=degree+4,diffusion=diffusion,tolerance=tolerance)
    tree=ET.parse(path);root=tree.getroot();c=root.find("CONDITIONS")
    names=["u","v"]+(["w"] if dimension==3 else [])+["k","eps" if standard else "eps_tilde","p"]
    wall="None" if standard else ("AWFOriginal" if model=="cls_original" else "AWFEquilibrium")
    info(root,"TurbulenceModel","StandardKEpsilon" if standard else "CLSHighRe");info(root,"RANSWallModel",wall)
    info(root,"GlobalSysSoln","IterativeStaticCond")
    parameters={"Kinvis":nu,"ResolvedNu":0 if diffusion=="HighRe" else nu,"Ustar":star,"Qstar":q,
        "Kappa":kappa,"MatchingVelocity":matching,"H":height,"Z":top,"L":length,"Span":width,"Uref":uref,
        "RANSNonlinearTolerance":tolerance,"RANSPressureTolerance":max(tolerance,1e-8) if dimension==3 else tolerance,
        "RANSMaxIterations":2000,"RANSStageAudit":0,"RANSReferenceMonitor":int(monitor),
        "IO_InfoSteps":1,"IO_CheckSteps":checkpoint,"RANSStripHeight":height,"RANSStripPointP":height/2,
        "RANSStripPointN":2*height,"RANSRoughnessHeight":roughness,"RANSAWFFixedWeight":1,
        "RANSAWFOn":.02,"RANSAWFOff":.1,"RANSAWFMinimumRe":100,"RANSAWFPointDiagnostics":0}
    for name,value in parameters.items():parameter(root,name,value)
    info(root,"RANSAWFWeight","Fixed" if model=="cls_fixed" else "Switched");info(root,"RANSAWFFailurePolicy","RequirePreservation")
    variables=c.find("VARIABLES");variables.clear()
    for i,name in enumerate(names):element(variables,"V",name,ID=i)
    root.remove(root.find("GEOMETRY"));xs=[length*i/nx for i in range(nx+1)]
    zs=[height*(top/height)**(i/layers) for i in range(layers+1)];zs[0]=height;zs[-1]=top
    if dimension==2:mesh(root,xs,zs)
    else:structured_hex(root,xs,[width*j/ny for j in range(ny+1)],zs)
    for node,fields,modes in zip(root.find("EXPANSIONS"),(",".join(names[:-1]),"p"),(degree+1,degree-1)):
        node.attrib.update(NUMMODES=",".join([str(modes)]*dimension),NUMPOINTS=",".join([str(degree+4)]*dimension),
            BASISTYPE=",".join(["Modified_A"]*dimension),POINTSTYPE=",".join(["GaussLobattoLegendre"]*dimension),FIELDS=fields)
    if backend=="matrixfree":info(root,"RANSLinearSolver","MatrixFree")
    if backend=="full":
        settings=element(c,"GLOBALSYSSOLNINFO")
        group=element(settings,"V",VAR=",".join(names[:-1]))
        element(group,"I",PROPERTY="GlobalSysSoln",VALUE="IterativeFull")
        element(group,"I",PROPERTY="Preconditioner",VALUE="Null")
    coord="y" if dimension==2 else "z"
    u=f"MatchingVelocity+Ustar/Kappa*log({coord}/H)"
    exact=[u,"0"]+(["0"] if dimension==3 else [])+["Qstar*Ustar^2",f"Ustar^3/(Kappa*{coord})","0"]
    if wrong_kappa:parameter(root,"Kappa",.41)
    regions=c.find("BOUNDARYREGIONS");regions.clear()
    boundaries=c.find("BOUNDARYCONDITIONS");boundaries.clear()
    for b in range(4 if dimension==2 else 6):
        element(regions,"B",f"C[{b+1}]",ID=b);region=element(boundaries,"REGION",REF=b)
        for i,name in enumerate(names):
            if b in (4,5) or (boundary=="periodic" and b in (1,3)):
                element(region,"P",VAR=name,VALUE=f"["+str({1:3,3:1,4:5,5:4}[b])+"]");continue
            if name=="p":element(region,"N",VAR=name,USERDEFINEDTYPE="RANSPressure",VALUE="0");continue
            if b==3 or (boundary=="matched" and b==1):
                if b==1 and i>=dimension:element(region,"N",VAR=name,VALUE="0")
                else:element(region,"D",VAR=name,VALUE=exact[i])
            elif b==1 and boundary=="outflow":
                if i<dimension:
                    traction="(-0.253333333333333333*Qstar*Ustar^2)" if i==0 and not standard else ("Ustar^2" if i==dimension-1 else "0")
                    element(region,"N",VAR=name,USERDEFINEDTYPE="RANSOutflow",VALUE=traction)
                else:element(region,"N",VAR=name,VALUE="0")
            elif b==0:
                if standard:
                    if i<dimension-1:element(region,"N",VAR=name,USERDEFINEDTYPE="RANSLogLaw",VALUE="0")
                    elif i==dimension-1:element(region,"D",VAR=name,VALUE="0")
                    elif name=="k":element(region,"N",VAR=name,VALUE="0")
                    else:element(region,"D",VAR=name,VALUE=exact[i])
                elif i<dimension-1 or name=="k":element(region,"R",VAR=name,USERDEFINEDTYPE=wall,VALUE="0",PRIMCOEFF="0")
                else:element(region,"D",VAR=name,VALUE="0" if i==dimension-1 else exact[i])
            elif i==dimension-1:element(region,"D",VAR=name,VALUE="0")
            elif name=="k":element(region,"N",VAR=name,VALUE="0")
            elif i<dimension-1:element(region,"N",VAR=name,VALUE="Ustar^2" if i==0 else "0")
            else:element(region,"N",VAR=name,VALUE="0" if wrong_top else "-Ustar^4/(1.3*Z)")
    for f in c.findall("FUNCTION"):c.remove(f)
    # An exactly divergence-free smooth perturbation, zero at H, N and Z.
    zz=f"({coord}-H)"; zn=f"({coord}-2*H)"; zt=f"(Z-{coord})"
    g=f"({zz}*{zn}*{zt})^2"
    gp=f"2*{zz}*{zn}*{zt}*({zn}*{zt}+{zz}*{zt}-{zz}*{zn})"
    sample=[height+(top-height)*i/10000 for i in range(10001)]
    gv=lambda z:((z-height)*(z-2*height)*(top-z))**2
    gpv=lambda z:2*(z-height)*(z-2*height)*(top-z)*((z-2*height)*(top-z)+(z-height)*(top-z)-(z-height)*(z-2*height))
    scale=max(abs(gpv(z)) for z in sample);gmax=max(gv(z) for z in sample)
    fx="16*x^2*(L-x)^2/L^4"; fpx="32*x*(L-x)*(L-2*x)/L^4"
    for label in ("InitialConditions","ExactSolution"):
        function=element(c,"FUNCTION",NAME=label)
        for name,value in zip(names,exact):
            if label=="InitialConditions" and perturbation:
                if name=="u":value=f"({value})+({perturbation}/({scale:.17g}))*Ustar*({fx})*({gp})"
                if name==names[dimension-1]:value=f"-({perturbation}/({scale:.17g}))*Ustar*({fpx})*({g})"
                if name=="k":value=f"({value})*(1+({perturbation}/({gmax:.17g}))*({fx})*({g}))"
            element(function,"E",VAR=name,VALUE=value)
    refs={"nut":f"Kappa*Ustar*{coord}","P_k":exact[-2],"eps_total":exact[-2],"Yap":"0",
          "Rxx":"(2/3)*Qstar*Ustar^2" if standard else ".92*Qstar*Ustar^2",
          "Ryy":"(2/3)*Qstar*Ustar^2" if standard else ".56*Qstar*Ustar^2",
          "Rzz":"(2/3)*Qstar*Ustar^2" if standard else ".52*Qstar*Ustar^2",
          "Rxy":"0","Rxz":"-Ustar^2","Ryz":"0","CLS_quadratic_norm":"0","CLS_cubic_norm":"0"}
    if not standard:
        r=k**1.5/(2.5*star**3/kappa)
        hy=.83*max(r-1,0)*r*r
        refs["Yap"]=f"({hy:.17g})*({exact[-2]})^2/({exact[-3]})"
        quadratic=math.sqrt(sum((a-2/3)**2 for a in (.92,.56,.52)))
        refs["CLS_quadratic_norm"]=f"({quadratic:.17g})*Qstar*Ustar^2"
    function=c.find("FUNCTION[@NAME='ExactSolution']")
    rename={"nut":"rans_nut","P_k":"rans_P","eps_total":"rans_eps_total","Yap":"rans_Yap",
            "Rxx":"rans_R_xx","Ryy":"rans_R_yy","Rzz":"rans_R_zz","Rxy":"rans_R_xy",
            "Rxz":"rans_R_xz","Ryz":"rans_R_yz","CLS_quadratic_norm":"rans_quadratic_norm","CLS_cubic_norm":"rans_cubic_norm"}
    for name,value in refs.items():element(function,"E",VAR=rename[name],VALUE=value)
    root.find("FILTERS/FILTER/PARAM[@NAME='Points']").text=f"{.25*length} {height+.1} 0\n{.5*length} {2*height} 0\n{.95*length} 10 0" if dimension==2 else f"{.25*length} {.5*width} {height+.1}\n{.5*length} {.5*width} {2*height}\n{.95*length} {.5*width} 10"
    ET.indent(root,space="  ");tree.write(path,encoding="utf-8",xml_declaration=True)
    return dict(model=model,dimension=dimension,degree=degree,layers=layers,dt=dt,steps=steps,physical_time=steps*dt,
                length=length,width=width,height=height,top=top,friction_velocity=star,uref=uref,flow_times=steps*dt*uref/length,
                q=q,kappa=kappa,matching_velocity=matching,log_z0=math.log(height)-kappa*matching/star)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("path",type=Path)
    for key,kind,default in (("model",str,"cls_switched"),("dimension",int,2),("degree",int,4),("layers",int,12),
                             ("dt",float,1.),("steps",int,4),("backend",str,"static"),("boundary",str,"periodic")):
        parser.add_argument("--"+key,type=kind,default=default)
    args=vars(parser.parse_args());path=args.pop("path");print(generate_abl(path,**args))
