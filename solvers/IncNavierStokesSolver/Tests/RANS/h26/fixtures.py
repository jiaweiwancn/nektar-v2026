#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H26 independent polynomial basis, physical MMS and native XML generators."""
import argparse,hashlib,json,math
from pathlib import Path
import xml.etree.ElementTree as ET


def modified(degree,x):
    # Independent Jacobi-polynomial definition, not native sampled basis data.
    import numpy as np
    from scipy.special import eval_jacobi
    x=np.asarray(x);v=[(1-x)/2,(1+x)/2];d=[np.full_like(x,-.5),np.full_like(x,.5)]
    for p in range(2,degree+1):
        j=eval_jacobi(p-2,1,1,x);dj=0*x if p==2 else .5*(p+1)*eval_jacobi(p-3,2,2,x)
        v.append((1-x*x)*j/4);d.append(-x*j/2+(1-x*x)*dj/4)
    return np.array(v).T,np.array(d).T


def symbolic():
    import sympy as s
    x,z,t,nu,a=s.symbols('x z t nu a',real=True)
    psi=s.exp(-t)*x*x*(1-x)*z*z*(1-z)
    u=[s.diff(psi,z),-s.diff(psi,x)];p=s.exp(-2*t)*(x*x+z*z-s.Rational(2,3))
    mu=nu+a*(1+x+z)/100;rem=[[x*z*s.exp(-t)/1000,(x+z)*s.exp(-t)/2000],[(x+z)*s.exp(-t)/2000,-x*z*s.exp(-t)/1000]]
    g=s.Matrix(u).jacobian([x,z]);stress=mu*(g+g.T)+s.Matrix(rem)
    f=[s.simplify(s.diff(u[i],t)+sum(u[j]*g[i,j]for j in range(2))+s.diff(p,[x,z][i])-sum(s.diff(stress[i,j],[x,z][j])for j in range(2)))for i in range(2)]
    assert s.simplify(g.trace())==0
    return dict(coordinates='physical Cartesian x,z; manufactured verification only',
        velocity=[str(v)for v in u],pressure=str(p),mu=str(mu),remainder=[[str(v)for v in row]for row in rem],
        continuous_forcing=[str(v)for v in f],divergence=str(s.simplify(g.trace())),
        c99_forcing=[s.ccode(v)for v in f],derivation='symbolic physical derivatives of continuous PDE, independent of native matrices')


def mesh(path,xs,rs,degree=3,quadrature=8,curve=0.,natural=False):
    n,m=len(xs)-1,len(rs)-1;x=ET.Element('NEKTAR');g=ET.SubElement(x,'GEOMETRY',DIM='2',SPACE='2')
    vs=ET.SubElement(g,'VERTEX');edges=ET.SubElement(g,'EDGE');elements=ET.SubElement(g,'ELEMENT');curved=ET.SubElement(g,'CURVED')
    vertex=lambda i,j:j*(n+1)+i
    ground=lambda a:curve*a*a*(1-a)*(1-a)
    for j,r in enumerate(rs):
        for i,a in enumerate(xs):ET.SubElement(vs,'V',ID=str(vertex(i,j))).text=f'{a:.17g} {r+ground(a):.17g} 0'
    hs={};es={};eid=0
    for j in range(m+1):
        for i in range(n):
            hs[i,j]=eid;ET.SubElement(edges,'E',ID=str(eid)).text=f'{vertex(i,j)} {vertex(i+1,j)}'
            if curve:
                pts=[]
                for k in range(5):
                    a=xs[i]+(xs[i+1]-xs[i])*k/4;pts.extend((a,rs[j]+ground(a),0.))
                ET.SubElement(curved,'E',ID=str(eid),EDGEID=str(eid),TYPE='PolyEvenlySpaced',NUMPOINTS='5').text=' '.join(f'{v:.17g}'for v in pts)
            eid+=1
    for j in range(m):
        for i in range(n+1):es[i,j]=eid;ET.SubElement(edges,'E',ID=str(eid)).text=f'{vertex(i,j)} {vertex(i,j+1)}';eid+=1
    for j in range(m):
        for i in range(n):ET.SubElement(elements,'Q',ID=str(j*n+i)).text=f'{hs[i,j]} {es[i+1,j]} {hs[i,j+1]} {es[i,j]}'
    cs=ET.SubElement(g,'COMPOSITE');ET.SubElement(cs,'C',ID='0').text=f'Q[0-{n*m-1}]'
    for bid,ids in enumerate(([es[0,j]for j in range(m)],[es[n,j]for j in range(m)],[hs[i,0]for i in range(n)],[hs[i,m]for i in range(n)]),1):ET.SubElement(cs,'C',ID=str(bid)).text='E['+','.join(map(str,ids))+']'
    ET.SubElement(g,'DOMAIN').text='C[0]'
    ex=ET.SubElement(x,'EXPANSIONS')
    for fields,modes in (('u,v',degree+1),('p',degree)):
        ET.SubElement(ex,'E',COMPOSITE='C[0]',NUMMODES=f'{modes},{modes}',BASISTYPE='Modified_A,Modified_A',FIELDS=fields,NUMPOINTS=f'{quadrature},{quadrature}',POINTSTYPE='GaussLobattoLegendre,GaussLobattoLegendre')
    c=ET.SubElement(x,'CONDITIONS');si=ET.SubElement(c,'SOLVERINFO')
    for k,v in dict(SolverType='VelocityCorrectionScheme',EQTYPE='UnsteadyNavierStokes',Projection='Continuous',GlobalSysSoln='DirectStaticCond',AdvectionType='Convective',SpectralVanishingViscosity='False').items():ET.SubElement(si,'I',PROPERTY=k,VALUE=v)
    ts=ET.SubElement(c,'TIMEINTEGRATIONSCHEME');ET.SubElement(ts,'METHOD').text='IMEX';ET.SubElement(ts,'ORDER').text='2'
    ps=ET.SubElement(c,'PARAMETERS')
    for k,v in dict(TimeStep=.005,NumSteps=420,Kinvis=6.423173803526447e-5,IO_CheckSteps=200,IO_InfoSteps=20).items():ET.SubElement(ps,'P').text=f'{k} = {v}'
    br=ET.SubElement(c,'BOUNDARYREGIONS');bc=ET.SubElement(c,'BOUNDARYCONDITIONS')
    for b in range(4):
        ET.SubElement(br,'B',ID=str(b)).text=f'C[{b+1}]';reg=ET.SubElement(bc,'REGION',REF=str(b))
        for name in ('u','v'):ET.SubElement(reg,'N'if natural and b==3 else'D',VAR=name,VALUE='0')
        ET.SubElement(reg,'N',VAR='p',VALUE='0',USERDEFINEDTYPE='H')
    variables=ET.SubElement(c,'VARIABLES')
    for i,name in enumerate(('u','v','p')):ET.SubElement(variables,'V',ID=str(i)).text=name
    ic=ET.SubElement(c,'FUNCTION',NAME='InitialConditions')
    for name in ('u','v','p'):ET.SubElement(ic,'E',VAR=name,VALUE='0')
    ET.indent(x);Path(path).write_text(ET.tostring(x,encoding='unicode')+'\n')


def generate(out):
    out.mkdir(exist_ok=False,parents=True)
    for curve in (0.,.1):
        for natural in (False,True):mesh(out/f'tiny-c{curve}-n{int(natural)}.xml',[0.,.5,1.],[0.,.15,.5,1.],quadrature=12,curve=curve,natural=natural)
    (out/'continuous-mms.json').write_text(json.dumps(symbolic(),indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps(dict(files=[dict(name=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in sorted(out.iterdir())]),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);generate(p.parse_args().out)
