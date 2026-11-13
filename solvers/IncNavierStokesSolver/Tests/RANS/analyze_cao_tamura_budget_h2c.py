#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent exact-radial mapped-seed witnesses and indexed covariance replay."""
import argparse,json,math,os
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import numpy as np
from validate_cls_epawf_m1 import indexed

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--preflight',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    pre=json.loads(a.preflight.read_text());progress=json.loads((a.root/'progress.json').read_text());results={}
    gx,gw=np.polynomial.legendre.leggauss(96)
    for name,record in progress['results'].items():
        if 'result_file' in record:record=json.loads(Path(record['result_file']).read_text())
        elif not name.endswith('-hill'):continue
        definition=json.loads((a.root/name/'definition.json').read_text());H=definition['H_A'];v=next(h for h in pre['heights'] if h['H']==H)
        adjusted=definition['initialization'].startswith(('PhysicalHeightInterfaceV1','PhysicalHeightInterfaceV2'))
        if 'seed_target_witnesses' in record:record=dict(record,target=record['seed_target_witnesses'])
        us=pre['u_tau'];kap=pre['kappa_CLS'];K0=v['K0'];nu=pre['nu'];tau=us*us;slope=pre['cMu_wall']*pre['cLength_wall']*math.sqrt(K0);delta=nu/slope-v['viscosity_intercept'];UH=v['U_H'];A=us**3/kap;B=.2;Z=.8;L=.1
        assert v['viscosity_intercept']<0 and delta>0
        def radial(r,d=0):
            if r<=H:
                if d<0:return tau/slope*((r+delta)*math.log1p(r/delta)-r)
                if d==0:return tau/slope*math.log1p(r/delta)
                if d==1:return tau/(slope*(r+delta))
                return -tau/(slope*(r+delta)**2)
            if d<0:return radial(H,-1)+UH*(r-H)+us/kap*(r*math.log(r/H)-r+H)
            if d==0:return UH+us/kap*math.log(r/H)
            return us/(kap*r) if d==1 else -us/(kap*r*r)
        def geom(s,r,foot=0,outer_side=False):
            inside=abs(s)<L or (abs(abs(s)-L)<1e-13 and ((s<0 and foot>0) or (s>0 and foot<0)))
            f=.04*math.cos(math.pi*s/(2*L))**2 if abs(s)<L else 0.
            fp=-.04*math.pi/(2*L)*math.sin(math.pi*s/L) if inside else 0.
            fpp=-.04*math.pi**2/(2*L**2)*math.cos(math.pi*s/L) if inside else 0.
            fppp=.04*math.pi**3/(2*L**3)*math.sin(math.pi*s/L) if inside else 0.
            q=math.hypot(1,fp);k=-fpp/q**3;h=f+H/q;x=s-H*fp/q;hs=fp*(1+H*k);xs=1+H*k
            if r<H or (r==H and not outer_side):
                J=1+r*k;pt=np.array([s-r*fp/q,f+r/q]);jac=np.array([[J,-fp/q],[fp*J,1/q]])
            else:
                eta=(r-H)/(Z-H);pt=np.array([x,h+eta*(Z-h)]);jac=np.array([[xs,0],[hs*(1-eta),(Z-h)/(Z-H)]])
            return f,fp,fpp,fppp,q,k,h,hs,xs,pt,jac
        def velocity(s,r,foot=0):
            f,fp,fpp,fppp,q,k,h,hs,dx,pt,jac=geom(s,r,foot);nz=1/q;nzs=fp*k
            amp=1.;amps=0.;ampss=0.
            if adjusted:
                nzss=-(fpp*fpp+fp*fppp)/q**3+3*fp*fp*fpp*fpp/q**5;hss=fpp+H*nzss
                amp=radial(h)/(UH*nz);amps=radial(h,1)*hs/(UH*nz)-radial(h)*nzs/(UH*nz*nz)
                ampss=(radial(h,2)*hs*hs+radial(h,1)*hss)/(UH*nz)-2*radial(h,1)*hs*nzs/(UH*nz*nz)-radial(h)*nzss/(UH*nz*nz)+2*radial(h)*nzs*nzs/(UH*nz**3)
            if r<=H:
                ps=amps*radial(r,-1);pr=amp*radial(r);inv=np.linalg.inv(jac)
                return np.array([ps*inv[0,1]+pr*inv[1,1],-ps*inv[0,0]-pr*inv[1,0]])
            z=pt[1]
            if z>=B:return np.array([radial(z),0.])
            FH=radial(H,-1);C=amp*FH-radial(h,-1);Cs=amps*FH-radial(h)*hs
            if adjusted:
                I=fp/(q*q*(1+H*k));ks=-fppp/q**3+3*fp*fpp*fpp/q**5;qs=fp*fpp/q
                Is=fpp/(q*q*(1+H*k))-fp*(2*q*qs*(1+H*k)+q*q*H*ks)/(q**4*(1+H*k)**2)
                D=amps*FH*I;Ds=ampss*FH*I+amps*FH*Is
            else:D=UH*nz-radial(h);Ds=UH*nzs-radial(h,1)*hs
            width=B-h;t=(z-h)/width;ts=-hs*(1-t)/width
            h00=1-3*t*t+2*t*t*t;h10=t-2*t*t+t*t*t;d00=-6*t+6*t*t;d10=1-4*t+3*t*t
            u=radial(z)+C*d00/width+D*d10
            ps=Cs*h00+C*d00*ts+(-hs*D+width*Ds)*h10+width*D*d10*ts
            return np.array([u,-ps/dx])
        def kref(r):
            if r>=H:return K0
            def integrand(z):
                ec=A*(z/v['cutoff'] if z<=v['cutoff'] else 1+math.log(z/v['cutoff']))
                prod=tau*tau/slope*(math.log1p(z/delta)+(nu/slope)*(1/(z+delta)-1/delta))
                return (ec-prod)/(slope*(z+delta))
            total=0
            knots=sorted(set([r,H]+([v['cutoff']] if r<v['cutoff']<H else [])))
            for lo,hi in zip(knots,knots[1:]):total+=sum(w*integrand(lo+(x+1)*(hi-lo)/2)*(hi-lo)/2 for x,w in zip(gx,gw))
            return K0-total
        assert abs(kref(0)-v['k_ground'])<1e-11
        from xml.etree import ElementTree as ET
        xml=ET.parse(a.root/name/'case.xml').getroot();verts={int(e.get('ID')):list(map(float,e.text.split())) for e in xml.findall('./GEOMETRY/VERTEX/V')};edges={int(e.get('ID')):list(map(int,e.text.split())) for e in xml.findall('./GEOMETRY/EDGE/E')};boxes={}
        for e in xml.findall('./GEOMETRY/ELEMENT/Q'):
            pts=[verts[i] for j in map(int,e.text.split()) for i in edges[j]]
            boxes[int(e.get('ID'))]=[(min(x[0] for x in pts),max(x[0] for x in pts)),(min(x[1] for x in pts),max(x[1] for x in pts))]
        def law(g,k,eps,physical,inner=False,nt=None):
            iso=2*k/3*np.eye(3);strain=g+g.T
            if inner:parts=dict(isotropic=iso,linear=-nt*strain,quadratic=np.zeros((3,3)),cubic=np.zeros((3,3)),regularization=np.zeros((3,3)))
            else:
                o=indexed(g,k,eps,'physical');sr=np.zeros((3,3));sr[0,2]=sr[2,0]=us/(kap*physical[1]);gamma=kap*us*physical[1]
                parts=dict(isotropic=iso,linear=-float(o['nut'].real)*strain,quadratic=o['quadratic'].real,cubic=.4*o['cubic'].real,regularization=-gamma*(strain-sr))
            stress=sum(parts.values());eig,vec=np.linalg.eigh(stress);d=vec[:,0]
            return dict(tensor=stress.tolist(),eigenvalues=eig.tolist(),critical_eigenvector=d.tolist(),critical_terms={k:float(d@x@d) for k,x in parts.items()},production=float(-np.sum(stress*g)))
        output=[]
        for region in ('inner','outer'):
            w=record['seed_constitutive']['minimum_'+region];s,r=w['coordinate'][0],w['coordinate'][2];bounds=boxes[int(w['geometry'])];foot=(-1 if bounds[0][1]<=s+1e-12 else 1) if abs(abs(s)-L)<1e-12 else 0
            def grad(scale):
                d=np.zeros((2,2));pt=velocity(s,r,foot)
                for axis in range(2):
                    step=scale*(1e-7 if axis==0 else 1e-8);c=[s,r]
                    def ev(v):
                        a=[s,r];a[axis]=v;return velocity(*a,foot)
                    if c[axis]-step<bounds[axis][0]:d[:,axis]=(-3*pt+4*ev(c[axis]+step)-ev(c[axis]+2*step))/(2*step)
                    elif c[axis]+step>bounds[axis][1]:d[:,axis]=(3*pt-4*ev(c[axis]-step)+ev(c[axis]-2*step))/(2*step)
                    else:d[:,axis]=(ev(c[axis]+step)-ev(c[axis]-step))/(2*step)
                cart=d@np.linalg.inv(geom(s,r,foot,region=='outer')[-1]);G=np.zeros((3,3));G[np.ix_([0,2],[0,2])]=cart;return G
            g=grad(.5);physical=geom(s,r,foot)[-2];actual=np.array(w['velocity_gradient']);nt=slope*(r-v['viscosity_intercept'])
            replay=law(actual,w['k'],w['epsilon'],physical,region=='inner',nt)
            error=np.max(np.abs(np.array(replay['tensor'])-np.array(w['selected']['tensor'])));assert error<1e-11,(name,region,error)
            ak=kref(r) if region=='inner' else K0;ae=A/max(r,w['cutoff']) if region=='inner' else A/r
            exact=law(g,ak,ae,physical,region=='inner',nt)
            one=next((t for t in record['target'] if int(t['geometry'])==int(w['geometry']) and bool(t['inner'])==(region=='inner')),None)
            output.append(dict(region=region,coordinate=w['coordinate'],physical_coordinate=physical.tolist(),geometry=w['geometry'],radial_side=(-1 if region=='inner' else 1) if abs(r-H)<1e-12 else 0,tangential_side=foot,
                native_tensor_replay_error=float(error),native_replay=replay,analytical_radial_target=exact,analytical_gradient=g.tolist(),
                gradient_step_halving_max=float(np.max(abs(g-grad(1)))),analytical_divergence=float(np.trace(g)),
                native_1D_target=one,actual_2D=w,analytical_k=ak,analytical_epsilon=ae,
                semantics='exact incoming AWF/log radial streamfunction; unprojected analytic scalars extended in chart, physical incoming gamma/strain; compared with native 1D target at fixed represented scalars and final projected native state'))
        results[name]=dict(witnesses=output,referenceK_independent_ground=kref(0),scalar_budget=record['scalar_inventory'])
    a.output.write_text(json.dumps(dict(status='passed',scope='independent indexed tensor and analytic radial construction, not an equilibrium residual or ellipticity certificate',results=results),indent=2)+'\n')
if __name__=='__main__':main()
