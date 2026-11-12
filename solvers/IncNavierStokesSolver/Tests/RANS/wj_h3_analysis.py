#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Passive H3 coefficient, principal and quadrature views; no governing solve.

The relocated NTFS mount is an explicit path alias, never a modification of a
producer report/checkpoint. Integration views do not change its signature.
"""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from run_wj_h3 import load,save,sha
from validate_cls_epawf_m3 import basis
from wj_h2_profiles import CD,D0,TOP,NU,constitutive,wall

OLD='/media/jay/新加卷/'
NEW='/media/jay/新加卷1/'
def actual(p):return Path(str(p).replace(OLD,NEW,1))
def identity(p):
    p=actual(p);return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)


def values(report,state,block,geometry,x,rho):
    """Independent modified-Jacobi reconstruction on an explicit owner."""
    e=next(e for e in report['fields'][block]['mesh'] if e['geometry']==geometry)
    x,rho=np.broadcast_arrays(x,rho);shape=x.shape;x=x.ravel();rho=rho.ravel()
    l,u=e['lower'],e['upper'];p=e['degree']
    px,dx=basis(p,2*(x-l[0])/(u[0]-l[0])-1)
    pr,dr=basis(p,2*(rho-l[2])/(u[2]-l[2])-1)
    coeff=np.asarray(state)[np.asarray(e['ids'],int)+int(report['offsets'][block])]*e['sign']
    tensor=lambda a,b:np.einsum('ni,nj->nij',a,b).reshape(len(x),-1)@coeff
    H=report['parameters'][1];J=1 if u[2]<=H else (TOP-D0-H)/(TOP-H)
    v=tensor(pr,px);g=np.zeros((len(x),3))
    g[:,0]=tensor(pr,dx)*2/(u[0]-l[0]);g[:,2]=tensor(dr,px)*2/(u[2]-l[2])/J
    return v.reshape(shape),g.reshape((*shape,3))


def algebra(k,e,G):
    """Vectorized independent full 3D WJ note algebra, physical cubic root."""
    k=np.asarray(k);e=np.asarray(e);G=np.asarray(G);I=np.eye(3)
    B=(k/e)[...,None,None]*(G-np.trace(G,axis1=-2,axis2=-1)[...,None,None]*I/3)
    S=(B+np.swapaxes(B,-1,-2))/2;O=(B-np.swapaxes(B,-1,-2))/2
    tr=lambda a:np.trace(a,axis1=-2,axis2=-1)
    s=tr(S@S);w=tr(O@O);i4=tr(S@O@O)
    N=1.8+np.sqrt(2.7*s)
    for _ in range(32):
        F=(N-1.8)*(N*N-2*w)-2.7*s*N
        N-=F/(3*N*N-3.6*N-2*w-2.7*s)
    Q=5/6*(N*N-2*w)*(2*N*N-w)
    m=lambda b:b[...,None,None]
    a=m(-N*(2*N*N-7*w)/Q)*S+m(-12*i4/N/Q)*(O@O-m(w/3)*I)
    a+=m(-2*(N*N-2*w)/Q)*(S@O-O@S)
    a+=m(-6*N/Q)*(S@O@O+O@O@S-m(2*i4/3)*I)+m(6/Q)*(O@S@O@O-O@O@S@O)
    R=m(k)*(a+2*I/3);nut=3*N/(5*(N*N-2*w))*k*k/e
    return R,nut,N


def point(report,state,geometry,x,rho):
    U=np.zeros(3);G=np.zeros((3,3));v={}
    for b,axis in ((0,0),(1,2)):
        U[axis],G[axis]=values(report,state,b,geometry,x,rho)
    k,gk=values(report,state,2,geometry,x,rho);ep,ge=values(report,state,3,geometry,x,rho)
    pi,_=values(report,state,5,geometry,x,rho)
    R,nut,N=constitutive(float(k),float(ep),G)
    return dict(k=float(k),epsilon=float(ep),U=U.tolist(),G=G.tolist(),pi=float(pi),R=R.tolist(),
        N=float(N),nu_eff=float(nut),nu_d=float(CD*k*k/ep),Praw=float(-np.sum(R*G)),
        Pdev=float(-np.sum(R*G)+2*k*np.trace(G)/3),covariance=np.linalg.eigvalsh(R).tolist(),
        divergence=float(np.trace(G)),geometry=geometry,chart=[x,0,rho],
        physical=[x,0,D0+report['parameters'][1]+(rho-report['parameters'][1])*(TOP-D0-report['parameters'][1])/(TOP-report['parameters'][1])])


def principal(k,e,G):
    """Continuous angular minimization of the in-plane quartic, transverse Q2.

    Differentiate the implicit physical root by complex continuation, not a
    frozen root. The three tensor directions determine the exact quartic.
    Floating point polynomial roots delimit this finite planar search scope.
    """
    R,nut,N=constitutive(k,e,G);L=np.zeros((9,9))
    for j,h in enumerate(np.eye(9)):
        L[:,j]=-constitutive(k,e,G.astype(complex)+1e-27j*h.reshape(3,3),N)[0].imag.ravel()/1e-27/nut
    A=np.zeros((3,3));A[2,0]=1
    B=np.diag([-1.,0,1]);C=np.zeros((3,3));C[0,2]=-1
    V=np.array([A.ravel(),B.ravel(),C.ravel()]);Q=V@(L+L.T)/2@V.T
    a,b,c,d,f=Q[0,0],2*Q[0,1],Q[1,1]+2*Q[0,2],2*Q[1,2],Q[2,2]
    coeff=[-d,4*f-2*c,3*d-3*b,2*c-4*a,b]
    roots=np.roots(np.trim_zeros(coeff,'f')) if np.any(coeff) else []
    angles=[0.,math.pi/2]+[math.atan(float(z.real))%math.pi for z in roots if abs(z.imag)<1e-8]
    def form(t):
        xi=np.array([math.cos(t),0,math.sin(t)]);a=np.array([-xi[2],0,xi[0]]);h=np.outer(a,xi).ravel()
        return float(h@L@h)
    t=min(angles,key=form);xi=np.array([math.cos(t),0,math.sin(t)]);pol=np.array([-xi[2],0,xi[0]])
    v=np.zeros((2,9));v[0,3]=1;v[1,5]=1
    et,vt=np.linalg.eigh(v@(L+L.T)/2@v.T)
    return dict(in_plane=form(t),angle=t,xi=xi.tolist(),polarization=pol.tolist(),a_dot_xi=float(pol@xi),
        transverse=float(et[0]),transverse_xi=[float(vt[0,0]),0,float(vt[1,0])],
        scope='all real stationary directions of the local planar quartic, floating point roots; transverse symmetric 2x2 action')


def outer_loads(report,state,history,alpha,order):
    """Full outer volume loads, same three native tests and actual physical J.

    Exterior loads are held fixed when comparing orders. With NormalOrder=16
    and q<=16 their production rule stays16. Explicit per-element signed
    arrays localize discrepancies; separately normed components are not added.
    """
    H=report['parameters'][1];J=(TOP-D0-H)/(TOP-H)
    support=min(e['upper'][2] for e in report['fields'][0]['mesh'] if e['lower'][2]==H)
    g,w=np.polynomial.legendre.leggauss(order);parts={};elements={}
    for el in report['fields'][3]['mesh']:
        lo,hi=el['lower'],el['upper'];geom=el['geometry']
        x,r=np.meshgrid((lo[0]+hi[0])/2+(hi[0]-lo[0])/2*g,(lo[2]+hi[2])/2+(hi[2]-lo[2])/2*g)
        x=x.ravel();r=r.ravel();weight=(np.outer(w,w)*(hi[0]-lo[0])*(hi[2]-lo[2])/4*J).ravel()
        U=np.zeros((len(x),3));G=np.zeros((len(x),3,3))
        for b,axis in ((0,0),(1,2)):U[:,axis],G[:,axis]=values(report,state,b,geom,x,r)
        k,gk=values(report,state,2,geom,x,r);e,ge=values(report,state,3,geom,x,r);pi,_=values(report,state,5,geom,x,r)
        hist=np.array([values(report,history,b,geom,x,r)[0] for b in range(4)]).T
        R,nt,N=algebra(k,e,G);P=-np.sum(R*G,axis=(1,2));dk=CD*k*k/e
        T=-R+(2*k/3-pi)[:,None,None]*np.eye(3);div=np.trace(G,axis1=1,axis2=2)
        phi=np.array([np.ones_like(x),x/5,np.maximum(0,(support-r)/(support-H))]);grad=np.zeros((3,len(x),3));grad[1,:,0]=.2
        grad[2,:,2]=np.where(r<support,-1/(support-H)/J,0)
        local={}
        def add(name,mode,block,v):local.setdefault(name,np.zeros((3,6)))[mode,block]+=float(weight@v)
        for m in range(3):
            adv=np.sum(grad[m]*U,axis=1)
            for b,axis in ((0,0),(1,2)):
                add('time_mass',m,b,alpha*phi[m]*(U[:,axis]-hist[:,b]))
                add('outer_momentum_stress',m,b,np.sum(grad[m]*T[:,axis],axis=1))
                add('conservative_advection',m,b,-adv*U[:,axis])
            for b,f,gf,diff in ((2,k,gk,dk),(3,e,ge,dk/1.3)):
                add('time_mass',m,b,alpha*phi[m]*(f-hist[:,b]))
                add('conservative_advection',m,b,-adv*f)
                add('k_diffusion' if b==2 else 'epsilon_diffusion',m,b,diff*np.sum(grad[m]*gf,axis=1))
            add('outer_k_production',m,2,-phi[m]*P);add('outer_k_destruction',m,2,phi[m]*e)
            add('epsilon_production',m,3,-phi[m]*1.44*P*e/k);add('epsilon_destruction',m,3,phi[m]*1.82*e*e/k)
            add('continuity',m,5,phi[m]*div)
        for name,v in local.items():parts[name]=parts.get(name,np.zeros((3,6)))+v
        elements[str(geom)]={n:v.tolist() for n,v in local.items()}
    return dict(order=order,components={n:v.tolist() for n,v in parts.items()},elements=elements,
        regional=sum(parts.values()).tolist(),scope='outer volume only, full epsilon domain; fixed exterior/inner loads not reinterpreted')


def quadrature(report,state,history,alpha,native=None):
    views={q:outer_loads(report,state,history,alpha,q) for q in (8,12,16,24,32)}
    def compare(a,b):
        changes={n:np.array(a['components'][n])-b['components'][n] for n in a['components']}
        tests={n:bool(np.all(abs(v)<=2e-9+2e-8*abs(np.asarray(b['components'][n])))) for n,v in changes.items()}
        regional=np.asarray(a['regional'])-b['regional'];tests['regional_sum']=bool(np.all(abs(regional)<=4e-9+2e-8*abs(np.asarray(b['regional']))))
        worst=max(((float(np.max(abs(v))),n) for n,v in changes.items()))
        element=max((float(np.max(abs(np.asarray(a['elements'][g][n])-b['elements'][g][n]))),g,n)
            for g in a['elements'] for n in a['elements'][g])
        return dict(status='PASS' if all(tests.values()) else 'UNRESOLVED',conditions=tests,
            component_errors={n:float(np.max(abs(v))) for n,v in changes.items()},regional_max=float(np.max(abs(regional))),
            worst_component=worst,worst_element=element)
    comparisons={f'{q}-32':compare(views[q],views[32]) for q in (8,12,16,24)}
    # At most one additional reference order, only if24/32 is unresolved.
    if comparisons['24-32']['status']!='PASS':
        views[48]=outer_loads(report,state,history,alpha,48)
        comparisons['32-48']=compare(views[32],views[48]);reference=48
    else:reference=32
    qualified=[q for q in (8,12,16) if compare(views[q],views[reference])['status']=='PASS']
    return dict(reference=reference,reference_status=comparisons.get('32-48',comparisons['24-32'])['status'],
        smallest_adequate=min(qualified) if qualified else None,comparisons=comparisons,views=views)


def state_context(report,step=None,rejected=False):
    if rejected:
        event=report['physical_range_events'][-1]
        return event['state'],event['history'],event['alpha'],event['time'],event['wj_packet']
    step=report['steps'] if step is None else step
    first=report['steps']+1-len(report['accepted_states']);i=step-first
    state=report['accepted_states'][i];packet=report['wj_audit_packets'][i]
    if step==0:return state,state,0.,0.,packet
    return state,packet['callback_history'],packet['alpha_mass'],packet['time'],packet


def pressure_complement(report,state):
    """Offline L2 projection into actual continuous regional Q1 pressure space.

    This diagnostic mass solve never changes velocity, histories or native solve
    counts. Its complement includes interelement continuity constraints as well
    as polynomial modes absent from broken Q1.
    """
    out={};H=report['parameters'][1]
    g,w=np.polynomial.legendre.leggauss(6)
    for block in (4,5):
        mesh=report['fields'][block]['mesh'];ids=sorted({i for el in mesh for i in el['ids']});idx={v:i for i,v in enumerate(ids)}
        M=np.zeros((len(ids),len(ids)));b=np.zeros(len(ids));norm=0;broken=0.;detail=[]
        for el in mesh:
            lo,hi=el['lower'],el['upper'];x,r=np.meshgrid((lo[0]+hi[0])/2+(hi[0]-lo[0])/2*g,(lo[2]+hi[2])/2+(hi[2]-lo[2])/2*g)
            x=x.ravel();r=r.ravel();J=1 if hi[2]<=H else (TOP-D0-H)/(TOP-H)
            weight=(np.outer(w,w)*(hi[0]-lo[0])*(hi[2]-lo[2])/4*J).ravel()
            div=values(report,state,0,el['geometry'],x,r)[1][:,0]+values(report,state,1,el['geometry'],x,r)[1][:,2]
            px=basis(1,2*(x-lo[0])/(hi[0]-lo[0])-1)[0];pr=basis(1,2*(r-lo[2])/(hi[2]-lo[2])-1)[0]
            B=np.einsum('ni,nj->nij',pr,px).reshape(len(x),4)*el['sign'];me=B.T@(weight[:,None]*B);be=B.T@(weight*div)
            ii=[idx[v] for v in el['ids']];M[np.ix_(ii,ii)]+=me;b[ii]+=be
            local=float(weight@(div*div));proj=float(be@np.linalg.solve(me,be));norm+=local;broken+=proj
            detail.append(dict(geometry=el['geometry'],divergence_L2=math.sqrt(local),broken_Q1_complement_L2=math.sqrt(max(0,local-proj))))
        proj=float(b@np.linalg.solve(M,b));out[str(block)]=dict(divergence_L2=math.sqrt(norm),continuous_Q1_projection_L2=math.sqrt(max(0,proj)),
            continuous_Q1_complement_L2=math.sqrt(max(0,norm-proj)),broken_Q1_complement_L2=math.sqrt(max(0,norm-broken)),
            maximum_pressure_test_moment=float(max(abs(b))),integrated_net_flux=float(sum(b)),elements=detail)
    return out


def interface_view(report,state,x):
    """Explicit common interface, both regional and horizontal face owners."""
    H=report['parameters'][1];h=report['parameters'][2];result=[]
    inner=[el for el in report['fields'][0]['mesh'] if el['upper'][2]==H and el['lower'][0]-1e-12<=x<=el['upper'][0]+1e-12]
    outer=[el for el in report['fields'][0]['mesh'] if el['lower'][2]==H and el['lower'][0]-1e-12<=x<=el['upper'][0]+1e-12]
    for ei in inner:
        eo=next(el for el in outer if el['lower'][0]==ei['lower'][0])
        ui=np.zeros(3);gi=np.zeros((3,3))
        for b,axis in ((0,0),(1,2)):ui[axis],gi[axis]=values(report,state,b,ei['geometry'],x,H)
        ki,gk=values(report,state,2,ei['geometry'],x,H);pi,gpi=values(report,state,4,ei['geometry'],x,H)
        o=point(report,state,eo['geometry'],x,H);K=o['k'];w=wall(K,H,h);nt=w['a']*max(H-w['zv'],0);D=NU+nt
        Tw=D*(gi+gi.T)-float(pi)*np.eye(3);To=-np.array(o['R'])+(2*K/3-o['pi'])*np.eye(3)
        gko=values(report,state,2,eo['geometry'],x,H)[1];epsH=w['A']/H
        # Initial reference traction mismatch remains visible. No strong row is
        # added to enforce these pointwise defects in the governing system.
        result.append(dict(x=x,geometry_inner_outer=[ei['geometry'],eo['geometry']],K=K,cutoff=w['rc'],
            U_inner=ui.tolist(),U_outer=o['U'],k_inner=float(ki),pi_inner=float(pi),pi_outer=o['pi'],
            traction_defect=(Tw[:,2]-To[:,2]).tolist(),normal_traction_defect=float(Tw[2,2]-To[2,2]),
            k_flux_defect=float(D*gk[2]-o['nu_d']*gko[2]),mass_flux_defect=float(ui[2]-o['U'][2]),
            epsilon_actual=o['epsilon'],epsilon_target=epsH,epsilon_defect=o['epsilon']-epsH,
            mechanical_jump=float(pi-2*ki/3-o['pi']+2*K/3),
            roughness_force_x=float(gpi[0]-2*gk[0]/3+2*ui[0]*gi[0,0]),
            inner_production=float(.5*nt*np.sum((gi+gi.T)**2)-2*ki*np.trace(gi)/3)))
    return result
