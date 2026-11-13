#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent arbitrary-precision CT6 continuum oracle (no native loads).

mpmath differentiation and a Lambert-W cutoff are independent of the C++
Taylor algebra and iterative cutoff implementation. Run with the owned guard.
"""
import argparse, json, os, math
from pathlib import Path
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key]='1'
import mpmath as mp


def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mat(a,b):return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def trace(a):return sum(a[i][i] for i in range(3))
def matrix(function):return [[function(i,j) for j in range(3)] for i in range(3)]
def plain(x):
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    if isinstance(x,(mp.mpf,mp.mpc)):return float(x)
    return x


class Continuous:
    H=mp.mpf('.5');nu=mp.mpf('1.5e-5');hr=mp.mpf('.05')
    def __init__(self,model,amplitude='.005'):
        self.model=model;self.amplitude=mp.mpf(amplitude)
        self.qm=((57+mp.sqrt(6249))/40 if model=='A93' else
                 (3+mp.sqrt(61 if model=='A95' else 41))/(2*mp.sqrt(2)))
        self.sigma=mp.mpf('1.2' if model=='E94E95' else '1.3')
        self.kappa=mp.sqrt(self.sigma*(mp.mpf('1.9')-mp.mpf('.43')*self.qm if model=='E94E95' else mp.mpf('.48'))/self.qm)
        self.K0=self.qm*mp.mpf('.16')

    def geometry(self,q,r):
        h=self.amplitude*mp.exp(-q*q);hp=-2*q*h;hpp=(4*q*q-2)*h
        a=mp.sqrt(1+hp*hp);gamma=-hpp/a**3;J=1+r*gamma
        return dict(a=a,gamma=gamma,J=J,t=[1/a,0,hp/a],n=[-hp/a,0,1/a],point=[q-r*hp/a,0,h+r/a])

    def grad(self,f,q,r,t):
        g=self.geometry(q,r);dq=mp.diff(lambda s:f(s,r,t),q);dr=mp.diff(lambda z:f(q,z,t),r)
        return [g['t'][j]*dq/(g['a']*g['J'])+g['n'][j]*dr for j in range(3)]

    def gradient(self,u,q,r,t):
        return [self.grad(lambda a,b,c:u(a,b,c)[i],q,r,t) for i in range(3)]

    def divergence(self,u,q,r,t):
        return trace(self.gradient(u,q,r,t))

    def K(self,q,t):return self.K0*(1+mp.mpf('.02')*q+mp.mpf('.01')*q*q)*(1+mp.mpf('.05')*t+mp.mpf('.1')*t*t)
    def A(self,k):return (k/self.qm)**mp.mpf('1.5')/self.kappa
    def b(self,t):return mp.mpf('.2')*(1+mp.mpf('.2')*t+mp.mpf('.5')*t*t)
    def nt(self,q,r,t):
        k=self.K(q,t);x=self.hr*mp.sqrt(k)/(70*self.nu)
        exponent=max(mp.mpf('.5')-mp.mpf('.4')*x**mp.mpf('.7'),1-mp.mpf('.79')*x**(-mp.mpf('.28')))
        B=mp.mpf('10.7')*(1-x**exponent);zv=self.nu*B/mp.sqrt(k)
        assert zv<0
        return self.kappa*mp.sqrt(k/self.qm)*(r-zv)
    def D(self,q,r,t):return self.nu+self.nt(q,r,t)

    def outer(self,k,e,G):
        divergence=trace(G);S=matrix(lambda i,j:(G[i][j]+G[j][i])/2-(divergence/3 if i==j else 0))
        O=matrix(lambda i,j:(G[i][j]-G[j][i])/2)
        ss=sum(v*v for row in S for v in row);oo=sum(v*v for row in O for v in row)
        tau=k/e;eta=tau*mp.sqrt(2*ss);xi=tau*mp.sqrt(2*oo);rad=None
        if self.model=='A93':cmu=(mp.mpf(2)/3)/(mp.mpf('1.25')+eta+mp.mpf('.9')*xi)
        else:
            largest=(S[0][0]+S[2][2])/2+mp.sqrt(((S[0][0]-S[2][2])/2)**2+S[0][2]**2)
            cmu=1/(mp.mpf('6.5' if self.model=='A95' else '4')+3*largest*tau*mp.sqrt((ss+oo)/ss))
            if self.model=='A95':rad=1-9*cmu**2*tau**2*ss
        nt=cmu*k*k/e;s2=mat(S,S);o2=mat(O,O);os=mat(O,S);so=mat(S,O)
        Q=matrix(lambda i,j:0)
        if self.model=='A93':
            Q=matrix(lambda i,j:k*tau*tau/(1000+eta**3)*(3*(s2[i][j]-(ss/3 if i==j else 0))-19*(o2[i][j]+(oo/3 if i==j else 0))+15*(os[i][j]-so[i][j])))
        elif self.model=='A95':
            if rad<0:raise ValueError('independent actual A95 domain failure')
            cq=mp.sqrt(rad)/(1+6*tau*tau*mp.sqrt(ss*oo))
            Q=matrix(lambda i,j:2*cq*k*tau*tau*(os[i][j]-so[i][j]))
        R=matrix(lambda i,j:(2*k/3 if i==j else 0)-2*nt*S[i][j]+Q[i][j])
        P=-sum(R[i][j]*G[i][j] for i in range(3) for j in range(3))
        ep=max(mp.mpf('.43'),eta/(eta+5))*mp.sqrt(2*ss)*e if self.model=='E94E95' else mp.mpf('1.44')*e*P/k
        ed=mp.mpf('1.9' if self.model=='E94E95' else '1.92')*e*e/k
        return dict(R=R,Q=Q,nt=nt,P=P,ep=ep,ed=ed,radicand=rad)

    def base(self,q,r,t):
        g=self.geometry(q,r);ut=self.b(t)*(1+mp.mpf('.02')*q)*r
        w=-self.b(t)*mp.mpf('.02')*r*r/(2*g['a']*g['J'])
        return [ut*g['t'][j]+w*g['n'][j] for j in range(3)]

    def ray(self,q,t):
        H=self.H;g=self.geometry(q,H);k=self.K(q,t);D=self.D(q,H,t)
        G=self.gradient(self.base,q,H,t);law=self.outer(k,self.A(k)/H,G)
        F=matrix(lambda i,j:-law['R'][i][j]+(2*k/3 if i==j else 0))
        traction=lambda left,right:sum(left[i]*F[i][j]*right[j] for i in range(3) for j in range(3))
        w=lambda s,r,t:-self.b(t)*mp.mpf('.02')*r*r/(2*self.geometry(s,r)['a']*self.geometry(s,r)['J'])
        ut=self.b(t)*(1+mp.mpf('.02')*q)*H
        c=traction(g['t'],g['n'])/D-(self.b(t)*(1+mp.mpf('.02')*q)+mp.diff(lambda s:w(s,H,t),q)/(g['a']*g['J'])-g['gamma']*ut/g['J'])
        jump=traction(g['n'],g['n'])-2*D*mp.diff(lambda r:w(q,r,t),H)
        l=law['nt']/D*(mp.mpf('.03')*self.K0/H)
        a=self.kappa*mp.sqrt(k/self.qm);d=self.D(q,0,t);V=d+a*H;lg=mp.log(V/d);t0=k/self.qm
        p0=t0*t0/a*(lg+self.nu/V-self.nu/d)
        p1=t0*t0/a**2*(a*H-(d+self.nu)*lg+d*self.nu*(1/d-1/V))
        C=1+mp.log(H)+g['gamma']*H-(p0+g['gamma']*p1)/self.A(k)
        gamma=g['gamma'];cut=mp.exp(C) if gamma==0 else 2/gamma*mp.lambertw(gamma*mp.exp(C)/2)
        return dict(c=c,l=l,jump=jump,cutoff=cut)

    def fields(self,q,r,t,outer):
        g=self.geometry(q,r);ray=self.ray(q,t);H=self.H
        phi=r*r*(r-H)**2/(2*H*H);phiR=r*(r-H)*(2*r-H)/(H*H)
        ut=self.b(t)*(1+mp.mpf('.02')*q)*r+(0 if outer else ray['c']*phiR)
        w=-(self.b(t)*mp.mpf('.02')*r*r/2+(0 if outer else mp.diff(lambda s:self.ray(s,t)['c'],q)*phi))/(g['a']*g['J'])
        u=[ut*g['t'][j]+w*g['n'][j] for j in range(3)]
        k=self.K(q,t)+(mp.mpf('.03')*self.K0/H*(r-H) if outer else ray['l']*(r*r-H*H)/(2*H))
        eps=self.A(self.K(q,t))/(r if outer else max(r,ray['cutoff']))
        pi=mp.mpf('.03')*(q+mp.mpf('.2')*r)+mp.mpf('.02')*t+mp.mpf('.01')*t*t+(ray['jump'] if outer else 0)
        return dict(u=u,k=k,epsilon=eps,Pi=pi,ray=ray)

    def sample(self,q,r,t,outer,sources=False):
        q,r,t=map(mp.mpf,(str(q),str(r),str(t)))
        fields=self.fields(q,r,t,outer);u=lambda a,b,c:self.fields(a,b,c,outer)['u']
        G=self.gradient(u,q,r,t);k=fields['k'];e=fields['epsilon'];pi=fields['Pi']
        g=self.geometry(q,r)
        if outer:
            law=self.outer(k,e,G);D=law['nt'];de=D/self.sigma
            T=matrix(lambda i,j:-law['R'][i][j]+(2*k/3-pi if i==j else 0))
        else:
            nt=self.nt(q,r,t);D=self.nu+nt;de=0
            R=matrix(lambda i,j:(2*k/3 if i==j else 0)-nt*(G[i][j]+G[j][i]))
            law=dict(R=R,Q=matrix(lambda i,j:0),P=-sum(R[i][j]*G[i][j] for i in range(3) for j in range(3)),ep=0,ed=0,radicand=None)
            T=matrix(lambda i,j:D*(G[i][j]+G[j][i])-(pi if i==j else 0))
        gradients=[self.grad(lambda a,b,c:self.fields(a,b,c,outer)[key],q,r,t) for key in ('k','epsilon','Pi')]
        dtU=[mp.diff(lambda tt:u(q,r,tt)[j],t) for j in range(3)]
        dt=[mp.diff(lambda tt:self.fields(q,r,tt,outer)[key],t) for key in ('k','epsilon')]
        result=dict(coordinate=[q,0,r],time=t,outer=outer,
            values=[fields['u'][0],fields['u'][2],k,e,pi],
            gradients=[G[0],G[2],*gradients],G=G,R=law['R'],Q=law['Q'],T=T,
            P=law['P'],Dk=D,Depsilon=de,epsilon_production=law['ep'],epsilon_destruction=law['ed'],
            A95_radicand=law['radicand'],divergence=trace(G),normal=g['n'],tangent=g['t'],
            K=self.K(q,t),cutoff=fields['ray']['cutoff'],velocity_time=dtU,k_time=dt[0],epsilon_time=dt[1])
        if sources:
            # These are continuous physical divergences, not FE residuals.
            def stress(a,b,c):return self.sample(a,b,c,outer,False)['T']
            body=[]
            for i in range(3):
                flux=lambda a,b,c:[u(a,b,c)[i]*u(a,b,c)[j]-stress(a,b,c)[i][j] for j in range(3)]
                body.append(dtU[i]+self.divergence(flux,q,r,t))
            def sf(a,b,c,key):
                v=self.fields(a,b,c,outer);grad=self.grad(lambda aa,bb,cc:self.fields(aa,bb,cc,outer)[key],a,b,c)
                sample=self.sample(a,b,c,outer,False);diff=sample['Dk' if key=='k' else 'Depsilon']
                return [v['u'][j]*v[key]-diff*grad[j] for j in range(3)]
            result.update(required_body=body,k_source=dt[0]+self.divergence(lambda a,b,c:sf(a,b,c,'k'),q,r,t)-law['P']+e,
                          epsilon_source=dt[1]+self.divergence(lambda a,b,c:sf(a,b,c,'epsilon'),q,r,t)-law['ep']+law['ed'] if outer else 0)
        return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--sources',action='store_true');a=p.parse_args();mp.mp.dps=35
    a.output.mkdir(parents=True,exist_ok=False);rows=[];interfaces=[]
    for model in ('A95','A93','E94E95'):
        c=Continuous(model)
        for q in ('-.9','0','.9'):
            for outer,r in ((False,'.23'),(True,'.72')):
                row=c.sample(q,r,'.0032',outer,a.sources)
                eigenvalues=mp.eigsy(mp.matrix(row['R']),eigvals_only=True)
                row['covariance_eigenvalues']=list(eigenvalues)
                rows.append(dict(model=model,**plain(row)))
            inner=c.sample(q,'.5','.0032',False);outer=c.sample(q,'.5','.0032',True)
            n=inner['normal'];tk=lambda v:[dot(row,n) for row in v['T']]
            defects=dict(velocity=max(abs(x-y) for x,y in zip(inner['values'][:2],outer['values'][:2],strict=True)),
                k=abs(inner['values'][2]-outer['values'][2]),
                traction=max(abs(x-y) for x,y in zip(tk(inner),tk(outer),strict=True)),
                k_flux=abs(inner['Dk']*dot(inner['gradients'][2],n)-outer['Dk']*dot(outer['gradients'][2],n)),
                epsilon_trace=abs(inner['values'][3]-outer['values'][3]))
            interfaces.append(dict(model=model,q=float(q),defects=plain(defects)))
            (a.output/'progress.json').write_text(json.dumps(dict(samples=rows,interfaces=interfaces),allow_nan=False,indent=2)+'\n')
    maximum=max([abs(row['divergence']) for row in rows]+[v for row in interfaces for v in row['defects'].values()])
    domain=dict(minimum_k=min(row['values'][2] for row in rows),
        minimum_epsilon=min(row['values'][3] for row in rows),
        minimum_covariance=min(min(row['covariance_eigenvalues']) for row in rows),
        minimum_A95_radicand=min(row['A95_radicand'] for row in rows if row['A95_radicand'] is not None))
    result=dict(schema='CT6IndependentContinuousOracleV1',status='pass' if maximum<=2e-10 and min(domain.values())>0 else 'fail',
                precision_digits=mp.mp.dps,maximum_primitive_defect=maximum,samples=rows,interfaces=interfaces,
                sampled_domain=domain,source_evaluation=a.sources,scope='independent physical derivatives; no native operator or evolution qualification')
    (a.output/'results.json').write_text(json.dumps(result,allow_nan=False,indent=2)+'\n')
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
