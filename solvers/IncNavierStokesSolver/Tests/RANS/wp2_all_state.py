# SPDX-License-Identifier: MIT
"""Read-only full-tensor diagnostics from immutable native CG rows.

The actual coefficient-derived gradients and actual retained epsilon are used.
Analytical reference gradients enter only the separate error/reference objects.
"""
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_h2g_constitutive import parts
from run_cls_epawf_picard import require,sha

US=.07426793792419734


def wall(K,r,nu,H,hs,q,kap):
    require(np.all(K>0),'nonpositive endpoint K in retrospective evaluation')
    root=np.sqrt(K);ratio=hs*root/(70*nu);exponent=np.maximum(.5-.4*ratio**.7,1-.79*ratio**(-.28))
    intercept=-10.7*np.expm1(exponent*np.log(ratio));zv=nu*intercept/root;slope=kap/np.sqrt(q)*root
    def moments(z):
        z=np.asarray(z);split=np.minimum(z,np.maximum(0.,zv));length=z-split;D=nu+slope*(split-zv)
        return split/nu+np.log1p(slope*length/D)/slope,split/(nu*nu)+length/(D*(D+slope*length))
    j0,j2=moments(H);A=(K/q)**1.5/kap;budget=(K/q)**2*(j0-nu*j2);rc=H*np.exp(1-budget/A)
    require(np.all(rc>0) and np.all(rc<=H),'independent flat cutoff applicability')
    nt=slope*np.maximum(0.,r-zv);E=A/np.maximum(r,rc)
    p0,p2=moments(r);integral_E=A*np.where(r<=rc,r/rc,1+np.log(np.maximum(r,rc)/rc))
    dk=(integral_E-(K/q)**2*(p0-nu*p2))/(nu+nt)
    return dict(nt=nt,E=E,rc=rc,reference_k_gradient=dk,j0=p0,intercept=zv)


class Rows:
    def __init__(self,metadata):
        self.path=Path(metadata);self.meta=load(self.path);self.a={}
        require(self.meta['version']=='WP2AffineRowsV1','row schema')
        for name,d in self.meta['arrays'].items():
            f=self.path.parent/d['file'];require(sha(f)==d['sha256'],'row binary changed')
            self.a[name]=np.memmap(f,dtype=d['dtype'],mode='r',shape=tuple(int(x) for x in d['shape']))
        self.size=int(sum(self.meta['sizes']));self.offset=[int(x) for x in self.meta['offsets']]
        self.inner=self.a['regions'][:,0]==0;self.weight=np.asarray(self.a['weights'][:,0]);self.physical=np.asarray(self.a['coordinates'])
    def evaluate(self,s,family,kind='phi',offset=0,mask=None):
        require(len(s)==self.size,'state/geometry dimensions differ')
        ids=self.a[family+'_ids'] if mask is None else self.a[family+'_ids'][mask]
        constant=self.a[family+'_constant'] if mask is None else self.a[family+'_constant'][mask]
        phi=self.a[family+'_'+kind] if mask is None else self.a[family+'_'+kind][mask]
        centre=s[ids[:,0]+offset]/constant[:,0];delta=s[ids+offset]-centre[:,None]*constant
        value=centre.astype(np.longdouble) if kind=='phi' else np.zeros(len(ids),dtype=np.longdouble)
        for j in range(ids.shape[1]):value+=phi[:,j].astype(np.longdouble)*delta[:,j].astype(np.longdouble)
        return np.asarray(value,dtype=np.float64)
    def tuple(self,state,report):
        s=np.asarray(state,dtype=np.float64);o=self.offset;U=self.evaluate(s,'v',offset=o[0]);W=self.evaluate(s,'v',offset=o[1]);k=self.evaluate(s,'v',offset=o[2])
        e=np.ones(len(k));e[~self.inner]=self.evaluate(s,'e',mask=~self.inner);pi=self.evaluate(s,'p')
        G=np.zeros((len(k),3,3));gk=np.zeros((len(k),3));ge=np.zeros_like(gk);gp=np.zeros_like(gk)
        for axis,kind in [(0,'dx'),(2,'dz')]:
            G[:,0,axis]=self.evaluate(s,'v',kind,o[0]);G[:,2,axis]=self.evaluate(s,'v',kind,o[1]);gk[:,axis]=self.evaluate(s,'v',kind,o[2]);ge[~self.inner,axis]=self.evaluate(s,'e',kind,mask=~self.inner)
            gp[:,axis]=self.evaluate(s,'p',kind)
        kh=self.evaluate(s,'h');require(len(k)%len(kh)==0,'wall/volume ordering');K=np.repeat(kh,len(k)//len(kh))
        nu,H,hs=report['parameters'][:3];q,kap=report['parameters'][6:8]
        inner_wall=wall(K,np.minimum(self.physical[:,2],H),nu,H,hs,q,kap)
        e[self.inner]=inner_wall['E'][self.inner]
        law='calibrated_standard' if report['model_specification']['name']=='StandardKEpsilon' else 'regularized_CLS'
        closure=parts(G,k,e,self.inner,inner_wall['nt'],self.physical,US,kap,law)
        R=closure['tensor'];T=-R+(2*k/3-pi)[:,None,None]*np.eye(3)
        T[self.inner]+=nu*(G[self.inner]+G[self.inner].transpose(0,2,1))
        values=np.column_stack((U,W,k,e,pi,pi-2*k/3,closure['nu_t'],T[:,0,2],T[:,2,2],closure['production'],G[:,0,0],G[:,0,2],G[:,2,0],G[:,2,2]))
        return dict(values=values,R=R,T=T,G=G,gk=gk,ge=ge,gp=gp,wall=inner_wall,endpoint_K=kh)
    def check_native(self,report):
        a=report['wp2_common_samples']
        require(a['physical_z_partition']==self.meta['physical_z_partition'] and
                np.array_equal(np.asarray(a['weights']),self.weight) and
                np.array_equal(np.asarray(a['regions']),self.a['regions'][:,0]) and
                [int(i) for i in report['offsets']]==self.offset,'native/geometry catalogue partition, measures, region or block identity differs')
        d=self.tuple(report['state'],report);native=np.asarray(report['wp2_common_samples']['values']);difference=np.max(abs(d['values']-native),axis=0)
        require(max(difference[6:10])<2e-10 and max(difference[:6])<2e-11 and max(difference[10:])<2e-8,'native/retrospective tuple differs')
        require(np.max(abs(d['endpoint_K']-np.asarray(report['physical_samples']['wall_K'])))<2e-12,'endpoint sampler differs')
        # Independent differentiation of the represented scalar polynomials
        # also checks the reused U/K row indexing and outer-epsilon rows.
        from numpy.polynomial.legendre import leggauss,legvander,legder
        old,w=leggauss(8);V=legvander(old,3);P=(np.arange(4)*2+1)[:,None]/2*V.T*w[None,:]
        D=legvander(old,2)@legder(np.eye(4),axis=0)@P;F=V@P
        zk=np.asarray(report['wp2_common_samples']['physical_z_partition']);nz=len(zk)-1
        xk=np.array([-.2,-.1,-.05,0,.05,.1,.6]);v=native.reshape(6,8,nz,8,14).transpose(0,2,1,3,4)[:,:,:,:,[2,3]]
        dx=np.einsum('ai,xyijf,bj->xyabf',D,v,F,optimize=True)*2/np.diff(xk)[:,None,None,None,None]
        dz=np.einsum('ai,xyijf,bj->xyabf',F,v,D,optimize=True)*2/np.diff(zk)[None,:,None,None,None]
        dx=dx.transpose(0,2,1,3,4).reshape(-1,2);dz=dz.transpose(0,2,1,3,4).reshape(-1,2)
        errors=[float(max(abs(dx[:,0]-d['gk'][:,0]))),float(max(abs(dz[:,0]-d['gk'][:,2]))),
                float(max(abs(dx[~self.inner,1]-d['ge'][~self.inner,0]))),float(max(abs(dz[~self.inner,1]-d['ge'][~self.inner,2])))]
        require(max(errors)<2e-8,'scalar derivative/native CG row replay differs')
        return dict(maximum_absolute_by_column=difference.tolist(),scalar_gradient_polynomial_replay_errors=errors,passed=True)
    def full_errors(self,state,report,tuple_data=None):
        d=self.tuple(state,report) if tuple_data is None else tuple_data
        ref=np.asarray(report['wp2_common_samples']['reference_values']);nu,H,hs=report['parameters'][:3];q,kap=report['parameters'][6:8];K0=report['parameters'][5]
        gr=np.zeros_like(d['G']);gr[:,0,2]=ref[:,11];kr=ref[:,2];er=ref[:,3];nt=ref[:,6]
        law='calibrated_standard' if report['model_specification']['name']=='StandardKEpsilon' else 'regularized_CLS'
        rr=parts(gr,kr,er,self.inner,nt,self.physical,US,kap,law)['tensor'];tr=-rr+(2*kr/3-ref[:,4])[:,None,None]*np.eye(3)
        tr[self.inner]+=nu*(gr[self.inner]+gr[self.inner].transpose(0,2,1))
        reference_wall=wall(np.full(len(kr),K0),np.minimum(self.physical[:,2],H),nu,H,hs,q,kap)
        gkr=np.zeros_like(d['gk']);gkr[self.inner,2]=reference_wall['reference_k_gradient'][self.inner]
        ger=np.zeros_like(d['ge']);ger[~self.inner,2]=-er[~self.inner]/self.physical[~self.inner,2]
        regions=[]
        for region in [0,1]:
            selected=self.inner if region==0 else ~self.inner;w=self.weight[selected];volume=sum(w);den=US**2*np.sqrt(volume)
            error_R=d['R'][selected]-rr[selected];error_T=d['T'][selected]-tr[selected]
            R=d['R'][selected];ev=np.minimum(R[:,1,1],.5*(R[:,0,0]+R[:,2,2]-np.hypot(R[:,0,0]-R[:,2,2],2*R[:,0,2])))
            where=np.flatnonzero(selected)[int(np.argmin(ev))]
            gradient_k_den=np.sqrt(np.sum(w[:,None]*gkr[selected]**2)) if region==0 else K0/H*np.sqrt(volume)
            regions.append(dict(reynolds_tensor_relative_RMS=float(np.sqrt(np.sum(w[:,None,None]*error_R**2))/den),
                mean_stress_tensor_relative_RMS=float(np.sqrt(np.sum(w[:,None,None]*error_T**2))/den),
                reynolds_component_relative_RMS=(np.sqrt(np.sum(w[:,None,None]*error_R**2,axis=0))/den).tolist(),
                mean_stress_component_relative_RMS=(np.sqrt(np.sum(w[:,None,None]*error_T**2,axis=0))/den).tolist(),
                covariance_minimum=float(min(ev)),covariance_minimum_point=self.physical[where].tolist(),
                full_k_gradient_error=float(np.sqrt(np.sum(w[:,None]*(d['gk'][selected]-gkr[selected])**2))/gradient_k_den),
                k_gradient_normalization='reference physical gradient L2' if region==0 else 'K0/H times sqrt(volume), since reference gradient is zero',
                full_outer_epsilon_gradient_relative_error=None if region==0 else float(np.sqrt(np.sum(w[:,None]*(d['ge'][selected]-ger[selected])**2)/np.sum(w[:,None]*ger[selected]**2))),
                modified_pressure_gradient_dimensional_RMS=float(np.sqrt(np.sum(w[:,None]*d['gp'][selected]**2)/volume)),
                mechanical_pressure_gradient_error_dimensional_RMS=float(np.sqrt(np.sum(w[:,None]*(d['gp'][selected]-2*d['gk'][selected]/3+2*gkr[selected]/3)**2)/volume)),
                integrated_actual_production=float(np.dot(w,d['values'][selected,9])),
                integrated_eddy_viscosity=float(np.dot(w,d['values'][selected,6]))))
        return dict(regions=regions,endpoint_K_range=[float(min(d['endpoint_K'])),float(max(d['endpoint_K']))],
                    independently_recomputed_cutoff_range=[float(min(d['wall']['rc'])),float(max(d['wall']['rc']))],
                    TKE_integral=float(np.dot(self.weight,d['values'][:,2])),
                    scope='actual native coefficient gradients and actual outer epsilon; full selected 3x3 tensors; stress errors normalized by u_tau^2; analytical reference evaluated separately; complements native sided hard inventory')

    def volume_audit(self,report,tuple_data,order=28):
        """Independent polynomial evaluation on the frozen union-panel GL28.

        Stresses are reevaluated, never interpolated. Native exact sided audits
        remain mandatory and complementary to this interior volume grid.
        """
        from wp2_common_quadrature import resample
        require(order==28,'WP2 retrospective denser grid is frozen at GL28')
        sample=dict(report['wp2_common_samples'],values=tuple_data['values'])
        v=resample(dict(report,wp2_common_samples=sample),order)
        R=v['R'];ev=np.minimum(R[:,1,1],.5*(R[:,0,0]+R[:,2,2]-np.hypot(R[:,0,0]-R[:,2,2],2*R[:,0,2])))
        nu,H,hs=report['parameters'][:3];K0,q,kap=report['parameters'][5:8]
        inner=v['regions']==0;E=v['values'][:,3].copy()
        E[inner]=wall(v['endpoint_K'][inner],v['physical'][inner,2],nu,H,hs,q,kap)['E']
        require(np.all(np.isfinite(R)) and np.all(np.isfinite(E)) and np.all(E>0) and np.all(v['values'][:,2]>0),'denser-grid scalar/tensor admissibility')
        witnesses=[];nz=len(sample['physical_z_partition'])-1
        geometries=self.a['geometry'][:,0].reshape(6,8,nz,8)[:,0,:,0];elements=self.a['elements'][:,0].reshape(6,8,nz,8)[:,0,:,0]
        for region in [0,1]:
            ids=np.flatnonzero(v['regions']==region);j=ids[int(np.argmin(ev[ids]))];ix,_,iz,_=np.unravel_index(j,(6,order,nz,order))
            eigen,vec=np.linalg.eigh(R[j]);require(abs(eigen[0]-ev[j])<1e-12,'independent denser-grid eigendecomposition')
            witnesses.append(dict(region=region,physical=v['physical'][j].tolist(),geometry=int(geometries[ix,iz]),element=int(elements[ix,iz]),radial_side=0,tangential_side=0,
                k=float(v['values'][j,2]),epsilon=float(E[j]),epsilon_semantics='dependent inner E_C' if region==0 else 'actual retained outer coefficients',
                gradient=v['G'][j].tolist(),tensor=R[j].tolist(),eigenvalues=eigen.tolist(),critical_eigenvector=vec[:,0].tolist()))
        divergence=v['G'][:,0,0]+v['G'][:,2,2]
        strong=float(np.sqrt(np.dot(v['weights'],divergence**2)/sum(v['weights']))/(US/H))
        require(min(ev)>=1e-4*K0 and strong<=1e-3,'denser-grid covariance/divergence budget')
        return dict(order=order,point_count=len(ev),minimum_covariance_inner_outer=[w['eigenvalues'][0] for w in witnesses],witnesses=witnesses,
            positive_scalars=True,minimum_k=float(min(v['values'][:,2])),minimum_epsilon=float(min(E)),physical_volume=float(sum(v['weights'])),
            strong_divergence_RMS_over_ustar_H=strong,polynomial_fit_error=v['fit_error'],passed=True,
            scope='fixed GL28 on the common union partition; actual affine native polynomials reconstructed from validated native rows, selected tensors freshly recomputed; finite sampling, complements exact native sided inventories; no freshness/selection claim')
