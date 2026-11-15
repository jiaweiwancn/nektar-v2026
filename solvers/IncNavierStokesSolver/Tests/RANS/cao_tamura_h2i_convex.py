# SPDX-License-Identifier: MIT
"""Fixed-k standard stress LMIs and physical convex caps, with LP separation.

The LP is only a necessary relaxation. Passing requires full true PSD and
quadratic scans; obstruction requires an actual-limit checked dual bound.
"""
from pathlib import Path
import time
import numpy as np
from scipy import sparse,linalg
from cao_tamura_h2d_state import load
from cao_tamura_h2h_scalar import Budget,BudgetExhausted,lp
from run_cls_epawf_picard import sha,require

def validate_direction_counts(version,counts):
 require(version in ('FixedKJointContextV1','FixedKJointContextV2','FixedKJointContextV3'),'wrong joint schema')
 nv=counts['velocity'];require(type(nv) is int and counts['k_anchor']==8 and counts['epsilon']==8,'joint scalar/count contract')
 require((version=='FixedKJointContextV1' and nv==41) or (version=='FixedKJointContextV2' and 41<=nv<=47) or (version=='FixedKJointContextV3' and 44<=nv<=128),'joint velocity count/schema mismatch')
 return nv

class Context:
 def __init__(self,path):
  self.path=Path(path);self.doc=load(path);self.meta=self.doc['context'];self.root=Path(self.doc['array_directory']);require(self.meta['version'] in ('FixedKJointContextV1','FixedKJointContextV2','FixedKJointContextV3'),'wrong H2i schema');r=self.doc['receipt'];require(sha(Path(r['output']))==r['output_sha256'],'joint native receipt mismatch');native=load(Path(r['output']));require(native['context']==self.meta,'joint context was changed after export')
  for n in ('geometry_identity','inventory_identity','reference_sha256','defect_sha256','patch_signature'):require(native[n]==self.doc[n],'joint identity mismatch '+n)
  spec=self.meta['model_specification'];require(spec['name']=='StandardKEpsilon' and spec['calibration']=='StandardNeutralCalibratedV1' and spec['coefficients']['Cmu']==.09 and not spec['physical_regularization_active'],'H2i requires calibrated standard with zero regularization')
  counts=self.meta['direction_counts'];nv=validate_direction_counts(self.meta['version'],counts);self.arrays={}
  for name,entry in self.meta['files'].items():self.arrays[name]=self.read(entry)
  self.algebra={name:self.read(entry) for name,entry in self.meta['algebra']['files'].items()};self.blocks=self.meta['observation_blocks'];self.N=len(self.arrays['points'])
  occupied=[]
  for b in self.blocks.values():occupied+=list(range(b['offset'],b['offset']+b['width']))
  require(sorted(occupied)==list(range(self.arrays['points'].shape[1])),'incomplete/overlapping joint blocks')
  for name,shape in (('velocity_basis',[nv,6]),('k_basis',[8]),('epsilon_basis',[8])):require(self.blocks[name]['shape']==shape and self.blocks[name]['width']==int(np.prod(shape)),'joint block/count mismatch '+name)
 def read(self,info):
  p=self.root/info['file'];require(sha(p)==info['sha256'],'joint array changed '+str(p));return np.memmap(p,dtype=info['dtype'],mode='r',shape=tuple(info['shape']))
 def block(self,name):
  b=self.blocks[name];v=self.arrays['points'][:,b['offset']:b['offset']+b['width']];return v[:,0] if b['width']==1 else v.reshape((len(v),)+tuple(b['shape']))
 def retention(self,lower=.501,upper=1.999):
  candidates=[]
  for label in ('k','epsilon'):
   ratio=self.block(label+'_B')/self.block(label+'_cold');d=ratio-1;limit=np.ones(len(d));positive=d>0;negative=d<0;limit[positive]=(upper-1)/d[positive];limit[negative]=(1-lower)/(-d[negative]);j=int(np.argmin(limit));candidates.append((float(limit[j]),label,j,float(ratio[j])))
  maximum=min(1.,min(v[0] for v in candidates));theta=max(0.,np.nextafter(maximum,0)-1e-12);require(0<=theta<=1,'invalid anchor retention');return theta,dict(maximum_retention=maximum,theta=theta,inward_guard=1e-12,limiting=min(candidates),scalar_band=[lower,upper])
 def problem(self,anchor,theta):
  problem=JointProblem(self.block('gradient_B').reshape(-1,4),self.block('velocity_B'),self.block('velocity_cold'),self.block('k_cold')+(theta if anchor=='A' else 0)*(self.block('k_B')-self.block('k_cold')),self.block('epsilon_cold')+theta*(self.block('epsilon_B')-self.block('epsilon_cold')),self.block('k_cold'),self.block('epsilon_cold'),self.block('inner_nu'),self.block('inner').astype(bool),self.block('covariance').astype(bool),self.block('velocity_basis'),self.block('epsilon_basis'),self.algebra['velocity_mass'],[self.algebra['divergence_GL'+str(o)] for o in (16,20,24)],self.meta['model_specification']['K0'])
  problem.geometry=np.empty(self.N,int);m=self.arrays['members'];problem.geometry[m[:,1]]=m[:,5];return problem

class JointProblem:
 def __init__(self,g,v,vc,k,e,kc,ec,nt,inner,cov,VB,EB,mass,div,K0):
  self.g=np.array(g);self.v=np.array(v);self.vcold=np.array(vc);self.k=np.array(k);self.e=np.array(e);self.kc=np.array(kc);self.ec=np.array(ec);self.nt=np.array(nt);self.inner=np.array(inner);self.cov=np.array(cov);self.n=len(k);self.nv=VB.shape[1];self.ne=EB.shape[1];self.dim=self.nv+self.ne;self.K0=K0
  self.VB=np.ascontiguousarray(VB);self.E=np.ascontiguousarray(EB);self.mass=np.array(mass);self.div=[np.array(d) for d in div];self.div_limits=np.square([.06030809942199109,.06030809942199109+1e-8,.06030809942199109+1e-8]);require(np.min(self.k)>0 and np.min(self.e)>0,'joint scalar anchor inadmissible');self.bounds,self.bound_report=self.coefficient_bounds()
 def coefficient_bounds(self):
  G=self.mass[1:,1:];h=self.mass[0,1:];ev=np.linalg.eigvalsh(G);require(ev[0]>1e-10*ev[-1],'velocity mass rank loss');inv=np.linalg.inv(G);centre=-inv@h;radius2=.12**2-self.mass[0,0]+h@inv@h;require(radius2>=0,'immutable B outside mass ellipsoid');rad=np.sqrt(max(0,radius2)*np.diag(inv));lower=centre-rad-1e-10;upper=centre+rad+1e-10
  outer=~self.inner;C=self.E[outer]/self.ec[outer,None];_,_,p=linalg.qr(C.T,mode='economic',pivoting=True);ids=p[:self.ne];Q=C[ids];require(np.linalg.matrix_rank(Q)==self.ne,'epsilon basis rank loss');N=np.linalg.inv(Q);eta=float(np.linalg.norm(np.eye(self.ne)-N@Q,np.inf));require(eta<.01,'epsilon bound inversion unstable');lo=.5-self.e[outer]/self.ec[outer];hi=2-self.e[outer]/self.ec[outer];box=float(np.max(abs(N)@np.maximum(abs(lo[ids]),abs(hi[ids])))/(1-eta))*(1+1e-10)+1e-10
  return list(zip(np.r_[lower,np.full(self.ne,-box)],np.r_[upper,np.full(self.ne,box)])),dict(velocity_rank=self.nv,velocity_mass_eigenvalues=ev.tolist(),velocity_mass_orthogonality_error=float(np.max(abs(G-np.eye(self.nv)))),velocity_ellipsoid_centre=centre.tolist(),velocity_ellipsoid_radius_squared=float(radius2),epsilon_rank=self.ne,epsilon_independent_rows=ids.tolist(),epsilon_condition=float(np.linalg.cond(Q)),epsilon_inverse_residual=eta,epsilon_box=box,method='consistent physical mass ellipsoid and independent actual scalar envelope rows; no trust boxes')
 def fields(self,x):
  d=np.einsum('njc,j->nc',self.VB,x[:self.nv],optimize=False);return self.g+d[:,:4],self.v-self.vcold+d[:,4:],self.e+self.E@x[self.nv:]
 def scan(self,x,delta,budget):
  budget.scan();g,dv,e=self.fields(x);Sx=2*g[:,0];Sz=2*g[:,3];Sxz=g[:,1]+g[:,2];smax=.5*(Sx+Sz+np.hypot(Sx-Sz,2*Sxz));c=2*self.k/3-delta;F=np.where(self.inner,self.nt/self.K0,.09*self.k*self.k/(self.K0*self.ec));iso=np.where(self.inner,c/self.K0,c*e/(self.K0*self.ec));xx=iso-F*Sx;zz=iso-F*Sz;xz=-F*Sxz;lam2=.5*(xx+zz-np.hypot(xx-zz,2*xz));lam=np.minimum(lam2,iso);angle=.5*np.arctan2(2*xz,xx-zz)+np.pi/2;v=np.c_[np.cos(angle),np.zeros(self.n),np.sin(angle)];v[iso<=lam2]=[0,1,0]
  nu=np.where(self.inner,self.nt,.09*self.k*self.k/e);raw=2*self.k/3-nu*np.maximum(0,smax);y=np.r_[1.,x[:self.nv]];mass=float(y@self.mass@y);div=np.array([y@D@y for D in self.div]);ratios=e/self.ec;cap=float(np.max(abs(dv)));scalar_violation=float(max(0,.5-ratios.min(),ratios.max()-2));hard=max(mass/.12**2-1,float(np.max(div/self.div_limits-1)),cap/.5-1,scalar_violation)
  return dict(gradient=g,velocity_increment=dv,epsilon=e,ratios=ratios,lambda_transformed=lam,eigenvectors=v,lambda_covariance=raw,slack=float(max(0,-np.min(lam[self.cov]))),inner_min=float(np.min(raw[self.cov&self.inner])),outer_min=float(np.min(raw[self.cov&~self.inner])),mass_squared=mass,strong_divergence_squared=div,velocity_max=cap,epsilon_ratio=[float(ratios.min()),float(ratios.max())],hard_violation=hard)
 def epsilon_constraints(self,band):
  outer=~self.inner;C=self.E[outer]/self.ec[outer,None];b=self.e[outer]/self.ec[outer];A=np.vstack((C,-C));rhs=np.r_[band[1]-b,b-band[0]];scale=np.maximum(1,np.maximum(np.max(abs(A),axis=1),abs(rhs)));A=sparse.csr_matrix(A/scale[:,None]);return sparse.hstack((sparse.csr_matrix((len(rhs),self.nv)),A,sparse.csr_matrix((len(rhs),1))),format='csr'),rhs/scale
 def rayleigh_cut(self,i,v,delta):
  v=np.array(v);norm=float(v@v);c=2*self.k[i]/3-delta;F=self.nt[i]/self.K0 if self.inner[i] else .09*self.k[i]**2/(self.K0*self.ec[i]);iso=c/self.K0 if self.inner[i] else c*self.e[i]/(self.K0*self.ec[i]);g=self.g[i];Svv=v[0]**2*2*g[0]+2*v[0]*v[2]*(g[1]+g[2])+v[2]**2*2*g[3];V=self.VB[i];d=-F*(2*v[0]**2*V[:,0]+2*v[0]*v[2]*(V[:,1]+V[:,2])+2*v[2]**2*V[:,3]);de=np.zeros(self.ne) if self.inner[i] else c/(self.K0*self.ec[i])*norm*self.E[i];return np.r_[-d,-de,-norm],float(iso*norm-F*Svv)
 def cap_cut(self,i,component,sign):return np.r_[sign*self.VB[i,:,4+component],np.zeros(self.ne+1)],float(.5-sign*(self.v[i,component]-self.vcold[i,component]))
 def quadratic_cut(self,Q,limit,x):
  a=x[:self.nv];value=float(Q[0,0]+2*Q[0,1:]@a+a@Q[1:,1:]@a);grad=2*(Q[0,1:]+Q[1:,1:]@a);return np.r_[grad,np.zeros(self.ne+1)]/limit,(limit-value+grad@a)/limit

def convex_segment(P,origin,target,start_scan,end_scan):
 """Stay in true cap/ellipsoid intersection; no new objective or projection."""
 alpha=1.;dv=start_scan['velocity_increment'];change=end_scan['velocity_increment']-dv
 for sign in (1.,-1.):
  d=sign*change;mask=d>0
  if np.any(mask):alpha=min(alpha,float(np.min((.5-sign*dv[mask])/d[mask])))
 a=origin[:P.nv];d=(target-origin)[:P.nv]
 for Q,limit in [(P.mass,.12**2)]+list(zip(P.div,P.div_limits)):
  y=np.r_[1.,a];z=np.r_[0.,d];A=float(z@Q@z);B=float(2*y@Q@z);C=float(y@Q@y-limit)
  if A+B+C>0:
   if A>1e-30:
    root=np.sqrt(max(0.,B*B-4*A*C));r=(-2*C/(B+root)) if B>=0 and B+root>0 else (-B+root)/(2*A)
   else:r=-C/B if B>0 else 0.
   alpha=min(alpha,max(0.,r))
 return max(0.,min(1.,alpha))*(1-1e-10)

def summary(s):return {k:s[k] if not isinstance(s[k],np.ndarray) else s[k].tolist() for k in ('slack','inner_min','outer_min','mass_squared','strong_divergence_squared','velocity_max','epsilon_ratio','hard_violation')}

def phase(problem,delta,band,budget,save,initial=None,lp_solver=None):
 P=problem;n=P.dim;x=np.zeros(n) if initial is None else np.array(initial);cuts=[];rhs=[];meta=[];rounds=[];best=None;best_summary=None;best_slack=float('inf');best_scan=None;terminated='inconclusive';
 try:initial_scan=P.scan(np.zeros(n),delta,budget)
 except BudgetExhausted as error:return dict(termination='budget_exhausted:'+str(error),best_coefficients=None,best_summary=None,budget=budget.report(),proof=None)
 tmax=initial_scan['slack']+1.;bounds=P.bounds+[(0,tmax)];c=np.r_[np.zeros(n),1.];hard,hb=P.epsilon_constraints(band);seen=set()
 def add(row,b,m):
  scale=max(1,float(np.max(abs(row))),abs(b));cuts.append(row/scale);rhs.append(b/scale);m['scale']=scale;meta.append(m)
 def separate(s,x,initial=False):
  nonlocal best,best_summary,best_slack,best_scan
  buffered_scalar=max(0,band[0]-s['ratios'].min(),s['ratios'].max()-band[1]);valid=s['hard_violation']<=1e-9 and buffered_scalar<=1e-10
  if valid and s['slack']<best_slack:
   best=x.copy();best_scan=s;best_slack=s['slack'];best_summary=summary(s);save(dict(event='best_hard_diagnostic',coefficients=x.tolist(),summary=best_summary))
  values=s['lambda_transformed'];eligible=np.flatnonzero(P.cov);ids=[]
  for inner in (True,False):
   group=eligible[P.inner[eligible]==inner];ids+=group[np.argsort(values[group])[:min(48,len(group))]].tolist()
  if hasattr(P,'geometry'):
   for geometry in np.unique(P.geometry):
    group=np.flatnonzero(P.cov&(P.geometry==geometry))
    if len(group):ids.append(int(group[np.argmin(values[group])]))
  for i in sorted(set(ids)):
   if initial or values[i]<-1e-10:
    v=s['eigenvectors'][i];key=(i,tuple(np.round(v,12)))
    if key not in seen:add(*P.rayleigh_cut(i,v,delta),dict(kind='PSD',point=int(i),inner=bool(P.inner[i]),vector=v.tolist()));seen.add(key)
  dv=s['velocity_increment'];flat=np.argsort(abs(dv).ravel())[-min(96,dv.size):]
  for q in flat:
   i,j=divmod(int(q),2)
   if abs(dv[i,j])>.5+1e-11 or (initial and abs(dv[i,j])>.44):
    sign=1 if dv[i,j]>=0 else -1;key=('cap',i,j,sign)
    if key not in seen:add(*P.cap_cut(i,j,sign),dict(kind='response_cap',point=i,component=j,sign=sign));seen.add(key)
  quads=[('mass',P.mass,.12**2,s['mass_squared'])]+[('divergence_GL'+str(o),D,limit,v) for o,D,limit,v in zip((16,20,24),P.div,P.div_limits,s['strong_divergence_squared'])]
  for name,Q,limit,value in quads:
   if initial or value>limit*(1+1e-10):add(*P.quadratic_cut(Q,limit,x),dict(kind=name,tangent=x[:P.nv].tolist()))
  return valid
 try:
  start=P.scan(x,delta,budget);valid=separate(start,x,True)
  if valid and start['slack']<=1e-9 and min(start['inner_min'],start['outer_min'])>=delta-1e-12:terminated='seen_feasible'
  while terminated!='seen_feasible':
   budget.check();A=sparse.vstack((hard,sparse.csr_matrix(np.asarray(cuts))),format='csr');b=np.r_[hb,rhs];result,receipt=(lp_solver or lp)(c,A,b,bounds,budget);rounds.append(receipt)
   if result.status!=0:terminated='numerical_or_work_limit';break
   x=result.x[:n];terminal=x.copy();cert=receipt['certificate'];y=np.minimum(0,np.array(result.ineqlin.marginals));box=np.array([max(abs(l),abs(u)) for l,u in bounds]);row_error=2e-9*(1+abs(b)+abs(A)@box);round_allowance=float(abs(y)@row_error);cert['native_data_rounding_allowance']=round_allowance;cert['limited_lower_bound']=cert['lower_bound']-round_allowance
   save(dict(event='LP',round=budget.rounds,receipt=receipt,coefficients=x.tolist()))
   if cert['verified'] and cert['limited_lower_bound']>1e-7:
    terminated='verified_relaxation_obstruction';dual=np.flatnonzero(abs(y)>1e-10);proof=dict(receipt=receipt,rows=A[dual].toarray().tolist(),rhs=b[dual].tolist(),multipliers=y[dual].tolist(),row_indices=dual.tolist(),hard_epsilon_rows=len(hb),cut_metadata=[meta[j-len(hb)] if j>=len(hb) else dict(kind='epsilon_bound',hard_row=int(j)) for j in dual]);break
   s=P.scan(x,delta,budget);valid=separate(s,x)
   if not valid and best is not None:
    alpha=convex_segment(P,best,x,best_scan,s)
    if alpha>1e-14:
     z=best+alpha*(x-best);zs=P.scan(z,delta,budget);zvalid=separate(zs,z,True);save(dict(event='convex_physical_segment',alpha=alpha,summary=summary(zs)))
     if zvalid and zs['slack']<=1e-9 and min(zs['inner_min'],zs['outer_min'])>=delta-1e-12:terminated='seen_feasible';break
   if valid and s['slack']<=1e-9 and min(s['inner_min'],s['outer_min'])>=delta-1e-12:terminated='seen_feasible';break
 except BudgetExhausted as e:terminated='budget_exhausted:'+str(e)
 finally:
  # No report-only scan is allowed to exceed the declared budget. Best and
  # terminal fields were evaluated/saved at the point they became available.
  out=dict(termination=terminated,delta=delta,epsilon_band=band,rounds=rounds,best_coefficients=best.tolist() if best is not None else None,best_summary=best_summary,terminal_coefficients=locals().get('terminal',x).tolist(),budget=budget.report(),cut_count=len(cuts),proof=locals().get('proof'))
  save(dict(event='phase_terminal',result=out))
 return out
