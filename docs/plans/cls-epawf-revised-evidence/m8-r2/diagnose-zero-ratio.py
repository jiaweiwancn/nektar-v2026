import argparse, hashlib, json, os, pathlib, sys
p=argparse.ArgumentParser();p.add_argument('--source',type=pathlib.Path,required=True);p.add_argument('--raw',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
sys.dont_write_bytecode=True;sys.path.insert(0,str(a.source/'solvers/IncNavierStokesSolver/Tests/RANS'))
import numpy as np
import mpmath as mp
mp.mp.dps=80
from cls_epawf_m7_energy_oracle import _fields,_velocity,_sample,Audit
from cls_epawf_m8_native_oracle import ray_source
from cls_epawf_m8_r2_preflight import selected_direct
from cls_epawf_m8_r2_oracle import indexed_inner
r=json.loads(a.raw.read_text());f=_fields(r);s=np.asarray(r['state']);off=r['offsets'];d=r['dimension'];v=r['constitutive_audit']['inner_stress_max_ratio_point'];point=np.array(v[:3]);report=dict(input=str(a.raw),input_sha256=hashlib.sha256(a.raw.read_bytes()).hexdigest(),native_point=v,native_stats=r['constitutive_audit']['inner_stress'],native_zero=r['constitutive_audit']['inner_zero_covariance'])
ind=selected_direct(r,f,s,Audit());report['independent_stats']=ind['inner_stress'];report['independent_point']=ind['inner_stress_max_ratio_point']
e=f[d].element(point);ids,phi,deriv=f[d].rows(e,[point]);coeff=s[off[d]+ids]
report['plain_k']=float(phi[0]@coeff);report['absolute_dot_sum']=float(np.sum(np.abs(phi[0]*coeff)));report['coefficient_l1']=float(np.sum(np.abs(coeff)))
def b(x):return [(1-x)/2,(1+x)/2,(1-x*x)/4,x*(1-x*x)/2]
def db(x):return [mp.mpf('-.5'),mp.mpf('.5'),-x/2,(1-3*x*x)/2]
def exact(block,point):
 field=f[block];el=field.element(point);axes=field.axes;n=el['degree']+1
 assert n==4
 bs=[];gs=[]
 for axis in axes:
  length=mp.mpf(el['upper'][axis])-mp.mpf(el['lower'][axis]);x=2*(mp.mpf(float(point[axis]))-mp.mpf(el['lower'][axis]))/length-1
  bs.append(b(x));gs.append([z*2/length for z in db(x)])
 value=mp.mpf(0);grad=[mp.mpf(0)]*3
 for j in range(n*n):
  x,z=j%n,j//n;c=mp.mpf(float(s[off[block]+el['ids'][j]]))*el['sign'][j]
  value+=c*bs[0][x]*bs[1][z];grad[0]+=c*gs[0][x]*bs[1][z];grad[2]+=c*bs[0][x]*gs[1][z]
 return value,grad
kmp,_=exact(d,point);Gmp=[[mp.mpf(0) for _ in range(3)] for _ in range(3)]
for block,axis in enumerate([0,2]):Gmp[axis]=exact(block,point)[1]
smp=mp.sqrt(sum((Gmp[i][j]+Gmp[j][i])**2 for i in range(3) for j in range(3)))
ray=ray_source(r,f,s,point);nuA=float(ray['profile'].a)*max(point[2]-float(ray['profile'].zv),0)
report.update(exact_polynomial_k_80digit=str(kmp),exact_polynomial_strain_80digit=str(smp),native_k_minus_exact=str(mp.mpf(v[3])-kmp),native_k_relative_error=float(abs(mp.mpf(v[3])-kmp)/abs(kmp)),native_basis_linf_error_sufficient=float(abs(mp.mpf(v[3])-kmp)/mp.mpf(report['coefficient_l1'])),plain_dot_condition=report['absolute_dot_sum']/abs(float(kmp)),independent_nuA=nuA)
G=np.array([[float(x) for x in row] for row in Gmp]);native_input=indexed_inner(G,v[3],v[5]);exact_input=indexed_inner(G,float(kmp),nuA)
report['native_point_recomputed']={k:float(np.real(native_input[k])) for k in ['s','a','c','t','radial']};report['native_point_recomputed']['ratio']=float(np.real(native_input['a']/native_input['c']))
report['exact_polynomial_ratio']=float(np.real(exact_input['a']/exact_input['c']))
report['scope']='Read-only finite exported-point conditioning diagnosis. 80-digit polynomial evaluation is a reference calculation, not a rigorous interval certificate. No native or solver calls.'
a.output.mkdir(parents=True,exist_ok=False);(a.output/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report,indent=2))
