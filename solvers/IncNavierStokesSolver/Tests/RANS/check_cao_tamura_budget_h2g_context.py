#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Common discretization, native model, smoothing and divergence comparison."""
import argparse,os
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import numpy as np
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,state_hash
from cao_tamura_h2g_constitutive import tensors_for_spec

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('context','old-context','preflight','output'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();c=NativeContext(a.context);o=NativeContext(a.old_context);pre=load(a.preflight);require(c.model['name']=='StandardKEpsilon' and c.model['calibration']=='StandardNeutralCalibratedV1','actual standard native model required');require(c.counts==o.counts==dict(velocity=41,k=8,epsilon=8),'common qualified pool changed')
 base=np.array(c.array('base')[0]);old=np.array(o.array('base')[0]);delta=float(np.max(abs(base[:1520]-old[:1520])));require(delta<1e-10,'cold native velocities differ materially')
 basis=np.array(c.array('basis'));v_equal=basis[:41].tobytes()==o.array('basis')[:41].tobytes();require(v_equal,'velocity pool coefficient bits changed')
 mass=float(np.max(abs(c.array('mass_gram')-o.array('mass_gram'))));require(mass<1e-12,'common physical mass/Gram differs');D=np.array(c.array('divergence'));Do=np.array(o.array('divergence'));dd=float(np.max(abs(D-Do)));require(dd<1e-8,'common physical divergence form differs')
 common=float(np.sqrt(Do[0,0])*(1+1e-6)+1e-9);native=float(np.sqrt(D[0,0])*(1+1e-6)+1e-9);require(abs(common-native)<1e-10,'common divergence bound inconsistent');strain=float(np.max(abs(c.array('strain')-o.array('strain'))));require(strain<1e-8,'standard Gamma=0 disabled or changed initializer smoothing')
 inner=c.block('inner')>0;dep=c.block('dependent_epsilon')>0;require(np.max(abs(c.block('k_basis')[inner]))==0 and np.max(abs(c.block('epsilon_basis')[inner|dep]))==0,'scalar basis changed inner/trace values');ids=np.flatnonzero((~inner)&(c.block('covariance')>0))[::997];g=c.gradient()[ids];k=c.block('k')[ids];e=c.block('epsilon')[ids];vec=np.array([np.linalg.eigh(R)[1][:,0] for R in tensors_for_spec(g,k,e,np.zeros(len(ids),bool),np.zeros(len(ids)),c.block('physical')[ids],pre['u_tau'],pre['kappa'],c.model)]);errors=[]
 for col,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):
  h=1e-5;gp=g.copy();gm=g.copy();gp[:,i,j]+=h;gm[:,i,j]-=h;R1=tensors_for_spec(gp,k,e,np.zeros(len(ids),bool),np.zeros(len(ids)),c.block('physical')[ids],pre['u_tau'],pre['kappa'],c.model);R0=tensors_for_spec(gm,k,e,np.zeros(len(ids),bool),np.zeros(len(ids)),c.block('physical')[ids],pre['u_tau'],pre['kappa'],c.model);fd=np.einsum('ni,nij,nj->n',vec,R1-R0,vec)/(2*h);exact=-2*.09*k*k/e*vec[:,i]*vec[:,j];errors.append(float(np.max(abs(fd-exact))))
 require(max(errors)<1e-9,'analytic standard Rayleigh derivative mismatch')
 relabelled=load(a.old_context);relabelled['context']['model_specification']=c.model;relabelled['expected_model']='StandardKEpsilon';bad=a.output.with_suffix('.relabelled-CLS.json');atomic(bad,relabelled)
 try:NativeContext(bad);raise RuntimeError('relabelled CLS context accepted')
 except ValueError:pass
 atomic(a.output,dict(native_export_binding_verified=True,relabelled_CLS_context_rejected=True,passed=True,model_specification=c.model,base_sha256=state_hash(base),retained_CLS_base_sha256=state_hash(old),cold_velocity_max_difference=delta,velocity_pool_bitwise=v_equal,velocity_directions=41,scalar_directions=[8,8],mass_Gram_max_difference=mass,divergence_form_max_difference=dd,initializer_strain_max_difference=strain,common_divergence_limit=common,standard_native_limit=native,scalar_inner_and_dependent_traces_exact=True,analytic_standard_Rayleigh_errors=errors,classification='common discretization and initializer controls; not matched optimization histories'))
if __name__=='__main__':main()
