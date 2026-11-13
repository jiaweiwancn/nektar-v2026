#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused native model identity, regression, restart and V2 rejection checks."""
import argparse,os,subprocess,time,copy
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_budget_h2g import case,canonical_bytes
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import THREADS,sha,state_hash,set_parameter,set_info,write,require
from compare_cls_epawf_pimple import same_history
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2g_constitutive import tensors_for_spec

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('solver','probe','h1','preflight','h2f','cold','flat','output'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();os.environ.update(THREADS);pre=load(a.preflight);a.output.mkdir(exist_ok=False);out=dict(version='H2gFocusedVerificationV1',tests={},native={},status='running');start=time.monotonic()
 def save():atomic(a.output/'verification.json',out)
 def native(name,model='standard',steps=0,cp=None,cert=None,mutate=None,expected=None,probe=False):
  xml=case(a.output/name,a.h1,pre,flat=True,model=model);tree=ET.parse(xml).getroot();set_parameter(tree,'NumSteps',steps);set_parameter(tree,'TimeStep',2e-6)
  if cp:set_info(tree,'RANSWallRestart',str(cp))
  if cert:set_info(tree,'RANSWallCertifiedStartupFile',str(cert));set_info(tree,'RANSWallCertifiedStartupSHA256',sha(cert));set_info(tree,'RANSWallCertifiedStartupVersion','CertifiedNeutralStartupV2')
  if mutate:mutate(tree,xml)
  write(tree,xml);env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str((a.probe if probe else a.solver).parent/'lib'));cmd=[str(a.solver),str(xml)]
  if probe:
   packet=xml.parent/'input.json';atomic(packet,dict(operation='h2c_seed',inner_stress_model='Linear',xml_sha256=sha(xml),expected_model='CLSHighRe'));cmd=[str(a.probe),str(xml),str(packet),str(xml.parent/'result.json')]
  with (xml.parent/'solver.log').open('w') as f:r=subprocess.run(cmd,cwd=xml.parent,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=300)
  text=(xml.parent/'solver.log').read_text();path=xml.with_suffix('.evolution.json');d=load(path) if path.exists() else None
  if expected:require(expected in text or (d and expected in str(d.get('failure'))),'unexpected rejection reason '+name)
  else:require(r.returncode==0 and d and not d['failure'],'regression failed '+name)
  out['native'][name]=dict(code=r.returncode,report=str(path) if d else None,failure=d.get('failure') if d else text[-1800:],expected=expected);save();return d
 try:
  native('reject-sigma',mutate=lambda t,x:set_parameter(t,'RANSSigmaEpsilon',1.3),expected='inconsistent effective calibrated coefficient')
  native('reject-Qstar',mutate=lambda t,x:set_parameter(t,'Qstar',pre['models'][0]['q']),expected='native coefficient/reference expression mismatch')
  native('reject-probe-model',probe=True,expected='probe/model-manifest mismatch')
  def manifest_mutation(t,x):
   path=x.parent/'model-manifest.json';m=load(path);m['coefficients']['Cmu']=.1;path.write_bytes(canonical_bytes(m)+b'\n');set_info(t,'RANSWallModelManifestSHA256',sha(path))
  native('reject-manifest',mutate=manifest_mutation,expected='native/model manifest coefficient mismatch')
  old=NativeContext(a.h2f/'enriched-context/context.json');corrupt=load(a.h2f/'enriched-context/context.json');corrupt['expected_model']='StandardKEpsilon';path=a.output/'reject-context.json';atomic(path,corrupt)
  try:NativeContext(path);raise RuntimeError('CLS context relabelled standard')
  except ValueError:out['tests']['CLS_context_as_standard_rejected']=True
  ids=np.flatnonzero((old.block('inner')==0)&(old.block('covariance')>0))[:4];spec=dict(name='StandardKEpsilon',calibration='DeployedCLSV1',physical_regularization_active=True)
  try:tensors_for_spec(old.gradient()[ids],old.block('k')[ids],old.block('epsilon')[ids],np.zeros(len(ids),bool),np.zeros(len(ids)),old.block('physical')[ids],pre['u_tau'],pre['kappa'],spec);raise RuntimeError('CLS evaluator relabelled standard')
  except ValueError:out['tests']['CLS_evaluator_as_standard_rejected']=True
  # Fresh unchanged CLS h0/two steps only; never rerun its search campaign.
  cls=native('CLS-h0',model='cls');prior=a.h2f.parent/'cls-cao-tamura-budget-h2e-20261113/verification-final/ordinary-h0/case.evolution.json';saved=load(prior);require(same_history(cls,saved),'touched CLS h0/history differs from retained regression');c2=native('CLS-two',model='cls',steps=2,cp=Path(cls['checkpoint']));old2=load(prior.parent.parent/'ordinary2/case.evolution.json');require(same_history(c2,old2),'touched CLS two-step/history differs');out['tests']['CLS_h0_two_steps_histories_bitwise']=True
  native('cross-model-restart',cp=Path(cls['checkpoint']),expected='checkpoint')
  cert=a.flat/'certificate.json';native('cross-model-certificate',model='cls',cert=cert,expected='certified startup physics/geometry/inventory/gate differs')
  standard=load(a.flat/'evolution/h0/case.evolution.json');native('ordinary-changed-dt',cp=Path(standard['checkpoint']),cert=cert,mutate=lambda t,x:set_parameter(t,'TimeStep',1e-6),expected='checkpoint')
  hillcls=load(a.cold/'cls-hill-probe/result.json');hillstd=load(a.cold/'standard-hill-probe/result.json');require(hillcls['state_sha256']=='83b1a4fd9a4fcda694f6e6f10bfb789700d1fc6fc83443b91f67bfcb045e79b5','cold CLS base changed');n=760;delta=float(np.max(abs(np.array(hillcls['state'])[:2*n]-np.array(hillstd['state'])[:2*n])));require(delta<1e-10,'model cold velocity differs materially');out['tests'].update(native_solver_probe_manifests_and_expression_rejections=True,cross_model_and_dt_restarts_rejected=True,CLS_cold_base_bitwise=True,cold_velocity_max_difference=delta)
  out.update(status='passed',seconds=time.monotonic()-start)
 except Exception as e:out.update(status='failed',failure=str(e));raise
 finally:save()
if __name__=='__main__':main()
