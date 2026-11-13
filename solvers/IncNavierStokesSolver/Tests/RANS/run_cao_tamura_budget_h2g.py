#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Isolated H2g native calibration, cold bases, certification and bounded pilots."""
import argparse,os,time,subprocess,resource
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_budget_h2g import case,prepare
from continue_cls_epawf_pimple import atomic,freeze,identity
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import sha,state_hash,require,THREADS
from run_cao_tamura_budget_h2e import certificate_v2
from run_cao_tamura_budget_h2f_startup import StartupRun
from compare_cls_epawf_pimple import same_history

class Study:
 def __init__(self,a):
  self.a=a;self.pre=load(a.preflight);self.solver=a.solver;self.probe=a.probe;self.root=a.output;self.root.mkdir(exist_ok=False);self.out=dict(version='H2gNativeWorkV1',phase=a.phase,results={},physical_hill_attempts=0,physical_hill_steps=0);self.save()
 def save(self):atomic(self.root/'progress.json',self.out)
 def read(self,name,flat,model,operation='h2c_seed',extra=None):
  xml=case(self.root/name,self.a.h1,self.pre,flat=flat,model=model);p=dict(operation=operation,inner_stress_model='Linear',xml_sha256=sha(xml),expected_model='StandardKEpsilon' if model=='standard' else 'CLSHighRe')
  if extra:p.update(extra)
  start=time.monotonic();d,r=run_probe(self.probe,xml,p,self.root/(name+'-probe'));d=load(Path(r['output']));self.out['results'][name]=dict(receipt=r,seconds=time.monotonic()-start,model=load(Path(r['output']+'.model.json')),state_sha256=d.get('state_sha256'),audit=d.get('audit'),output=r['output'],xml=str(xml),peak_RSS_KiB_cumulative=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);self.save();return d,r,xml
 def cold(self):
  if self.a.reuse_cold:
   previous=load(self.a.reuse_cold/'progress.json')
   for name in ('cls-flat','cls-hill'):
    self.out['results'][name]=previous['results'][name];(self.root/(name+'-probe')).symlink_to(self.a.reuse_cold/(name+'-probe'))
  for model in (('standard',) if self.a.reuse_cold else ('cls','standard')):
   for flat in (True,False):self.read(model+('-flat' if flat else '-hill'),flat,model)
  self.out['status']='completed_cold_audits';self.save()
 def frozen(self):
  for label,path in [('base',None),('A',self.a.h2f/'stage-A-final/best-candidate.json'),('B',self.a.h2f/'stage-B/best-candidate.json')]:
   folder=self.root/(label+'-arrays');folder.mkdir();extra=dict(context_directory=str(folder))
   if path:
    c=load(path);extra.update(state=c['state'],state_sha256=c['state_sha256'])
   d,r,x=self.read(label,False,'cls','h2g_frozen',extra);self.out['results'][label].update(arrays=str(folder),state_sha256=d['state_sha256']);self.save()
  self.out['status']='completed_zero_governing_frozen_fields';self.save()
 def context(self):
  cold=load(self.a.cold/'standard-hill-probe/result.json');basis=self.a.h2f/'enrichment/velocity.f64';folder=self.root/'arrays';folder.mkdir()
  d,r,xml=self.read('standard-context',False,'standard','h2f_context',dict(base_sha256=cold['state_sha256'],basis_file=str(basis),basis_sha256=sha(basis),velocity_direction_count=41,context_directory=str(folder)))
  d.update(array_directory=str(folder),xml=str(xml),xml_sha256=sha(xml),runtime=identity(self.probe),preflight=str(self.a.preflight),expected_model='StandardKEpsilon',receipt=r)
  atomic(self.root/'context.json',d);base=d['context']['base_sha256'];c=dict(version='H2gColdCandidateV1',stage='A',state=cold['state'],state_sha256=base,base_sha256=base,conditioned_velocity_coefficients=[0.]*41,conditioned_scalar_coefficients=[]);atomic(self.root/'cold-candidate.json',c);self.out['status']='completed_standard_context';self.save()
 def certification(self,candidate,flat=False):
  c=load(candidate);d,r,x=self.read('candidate',flat,'standard','h2e_candidate',dict(state=c['state'],state_sha256=c['state_sha256'],base_sha256=c['base_sha256'],coupled=True));K0=self.pre['models'][1]['K0'];au=d['audit'];ex=d['extra_sampling'];limit=c.get('divergence_limit',1e-8)
  gates=dict(native_admissible=au['accepted_awf_hard_admissible']==[1] and au['admissible']==[1],budget=au['accepted_awf_scalar_domain'][8]<=1e-10,algebraic=max(au['coupled_norms'][4:])<1e-9,extra_covariance=ex['minimum_covariance']>=1e-4*K0-1e-12,strong_divergence=ex['strong_divergence_L2']<=limit+1e-8,velocity=ex['velocity_change_L2']<=.12*(1+1e-10) and ex['maximum_sampled_velocity_increment']<=.5+1e-10,optimizer_interior_margin=min(au['accepted_awf_covariance'])>=1e-4*K0-1e-12)
  receipt=self.root/'native-receipt.json';atomic(receipt,dict(passed=all(gates.values()),gates=gates,K0=K0,state_sha256=d['state_sha256'],base_sha256=c['base_sha256'],probe_output_file=r['output'],probe_output_sha256=sha(Path(r['output'])),audit=au,extra_sampling=ex))
  if not all(gates.values()):self.out.update(status='candidate_rejected',gates=gates);self.save();return None
  base=d['fixed_reference'];source=Path(__file__).resolve().parents[4];names=['solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallModelSpec.h','library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp'];provenance=dict(K0=K0,divergence_limit=limit,basis_sha256=sha(candidate),context_sha256=sha(candidate),native_receipt_file=str(receipt),source_revision='H2g exact source manifest',source_hashes={n:sha(source/n) for n in names});cert=certificate_v2(d,load(receipt),base,identity(self.solver),identity(self.probe),provenance);cert['model_specification']=load(Path(r['output']+'.model.json'));path=self.root/'certificate.json';atomic(path,cert);self.out.update(status='candidate_certified',gates=gates,certificate=str(path));self.save();return path
 def flat(self):
  cold=load(self.a.cold/'standard-flat-probe/result.json');base=np.array(cold['state']);xml=case(self.root/'contract-case',self.a.h1,self.pre,flat=True)
  # Obtain the native permitted-row contract on the actual standard base.
  temp=self.root/'base-candidate.json';atomic(temp,dict(state=base.tolist(),state_sha256=state_hash(base),base_sha256=state_hash(base),divergence_limit=1e-8))
  if self.a.reuse_flat:
   source=self.a.reuse_flat/'base-contract';d=load(source/'result.json');r=load(source/'receipt.json');require(d['state_sha256']==state_hash(base),'retained flat contract base changed');self.out['base_contract_reused']=str(source);self.save()
  else:
   d,r=run_probe(self.probe,xml,dict(operation='h2e_candidate',inner_stress_model='Linear',expected_model='StandardKEpsilon',xml_sha256=sha(xml),state=base.tolist(),state_sha256=state_hash(base),base_sha256=state_hash(base),coupled=True),self.root/'base-contract');d=load(Path(r['output']))
  rows=d['outer_initializer_contract'];state=base.copy();kr=rows['outer_k_rows'][len(rows['outer_k_rows'])//2];er=rows['outer_epsilon_rows'][len(rows['outer_epsilon_rows'])//2];state[kr]+=1e-4*self.pre['models'][1]['K0'];state[er]+=1e-4*self.pre['models'][1]['E_C_H'];path=self.root/'changed-candidate.json';atomic(path,dict(stage='B',state=state.tolist(),state_sha256=state_hash(state),base_sha256=state_hash(base),divergence_limit=1e-8,permitted_nonzero_rows=[kr,er]));cert=self.certification(path,True);require(cert,'standard changed affine fixture failed')
  builder=lambda folder,h1,pre,flat=False:case(folder,h1,pre,flat=flat)
  run=StartupRun(self.solver,cert,self.a.preflight,self.a.h1,self.root/'evolution',flat=True,builder=builder)
  try:
   h=run.native('h0',2e-6);run.accepted(h);c=load(cert);require(state_hash(h['state'])==c['state_sha256'] and c['state_sha256']!=c['immutable_base_sha256'],'V2 candidate distinction/import');require(h['initialization_lift_calls_iterations']==[0,0],'import Lift/projection');cp=Path(h['checkpoint']);parent=sha(cp/'State.bin');q=run.native('direct2',2e-6,2,cp);run.accepted(q);f=run.native('first1',2e-6,1,cp);run.accepted(f);s=run.native('split1',2e-6,1,Path(f['checkpoint']));run.accepted(s);require(same_history(q,s),'standard direct/split differ');run.rollback(h,2e-6);require(sha(cp/'State.bin')==parent,'parent changed')
   for v in (h,q,f,s):require(state_hash(v['fixed_reference'])==c['fixed_reference_sha256'] and state_hash(v['fixed_defect'])==c['fixed_defect_sha256'],'reference/defect changed')
   run.out.update(status='passed_standard_nonzero_V2',tests=dict(candidate_distinct=True,exact_import=True,no_Lift_projection=True,direct_split_bitwise=True,rollback_parent_preserved=True,one_sweep=True,solved_epsilon_retained=True))
  finally:run.save();self.out['evolution']=run.out;self.save()
 def unperturbed(self):
  from run_cls_epawf_picard import set_parameter,set_info,write
  import xml.etree.ElementTree as ET
  parent=None
  for label,steps in (('h0',0),('two',2)):
   xml=case(self.root/label,self.a.h1,self.pre,flat=True);tree=ET.parse(xml).getroot();set_parameter(tree,'NumSteps',steps);set_parameter(tree,'TimeStep',2e-6)
   if parent:set_info(tree,'RANSWallRestart',str(parent))
   write(tree,xml);env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(self.solver.parent/'lib'));started=time.monotonic()
   with (xml.parent/'solver.log').open('w') as log:run=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(self.solver),str(xml)],cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300)
   d=load(xml.with_suffix('.evolution.json'));require(run.returncode==0 and not d['failure'],'unperturbed standard flat rejected')
   for v in d['pimple_attempts']:
    for key in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(v[key]==1,'wrong standard schedule')
    require(v['solved_epsilon_retained_bitwise'],'epsilon overwritten')
   self.out['results'][label]=dict(report=str(xml.with_suffix('.evolution.json')),report_sha256=sha(xml.with_suffix('.evolution.json')),seconds=time.monotonic()-started,usage=(xml.parent/'usage.txt').read_text(),accepted_steps=d['steps']);self.save();parent=Path(d['checkpoint'])
  self.out['status']='passed_unperturbed_standard_two_steps';self.save()
 def pilot(self):
  cert=self.certification(self.a.candidate)
  if cert:
   run=StartupRun(self.solver,cert,self.a.preflight,self.a.h1,self.root/'evolution',builder=lambda folder,h1,pre,flat=False:case(folder,h1,pre,flat=flat));run.pair();self.out['evolution']=run.out;self.save()

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['prepare','cold','frozen','context','flat','certify','pilot','unperturbed'])
 for n in ('preflight','solver','probe','h1','h2f','cold','candidate','calibration','cls-preflight','mesh-definition','reuse-cold','reuse-flat'):p.add_argument('--'+n,type=Path)
 p.add_argument('--output',type=Path,required=True);a=p.parse_args();os.environ.update(THREADS)
 if a.phase=='prepare':prepare(a.calibration,a.cls_preflight,a.mesh_definition,a.output);return
 s=Study(a)
 if a.phase=='certify':s.certification(a.candidate)
 else:getattr(s,a.phase)()
if __name__=='__main__':main()
