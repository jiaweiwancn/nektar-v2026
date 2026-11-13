#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused conditioning, coupled-row and V2 startup/rollback verification."""
import argparse,json,os,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from cao_tamura_h2d_state import load
from run_cao_tamura_budget_h2d import case
from continue_cls_epawf_pimple import freeze,atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import sha,state_hash,set_parameter,set_info,write,require
from compare_cls_epawf_pimple import same_history
from search_cao_tamura_budget_h2e import condition
from search_cao_tamura_budget_h2d import tensors

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','probe','context','h1','preflight','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(exist_ok=False);solver,image=freeze(a.solver,a.output/'runtime');probe,pimage=freeze(a.probe,a.output/'probe-runtime');pre=load(a.preflight);ctx=load(a.context);f=ctx['context'];root=Path(ctx['array_directory']);out=dict(version='H2eFocusedVerificationV1',status='running',runtime=image,probe=pimage,tests={},native={})
    def save():atomic(a.output/'verification.json',out)
    def arr(n):d=f['files'][n];require(sha(root/d['file'])==d['sha256'],'binary context corruption');return np.memmap(root/d['file'],dtype=d['dtype'],mode='r',shape=tuple(d['shape']))
    Z=np.array(arr('basis'));G=np.array(arr('mass_gram'));T,conditioning=condition(G);old=np.cos(np.arange(32)+.1)*1e-4;y=np.linalg.solve(T,old);physical_error=float(np.max(abs(old@Z[:32]-y@(T.T@Z[:32]))));obs=arr('observations');raw=obs[:,18:210].reshape(len(obs),32,6);idx=np.arange(0,len(obs),max(1,len(obs)//200));r=np.array(raw[idx]);original=np.einsum('nbc,b->nc',r,old);conditioned=np.einsum('ncj,j->nc',np.einsum('nbc,bj->ncj',r,T),y);gradient_error=float(np.max(abs(original-conditioned)));distance_error=abs(float(old@G@old-y@y));require(physical_error<1e-11 and gradient_error<1e-9 and distance_error<1e-12,'conditioning changed physical span')
    # Stress arithmetic and fresh scalar consumption are checked on actual
    # native observation gradients, including both strain/rotation branches.
    g=np.zeros((len(idx),3,3));g[:,0,0]=obs[idx,6];g[:,0,2]=obs[idx,7];g[:,2,0]=obs[idx,8];g[:,2,2]=obs[idx,9];g2=g.copy()
    for a1,b1,c in ((0,0,0),(0,2,1),(2,0,2),(2,2,3)):g[:,a1,b1]+=original[:,c];g2[:,a1,b1]+=conditioned[:,c]
    law=lambda g,k,e:tensors(g,k,e,obs[idx,15].astype(bool),obs[idx,14],obs[idx,3:6],pre['u_tau'],pre['kappa_CLS'])
    stress_error=float(np.max(abs(law(g,obs[idx,12],obs[idx,13])-law(g2,obs[idx,12],obs[idx,13]))));require(stress_error<1e-10,'conditioning stress differs')
    outer=~obs[idx,15].astype(bool);fresh_change=float(np.max(abs(law(g,obs[idx,12]*1.1,obs[idx,13]*1.05)[outer]-law(g,obs[idx,12],obs[idx,13])[outer])));require(fresh_change>1e-8,'trial k/epsilon not consumed')
    out['tests']['conditioning']=dict(**conditioning,physical_coefficient_error=physical_error,gradient_error=gradient_error,mass_distance_error=distance_error,stress_error=stress_error,fresh_scalar_stress_change=fresh_change);save()
    # One native process verifies positive scalar support and all prohibited
    # updates. This is a contract test, not covariance acceptance of the hill.
    base=np.array(arr('base')[0]);candidate=base.copy();row=Z[32]+Z[40];amplitude=1e-3;allowed=np.array(f['rows']['k']+f['rows']['epsilon'],int);candidate[allowed]+=amplitude*row[allowed]
    cxml=Path(ctx['xml']);d,receipt=run_probe(probe,cxml,dict(operation='h2e_contract_checks',inner_stress_model='Linear',xml_sha256=sha(cxml),state=candidate.tolist()),a.output/'contract-fixtures');out['tests']['native_contract']=d;out['contract_receipt']=receipt;save()
    def native(name,steps=0,cp=None,cert=None,version=None,dt=2e-6,mutate=None):
        xml=case(a.output/name,a.h1,pre,flat=True);tree=ET.parse(xml).getroot();set_parameter(tree,'NumSteps',steps);set_parameter(tree,'TimeStep',dt)
        if cp:set_info(tree,'RANSWallRestart',str(cp))
        if cert:set_info(tree,'RANSWallCertifiedStartupFile',str(cert));set_info(tree,'RANSWallCertifiedStartupSHA256',sha(cert))
        if version:set_info(tree,'RANSWallCertifiedStartupVersion',version)
        if mutate:mutate(tree)
        write(tree,xml);env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),NEKTAR_RANS_TIMING_FILE=str(xml.parent/'timing.jsonl'));started=time.monotonic()
        with (xml.parent/'solver.log').open('w') as log:r=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(solver),str(xml)],cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        report=xml.with_suffix('.evolution.json');d=load(report) if report.exists() else {};out['native'][name]=dict(code=r.returncode,failure=d.get('failure'),report=str(report),seconds=time.monotonic()-started,usage=(xml.parent/'usage.txt').read_text());save();return d,xml
    parent,xml=native('ordinary-h0');require(not parent['failure'],'flat fixture h0 failed')
    packet=dict(operation='h2e_candidate',inner_stress_model='Linear',xml_sha256=sha(xml),state=parent['state'],state_sha256=state_hash(parent['state']),base_sha256=state_hash(parent['state']),coupled=True)
    d,receipt=run_probe(probe,xml,packet,a.output/'flat-certification');d=load(Path(receipt['output']));require(d['audit']['accepted_awf_hard_admissible']==[1] and d['extra_sampling']['minimum_covariance']>0,'positive flat certificate fixture failed')
    native_receipt=a.output/'flat-native-receipt.json';atomic(native_receipt,dict(passed=True,state_sha256=d['state_sha256'],base_sha256=state_hash(parent['state']),probe=pimage,probe_output_file=receipt['output'],probe_output_sha256=sha(Path(receipt['output'])),gates=dict(extra_covariance=True,strong_divergence=True)))
    source=Path(__file__).resolve().parents[4];paths=['solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionComposite.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp','library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp'];hashes={p:sha(source/p) for p in paths}
    from run_cao_tamura_budget_h2e import certificate_v2
    fixture_receipt=load(native_receipt);fixture_receipt['passed']=d['audit']['accepted_awf_hard_admissible']==[1] and d['audit']['accepted_awf_scalar_domain'][8]<=1e-10 and max(d['audit']['coupled_norms'][4:])<1e-9 and d['extra_sampling']['strong_divergence_L2']<1e-8
    atomic(native_receipt,fixture_receipt)
    cert=certificate_v2(d,fixture_receipt,parent['state'],image,pimage,dict(K0=pre['heights'][-1]['K0'],divergence_limit=1e-8,basis_sha256=f['files']['basis']['sha256'],context_sha256=sha(a.context),native_receipt_file=str(native_receipt),source_revision='29b6da542cd4a38c9facb1417b59e7590376e5ce plus H2e recorded hashes',source_hashes=hashes))
    # The real rejected hill receipts must never yield an import certificate.
    root=a.context.parent.parent
    for stage in ('A','B'):
        receipt_doc=load(root/('certify-'+stage)/'certification.json');report=load(Path(receipt_doc['receipt']['output']))
        try:certificate_v2(report,receipt_doc,base.tolist(),image,pimage,{});raise RuntimeError('rejected hill certified')
        except ValueError as e:require('cannot certify rejected native/extra initializer' in str(e),'unexpected rejected certificate failure')
    out['tests']['rejected_hill_certificates_refused']=True
    certfile=a.output/'startup-certificate.json';atomic(certfile,cert)
    seed,_=native('v2-import-h0',cert=certfile,version=cert['version']);require(not seed['failure'] and same_history(seed,parent),'V2 exact seed/history import differs');require(seed['initialization_lift_calls_iterations']==[0,0],'V2 performed forbidden Lift')
    reference,_=native('ordinary2',2,Path(parent['checkpoint']));direct,_=native('v2-direct2',2,Path(seed['checkpoint']),certfile,cert['version']);require(not direct['failure'] and same_history(reference,direct),'V2 physical path/history differs')
    for key in ('state','fixed_reference','fixed_defect'):require(state_hash(reference[key])==state_hash(direct[key]),'V2 changed '+key)
    count=0
    for x,y in zip(reference['pimple_attempts'],direct['pimple_attempts']):
        for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(x[k]==y[k]==1,'one-sweep schedule changed')
        require(x['post_flow_sha256']==y['post_flow_sha256'] and x['pre_k_sha256']==y['pre_k_sha256'],'governing stages differ');require(x['solved_epsilon_retained_bitwise'] and y['solved_epsilon_retained_bitwise'],'solved epsilon replaced')
        if x.get('k_snapshot') and Path(x['k_snapshot']).exists():
            px,py=load(Path(x['k_snapshot'])),load(Path(y['k_snapshot']))
            for key in ('committed','post_flow','pre_k','raw','corrected','history','fixed_reference'):require(state_hash(px[key])==state_hash(py[key]),'snapshot field differs '+key)
            count+=1
    require(count,'no retained verification packet')
    first,_=native('v2-first',1,Path(seed['checkpoint']),certfile,cert['version']);split,_=native('v2-split',1,Path(first['checkpoint']),certfile,cert['version']);require(same_history(direct,split),'V2 split restart differs')
    cp_hash=sha(Path(seed['checkpoint'])/'State.bin');bad,_=native('postassign-reject',1,Path(seed['checkpoint']),certfile,cert['version'],mutate=lambda t:set_parameter(t,'RANSWallPimplePostAssignRejectStep',1));require(bad['failure'] and same_history(seed,bad),'V2 rollback history/state changed');require(sha(Path(seed['checkpoint'])/'State.bin')==cp_hash,'rollback changed parent bytes')
    for key in ('fixed_reference','fixed_defect'):require(state_hash(seed[key])==state_hash(bad[key]),'rollback changed '+key)
    require(seed['physical_samples']==bad['physical_samples'],'rollback changed accepted physical samples')
    bad,_=native('changed-dt',1,Path(seed['checkpoint']),certfile,cert['version'],dt=1e-6);require(bad['failure'],'ordinary changed-dt restart accepted')
    for name,edit in [('runtime',lambda c:c['runtime'].update(sha256='0'*64)),('receipt',lambda c:c.update(native_receipt_sha256='0'*64)),('boundary',lambda c:c.update(prescribed_boundary_sha256='0'*64))]:
        c=load(certfile);edit(c);path=a.output/('corrupt-'+name+'.json');atomic(path,c);bad,_=native('reject-'+name,cert=path,version=cert['version']);require(bad.get('failure'),'corrupt V2 '+name+' accepted')
    # V1 remains a supported regression contract, with its original signature.
    from validate_cao_tamura_budget_h2d import certificate
    old=a.output/'v1-certificate.json';atomic(old,certificate(d,Path(receipt['output'])));v1,_=native('v1-import',cert=old);require(not v1['failure'] and same_history(parent,v1),'V1 regression changed')
    out['tests'].update(v2_seed_reference_history_bitwise=True,v2_same_dt_stage_packets_bitwise=count,v2_split_restart_bitwise=True,post_assignment_rollback=True,ordinary_changed_dt_rejected=True,corrupt_runtime_receipt_boundary_rejected=True,v1_regression_bitwise=True)
    out['status']='passed';save()
if __name__=='__main__':main()
