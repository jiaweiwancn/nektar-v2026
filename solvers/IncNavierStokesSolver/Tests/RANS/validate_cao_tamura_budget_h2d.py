#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact certified-startup import, restart and rollback; no seed regeneration."""
import argparse,os,json,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cao_tamura_budget_h2d import case
from continue_cls_epawf_pimple import freeze,atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import THREADS,load,sha,write,set_info,set_parameter,require,state_hash
from compare_cls_epawf_pimple import same_history

from cao_tamura_h2d_state import load

def certificate(d,source):
    require(d['audit']['accepted_awf_hard_admissible']==[1] and d['audit']['admissible']==[1] and d['audit']['accepted_awf_scalar_domain'][8]<=1e-10,'cannot certify failed native inventory')
    require(max(d['audit']['coupled_norms'][4:])<1e-9,'native algebraic/boundary initializer constraints failed')
    if d.get('extra_sampling'):require(d['extra_sampling']['minimum_covariance']>=0,'additional covariance sampling failed')
    return dict(version='CertifiedNeutralStartupV1',physical_t0=0,complete_native_admissible=True,
        state=d['state'],state_sha256=d['state_sha256'],patch_signature=d['patch_signature'],geometry_identity=d['geometry_identity'],
        inventory_identity=d['inventory_identity'],validation_report=str(source),validation_sha256=sha(source),
        complete_audit=d['audit'],extra_sampling=d.get('extra_sampling'),immutable_seed=True,
        semantics='fresh each-dt IMEX1 histories; exact coupled seed; no optimizer/Lift/projection on import')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('solver','probe','preflight','h1','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(exist_ok=False);solver,runtime=freeze(a.solver,a.output/'runtime');probe,pimage=freeze(a.probe,a.output/'probe-runtime');pre=load(a.preflight);out=dict(status='running',runtime=runtime,probe=pimage,results={},tests={})
    def save():atomic(a.output/'progress.json',out)
    def native(name,steps=0,cp=None,cert=None,dt=2e-6,mutate=None):
        xml=case(a.output/name,a.h1,pre,flat=True);t=ET.parse(xml).getroot();set_parameter(t,'NumSteps',steps);set_parameter(t,'TimeStep',dt)
        if cp:set_info(t,'RANSWallRestart',str(cp))
        if cert:set_info(t,'RANSWallCertifiedStartupFile',str(cert));set_info(t,'RANSWallCertifiedStartupSHA256',sha(cert))
        if mutate:mutate(t)
        write(t,xml);env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),NEKTAR_RANS_TIMING_FILE=str(xml.parent/'timing.jsonl'))
        cmd=['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(solver),str(xml)]
        with (xml.parent/'solver.log').open('w') as f:r=subprocess.run(cmd,cwd=xml.parent,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
        report=xml.with_suffix('.evolution.json');d=load(report) if report.exists() else {};out['results'][name]=dict(code=r.returncode,report=str(report),report_sha256=sha(report) if report.exists() else None,usage=(xml.parent/'usage.txt').read_text(),failure=d.get('failure'));save();return d
    save()
    try:
        xml=case(a.output/'flat-certification',a.h1,pre,flat=True)
        d,receipt=run_probe(probe,xml,dict(operation='h2c_seed',inner_stress_model='Linear',xml_sha256=sha(xml)),a.output/'flat-certification-probe')
        d=load(Path(receipt['output']))
        cert=a.output/'certificate.json';atomic(cert,certificate(d,Path(receipt['output'])))
        q=native('ordinary-h0');i=native('import-h0',cert=cert)
        require(not q['failure'] and not i['failure'] and same_history(q,i),'exact imported zero seed/history differs')
        require(state_hash(i['state'])==d['state_sha256'],'import silently regenerated/changed seed')
        ordinary=native('ordinary2',steps=2,cp=Path(q['checkpoint']));direct=native('import2',steps=2,cp=Path(i['checkpoint']),cert=cert)
        require(same_history(ordinary,direct),'import differs from ordinary same-dt first/second step')
        packet_comparisons=0
        # Stage vectors are separately checked; metadata identity may differ.
        for x,y in zip(ordinary['pimple_attempts'],direct['pimple_attempts']):
            for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(x[k]==y[k]==1,'wrong governing schedule')
            require(x['solved_epsilon_retained_bitwise'] and y['solved_epsilon_retained_bitwise'],'epsilon overwritten')
            if x.get('k_snapshot') and Path(x['k_snapshot']).exists():
                px,py=load(Path(x['k_snapshot'])),load(Path(y['k_snapshot']))
                for k in ('committed','post_flow','pre_k','raw','corrected','history','fixed_reference'):require(state_hash(px[k])==state_hash(py[k]),'certified native stage differs '+k)
                packet_comparisons+=1
        require(packet_comparisons>0,'no complete verification attempt packet retained')
        first=native('import1',steps=1,cp=Path(i['checkpoint']),cert=cert);split=native('import1-restart',steps=1,cp=Path(first['checkpoint']),cert=cert)
        require(same_history(direct,split),'certified split restart differs')
        parent_sha=sha(Path(i['checkpoint'])/'State.bin');bad=native('post-assignment-reject',steps=1,cp=Path(i['checkpoint']),cert=cert,mutate=lambda t:set_parameter(t,'RANSWallPimplePostAssignRejectStep',1))
        require(bad['failure'] and same_history(bad,i) and sha(Path(i['checkpoint'])/'State.bin')==parent_sha,'rollback changed parent or histories')
        for name,mutation,dt in [('fit',lambda t:set_info(t,'RANSWallNeutralFitSHA256','0'*64),2e-6),('ordinary-dt',None,1e-6),('mesh',lambda t:set_parameter(t,'RANSWallHeight',.006),2e-6)]:
            bad=native('reject-'+name,cp=Path(i['checkpoint']),cert=cert,dt=dt,mutate=mutation);require(not bad or bad['failure'],'incompatible restart/import accepted')
        corrupt=load(cert);corrupt['state'][10]+=1e-4;badcert=a.output/'corrupt-certificate.json';atomic(badcert,corrupt)
        bad=native('reject-corrupt-state',cert=badcert);require(not bad or bad['failure'],'corrupt coupled vector accepted')
        out.update(status='passed',tests=dict(seed_reference_history_bitwise=True,two_steps_bitwise=True,stage_packets_bitwise=True,split_restart_bitwise=True,post_assignment_rollback=True,fit_mesh_dt_and_corruption_rejected=True))
    except Exception as e:out.update(status='failed',failure=str(e));raise
    finally:save()
if __name__=='__main__':main()
