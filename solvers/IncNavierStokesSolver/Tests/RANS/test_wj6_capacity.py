#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compiled capacity/preflight and exact physical off/on transaction controls."""
import argparse
import json
import os
from pathlib import Path
import time
import uuid
import xml.etree.ElementTree as ET

from wj2_cases import case_xml, normal_knots, local_normal_knots
from cls_epawf_m3_cases import set_info, write
from wj_immutable import Ledger, artifacts, digest, file_hash, read, write_once
from wj6_resources import DEFAULT, run_limited
from run_cls_epawf_picard import THREADS, require, state_hash
from test_wj2_native import digest_object
from compare_cls_epawf_pimple import same_history


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--policy',choices=('EP','Budget'),required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS)
    root=a.output.resolve()/a.policy
    root.mkdir(parents=True,exist_ok=False)
    ledger=Ledger(Path(os.environ['WJ4_PACKAGE_ROOT']),a.policy,cap=800)
    session='capacity-'+uuid.uuid4().hex;ledger.begin_session(session)
    design=dict(schema='WJ6CapacityControlResourceV1',limits=DEFAULT,policy=a.policy,
        native_capacity=4331,scope='finite serial WJ capacity control, not scientific qualification')
    write_once(root/'resource-design.json',design);resource_hash=file_hash(root/'resource-design.json')
    records=[]

    def case(name,steps=0,capacity=None,cp=None,mesh='small',extras=None):
        folder=root/name
        z=normal_knots(2) if mesh in ('L2Nx2','L2Nx1') else local_normal_knots('M2') if mesh=='M2' else normal_knots(0)
        x=[0.,1.] if mesh in ('small','L2Nx1') else [0.,.5,1.]
        case_xml(folder,a.policy,steps,cp,dt=3.2e-5,xknots=x,zknots=z)
        xml=folder/'case.xml';tree=ET.parse(xml).getroot()
        if capacity is not None:
            set_info(tree,'RANSWallDenseEvolutionCapacity',str(capacity))
            set_info(tree,'RANSWallDenseEvolutionResourceDesign',resource_hash)
        for key,value in (extras or {}).items():set_info(tree,key,value)
        write(tree,xml)
        return folder,xml

    def execute(folder,argv,steps=0,reason=None):
        label=str(folder)+'/'+str(len(records));intent=dict(argv=list(map(str,argv)),inputs=artifacts([folder/'case.xml']),code=artifacts([Path(__file__)]),resource_design_sha256=resource_hash)
        write_once(folder/('intent-'+str(len(records))+'.json'),intent)
        ledger.reserve(label,intent,steps,600)
        start=time.time()
        try:
            metrics=run_limited(argv,folder,folder/('command-'+str(len(records))+'.log'),DEFAULT)
            report=folder/'case.evolution.json'
            data=read(report) if report.exists() and metrics['exit_code']==0 else None
            if reason:
                require(metrics['exit_code']!=0 and reason in (folder/('command-'+str(len(records))+'.log')).read_text(),
                        'capacity rejection did not retain intended reason: '+str(reason))
            else:require(metrics['exit_code']==0,'compiled capacity command failed: '+str(folder))
            record=dict(folder=str(folder),intent=intent,resource=metrics,expected_failure=reason,artifacts=artifacts([f for f in folder.rglob('*') if f.is_file()]))
            records.append(record);write_once(folder/'capacity-receipt.json',record)
            return data
        finally:ledger.finish(label,steps,time.time()-start,str(folder/'capacity-receipt.json'))

    try:
        for token in ('2999','4332','4330.5','nan','inf','-inf','1e309','4331junk'):
            folder,xml=case('invalid-'+token.replace('.','p').replace('-','m'),capacity=token)
            execute(folder,[a.solver,xml],reason='dense_capacity_invalid')
        folder,xml=case('out-of-route',capacity=4331,extras={'TurbulenceModel':'CLSHighRe','RANSWallFormulation':'EPAWF','RANSWallKSampling':'Point','RANSWallTransition':'Fixed'})
        execute(folder,[a.solver,xml],reason='dense_capacity_route')
        folder,xml=case('missing-resource',capacity=4331,extras={'RANSWallDenseEvolutionResourceDesign':'None'})
        execute(folder,[a.solver,xml],reason='dense_capacity_resource_design')
        large,xml=case('large-default-audit',mesh='L2Nx2')
        large_data=execute(large,[a.solver,xml]);require(len(large_data['state'])==4331,'large default audit blocked or wrong size')
        folder,xml=case('large-default-reject',steps=1,cp=large_data['checkpoint'],mesh='L2Nx2')
        execute(folder,[a.solver,xml],1,reason='dense_evolution_capacity_exceeded')
        folder,xml=case('large-default-readback',cp=large_data['checkpoint'],mesh='L2Nx2')
        back=execute(folder,[a.solver,xml]);same_history(large_data,back)
        dimensions={}
        for mesh,expected in [('L2Nx2',[4331,1790,187,1977,900,242]),('M2',[2975,1250,133,1383,630,122]),('L2Nx1',[2493,716,125,841,360,98])]:
            folder,xml=case('preflight-'+mesh,capacity=4331,mesh=mesh)
            packet=folder/'input.json';write_once(packet,dict(operation='dense_capacity_preflight',inner_stress_model='Linear'))
            output=folder/'result.json';execute(folder,[a.probe,xml,packet,output]);q=read(output)
            actual=[q[k] for k in ('state_size','velocity_size','pressure_gauge_size','assembled_flow_size','k_solve_size','epsilon_solve_size')]
            require(actual==expected and q['governing_solves']==0 and q['state_unchanged'],'capacity actual frozen dimensions differ')
            dimensions[mesh]=q
        zero,xml=case('small-cp0');initial=execute(zero,[a.solver,xml])
        off,xml=case('small-off',2,cp=initial['checkpoint']);d0=execute(off,[a.solver,xml],2)
        on,xml=case('small-on',2,4331,initial['checkpoint']);d1=execute(on,[a.solver,xml],2)
        require(d0['signature']==d1['signature'],'execution capacity changed physical checkpoint identity')
        same_history(d0,d1)
        for key in ('accepted_states','audit_times','state_audits','wj_audit_packets','physical_samples'):
            require(digest_object(d0[key])==digest_object(d1[key]),'capacity changed native accepted packet: '+key)
        for x,y in zip(d0['pimple_attempts'],d1['pimple_attempts']):
            for key in ('history_sha256','incoming_sha256','momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls','bound_activations'):
                require(x[key]==y[key],'capacity changed scheduled action/count: '+key)
            for key in ('raw','corrected'):
                require(state_hash(read(x['k_snapshot'])[key])==state_hash(read(y['k_snapshot'])[key]),'capacity changed raw/corrected k')
        folder,xml=case('small-on-to-off-readback',cp=d1['checkpoint']);back=execute(folder,[a.solver,xml]);same_history(d1,back)
        require(d1['dense_evolution_capacity']['effective']==4331 and d0['dense_evolution_capacity']['effective']==3000,'capacity report is missing')
        write_once(root/'result.json',dict(schema='WJ6CapacityControlsV1',status='PASS',dimensions=dimensions,
            exact_signatures_accepted_states_gauge_histories_audits_samples=True,unchanged_counts_raw_corrected_k=True,
            inverse_checkpoint_readback=True,records=records,resource_design=artifacts([root/'resource-design.json'])))
        print(root/'result.json')
    finally:ledger.finish_session(session)


if __name__=='__main__':main()
