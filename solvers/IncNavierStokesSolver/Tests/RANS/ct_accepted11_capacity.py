#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Required observer capacity on an actual target frozen field, not evolution."""
import argparse,os,resource,subprocess,time,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted_pilot import replay,budget_classification
from ct_accepted_open import run,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha
from run_cls_epawf_picard import set_parameter,write
from ct_accepted11_gate import observer,capacity as capacity_gate
from ct_accepted11_equivalence import POLICY


def bounded_probe(command,folder):
    """Observe the native PID directly, including libraries and peak memory."""
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);mapped={};peak=virtual=0;expired=False;stopped=None
    with (folder/'probe.log').open('x') as log:
        p=subprocess.Popen([str(v) for v in command],cwd=folder,env=dict(os.environ,**THREADS),stdout=log,stderr=subprocess.STDOUT)
        while p.poll() is None:
            try:
                for line in (Path('/proc')/str(p.pid)/'maps').read_text().splitlines():
                    word=line.split()[-1];path=Path(word)
                    if word.startswith('/') and path.is_file() and word not in mapped:mapped[word]=sha(path)
                for line in (Path('/proc')/str(p.pid)/'status').read_text().splitlines():
                    if line.startswith('VmHWM:'):peak=max(peak,int(line.split()[1]))
                    if line.startswith('VmPeak:'):virtual=max(virtual,int(line.split()[1]))
            except (FileNotFoundError,ProcessLookupError):pass
            elapsed=time.monotonic()-start
            if elapsed>1200 and stopped is None:expired=True;stopped=elapsed;p.terminate()
            if stopped is not None and elapsed-stopped>30:p.kill()
            time.sleep(.1)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(command=[str(v) for v in command],pid=p.pid,exit_code=p.returncode,wall_seconds=time.monotonic()-start,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,peak_RSS_KiB=peak,peak_virtual_KiB=virtual,control_wall_limit_reached=expired,actual_mapped_files=mapped)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct9','runtime','preflight','contract','equivalence','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    equiv=sealed(a.equivalence);require(observer(equiv)['status']=='pass','observer equivalence prerequisites failed');eq_req=n.read(a.equivalence/'request-manifest.json')['request'];require(eq_req['runtime']['sha256']==sha(a.runtime/'RANSPreservationProbe'),'foreign observer capacity runtime');pre=sealed(a.preflight);capacity=sealed(a.ct9/'capacity-v2');c=capacity['cases']['target-p4']['case'];old=a.ct9/'capacity-v2/target-p4';initial=n.read(old/'frozen.json');state=initial['observation']['state'];require(exact_vector_hash(state)==pre['actual_initial_sha256'],'foreign target initial coefficients')
    m=RequestManifest(a.output,dict(schema='CT11AuditCapacityV1',case=c,equivalence=sha(a.equivalence/'results.json'),policy=POLICY,contract=sha(a.contract),preflight=sha(a.preflight/'results.json'),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source=sha(Path(__file__)),actual_initial_sha256=exact_vector_hash(state),controls=n.CONTROLS,nominal_dt=.001,scope='frozen actual initialized field only; no accepted numerical history or physical budget admission'))
    out=dict(schema='CT11AuditCapacityV1',status='incomplete',actual_initial_sha256=exact_vector_hash(state),policy=POLICY,controls={},physical_budget_qualification='not_supplied_by_this_static_fixture',native_trajectory='not_run',first_failure=None)
    for label,orders,face in n.CONTROLS:
        folder=a.output/label;folder.mkdir();root=ET.parse(old/'case.xml').getroot()
        for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
        write(root,folder/'case.xml');req=replay(initial,state,.001,[],export_physical_scalar_budget=True,physical_budget_volume_policy=POLICY,export_frozen_scalar_reconciliation=True,previous_state=state,stage_history=state,actual_flow_stage=state,actual_pre_k_stage=state,physical_face_order=face)
        write_json(folder/'request.json',req);e=bounded_probe([a.runtime/'RANSPreservationProbe','case.xml',folder/'request.json',folder/'packet.json'],folder);row=dict(orders=orders,face_order=face,execution=e,request_sha256=sha(folder/'request.json'),configuration_sha256=sha(folder/'case.xml'))
        packet=n.read(folder/'packet.json') if (folder/'packet.json').exists() else {};row['packet_sha256']=sha(folder/'packet.json') if packet else None
        if e['exit_code']!=0 or packet.get('failure'):
            row.update(status='failed_required_audit_capacity',reason=packet.get('failure','process terminated or exceeded declared bound'));out['first_failure']=dict(control=label,reason=row['reason'],state_sha256=exact_vector_hash(state),directory=str(folder));out['controls'][label]=row;out['status']='fail';write_json(a.output/'progress.json',out);break
        require(packet['observation']['state']==state and packet['observation']['identity']==initial['observation']['identity'],'capacity observer state or identity changed');b=budget_classification(packet,c);row['actual_initial_sha256']=exact_vector_hash(packet['observation']['state']);row['equations']=b['equations'];row['volume_observer']=packet['physical_scalar_budget']['volume_observer'];row['output_bytes']=(folder/'packet.json').stat().st_size;row['frozen_load_reconstruction']=b['reconstruction'];row['equation_names']=sorted(b['equations']);row['status']='pass' if b['reconstruction_status']=='pass' and set(b['equations'])=={'inner_k','outer_k','outer_epsilon','whole_k'} else 'failed_required_reconstruction';out['controls'][label]=row
        if row['status']!='pass':out['status']='fail';out['first_failure']=dict(control=label,reason=row['status'],directory=str(folder));write_json(a.output/'progress.json',out);break
        write_json(a.output/'progress.json',out);print('CT11 required audit capacity',label,row['status'],flush=True)
    if len(out['controls'])==4 and all(q['status']=='pass' for q in out['controls'].values()):out['status']='pass'
    out['qualification']=capacity_gate(out);out['status']=out['qualification']['status'];write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT11 required audit capacity',out['status'],out['first_failure'],flush=True);return 0 if out['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
