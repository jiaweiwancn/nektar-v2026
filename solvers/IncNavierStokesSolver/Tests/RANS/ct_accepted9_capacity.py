#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Increasing actual flat layouts: count, fresh startup, frozen CSR actions.

Zero time only. A pressure-study failure forbids physical full-domain evolution
but does not forbid this separately authorized capacity measurement.
"""
import argparse,copy,os,resource
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_controls import enable
from ct_accepted9_pressure_run import sealed
from ct_accepted7_domain import X
from ct_accepted_open import run,probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

def samples(c):
    g=[];s=[]
    for x in sorted(set([c['x'][0],c['x'][-1],0.,-2.5,2.5])):
        if not c['x'][0]<=x<=c['x'][-1]:continue
        for z in sorted(set(c['knots']+[c['roughness_input'],n.col(c).cutoff])):
            if z==c['H']:
                for side in (-1,1):g.append([x,0.,z]);s.append(side)
            else:g.append([x,0.,z]);s.append(0)
    return g,s

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','ct8','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    pre=sealed(a.ct8/'preflight-v6');old=sealed(a.ct8/'assessment-v3');inventory=n.read(a.ct8/'assessment-v3/domain.json');c=copy.deepcopy(pre['case']);radial=inventory['radial_grid'];xs=sorted(set(X[::4]+[X[-1]]))
    m=RequestManifest(a.output,dict(schema='CT9FlatCapacityV1',case=c,source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.runtime/'RANSPreservationProbe'),ct8_domain=sha(a.ct8/'assessment-v3/domain.json'),radial_grid=radial,layouts=[dict(name='intermediate-p3',x=xs,degree=3),dict(name='target-p3',x=X,degree=3),dict(name='target-p4',x=X,degree=4)],maximum_address_space_bytes=16*1024**3,physical_evolution=False,primary_prospective_layout='target-p4, contingent on static representation and all independent admissions'))
    result=dict(schema='CT9FlatCapacityV1',cases={},physical_evolution='not_admitted_pressure_uncertainty_failed',terrain_field_count='not admitted by affine global-height region classifier',status='incomplete')
    for label,x,degree in [('intermediate-p3',xs,3),('target-p3',X,3),('target-p4',X,4)]:
        cc=copy.deepcopy(c);cc.update(x=x,knots=radial,Z=20.,degree=degree,pressure_degree=degree-2);folder=a.output/label;n.make(folder,a.runtime,cc,.001,0);enable(folder);row=dict(case=cc,configuration_sha256=sha(folder/'case.xml'))
        try:
            count,e=probe(a.runtime,folder,dict(operation='native_layout_inventory'),'count');require(e['exit_code']==0,'actual native count failed');row.update(native_count=count,count_execution=e);require(count['stored_slots']<=70000,'new backend count guard')
            if label.startswith('target'):require(count['stored_slots']==(33267 if degree==3 else 60907),'target changed native count')
            os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'initial-timing.jsonl');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE',None);row['startup_execution']=e
            require((folder/'case.evolution.json').exists(),'fresh zero-step report unavailable');native=n.read(folder/'case.evolution.json');row.update(zero_time=native['time'],zero_state_sha256=exact_vector_hash(native['state']),history=native['checkpoint'],native_report_sha256=sha(folder/'case.evolution.json'),initialization=native['initialization_receipt'],startup_physical_failure=native['failure'],initial_raw_spatial_norms=native.get('raw_reference_norms'))
            require(e['exit_code']==0 and not native['failure'] and native['time']==0.,'zero-step startup/domain failed:'+str(native['failure']))
            g,s=samples(cc);req=dict(operation='functions',grid=g,grid_sides=s,export_test_basis=True,export_constitutive=True,export_neutral_composite_reference=True,export_sparse_backend=True,sparse_alpha=500.)
            os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'frozen-timing.jsonl');q,e=probe(a.runtime,folder,req,'frozen');os.environ.pop('NEKTAR_RANS_TIMING_FILE',None);row['frozen_execution']=e;require(e['exit_code']==0,'frozen native sparse operation failed:'+str(q.get('failure')))
            n.binding(q,cc);require(q['observation']['state']==native['state'],'fresh provider/replay projected state mismatch');row.update(frozen_packet_sha256=sha(folder/'frozen.json'),sparse_backend=q['sparse_backend'],reference_agreement=n.oracle_agreement(q,cc),native_domain_and_spatial=q['observation']['audit_alpha_zero'],actual_quadrature=q['observation']['actual_quadrature'])
            require(row['reference_agreement']['status']=='pass','continuous physical reference differs')
            row['integration_breakpoints']=dict(mesh_nodes=radial,analytic_cutoff=n.col(cc).cutoff,roughness=cc['roughness_input'],viscosity_intercept=n.col(cc).zv,cutoff_is_mesh_node=n.col(cc).cutoff in radial,roughness_is_mesh_node=cc['roughness_input'] in radial,scope='governing wall integration retains actual analytic/current breakpoints; finite-element projection can cross non-mesh breaks, and is assessed separately')
            row['status']='static_capacity_pass';row['physical_preservation_qualification']='not_run'
        except Exception as error:row.update(status='failed_capacity_or_reference_preflight',failure=repr(error));os.environ.pop('NEKTAR_RANS_TIMING_FILE',None)
        row['output_bytes']=sum(f.stat().st_size for f in folder.rglob('*') if f.is_file());result['cases'][label]=row;write_json(a.output/'progress.json',result);print('CT9 capacity',label,row['status'],flush=True)
    result['status']='pass' if len(result['cases'])==3 and all(q['status']=='static_capacity_pass' for q in result['cases'].values()) else 'fail';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
