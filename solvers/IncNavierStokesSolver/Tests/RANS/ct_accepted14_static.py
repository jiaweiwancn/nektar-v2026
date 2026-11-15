#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete selected-space static admission; no advancement of physical histories."""
import argparse,gc,json,os,resource,shutil,xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.sparse import csr_matrix
import ct_accepted8_native as n
from ct_accepted14_component import ordering,COMPONENT
from ct_accepted9_controls import enable
from ct_accepted9_pressure_run import sealed
from ct_accepted11_equivalence import POLICY
from ct_accepted11_capacity import bounded_probe
from ct_accepted3_metrics import inventory
from ct_accepted10_compare import gradients,outlet,regional_pressure
from ct_accepted_pilot import replay,budget_classification
from ct_accepted_open import probe,run,require
from ct_accepted_feasibility import write_json
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct12','ct13','runtime','component','capacity','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1';n.reference=lru_cache(None)(n.reference)
    component=sealed(a.component);capacity=sealed(a.capacity);require(component['status']==capacity['status']=='pass','component/capacity prerequisite')
    previous=sealed(a.ct12/'static-v1');c=previous['case'];frozen=sealed(a.ct13/'exports-v2')
    m=RequestManifest(a.output,dict(schema='CT14CompleteSelectedStaticV1',source=sha(Path(__file__)),case=c,contract=sha(a.contract),component=sha(a.component/'results.json'),capacity=sha(a.capacity/'results.json'),CT13_exports=sha(a.ct13/'exports-v2/results.json'),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.runtime/'RANSPreservationProbe'),physical_advance=False))
    out=dict(schema='CT14CompleteSelectedStaticV1',status='incomplete',case=c,controls={},physical_scalar_controls={},canonical_exports={},physical_time=0.)
    folder=a.output/'target';n.make(folder,a.runtime,c,.001,0);enable(folder);ordering(folder)
    os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(folder/'sparse-receipts.jsonl')
    q,e=probe(a.runtime,folder,dict(operation='native_layout_inventory'),'count');require(e['exit_code']==0 and q['stored_slots']==62533 and q['elements']==961,'selected native inventory')
    out['count']=q;out['count_execution']=e
    e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');native=n.read(folder/'case.evolution.json');require(e['exit_code']==0 and not native['failure'] and native['time']==0.,'selected actual initialization')
    original=n.read(a.ct12/'static-v1/target/case.evolution.json');require(native['state']==original['state'] and COMPONENT in native['signature'],'initial representation/component changed')
    out['initial']=dict(execution=e,state_sha256=exact_vector_hash(native['state']),configuration_sha256=sha(folder/'case.xml'),native_report_sha256=sha(folder/'case.evolution.json'),history=str(native['checkpoint']),same_CT12_actual_coefficients=True,initialization_receipt=native['initialization_receipt'])
    grid,sides,weights,count=n.inventory(c);req=dict(operation='functions',grid=grid,grid_sides=sides,export_faces_only=True,export_test_basis=True,export_neutral_composite_reference=True,export_ct10_momentum=True,ct10_export_points=False,export_sparse_backend=True,sparse_alpha=1000.)
    q,e=probe(a.runtime,folder,req,'frozen');require(e['exit_code']==0,'native selected static:'+str(q.get('failure')));n.binding(q,c);require(q['observation']['state']==native['state'],'selected static changed coefficients')
    dom=inventory(q['test_basis'],q['observation']['grid'],np.asarray(native['state']),q['observation']['offsets'],'A95',q['observation']['grid_sides']);metrics=n.physical_metrics(q,native['state'],c,weights,count);reference=n.oracle_agreement(q,c)
    require(metrics['status']==reference['status']=='pass' and dom['domain_failure_count']==0 and dom['minimum_inner_k']['value']>0 and dom['minima']['covariance']['value']>=-1e-10*c['K0'] and dom['minimum_radicand']['radicand']>=0,'selected initial physical admission')
    require(all(v['scaled']<=(1e-4 if k=='w' else 1e-3) for k,v in metrics['absolute_reference_error'].items()),'selected physical reference projection error')
    mom=q['ct10_momentum'];require(mom['maximum_reconstruction_error']<=2e-10 and mom['maximum_input_lift_difference']==0.,'all native signed load rows/lift')
    out.update(initial_metrics=metrics,initial_gradients=gradients(q,native['state'],c,weights,count),outlet=outlet(q,native['state'],c),reference_agreement=reference,physical_domain=dom,native_audit=q['observation']['audit_alpha_zero'],all_mode_momentum_reconstruction=mom['maximum_reconstruction_error'],native_packet_sha256=sha(folder/'frozen.json'))
    off=q['observation']['offsets'];npres=off[6]-off[4];rr=[];cc=[];vv=[]
    for i,(row,region) in enumerate(zip(q['test_basis'][:count],q['observation']['grid_regions'][:count],strict=True)):
        for j,v in zip(row['pressure']['ids'],row['pressure']['phi'],strict=True):rr.append(i);cc.append(off[4+region]+j-off[4]);vv.append(v)
    B=csr_matrix((vv,(rr,cc)),shape=(count,npres));base=np.asarray(q['sparse_backend']['pressure_increment'])
    def record(packet,execution,alpha):
        b=packet['sparse_backend'];require(packet['observation']['state']==native['state'] and b['pressure_C_nnz']==0 and b['pressure_rows']==npres,'canonical state/pressure convention')
        for kind in ('predictor','pressure'):
            key='true_pressure_residual' if kind=='pressure' else 'predictor_residual';require(b[key]<=2e-12*b[kind+'_rhs_maximum']+1e-25,'original true residual:'+kind)
        for v in b['scalars']:require(v['maximum_true_residual']<=2e-12*v['rhs_maximum']+1e-25 and v['factor_nnz']<=40000000,'scalar action/residual/fill')
        require(max(b[k] for k in ('predictor_factor_nnz','consistent_mass_factor_nnz','pressure_preconditioner_factor_nnz'))<=40000000,'static factor fill')
        require(max(b[k] for k in ('full_native_G_action_defect','consistent_mass_correction_true_residual','splitting_identity'))<=2e-10,'full G/consistent-mass/splitting actions')
        delta=B@(np.asarray(b['pressure_increment'])-base);row={k:v for k,v in b.items() if k not in ('predictor','pressure_increment','consistent_mass_correction','scalars','identity')};row['scalars']=[{k:v for k,v in x.items() if k!='diagnostic_increment'} for x in b['scalars']]
        row.update(execution=execution,alpha=alpha,pressure_change_scaled_maximum=float(np.max(abs(delta),initial=0))/c['ustar']**2,pressure_change_scaled_RMS=float(np.sqrt(np.asarray(weights)@(delta*delta)/sum(weights)))/c['ustar']**2)
        return row
    out['controls']['native']=record(q,e,1000.)
    for label,orders,alpha in [('alpha2000',None,2000.),('q48',(48,12,12,16,64),1000.),('q64',(64,16,12,20,64),1000.)]:
        dest=a.output/label;dest.mkdir();root=ET.parse(folder/'case.xml').getroot()
        if orders:
            for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
        write(root,dest/'case.xml');os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(dest/'sparse-receipts.jsonl')
        request=replay(q,native['state'],0.,[],export_basis=False,export_faces_only=True,export_sparse_backend=True,sparse_alpha=alpha)
        packet,execution=probe(a.runtime,dest,request,'frozen');require(execution['exit_code']==0,'selected '+label+':'+str(packet.get('failure')))
        row=record(packet,execution,alpha);row['packet_sha256']=sha(dest/'frozen.json');out['controls'][label]=row
        if orders:require(row['pressure_change_scaled_maximum']<1e-8,'static governing-quadrature pressure uncertainty')
        write_json(a.output/'progress.json',out);del packet;gc.collect();print('CT14 static',label,'pass',flush=True)
    # Canonical native assembly and restriction must still equal the four-bound
    # capacity inputs; the export does not solve or advance a history.
    for alpha in (1000.,2000.):
        label='candidate-alpha'+str(int(alpha));dest=a.output/('canonical-'+str(int(alpha)));dest.mkdir();shutil.copy2(folder/'case.xml',dest/'case.xml')
        request=replay(q,native['state'],0.,[],stored_reference_sha256=q['initialization_receipt']['fixed_reference_sha256'],export_pre_factor_frozen=True,actual_history=native['state'],sparse_alpha=alpha,frozen_request_kind='InitializedFrozenV1',pre_factor_prefix=str(dest/'frozen'))
        packet,e=probe(a.runtime,dest,request,'prefactor');require(e['exit_code']==0,'new component canonical export:'+str(packet.get('failure')));f=packet['pre_factor_frozen'];prior=frozen['systems'][label];old=n.read(Path(prior['packet']))['pre_factor_frozen']
        for key in ('A','G','D','C','free_velocity_mass','frozen_flow_matrix'):require(f[key]['sha256']==old[key]['sha256'],'canonical '+key+' changed')
        for key in ('rows','frozen_flow_residual','history','input','effective_input','lifting'):require(f[key]==old[key],'canonical '+key+' changed')
        require(f['whole_command_sparse_factorizations']==0,'factor-free runtime integration')
        out['canonical_exports'][str(int(alpha))]=dict(status='pass',packet_sha256=sha(dest/'prefactor.json'),A=f['A'],actual_history_sha256=f['history_sha256'],complete_canonical_equality=True,whole_command_sparse_factorizations=0,execution=e)
        del packet,old;gc.collect()
    os.environ.pop('NEKTAR_RANS_SPARSE_RECEIPTS',None)
    # These explicitly initialized frozen fixtures qualify observer capacity,
    # complete scalar source/face/reconstruction modes, not physical evolution.
    for label,orders,face in n.CONTROLS:
        dest=a.output/('scalar-'+label);dest.mkdir();root=ET.parse(folder/'case.xml').getroot()
        for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
        write(root,dest/'case.xml');request=replay(q,native['state'],.001,[],export_basis=False,export_physical_scalar_budget=True,physical_budget_volume_policy=POLICY,export_frozen_scalar_reconciliation=True,previous_state=native['state'],stage_history=native['state'],actual_flow_stage=native['state'],actual_pre_k_stage=native['state'],physical_face_order=face)
        write_json(dest/'request.json',request);execution=bounded_probe([a.runtime/'RANSPreservationProbe','case.xml',dest/'request.json',dest/'packet.json'],dest);packet=n.read(dest/'packet.json') if (dest/'packet.json').exists() else {}
        require(execution['exit_code']==0 and not packet.get('failure'),'selected scalar capacity:'+label+':'+str(packet.get('failure')));b=budget_classification(packet,c)
        require(b['reconstruction_status']=='pass' and set(b['equations'])=={'inner_k','outer_k','outer_epsilon','whole_k'},'complete static physical scalar modes')
        require(packet['observation']['state']==native['state'] and packet['observation']['identity']==q['observation']['identity'],'scalar fixture actual state/component')
        out['physical_scalar_controls'][label]=dict(status='pass',classification=b,volume_observer=packet['physical_scalar_budget']['volume_observer'],execution=execution,packet_sha256=sha(dest/'packet.json'),configuration_sha256=sha(dest/'case.xml'),actual_initial_sha256=exact_vector_hash(native['state']),role='InitializedFrozenV1; capacity/reconstruction,not an evolved consecutive-state budget')
        write_json(a.output/'progress.json',out);del packet;gc.collect();print('CT14 static scalar',label,'pass',flush=True)
    out['status']='pass';out['scope']='complete selected-space static/operator/diagnostic admission only; fresh physical pair separately contracted and gated';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT14 selected static pass',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
