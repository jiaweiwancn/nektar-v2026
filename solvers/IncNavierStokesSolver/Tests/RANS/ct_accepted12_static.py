#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""The single target-specific static inner-slab contrast; no time advancement."""
import argparse,copy,os,resource,time,xml.etree.ElementTree as ET
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.sparse import csr_matrix
import ct_accepted8_native as n
from ct_accepted9_controls import enable
from ct_accepted9_pressure_run import sealed
from ct_accepted10_compare import gradients
from ct_accepted3_metrics import inventory
from ct_accepted_open import probe,run,require
from ct_accepted_pilot import replay,artifact_hash
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha
from run_cls_epawf_picard import set_parameter,write


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct10','ct11','runtime','probe','pressure','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    start=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    def spent():
        now=resource.getrusage(resource.RUSAGE_CHILDREN)
        return time.process_time()-start+now.ru_utime+now.ru_stime-before.ru_utime-before.ru_stime
    contract=n.read(a.contract);pressure=sealed(a.pressure);old=sealed(a.ct11/'target-v1');pre=sealed(a.ct10/'target-preflight-v1')
    require(sha(a.pressure/'results.json')==contract['pressure_receipt']['sha256'],'foreign target attribution')
    band=max(pressure['bands']['rows'],key=lambda r:abs(r['signed_alignment']))
    require([band['lower'],band['upper']]==contract['selection']['interval'],'unsupported static selection')
    c=copy.deepcopy(old['case']);original=copy.deepcopy(c);c['knots']=sorted(set(c['knots']+[contract['selection']['added_point']]))
    require(c['degree']==4 and c['pressure_degree']==2 and c['x'][0]==-25. and c['x'][-1]==30. and c['Z']==20.,'foreign target space/domain')
    require(all(c[k]==v for k,v in original.items() if k!='knots'),'physical input changed')
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');gov=runtime_identity(a.runtime/'IncNavierStokesSolver')
    os.environ['LD_LIBRARY_PATH']=str(a.probe/'lib');pr=runtime_identity(a.probe/'RANSPreservationProbe')
    require(gov['sha256']=='5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7' and pr['sha256']=='20236b14eb4a90dea0f2a9a38b92e41135a83e9241c978fb265e185098f586e8','foreign runtime')
    manifest=RequestManifest(a.output,dict(schema='CT12StaticInnerContrastV1',source=sha(Path(__file__)),contract=sha(a.contract),pressure=sha(a.pressure/'results.json'),target=sha(a.ct11/'target-v1/results.json'),case=c,runtime=gov,probe=pr,physical_evolution=False))
    write_json(a.output/'paired-definition.json',dict(x=c['x'],radial=c['knots'],roof=c['Z'],H=c['H'],flat_amplitude=0.,terrain_definition='same CT8 joined collar/roof map with this identical logical radial amendment; no admitted mapped native field inventory',physical_manifest_sha256=contract['physical_manifest_sha256']))
    out=dict(schema='CT12StaticInnerContrastV1',status='incomplete',case=c,controls={},physical_evolution=False,physical_horizon=0.)
    folder=a.output/'target';os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');n.make(folder,a.runtime,c,.001,0);enable(folder)
    os.environ['LD_LIBRARY_PATH']=str(a.probe/'lib')
    try:
        q,e=probe(a.probe,folder,dict(operation='native_layout_inventory'),'count');out['count']=q;out['count_execution']=e
        require(e['exit_code']==0 and q['stored_slots']==62533 and q['elements']==961,'native count differs from prospective prediction')
        require(spent()<1400,'static process CPU budget before startup')
        os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');out['startup_execution']=e
        r=n.read(folder/'case.evolution.json');out['initial']=dict(time=r['time'],failure=r['failure'],state_sha256=exact_vector_hash(r['state']),history_sha256=artifact_hash(r['checkpoint']) if r['checkpoint'] else None,native_report_sha256=sha(folder/'case.evolution.json'),configuration_sha256=sha(folder/'case.xml'),raw_spatial_norms=r.get('raw_reference_norms'))
        require(e['exit_code']==0 and not r['failure'] and r['time']==0.,'static zero-step reference/domain rejection')
        g,s,w,count=n.inventory(c);n.reference=lru_cache(None)(n.reference)
        os.environ['LD_LIBRARY_PATH']=str(a.probe/'lib')
        req=dict(operation='functions',grid=g,grid_sides=s,export_faces_only=True,export_test_basis=True,export_neutral_composite_reference=True,export_ct10_momentum=True,ct10_export_points=False,export_sparse_backend=True,sparse_alpha=1000.)
        q,e=probe(a.probe,folder,req,'frozen');out['frozen_execution']=e
        require(e['exit_code']==0,'static frozen operation:'+str(q.get('failure')))
        n.binding(q,c);require(q['observation']['state']==r['state'],'static observer changed actual projection')
        mom=q['ct10_momentum'];off=q['observation']['offsets'];npr=off[6]-off[4]
        domain=inventory(q['test_basis'],q['observation']['grid'],np.asarray(r['state']),off,'A95',q['observation']['grid_sides'])
        out.update(initial_metrics=n.physical_metrics(q,r['state'],c,w,count),initial_gradients=gradients(q,r['state'],c,w,count),reference_agreement=n.oracle_agreement(q,c),physical_domain=domain,
            native_audit=q['observation']['audit_alpha_zero'],all_row_load_reconstruction=mom['maximum_reconstruction_error'],effective_lift=mom['maximum_input_lift_difference'],packet_sha256=sha(folder/'frozen.json'))
        rr=[];cc=[];vv=[]
        for i,(row,region) in enumerate(zip(q['test_basis'][:count],q['observation']['grid_regions'][:count],strict=True)):
            for j,value in zip(row['pressure']['ids'],row['pressure']['phi'],strict=True):rr.append(i);cc.append(off[4+region]+j-off[4]);vv.append(value)
        B=csr_matrix((vv,(rr,cc)),shape=(count,npr));base=np.asarray(q['sparse_backend']['pressure_increment']);us2=c['ustar']**2
        def record(packet,execution):
            b=packet['sparse_backend'];dp=B@(np.asarray(b['pressure_increment'])-base)
            return dict(execution=execution,pressure_C_nnz=b['pressure_C_nnz'],predictor_factor_nnz=b['predictor_factor_nnz'],flow_nnz=b['flow_nnz'],pressure_rows=b['pressure_rows'],pressure_iterations=b['pressure_iterations'],
                true_pressure_residual=b['true_pressure_residual'],pressure_rhs_maximum=b['pressure_rhs_maximum'],true_continuity=b['true_continuity'],
                predictor_residual=b['predictor_residual'],predictor_rhs_maximum=b['predictor_rhs_maximum'],full_native_G_action_defect=b['full_native_G_action_defect'],consistent_mass_correction_true_residual=b['consistent_mass_correction_true_residual'],splitting_identity=b['splitting_identity'],
                pressure_change_scaled_maximum=float(np.max(abs(dp),initial=0.))/us2,pressure_change_scaled_RMS=float(np.sqrt(np.asarray(w)@(dp*dp)/sum(w)))/us2,
                state_sha256=exact_vector_hash(packet['observation']['state']),input_unchanged=packet['observation']['state']==r['state'])
        # Field names are checked against the actual established native schema.
        out['controls']['native']=record(q,e)
        require(out['controls']['native']['pressure_C_nnz']==0 and out['controls']['native']['pressure_rows']==npr and out['controls']['native']['predictor_factor_nnz']<=40000000,'static native pressure/fill guard')
        require(domain['minima']['covariance']['value']>=-1e-10*c['K0'] and domain['minimum_radicand']['radicand']>=0 and domain['minimum_inner_k']['value']>0,'static sampled physical-domain margins')
        metrics=out['initial_metrics'];require(metrics['status']=='pass' and out['reference_agreement']['status']=='pass' and domain['domain_failure_count']==0 and out['all_row_load_reconstruction']<=2e-10 and out['effective_lift']==0.,'static initial representation/domain/load gate')
        require(all(v['scaled']<=(1e-4 if k=='w' else 1e-3) for k,v in metrics['absolute_reference_error'].items()),'static absolute reference accuracy')
        for label,orders,alpha in [('alpha2000',None,2000.),('q48',(48,12,12,16,64),1000.),('q64',(64,16,12,20,64),1000.)]:
            require(spent()<1400,'static process CPU budget before '+label)
            dest=a.output/label;dest.mkdir();root=ET.parse(folder/'case.xml').getroot()
            if orders:
                for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
            write(root,dest/'case.xml');req=replay(q,r['state'],0.,[],export_neutral_composite_reference=True,export_sparse_backend=True,sparse_alpha=alpha)
            packet,execution=probe(a.probe,dest,req,'frozen');out['controls'][label]=dict(execution=execution,packet_sha256=sha(dest/'frozen.json'),failure=packet.get('failure'))
            require(execution['exit_code']==0,'static '+label+':'+str(packet.get('failure')))
            n.binding(packet,c);row=record(packet,execution);row['packet_sha256']=sha(dest/'frozen.json');out['controls'][label]=row
            require(row['input_unchanged'] and row['pressure_C_nnz']==0 and row['predictor_factor_nnz']<=40000000 and row['true_pressure_residual']<=2e-10*max(1.,row['pressure_rhs_maximum']) and row['full_native_G_action_defect']<=2e-10 and row['consistent_mass_correction_true_residual']<=2e-10 and row['splitting_identity']<=2e-10,'static pressure/operator/capacity gate')
            if orders:require(row['pressure_change_scaled_maximum']<1e-8,'static governing-quadrature pressure uncertainty')
            write_json(a.output/'progress.json',out)
        require(spent()<1400,'static process CPU budget exhausted');out['status']='static_pass_only'
    except Exception as error:
        out['status']='blocked_or_failed_static_preflight';out['failure']=repr(error)
    out.update(process_CPU=spent(),scope='fresh zero-time flat representation only; no new timestep trajectory, no pressure-uncertainty reclassification or terrain admission')
    write_json(a.output/'results.json',out);manifest.seal(a.output,['results.json','paired-definition.json']);print(out['status'],out.get('failure'),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
