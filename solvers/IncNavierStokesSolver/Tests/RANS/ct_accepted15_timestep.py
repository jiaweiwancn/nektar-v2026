#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Admit only the changed alpha4000/8000 actions before the fresh CT15 retry."""
import argparse,gc,json,os,resource,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.sparse import csr_matrix
import ct_accepted8_native as n
from ct_accepted14_pair import numerical_receipts
from ct_accepted14_gate import equivalence_controls
from ct_accepted15_exec import native_run
from ct_accepted10_diagnose import matrix
from ct_accepted9_pressure_run import sealed
from ct_accepted_pilot import replay
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,sha,runtime_identity

def maximum(v):return float(np.max(abs(v),initial=0.))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('runtime','static','previous-pair','contract','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1'
    s=sealed(a.static);old=sealed(a.previous_pair);contract=n.read(a.contract)
    require(s['status']=='pass' and old['qualification']['reasons']==['temporal:Pi'] and contract['dt']==[.00025,.000125] and contract['T']==.002 and contract['previous_pair_sha256']==sha(a.previous_pair/'results.json') and contract['static_sha256']==sha(a.static/'results.json'),'foreign retry rationale/cohort')
    rt=runtime_identity(a.runtime/'IncNavierStokesSolver');pr=runtime_identity(a.runtime/'RANSPreservationProbe');sr=n.read(a.static/'request-manifest.json')['request'];require(rt==sr['runtime'] and pr==sr['probe'],'foreign affected-gate runtime')
    q=n.read(a.static/'target/frozen.json');require(sha(a.static/'target/frozen.json')==s['native_packet_sha256'],'changed initial static packet');state=q['observation']['state'];off=q['observation']['offsets'];count=len(n.inventory(s['case'])[2]);rr=[];cc=[];vv=[]
    for i,(row,region) in enumerate(zip(q['test_basis'][:count],q['observation']['grid_regions'][:count],strict=True)):
        for j,v in zip(row['pressure']['ids'],row['pressure']['phi'],strict=True):rr.append(i);cc.append(off[4+region]+j-off[4]);vv.append(v)
    B=csr_matrix((vv,(rr,cc)),shape=(count,off[6]-off[4]));prior=n.read(a.static/'canonical-1000/prefactor.json')['pre_factor_frozen'];A0=matrix(prior['A']);M=matrix(prior['free_velocity_mass'])
    m=RequestManifest(a.output,dict(schema='CT15ChangedTimestepAdmissionV1',source=sha(Path(__file__)),contract=sha(a.contract),static=sha(a.static/'results.json'),previous_pair=sha(a.previous_pair/'results.json'),runtime=rt,probe=pr,dt=contract['dt'],alpha=[4000.,8000.],physical_advance=False,classifier_template_source=sha(Path(__file__).with_name('ct_accepted11_target.py')),request_specialization_source=sha(Path(__file__).with_name('ct_accepted14_gate.py'))))
    out=dict(schema='CT15ChangedTimestepAdmissionV1',status='incomplete',dt=contract['dt'],controls={},canonical={},classifier_controls=equivalence_controls(old))
    def execute(folder,request,tag):
        req=folder/(tag+'-request.json');packet=folder/(tag+'.json');write_json(req,request)
        e=native_run([a.runtime/'RANSPreservationProbe',folder/'case.xml',req,packet],folder,tag+'.log',600,900);r=n.read(packet) if packet.exists() else {};require(e['exit_code']==0 and not r.get('failure'),'changed-alpha control:'+str(r.get('failure')));return r,e
    for alpha in (4000.,8000.):
        os.environ.pop('NEKTAR_RANS_SPARSE_RECEIPTS',None)
        dest=a.output/('canonical-'+str(int(alpha)));dest.mkdir();root=ET.parse(a.static/'target/case.xml').getroot();set_parameter(root,'TimeStep',1/alpha);write(root,dest/'case.xml')
        packet,e=execute(dest,replay(q,state,0.,[],stored_reference_sha256=q['initialization_receipt']['fixed_reference_sha256'],export_pre_factor_frozen=True,actual_history=state,sparse_alpha=alpha,frozen_request_kind='InitializedFrozenV1',pre_factor_prefix=str(dest/'frozen')),'prefactor');f=packet['pre_factor_frozen']
        require(f['whole_command_sparse_factorizations']==0 and f['input']==state and f['history']==state and f['rows']==prior['rows'],'changed-alpha canonical identity/history')
        for key in ('G','D','C','free_velocity_mass'):require(f[key]['sha256']==prior[key]['sha256'],'alpha changed canonical '+key)
        defect=matrix(f['A'])-A0-(alpha-1000.)*M;scaled=maximum(defect.data)/max(1.,maximum(A0.data),maximum((alpha-1000.)*M.data));require(scaled<=2e-10,'independent canonical alpha*M relation')
        rhs=dest/'rhs';rhs.write_text('\n'.join(format(-x,'.17g') for x in f['frozen_flow_residual'][:f['velocity_size']])+'\n');groups=dest/'groups';groups.write_text('\n'.join(str(i%f['block_sizes'][0]) for i in f['rows'][:f['velocity_size']])+'\n');result=dest/'ordering.json'
        oe=native_run([a.runtime/'CT13SparseOrderingProbe',f['A']['file'],rhs,'VectorCOLAMDV1',result,groups],dest,'ordering.log',600,900);r=n.read(result) if result.exists() else {};require(oe['exit_code']==0 and r['factor_nnz']<=40000000,'changed-alpha frozen capacity')
        require(sorted(r['old_to_new'])==list(range(f['velocity_size'])),'changed-alpha inverse permutation')
        for prefix in ('','tiny_','cancelling_'):require(r[prefix+'true_residual']<=2e-12*r[prefix+'rhs_maximum']+1e-25,'changed-alpha actual/tiny/cancelling residual')
        out['canonical'][str(int(alpha))]=dict(status='pass',packet_sha256=sha(dest/'prefactor.json'),A=f['A'],mass_relation_scaled_defect=scaled,execution=e,ordering={k:v for k,v in r.items() if k not in ('solution','tiny_solution','old_to_new')},ordering_execution=oe,ordering_result_sha256=sha(result));del packet,f,r;gc.collect()
        base=None
        for label,orders in (('native',None),('q48',(48,12,12,16,64)),('q64',(64,16,12,20,64))):
            folder=a.output/(str(int(alpha))+'-'+label);folder.mkdir();tree=ET.parse(dest/'case.xml').getroot()
            if orders:
                for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(tree,key,value)
            write(tree,folder/'case.xml');os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(folder/'sparse-receipts.jsonl')
            packet,e=execute(folder,replay(q,state,0.,[],export_faces_only=True,export_sparse_backend=True,sparse_alpha=alpha),'frozen');b=packet['sparse_backend'];require(packet['observation']['state']==state and b['pressure_C_nnz']==0,'changed-alpha state/absolute-pressure convention')
            for name in ('predictor','pressure'):
                key='true_pressure_residual' if name=='pressure' else 'predictor_residual';require(b[key]<=2e-12*b[name+'_rhs_maximum']+1e-25,'changed-alpha original residual')
            for v in b['scalars']:require(v['maximum_true_residual']<=2e-12*v['rhs_maximum']+1e-25 and v['factor_nnz']<=40000000,'changed-alpha scalar action')
            require(max(b[k] for k in ('predictor_factor_nnz','consistent_mass_factor_nnz','pressure_preconditioner_factor_nnz'))<=40000000,'changed-alpha produced fill')
            require(max(b[k] for k in ('full_native_G_action_defect','consistent_mass_correction_true_residual','splitting_identity'))<=2e-10,'changed-alpha full G/mass/splitting action')
            if base is None:base=np.asarray(b['pressure_increment'])
            change=maximum(B@(np.asarray(b['pressure_increment'])-base))/s['case']['ustar']**2;require(change<1e-8,'changed-alpha governing quadrature uncertainty')
            row={k:v for k,v in b.items() if k not in ('predictor','pressure_increment','consistent_mass_correction','scalars','identity')};row.update(scalars=[{k:v for k,v in x.items() if k!='diagnostic_increment'} for x in b['scalars']],status='pass',pressure_change_scaled_maximum=change,execution=e,packet_sha256=sha(folder/'frozen.json'),numerical_receipts=numerical_receipts(folder/'sparse-receipts.jsonl'));out['controls'][str(int(alpha))+'-'+label]=row
            os.environ.pop('NEKTAR_RANS_SPARSE_RECEIPTS',None)
            write_json(a.output/'progress.json',out);del packet;gc.collect();print('CT15 changed alpha',alpha,label,'pass',flush=True)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT15 changed timestep admission pass',flush=True)

if __name__=='__main__':main()
