#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only reconciliation of sealed conservative scalar stages and physical fluxes."""
import argparse,os,shutil,xml.etree.ElementTree as ET
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from qualify_rans_preservation import load,replay_request,Q
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import write_json,exact_vector_hash
from rans_duration_contract import MEMBERS,require
from diagnose_rans_balanced_spatial import run_probe
from run_cls_epawf_picard import set_parameter,write


def classify(budget):
    r=budget['reconciliation'];require(r['schema']=='ConservativeFrozenScalarReconciliationV1','scalar_reconciliation_schema')
    checks={}
    for name,row in r['assembled'].items():
        defect=float(np.max(abs(np.array(row['reconstruction_defect'])),initial=0.))
        raw=np.array(row['native_raw_frozen_load']);free=np.array(row['free_rows'],int);fixed=np.array(row['prescribed_or_dependent_rows'],int)
        history=float(np.max(abs(np.array(row['previous_minus_stage_history'])),initial=0.))
        require(defect<=2e-11 and history==0.,'frozen_scalar_load_or_actual_storage_mismatch:'+name)
        maximum=float(np.max(abs(raw[free]),initial=0.));require(maximum<=1e-10,'actual_frozen_free_row_not_solved:'+name)
        checks[name]=dict(term_sum_maximum_defect=defect,free_row_maximum=maximum,
            prescribed_row_maximum=float(np.max(abs(raw[fixed]),initial=0.)),previous_stage_history_difference=history,
            free_count=len(free),prescribed_count=len(fixed),interpretation='prescribed reactions are algebraic, not independent physical fluxes')
    for name,row in r['regions'].items():
        adv=float(np.max(abs(np.array(row['advective_IBP_defect'])),initial=0.))
        require(adv<=2e-11,'conservative_advection_IBP_failed:'+name)
        diff=None if row['diffusion_IBP_defect'] is None else float(np.max(abs(np.array(row['diffusion_IBP_defect'])),initial=0.))
        if diff is not None:require(diff<=2e-10,'physical_diffusion_IBP_failed:'+name)
        checks[name+'-IBP']=dict(advective_maximum=adv,advective_limit=2e-11,diffusive_maximum=diff,diffusive_limit=2e-10)
    return checks


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--ownership-root',type=Path,required=True)
    p.add_argument('--runtime',type=Path,required=True);p.add_argument('--output-label',default='frozen-scalar-reconciliation-v1');a=p.parse_args()
    root=a.root.resolve();old=a.ownership_root.resolve();runtime=a.runtime.resolve();out=root/a.output_label
    os.environ.update(LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(root/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
    sources={}
    for kind in ('baseline','perturbed'):
        cohort='A95-N-Accepted-'+('duration-0.002' if kind=='baseline' else 'response-2e-05-dt-0.002')
        for member in MEMBERS:
            folder=old/'campaigns'/(cohort+'-'+member);request=load(folder/'request-manifest.json')
            endpoint=folder/'T-0.0160';require(RequestManifest(folder,request['request'],True).completed(endpoint),'reconciliation_unsealed_source')
            sources[kind+'-'+member]=dict(folder=str(folder),request_sha256=sha(folder/'request-manifest.json'),report_sha256=sha(endpoint/'case.evolution.json'))
    intent=dict(schema='FrozenScalarReconciliationRequestV1',configuration=dict(model='A95',policy='Natural',wall='Accepted'),
        time=.016,sources=sources,runtime=runtime_identity(runtime/'RANSPreservationProbe'),volume_orders=[48,64],face_orders=[16,32],
        quadrature_reconciliation_members=['baseline-B5','perturbed-B5'],
        independent_diffusion_IBP=dict(element_step_fractions=[2**-11,2**-12],members=['baseline-B5','perturbed-B5'],
            method='fourth-order physical flux differences within the generating owner; h/h2 uncertainty, no constitutive tangent'),
        limits=dict(native_row_reconstruction=2e-11,actual_frozen_free_rows=1e-10,advective_IBP=2e-11,diffusive_IBP=2e-10),
        scope='reconciliation only; unchanged historical physical-budget accuracy, no scalar*div(U) correction or flux replacement')
    seal=RequestManifest(out,intent);checks={};files=[]
    try:
        for key,source in sources.items():
            folder=Path(source['folder']);n=load(folder/'T-0.0160/case.evolution.json');req=load(folder/'request-manifest.json')['request']
            previous=n['accepted_states'][-2];require(exact_vector_hash(previous)==n['pimple_attempts'][-1]['incoming_sha256'],'reconciliation_incoming_state_differs')
            history=n['last_stage_history'];require(len(history)==len(n['state']),'missing_actual_coupled_stage_history')
            trials=[(2,32,0.)]
            if key in ('baseline-B5','perturbed-B5'):trials += [(1,32,0.),(2,16,0.),(2,32,2**-11),(2,32,2**-12)]
            for q,face,h in trials:
                name=key+'-Q'+str(q)+'-F'+str(face)+'-h'+str(h);dest=out/name;dest.mkdir()
                shutil.copy2(folder/'T-0.0160/case.xml',dest/'case.xml');xml=ET.parse(dest/'case.xml').getroot()
                for setting,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),Q[q],strict=True):set_parameter(xml,setting,value)
                write(xml,dest/'case.xml');seed=load(folder/'initial-packet.json')
                request=replay_request(seed,[],state=n['state'],time=n['time'])
                request.update(fixed_boundary=n['initialization_receipt']['live_fixed_boundary'],export_physical_scalar_budget=True,
                    export_frozen_scalar_reconciliation=True,physical_face_order=face,previous_state=previous,
                    stage_history=history,diffusion_IBP_element_step_fraction=h)
                packet,_=run_probe(runtime,dest,request);budget=packet['physical_scalar_budget']
                require(exact_vector_hash(packet['observation']['state'])==exact_vector_hash(n['state']),'reconciliation_modified_state')
                c=classify(budget);checks[name]=dict(status='pass',checks=c,source=source,budget=budget['equations'],packet_sha256=sha(dest/'packet.json.gz'))
                files += [dest/'case.xml',dest/'request.json',dest/'packet.json.gz',dest/'probe.log']
                write_json(out/'progress.json',dict(status='running',completed=list(checks)))
        # The h/h2 test is separately reported; it never modifies budget values.
        diffusion=[]
        for key in ('baseline-B5','perturbed-B5'):
            rows=[checks[key+'-Q2-F32-h'+str(h)]['checks'] for h in (2**-11,2**-12)]
            for region in ('inner_k','outer_k','outer_epsilon'):
                v=[r[region+'-IBP']['diffusive_maximum'] for r in rows]
                diffusion.append(dict(member=key,region=region,maximum_defects=v,limit=2e-10,status='pass' if max(v)<=2e-10 else 'unresolved'))
        result=dict(status='pass',reconstruction_status='pass',checks=checks,diffusion_IBP=diffusion,
            historical_budget_accuracy_status='failed_or_unresolved',accuracy_reclassified=False,
            scope='actual frozen matrix and conservative physical flux reconciliation; not a physical-budget accuracy qualification')
    except Exception as e:result=dict(status='fail',failure=repr(e),checks=checks,historical_budget_accuracy_status='failed_or_unresolved',accuracy_reclassified=False)
    write_json(out/'results.json',result)
    files=[p for p in out.rglob('*') if p.is_file() and p.name not in ('completed-artifacts.json','request-manifest.json')]
    seal.seal(out,files);print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
