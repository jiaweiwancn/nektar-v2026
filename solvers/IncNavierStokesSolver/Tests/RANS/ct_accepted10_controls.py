#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused CT10 observational parity and deliberately synthetic gate negatives."""
import argparse,copy,os,shutil
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted9_gate import pair
from ct_accepted10_target_run import target_pair
from ct_accepted_pilot import replay
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity


def fixture():
    """Schema unit fixture only, not native scientific evidence."""
    out=dict(schema='CT10TargetPhysicalPairV1',target=.002,required_dt=[.001,.0005],same_actual_initial=True,members={},temporal={},traction_temporal={})
    fields={f:dict(weighted_RMS=0.,scale=1.,scaled=0.) for f in ('U','w','k','outer_epsilon','Pi')}
    for dt in (.001,.0005):
        name=str(dt);binding=dict(state_sha256='synthetic_state',configuration_sha256='synthetic_config',native_report_sha256='synthetic_report',dt=dt,time=.002)
        budget=dict(budget_binding=binding,budgets={},captured_k_matrix=dict(maximum_defect=0.,scale=1.))
        for control in ('native','q48f16','q64f16','q64f32'):
            budget['budgets'][control]=dict(binding=binding,equations={f:dict(physical_scale=1.,storage_backward_difference=0.,stage_source=0.,stage_physical_boundary_flux=0.,frozen_stage_physical_defect=0.) for f in ('inner_k','outer_k','outer_epsilon','whole_k')},reconstruction={f:dict(normalization=1.,maximum_defect=0.,scaled_defect=0.) for f in ('k','outer_epsilon')})
        out['members'][name]=dict(dt=dt,time=.002,failure='',state_sha256='synthetic_state',configuration_sha256='synthetic_config',native_report_sha256='synthetic_report',signature=':resolved=A95:OriginalHighRe-v1 SamePhysicalWallDisplacedAcceptedV1 DisplacedNeutralCompositeV1 AcceptedRichardsNorrisFlatV1 EigenLU-LDLT-IR3-F40M',accepted_metrics=[dict(time=round(i*dt,12),drift=copy.deepcopy(fields),mass=dict(outward_net=0.,physical_scale=1.),traction=dict(normalized_wall_traction_error=0.,normalized_interface_traction_RMS=0.)) for i in range(round(.002/dt)+1)],budgets=budget,initialization_receipt=dict(open_explicit_implicit_agreement=True,explicit_boundary_sha256='synthetic_boundary',implicit_boundary_sha256='synthetic_boundary'),domain=dict(status='pass'),constitutive=dict(status='pass'),corrections=0,extra_solves=0)
    out['temporal']={f:dict(difference=v,measured_drift=0.) for f,v in fields.items()};out['traction_temporal']={f:dict(spread=0.,measured_drift=0.) for f in ('ground','interface','signed_mean')};return out


def gate_controls():
    good=fixture();require(target_pair(good)['status']=='pass','synthetic target schema positive fixture failed');rows={}
    def fault(label,word,fn):
        q=copy.deepcopy(good);fn(q);q['qualification']=dict(status='pass');got=target_pair(q);require(got['status']=='fail' and any(word in s for s in got['reasons']),'wrong negative reason:'+label+str(got));rows[label]=dict(status='pass',expected_rejection=word,observed=got)
    fault('missing_member','missing_member',lambda q:q['members'].pop('0.0005'))
    fault('wrong_horizon','foreign_target_request',lambda q:q.update(target=.008))
    fault('incomplete_accepted_endpoint','accepted_coverage',lambda q:q['members']['0.001']['accepted_metrics'].pop())
    fault('wrong_model','foreign_physics',lambda q:q['members']['0.001'].update(signature=q['members']['0.001']['signature'].replace('A95:','A93:')))
    fault('missing_outer_epsilon','budget_coverage',lambda q:q['members']['0.001']['budgets']['budgets']['native']['equations'].pop('outer_epsilon'))
    fault('missing_control','budget_coverage',lambda q:q['members']['0.001']['budgets']['budgets'].pop('q64f32'))
    fault('foreign_budget_state','budget_state_configuration_identity',lambda q:q['members']['0.001'].update(state_sha256='foreign'))
    fault('nonfinite_primitive','budget_coverage',lambda q:q['members']['0.001']['budgets']['budgets']['native']['equations']['outer_k'].update(stage_source=float('nan')))
    fault('stale_boundary','boundary_identity',lambda q:q['members']['0.001']['initialization_receipt'].update(implicit_boundary_sha256='stale'))
    fault('pass_label_failed_pressure','drift:Pi',lambda q:q['members']['0.001']['accepted_metrics'][-1]['drift']['Pi'].update(weighted_RMS=.002,scaled=0.))
    fault('pass_label_relaxed_temporal_limit','temporal:Pi',lambda q:q['temporal']['Pi'].update(difference=dict(weighted_RMS=2e-7,scale=1.,scaled=0.),limit=1.))
    return dict(scope='synthetic classifier fixtures only; no native or physical qualification claim',positive='pass',negative=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct9','runtime','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');old=a.ct9/'equivalence-v2/zero-0.001';initial=n.read(old/'initial.json')
    m=RequestManifest(a.output,dict(schema='CT10ObservationalControlsV1',initial=sha(old/'initial.json'),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source=sha(Path(__file__)),scope='same native state/reference/physical faces; original dense-export guard preserved; classifier negatives are synthetic'))
    folder=a.output/'faces';folder.mkdir();shutil.copy2(old/'case.xml',folder/'case.xml');req=replay(initial,initial['observation']['state'],0.,[],export_faces_only=True)
    q,e=probe(a.runtime,folder,req,'faces');require(e['exit_code']==0,'face-only export failed');require(q['observation']['state']==initial['observation']['state'] and q['reference']==initial['reference'] and q['fixed_boundary']==initial['fixed_boundary'],'observational probe changed state/reference/boundary')
    require(q['observation']['faces']==initial['observation']['faces'] and q['observation']['volume']==[],'face-only rows differ or volume catalog was serialized')
    differences={}
    for k,v in initial['observation']['audit_alpha_zero'].items():
        other=q['observation']['audit_alpha_zero'].get(k)
        if isinstance(v,list) and all(isinstance(x,(int,float)) for x in v):
            require(len(other)==len(v),'audit length changed');differences[k]=float(max(abs(np.asarray(v)-np.asarray(other)),default=0.));require(differences[k]<=2e-10*max(1.,max(abs(np.asarray(v)),default=0.)),'old/new native audit differs:'+k)
    guard=a.output/'dense-guard';guard.mkdir();shutil.copy2(old/'case.xml',guard/'case.xml');req=replay(initial,initial['observation']['state'],0.,[],export_resolution=True,resolution_localize=True);wrong,we=probe(a.runtime,guard,req,'negative')
    require(we['exit_code']!=0 and '1000-slot' in wrong.get('failure',''),'original dense per-element guard was bypassed')
    gates=gate_controls();result=dict(schema='CT10ObservationalControlsV1',status='pass',face_count=len(q['observation']['faces']),state_sha256=exact_vector_hash(q['observation']['state']),face_parity='bitwise JSON row/weight/owner/side equality with sealed CT9',audit_defects=differences,execution=e,packet_sha256=sha(folder/'faces.json'),unchanged_dense_export_guard=dict(status='pass',reason=wrong['failure'],execution=we),gate_controls=gates)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT10 observational and gate controls pass',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
