#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Required CT11 primitive negatives plus native unsupported-observer rejection."""
import argparse,copy,os,shutil
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted11_gate import observer,capacity,prerequisites
from ct_accepted11_target import target_pair
from ct_accepted10_controls import fixture
from ct_accepted_pilot import replay
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct9','runtime','equivalence','capacity','preflight','observer-resources','capacity-resources','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');e=sealed(a.equivalence);c=sealed(a.capacity);pre=sealed(a.preflight);er=n.read(a.equivalence/'request-manifest.json')['request'];cr=n.read(a.capacity/'request-manifest.json')['request'];resources={k:n.read(v/'execution.json') for k,v in [('equivalence',a.observer_resources),('capacity',a.capacity_resources)]}
    required=dict(probe_sha=sha(a.runtime/'RANSPreservationProbe'),case=pre['mesh'],state=pre['actual_initial_sha256'],equivalence_sha=sha(a.equivalence/'results.json'))
    good=[e,c,er,cr,resources]
    def admit(v):return prerequisites(*v[:4],**required,resources=v[4])
    require(admit(good)['status']=='pass','real prerequisite positive failed')
    m=RequestManifest(a.output,dict(schema='CT11NegativeControlsV1',equivalence=sha(a.equivalence/'results.json'),capacity=sha(a.capacity/'results.json'),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source=sha(Path(__file__)),scope='deliberately mutated classifier fixtures; not trajectories; one native unsupported-policy rejection'))
    out=dict(schema='CT11NegativeControlsV1',status='incomplete',positive=admit(good),negative={})
    def fault(name,word,fn):
        v=copy.deepcopy(good);fn(v);v[0]['status']=v[1]['status']='pass';r=admit(v);require(r['status']=='fail' and any(word in s for s in r['reasons']),'nondiscriminating fixture:'+name+str(r));out['negative'][name]=dict(status='pass',expected=word,observed=r)
    fault('missing_observer_control','observer_equivalence',lambda v:v[0]['checks']['d0-p3'].pop('q64f32'))
    fault('missing_fixture','observer_equivalence',lambda v:v[0]['checks'].pop('displaced-p4'))
    fault('missing_equation','observer_equivalence',lambda v:v[0]['checks']['d0-p3']['native']['classification']['equations'].pop('outer_epsilon'))
    fault('nonfinite_primitive','observer_equivalence',lambda v:v[0]['checks']['d0-p3']['native']['all_primitive_difference'].update(maximum_scaled_difference=float('nan')))
    fault('pass_label_bad_row','observer_equivalence',lambda v:v[0]['checks']['d0-p3']['native']['all_primitive_difference'].update(maximum_scaled_difference=1e-7))
    fault('wrong_owner_digest','observer_equivalence',lambda v:v[0]['checks']['d0-p3']['native']['native_provenance']['stream'].update(native_rows_provenance_fnv1a64='foreign'))
    fault('missing_capacity_control','observer_capacity',lambda v:v[1]['controls'].pop('q48f16'))
    fault('missing_reconstruction','observer_capacity',lambda v:v[1]['controls']['native']['frozen_load_reconstruction'].pop('outer_epsilon'))
    fault('nonfinite_capacity_primitive','observer_capacity',lambda v:v[1]['controls']['native']['equations']['outer_k'].update(stage_source=float('inf')))
    fault('pass_label_bad_reconstruction','observer_capacity',lambda v:v[1]['controls']['native']['frozen_load_reconstruction']['k'].update(maximum_defect=1.,normalization=1.,status='pass'))
    fault('wrong_capacity_state','foreign_capacity_input',lambda v:v[3].update(actual_initial_sha256='foreign'))
    fault('foreign_control_state','observer_capacity',lambda v:v[1]['controls']['native'].update(actual_initial_sha256='foreign'))
    fault('foreign_probe','foreign_observer_runtime',lambda v:v[3]['runtime'].update(sha256='foreign'))
    fault('stale_equivalence','stale_equivalence_receipt',lambda v:v[3].update(equivalence='stale'))
    fault('failed_shared_resources','closed_resource_failure',lambda v:v[4]['capacity']['resource'].update(project_status='fail'))
    q=fixture();q['schema']='CT11TargetPhysicalPairV1'
    for r in q['members'].values():
        r['actual_stage_checks']=[dict(predictor_factor_nonzeros=100,mass_factor_nonzeros=100,pressure_preconditioner_factor_nonzeros=100)]
        for b in r['budgets']['budgets'].values():b['volume_observer']=dict(policy='NativeOwnerVolumeStreamV1',materialized_global_volume_rows=False)
    require(target_pair(q)['status']=='pass','target synthetic positive failed');neg={}
    def bad(name,word,fn):
        z=copy.deepcopy(q);fn(z);z['qualification']=dict(status='pass');r=target_pair(z);require(r['status']=='fail' and any(word in s for s in r['reasons']),'target mutation reason:'+name);neg[name]=r
    bad('missing_dt','missing_member',lambda z:z['members'].pop('0.0005'))
    bad('wrong_horizon','foreign_target_request',lambda z:z.update(target=.001))
    bad('incomplete_endpoint','accepted_coverage',lambda z:z['members']['0.001']['accepted_metrics'].pop())
    bad('foreign_budget_identity','budget_state_configuration_identity',lambda z:z['members']['0.001'].update(state_sha256='foreign'))
    bad('wrong_model','foreign_physics',lambda z:z['members']['0.001'].update(signature='A93'))
    bad('pass_label_pressure_failure','temporal:Pi',lambda z:z['temporal']['Pi'].update(difference=dict(weighted_RMS=2e-7,scale=1.)))
    bad('stale_observer_policy','observer_policy',lambda z:z['members']['0.001']['budgets']['budgets']['native']['volume_observer'].update(policy='stale'))
    out['target_negatives']=dict(scope='schema fixtures only',results=neg)
    post={};trace_storage_signal=0.
    for name,rows in e['checks'].items():
        post[name]={}
        for label,row in rows.items():
            file=a.equivalence/name/label/'stream.json';require(sha(file)==row['packet_sha256'],'changed streamed audit packet');packet=n.read(file);rec=packet['physical_scalar_budget']['reconciliation'];regions={}
            for region,v in rec['regions'].items():
                maximum=max(abs(x) for x in v['advective_IBP_defect']);normal=max(1.,*(abs(x) for k in ('strong_conservative_advection','conservative_advection','all_element_advective_faces') for x in v[k]));require(maximum/normal<=2e-10,'conservative IBP failed:'+name+':'+label+':'+region);regions[region]=dict(maximum_defect=maximum,normalization=normal,scaled_defect=maximum/normal)
            prescribed=rec['assembled']['outer_epsilon']['prescribed_or_dependent_rows'];storage=rec['regions']['outer_epsilon']['backward_storage'];signal=max((abs(storage[j]) for j in prescribed),default=0.)
            if name.startswith(('baseline','candidate')):trace_storage_signal=max(trace_storage_signal,signal)
            post[name][label]=dict(advective_IBP=regions,prescribed_epsilon_backward_storage_maximum=signal,coefficient_policy=rec['coefficient_policy'],old_boundary_lift=rec['current_time_essential_lift_difference'],actual_flow_stage_difference=rec['actual_flow_stage_maximum_difference'],actual_pre_k_stage_difference=rec['actual_pre_k_stage_maximum_difference'],packet_sha256=sha(file))
    require(trace_storage_signal>0.,'actual strip dependent-boundary storage control has no signal');out['postchecks']=dict(status='pass',cases=post,actual_prescribed_epsilon_storage_signal=trace_storage_signal,scope='actual native weak conservative IBP and retained prescribed/dependent epsilon storage; no scalar*div(U) correction, physical flux replacement or extra native solve')
    f=a.ct9/'controls-v4/d0-p3-sparse';initial=n.read(f/'operator.json');dest=a.output/'unsupported-policy';dest.mkdir();shutil.copy2(f/'case.xml',dest/'case.xml');state=initial['observation']['state'];req=replay(initial,state,.002,[],export_physical_scalar_budget=True,previous_state=state,physical_budget_volume_policy='UnsupportedObserverV0');r,x=probe(a.runtime,dest,req,'negative');require(x['exit_code']!=0 and 'unsupported physical scalar observer policy' in r.get('failure',''),'wrong native policy rejection')
    out['native_negative']=dict(status='pass',reason=r['failure'],execution=x,request_sha256=sha(dest/'negative-request.json'));out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT11 admission negatives pass',len(out['negative'])+len(neg)+1,flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
