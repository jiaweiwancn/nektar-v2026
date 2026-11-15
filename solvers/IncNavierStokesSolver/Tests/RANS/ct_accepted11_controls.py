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
    for r in q['members'].values():r['actual_stage_checks']=[dict(predictor_factor_nonzeros=100,mass_factor_nonzeros=100,pressure_preconditioner_factor_nonzeros=100)]
    require(target_pair(q)['status']=='pass','target synthetic positive failed');neg={}
    def bad(name,word,fn):
        z=copy.deepcopy(q);fn(z);z['qualification']=dict(status='pass');r=target_pair(z);require(r['status']=='fail' and any(word in s for s in r['reasons']),'target mutation reason:'+name);neg[name]=r
    bad('missing_dt','missing_member',lambda z:z['members'].pop('0.0005'))
    bad('wrong_horizon','foreign_target_request',lambda z:z.update(target=.001))
    bad('incomplete_endpoint','accepted_coverage',lambda z:z['members']['0.001']['accepted_metrics'].pop())
    bad('foreign_budget_identity','budget_state_configuration_identity',lambda z:z['members']['0.001'].update(state_sha256='foreign'))
    bad('wrong_model','foreign_physics',lambda z:z['members']['0.001'].update(signature='A93'))
    bad('pass_label_pressure_failure','temporal:Pi',lambda z:z['temporal']['Pi'].update(difference=dict(weighted_RMS=2e-7,scale=1.)))
    out['target_negatives']=dict(scope='schema fixtures only',results=neg)
    f=a.ct9/'controls-v4/d0-p3-sparse';initial=n.read(f/'operator.json');dest=a.output/'unsupported-policy';dest.mkdir();shutil.copy2(f/'case.xml',dest/'case.xml');state=initial['observation']['state'];req=replay(initial,state,.002,[],export_physical_scalar_budget=True,previous_state=state,physical_budget_volume_policy='UnsupportedObserverV0');r,x=probe(a.runtime,dest,req,'negative');require(x['exit_code']!=0 and 'unsupported physical scalar observer policy' in r.get('failure',''),'wrong native policy rejection')
    out['native_negative']=dict(status='pass',reason=r['failure'],execution=x,request_sha256=sha(dest/'negative-request.json'));out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT11 admission negatives pass',len(out['negative'])+len(neg)+1,flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
