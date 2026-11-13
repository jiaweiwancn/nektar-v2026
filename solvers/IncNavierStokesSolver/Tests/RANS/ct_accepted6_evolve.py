#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One immutable CT6 member, fresh history, actual accepted-state MMS gates.

Run only through the owned command launcher. A failed trajectory is sealed;
analysis may continue read-only, but no failure is reset or resumed here.
"""
import argparse,json,os,traceback
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
from ct_accepted6_cases import make,request
from ct_accepted6_metrics import metrics
from ct_accepted6_operator import reconstruct,continuous_weak,represented_consistency
from ct_accepted_open import run,probe,require
from ct_accepted_pilot import replay,artifact_hash
from rans_request_manifest import RequestManifest,runtime_identity,sha


def report_summary(path,execution):
    r=json.loads(path.read_text())
    keys=('failure','time','dt','steps','audit_times','state_audits','physical_range_events',
          'stage_audits','step_audits','explicit_audits','initialization_receipt',
          'native_history_times','last_stage_time','last_stage_alpha','pimple_attempts',
          'roughness_budget','energy_budget','signature','pilot_accuracy_semantics')
    result={k:r[k] for k in keys if k in r}
    result.update(execution=execution,report_sha256=sha(path),state_sha256=exact_vector_hash(r['state']),
        initial_sha256=exact_vector_hash(r['initial']),checkpoint=r['checkpoint'],
        history_sha256=artifact_hash(r['checkpoint']) if r['checkpoint'] else None,stored_slots=len(r['state']))
    return r,result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('runtime','contract','admission','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--model',choices=['A95','A93','E94E95'],default='A95')
    p.add_argument('--degree',type=int,choices=[3,4],default=4)
    # CT7 authorizes one separately contracted finer member. Old CT6
    # requests, histories and failed consistency receipts remain immutable.
    p.add_argument('--dt',type=float,choices=[.004,.002,.001,.0005],default=.002)
    p.add_argument('--target',type=float,default=.016)
    p.add_argument('--prepared',type=Path,help='explicit reuse of a sealed zero-step preparation, never an evolved state')
    a=p.parse_args();require(0<a.target<=.016,'CT6 bounded horizon exceeded')
    admission=json.loads(a.admission.read_text())
    require(admission['status']=='pass','required native operator/geometry/zero-flat admission failed')
    steps=round(a.target/a.dt);require(abs(steps*a.dt-a.target)<1e-14,'nonintegral time horizon')
    prepared=None
    if a.prepared:
        old=json.loads((a.prepared/'results.json').read_text());seal=json.loads((a.prepared/'completed-artifacts.json').read_text())
        require(all(sha(a.prepared/n)==v for n,v in seal['files'].items()),'changed sealed preparation')
        require(old['model']==a.model and old['degree']==a.degree and old['dt']==a.dt,'foreign zero-step preparation')
        require('trajectory' not in old,'preparation already advanced a trajectory')
        zero=a.prepared/'zero';zr=json.loads((zero/'case.evolution.json').read_text())
        require(zr['time']==0. and zr['steps']==0 and not zr['failure'],'prepared history is not a successful zero-step state')
        require(sha(zero/'case.evolution.json')==old['zero']['report_sha256'] and artifact_hash(zr['checkpoint'])==old['zero']['history_sha256'],'changed zero-step state/history')
        require(old['initial_metrics']['status']=='pass' and sha(zero/'initial.json')==old['initial_packet']['sha256'],'changed initial packet or failed initial accuracy')
        prior_request=json.loads((a.prepared/'request-manifest.json').read_text())['request']
        require(prior_request['runtime']['IncNavierStokesSolver']==runtime_identity(a.runtime/'IncNavierStokesSolver'),'foreign prepared runtime')
        prepared=dict(directory=str(a.prepared),receipt_sha256=sha(a.prepared/'results.json'),state=old['zero']['state_sha256'],history=old['zero']['history_sha256'])
    m=RequestManifest(a.output,dict(schema='CT6FreshManufacturedMemberV1',model=a.model,degree=a.degree,
        dt=a.dt,target=a.target,contract_sha256=sha(a.contract),admission_sha256=sha(a.admission),
        runtime={n:runtime_identity(a.runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},
        scripts={n:sha(Path(__file__).with_name(n)) for n in ('ct_accepted6_evolve.py','ct_accepted6_cases.py','ct_accepted6_metrics.py','ct_accepted6_operator.py')},
        actual_initialization='fresh independently projected at zero, zero-step native history; same-space dt vectors checked separately',
        pressure='absolute Pi, no fitted gauge',forcing='continuous CT6 MMS only',explicit_preparation_reuse=prepared))
    result=dict(schema='CT6FreshManufacturedMemberV1',model=a.model,degree=a.degree,dt=a.dt,target=a.target,
                execution_status='incomplete',physical_completion='unattempted',diagnostic_status='unattempted')
    try:
        if a.prepared:
            r=zr;result['zero']=old['zero'];initial=json.loads((zero/'initial.json').read_text())
        else:
            zero=a.output/'zero';make(zero,a.runtime,model=a.model,degree=a.degree,dt=a.dt)
            os.environ['NEKTAR_RANS_TIMING_FILE']=str(zero/'timing.jsonl')
            e=run([a.runtime/'IncNavierStokesSolver','case.xml'],zero,'solver.log')
            r,z=report_summary(zero/'case.evolution.json',e);result['zero']=z
            write_json(a.output/'progress.json',result)
            require(e['exit_code']==0 and not r['failure'] and r['time']==0.,'fresh zero-step initialization rejected')
            require(len(r['state'])==(1273 if a.degree==3 else 2283),'actual native slots differ')
            req=request();req.update(export_tensor_validation=False,export_constitutive=False,export_resolution=False)
            initial,pe=probe(a.runtime,zero,req,'initial')
            result['initial_probe_execution']=pe
            require(pe['exit_code']==0 and initial['observation']['state']==r['state'],'initial effective driver/probe mismatch')
        result['initial_packet']=dict(path=str(zero/'initial.json'),sha256=sha(zero/'initial.json'))
        result['initial_metrics']=metrics(initial,a.model)
        require(result['initial_metrics']['status']=='pass','initial physical accuracy gate failed')
        result['zero_directory']=str(zero)
        folder=a.output/'trajectory';make(folder,a.runtime,model=a.model,degree=a.degree,dt=a.dt,steps=steps,restart=Path(r['checkpoint']))
        os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl')
        print(a.model,a.degree,a.dt,'starting fresh trajectory',flush=True)
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log')
        native,summary=report_summary(folder/'case.evolution.json',e);result['trajectory']=summary
        result['execution_status']='pass' if e['exit_code']==0 else 'fail'
        completed=e['exit_code']==0 and not native['failure'] and abs(native['time']-a.target)<1e-13
        if completed:require(bool(native['checkpoint']),'completed evolution did not publish a native history')
        if not native['checkpoint']:
            require(native['time']==0. and native['state']==r['state'],'unpublished accepted history requires a separate failure witness')
            summary['last_accepted_checkpoint']=r['checkpoint']
            summary['last_accepted_history_sha256']=artifact_hash(r['checkpoint'])
        result['physical_completion']='pass' if completed else 'fail'
        write_json(a.output/'progress.json',result)
        # A last accepted state can always be audited without erasing a failed candidate.
        req=replay(initial,native['state'],native['time'],initial['observation']['grid'],
            grid_sides=initial['observation']['grid_sides'],export_basis=True,export_ct6_primitives=True,
            export_manufactured=True,export_resolution=True,export_constitutive=True)
        end,ee=probe(a.runtime,folder,req,'endpoint');result['endpoint_probe_execution']=ee
        require(ee['exit_code']==0,'last accepted state diagnostic failed')
        result.update(endpoint_packet=dict(path=str(folder/'endpoint.json'),sha256=sha(folder/'endpoint.json')),
            endpoint_metrics=metrics(end,a.model,initial['observation']['state']),
            all_row=reconstruct(end,a.model),continuous_weak=continuous_weak(end,a.model),
            represented_consistency=represented_consistency(end,a.model))
        result['diagnostic_status']='pass' if all(result[k]['status']=='pass' for k in ('all_row','continuous_weak')) else 'fail'
    except Exception as error:
        result.update(error=repr(error));traceback.print_exc()
    finally:os.environ.pop('NEKTAR_RANS_TIMING_FILE',None)
    result['status']='pass' if all(result[k]=='pass' for k in ('execution_status','physical_completion','diagnostic_status')) else 'fail'
    write_json(a.output/'results.json',result)
    m.seal(a.output,['results.json']+[str(x.relative_to(a.output)) for x in a.output.rglob('case.xml')])
    print('CT6 member',result['status'],result.get('error',''),flush=True)
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
