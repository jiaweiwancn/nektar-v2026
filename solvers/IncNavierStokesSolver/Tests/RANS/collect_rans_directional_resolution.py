#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete compact export; preserve scientific failures and return nonzero."""
import argparse,gzip,hashlib,json,math,os,shutil,subprocess,time
from pathlib import Path
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest,sha
from rans_qualification import write_json,exact_vector_hash
from rans_directional_scope import classify,SCHEMA
from rans_ownership_contract import budget_gate
from rans_duration_contract import MEMBERS,require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--old-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();old=a.old_root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();bodies={};index={};native=[];costs=[];manifests={}
    excluded={'source','build','dependency','runtime','runtime-resolution-v3','runtime-budget-v5','runtime-export-check-v2','tmp'}
    def file_sha(path):
        h=hashlib.sha256()
        with path.open('rb') as f:
            for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
        return h.hexdigest()
    for base,dirs,files in os.walk(r,followlinks=False):
        folder=Path(base);dirs[:]=[d for d in dirs if d not in excluded and not d.startswith(('a93-python-support','runtime-')) and not (folder/d).is_symlink() and folder/d!=out]
        for name in files:
            path=folder/name;relative=str(path.relative_to(r))
            if path.is_symlink():continue
            if name.endswith(('.json','.json.gz','.jsonl','.xml','.opt','.bin','.npz','.py','.log')) or name in ('SHA256SUMS',):
                index[relative]=dict(path=str(path),sha256=file_sha(path),bytes=path.stat().st_size)
            if name=='results.json' or name.endswith('-execution.json') or name in ('metrics.json','initial-metrics.json','initialization-only.json','screen.json','runtime-identity.json','export-check-runtime-identity.json','neutral-comparison.json'):
                bodies[relative]=load(path)
            if name in ('request-manifest.json','completed-artifacts.json','admission.json','hypothesis.json') or name.endswith(('-hypothesis.json','-mesh.json')):
                manifests[relative]=load(path)
            if name=='case.evolution.json':
                d=load(path);states=d.get('accepted_states',[]);times=d.get('audit_times',[])
                row={k:v for k,v in d.items() if k not in ('state','accepted_states','native_history','last_stage_history','rejected_state','initial','raw_reference','reference')}
                row.update(source=str(path),source_sha256=index[relative]['sha256'],
                    accepted_state_identities=[dict(time=t,sha256=exact_vector_hash(v)) for t,v in zip(times,states,strict=True)])
                native.append(row)
            if name=='probe-process-cost.json' or name.endswith('-execution.json'):costs.append(dict(source=str(path),value=load(path)))
    def result(name):return bodies.get(name,dict(status='not_attempted'))
    def gate(name,key='status'):
        d=result(name);status=d.get(key,'not_attempted')
        return dict(status=status if status in ('pass','fail','failed_or_unresolved','unresolved','not_attempted','blocked') else 'failed_or_unresolved',
            source=name,sha256=index.get(name,{}).get('sha256'))
    capability=['directional-layout-v2/results.json','directional-export-layout-v1/results.json']
    gates=dict(directional_capability=dict(status='pass' if all(result(x)['status']=='pass' for x in capability) else 'failed_or_unresolved',sources=capability),
        isotropic_compatibility=gate('isotropic-compatibility-v1/results.json'),
        frozen_scalar_reconciliation=gate('frozen-scalar-reconciliation-v1/results.json'),geometry=gate('directional-geometry-v1/results.json'),
        history_cache_CLS=gate('directional-history-v3/results.json'),required_short_verification=gate('scope-negative-execution.json'),
        A93_mesh_candidate=gate('a93-redistribution-v2/results.json'))
    regression=result('compatibility-v1-execution.json')
    if len(regression.get('events',[]))!=4 or any(regression['events'][i]['exit_code'] for i in (0,2,3)):
        gates['required_short_verification']['status']='failed_or_unresolved'
    scope_ctest=result('scope-ctest-execution.json')
    if scope_ctest.get('status')!='pass':gates['required_short_verification']['status']='failed_or_unresolved'
    cohort='A95-N-Accepted-directional-dt-0.002';experiment=r/'directional-resolution-v1'
    analysis_root=r/'directional-resolution-supplement-v1' if (r/'directional-resolution-supplement-v1').exists() else experiment
    analysis_label=str(analysis_root.relative_to(r))
    gates['A95_neutral_control']=dict(status='pass' if result('directional-resolution-v1/neutral-comparison.json')['status']=='pass' and
        result('response-audits/'+cohort+'-neutral-0.016-canonical-v1/results.json')['status']=='pass' else 'failed_or_unresolved')
    gates['A95_response_fields']=gate(analysis_label+'/response-comparison/results.json')
    response=result(analysis_label+'/response-comparison/results.json')
    if 'fields' in response:
        expected_fields={'u','w','k','epsilon_outer','Pi','pi_mechanical','Tzz','Txz','Rxx','Ryy','Rzz','Rxz','Qxx','Qyy','Qzz','Qxz','P'}
        require(set(response['fields'])==expected_fields,'collector_missing_response_observable')
        for name,row in response['fields'].items():
            require(len(row['response_RMS'])==4 and len(row['timestep_differences'])==2,'collector_response_companion_count')
            values=[*row['response_RMS'],row['spatial_difference'],*row['timestep_differences']]
            require(all(math.isfinite(v) and v>=0 for v in values),'collector_nonfinite_response')
            amplitude=max(row['response_RMS']);h=max(5e-7,.1*amplitude);t=max(1e-7,.05*amplitude)
            require(row['spatial_limit']==h and row['temporal_limit']==t,'collector_response_limits_changed')
            require(row['spatial_status']==('pass' if row['spatial_difference']<=h else 'unresolved') and
                    row['temporal_status']==('pass' if max(row['timestep_differences'])<=t else 'unresolved'),
                    'collector_response_pass_label_disagrees_with_primitives:'+name)
    if gates['A95_response_fields']['status']=='pass' and not response.get('core_pressure_velocity_nonlinear_production_response_resolved',False):
        gates['A95_response_fields']['status']='unresolved'
    gates['A95_response_native_weak']=gate('response-audits/'+cohort+'-response-0.016/results.json')
    binding=result(analysis_label+'/budget-binding/results.json')
    if 'binding' in binding:
        b=binding['binding'];analysis=result(analysis_label+'/budget-analysis/results.json');accuracy=result(analysis_label+'/budget-accuracy/results.json')
        members={m:{} for m in MEMBERS}
        for member in MEMBERS:
            for kind,suffix in (('baseline','neutral'),('perturbed','response')):
                folder=r/'campaigns'/(cohort+'-'+suffix+'-'+member);payload=load(folder/'request-manifest.json')
                endpoint=folder/'T-0.0160';require(RequestManifest(folder,payload['request'],True).completed(endpoint),'collector_unsealed_budget_native')
                members[member][kind]=dict(trajectory_request_sha256=payload['sha256'],source_report_sha256=sha(endpoint/'case.evolution.json'))
                require(b['member_spaces'][member][kind]==payload['request']['spatial_control'],'collector_foreign_budget_space')
        expected=dict(configuration=load(experiment/'admission.json')['configuration'],cohort=cohort+'-response',time=.016,members=members,
            analysis_sha256=sha(analysis_root/'budget-analysis/results.json'),accuracy_sha256=sha(analysis_root/'budget-accuracy/results.json'))
        recalculated=budget_gate(analysis,accuracy,b,expected)
        if recalculated['checks']!=binding['checks']:raise ValueError('exported_budget_binding_primitive_mismatch')
        gates['physical_budget_reconstruction']=dict(status=recalculated['reconstruction_status'])
        gates['physical_budget_accuracy']=dict(status=recalculated['accuracy_status'],primitive_status=recalculated['accuracy_status'],
            binding_sha256=sha(analysis_root/'budget-binding/results.json'))
    else:
        gates['physical_budget_reconstruction']=dict(status='not_attempted');gates['physical_budget_accuracy']=dict(status='not_attempted')
    scope=dict(schema=SCHEMA,numerical_operator='GeneratedVolumeOwnerV1',layout=dict(fields=[6,5],pressure=[4,3]),unknown_limit=3000,
        response_max_time=.016,A93_response_out_of_scope=True,gates=gates,historical_status_unchanged=True)
    decision=classify(scope);write_json(out/'required-scope-input.json',scope);write_json(out/'decision.json',decision)
    with gzip.open(out/'complete-results.json.gz','wt') as f:json.dump(bodies,f,allow_nan=False)
    with gzip.open(out/'native-state-stage-receipts.json.gz','wt') as f:json.dump(native,f,allow_nan=False)
    with gzip.open(out/'request-and-completion-manifests.json.gz','wt') as f:json.dump(manifests,f,allow_nan=False)
    write_json(out/'artifact-index.json',index);write_json(out/'process-cost-receipts.json',costs)
    for name in ('runtime-identity.json','export-check-runtime-identity.json','face-probe-runtime-identity.json',
        'worker-guard-status.json','worker-bookings.jsonl','worker-guard-events.jsonl',
        'actual-owned-runtime-maps.json','actual-owned-runtime-maps-supplement.json','actual-owned-runtime-maps-faces.json','scope-negative-execution.json'):
        path=r/name
        if path.exists():
            with gzip.open(out/(name+'.gz'),'wb') as f:f.write(path.read_bytes())
    controllers=out/'controllers';controllers.mkdir()
    for path in r.glob('*controller*.py'):shutil.copy2(path,controllers/path.name)
    for path in r.glob('*monitor*.py'):shutil.copy2(path,controllers/path.name)
    # Minimal failed inputs and actual unchanged native histories. Large volume
    # packets and all trajectories remain at their hashed Linux destinations.
    witness=out/'failure-reproducers';witness.mkdir()
    for label in ('directional-history-v1','directional-history-v2','directional-history-v3'):
        source=r/label/'changed-axis-points-1'
        if source.exists():
            dest=witness/label;dest.mkdir()
            for name in ('case.xml','case.opt','case-m4-native-rank0.opt','solver.log'):
                if (source/name).exists():shutil.copy2(source/name,dest/name)
    for kind in ('baseline','perturbed'):
        for member in ('B5','B11','B5-dt2','B11-dt2'):
            leaf=r/'physical-budgets'/(cohort+'-'+kind+'-'+member+'-0.016-Q2-F32')
            if leaf.exists():
                dest=witness/(kind+'-'+member);dest.mkdir();shutil.copy2(leaf/'case.xml',dest/'case.xml')
                with gzip.open(dest/'request.json.gz','wb') as f:f.write((leaf/'request.json').read_bytes())
    for path in (r/'a93-redistribution-v2/campaigns').glob('*B11*/T-0.0160/case.evolution.json'):
        report=load(path);dest=witness/path.parent.parent.name;dest.mkdir()
        shutil.copy2(path.parent/'case.xml',dest/'case.xml');checkpoint=Path(report['checkpoint'])/'State.bin';shutil.copy2(checkpoint,dest/'State.bin')
        write_json(dest/'identity.json',dict(source_report_sha256=sha(path),source_report=str(path),checkpoint_sha256=sha(checkpoint),
            failure='unchanged fieldwise epsilon spatial reduction; not a native physical step rejection'))
    reused=dict(old_root=str(old),historical_decision_sha256=sha(old/'source/docs/plans/a93-a95-e94e95-awf-native-ownership-requalification-evidence/decision.json'),
        old_operator_evidence_reused_as_changed_operator=False,
        scope='unchanged kernel source identities, pinned A95 .5 admission and exact sealed A93 B5 trajectories; all A95 directional experiment histories are new')
    write_json(out/'reused-evidence.json',reused)
    write_json(out/'collection.json',dict(status='completed',scientific_status=decision['status'],complete_result_bodies=len(bodies),
        native_reports=len(native),indexed_artifacts=len(index),wall_seconds=time.monotonic()-start,
        report_source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=r/'source',text=True).strip(),
        scope='complete numerical bodies and accepted/stage receipts; heavy native fields and volume matrices retained by hash, no retrospective threshold change'))
    sums=[]
    for path in sorted(out.rglob('*')):
        if path.is_file():sums.append(file_sha(path)+'  '+str(path.relative_to(out)))
    (out/'SHA256SUMS').write_text('\n'.join(sums)+'\n');print(decision['status'])
    return 0 if decision['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
