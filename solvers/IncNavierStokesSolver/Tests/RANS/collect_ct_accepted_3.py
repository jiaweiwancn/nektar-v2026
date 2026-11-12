#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT3 primitive decisions. Heavy states/matrices remain on Linux."""
import argparse, copy, hashlib, json, os, subprocess, xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS, write_json, exact_vector_hash
os.environ.update(THREADS)
from ct_accepted3_run import CASES
from ct_accepted_open import require
from ct_accepted_pilot import budget_classification, artifact_hash
from ct_accepted_budget_gate import physical_budgets, paired_coverage, required_times
from ct_accepted2_diagnosis import function_identity
from rans_request_manifest import RequestManifest, sha, runtime_identity, canonical
from ct_accepted2_storage import canonical_location


def small_initial(data):
    data=copy.deepcopy(data)
    for key in ('grid','sides','global_weights'):
        value=data.pop(key,None)
        if value is not None:data[key+'_sha256']=hashlib.sha256(canonical(value)).hexdigest();data[key+'_count']=len(value)
    return data


def reclassify(root,receipt,initial,ct1,runtime):
    result=copy.deepcopy(receipt);c=initial['construction'];members={};expected={};bindings={}
    for dt,row in receipt['members'].items():
        directory=Path(row['directory']);report=json.loads((directory/'case.evolution.json').read_text())
        require(sha(directory/'case.evolution.json')==row['native_report_sha256'],'changed_native_report')
        require(exact_vector_hash(report['state'])==row['state_sha256'],'changed_actual_state')
        require(artifact_hash(report['checkpoint'])==row['history_sha256'],'changed_native_history')
        require(sha(directory/'case.xml')==row['identity']['configuration_sha256'],'changed_actual_configuration')
        root=ET.parse(directory/'case.xml').getroot()
        params={n.text.split('=',1)[0].strip():n.text.split('=',1)[1].strip() for n in root.findall('./CONDITIONS/PARAMETERS/P') if n.text and '=' in n.text}
        require(abs(float(params['TimeStep'])-report['dt'])<1e-15,'actual_XML_timestep_differs')
        require(row['identity']['runtime']==runtime_identity(runtime/'IncNavierStokesSolver'),'foreign_evolution_runtime')
        model=c['model'];frozen=ct1/'pilot-v3'/(model+'-dt-'+('0.002' if model=='A95' else '0.001')+'-zero')/'case.xml'
        require(function_identity(ET.parse(directory/'case.xml').getroot())==function_identity(ET.parse(frozen).getroot()),'changed_frozen_physical_functions')
        budgets={}
        for control,b in row['budgets'].items():
            path=Path(b['execution']['command'][-1]);require(sha(path)==b['packet_sha256'],'changed_physical_budget_packet')
            packet=json.loads(path.read_text())
            require(packet['observation']['state']==report['state'],'budget_actual_state_identity')
            budgets[control]=budget_classification(packet,c)
        identity=dict(state_sha256=exact_vector_hash(report['state']),configuration_sha256=sha(directory/'case.xml'),
            history_sha256=row['history_sha256'],model=model,dt=report['dt'],time=report['time'],runtime_sha256=sha(runtime/'IncNavierStokesSolver'))
        expected[dt]=identity
        members[dt]=dict(budgets=budgets,identity=identity,time=report['time'],failure=report['failure'],initial_state_sha256=initial['members'][dt]['state_sha256'])
        maxima={}
        for t,audit in zip(report['audit_times'],report['state_audits'],strict=True):
            for key,witness in [('native_raw_divergence_range','native_raw_divergence_maximum_witness'),
                                ('native_outer_isotropic_production_work_maximum','native_outer_isotropic_production_work_witness')]:
                value=max((abs(v) for v in audit.get(key,[])),default=0.)
                if key not in maxima or value>maxima[key]['maximum_absolute']:
                    maxima[key]=dict(maximum_absolute=value,time=t,witness=audit.get(witness),stage='accepted_native_inventory')
        bindings[dt]=dict(identity=identity,budget_gate=physical_budgets(members[dt]),budgets=budgets,native_divergence_and_isotropic_work=maxima,
            actual_execution_metadata=dict(native_dt=report['dt'],XML_dt=float(params['TimeStep']),XML_segment_steps=int(params['NumSteps']),
                native_absolute_time=report['time'],native_total_accepted_steps=report['steps'],
                original_construction_metadata={k:json.loads((directory/'construction.json').read_text()).get(k) for k in ('dt','steps','target_time','xml_sha256')},
                scope='actual consumed XML/report/history are authoritative; pre-fix CT3 construction JSON copied transient CT1 template fields; original artifacts preserved'),
            functions_sha256=hashlib.sha256(canonical(function_identity(ET.parse(frozen).getroot()))).hexdigest())
    coverage=paired_coverage(members,receipt['target'],expected)
    reasons=coverage['reasons'][:]
    for dt,row in receipt['members'].items():
        if row['domain']['status']!='pass':reasons.append(dt+':independent_domain')
        if row['mass']['scaled_defect']>1e-8:reasons.append(dt+':mass')
        if row['corrections'] or row['additional_governing_solves']:reasons.append(dt+':extra_correction_or_solve')
    comparisons=receipt['temporal']['comparisons']
    if not any(abs(q['time']-receipt['target'])<1e-13 for q in comparisons):reasons.append('missing_required_endpoint_comparison')
    if any(q['status']!='pass' for q in comparisons):reasons.append('temporal_difference')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,coverage=coverage,actual_bindings=bindings,
        physical_completion=receipt['qualification']['physical_completion'],spatial_uncertainty='separate_common_physical_contrast',resource_compliance='separate_fresh_resource_audit')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','ct1','ct2','comparisons','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();runtime=a.ct2/'runtime';os.environ['LD_LIBRARY_PATH']=str(runtime/'lib')
    m=RequestManifest(a.output,dict(schema='CT3CompactRequiredScopeV1',source_sha256=sha(Path(__file__)),
        comparisons_sha256=sha(a.comparisons/'results.json'),runtime=runtime_identity(runtime/'IncNavierStokesSolver'),
        source_revision=subprocess.check_output(['git','-c','index.threads=1','-c','core.preloadIndex=false','rev-parse','HEAD'],cwd=a.root/'source',text=True).strip()))
    def save(name,value):write_json(a.output/name,value)
    decision=dict(schema='CT-Accepted-3-DecisionV1',diagnostic_execution='pass',cases={},resource_compliance={},
        physical_completion='per_case_only',spatial_uncertainty='unresolved_until_common_contrast',precursor_admission=False,benchmark_admission=False)
    raw={};exports={};requests={}
    for case in CASES:
        files=list((a.root/case).rglob('results.json'));initials=[];segments={}
        for f in files:
            q=json.loads(f.read_text())
            if q.get('schema')=='CT3InitialV1':initials.append((f,q))
            if q.get('schema')=='CT3SegmentV1':
                require(q['target'] not in segments,'duplicate_required_milestone_receipt')
                segments[q['target']]=(f,q)
        row=dict(required_milestones=[.008,.02,.2],milestones={},initialization='failed_or_unattempted')
        decision['cases'][case]=row
        if not initials:
            row['diagnostic_execution']='missing_initial_result';continue
        require(len(initials)==1,'ambiguous_initialization_receipt');f,ini=initials[0]
        save(case+'-initial.json',small_initial(ini));exports[case+'-initial.json']=dict(source=str(f),sha256=sha(f))
        row['initialization']=ini['status'];row['native_slots']={dt:r['slots'] for dt,r in ini['members'].items()}
        for t in (.008,.02,.2):
            if t not in segments:
                row['milestones'][str(t)]=dict(status='unattempted_or_incomplete',reason='required prior gate did not admit or diagnostic execution incomplete')
                continue
            f,q=segments[t];derived=reclassify(a.root,q,ini,a.ct1,runtime)
            name=case+'-T'+str(t)+'.json';save(name,q);exports[name]=dict(source=str(f),sha256=sha(f))
            save(case+'-T'+str(t)+'-derived.json',derived)
            row['milestones'][str(t)]=dict(status=derived['status'],physical_completion=derived['physical_completion'],reasons=derived['reasons'],complete_body=name,
                times={dt:r['time'] for dt,r in q['members'].items()},failures={dt:r['failure'] for dt,r in q['members'].items()})
        row['required_time_coverage']=required_times(list(segments))
        row['pilot_completion']='pass_T0p2_bounded_flat_only' if row['required_time_coverage']['status']=='pass' and all(q['status']=='pass' for q in row['milestones'].values()) else 'not_qualified_T0p2'
    for f in sorted(a.root.rglob('request-manifest.json')):
        requests[str(f.relative_to(a.root))]=dict(sha256=sha(f),body=json.loads(f.read_text()))
    for f in sorted(a.root.rglob('completed-artifacts.json')):
        requests[str(f.relative_to(a.root))]=dict(sha256=sha(f),body=json.loads(f.read_text()))
    save('requests-and-seals.json',requests)
    resources={};highwater=[];events=[]
    for directory in sorted(a.root.glob('launch-*')):
        file=directory/'resource-final.json'
        if file.exists():resources[directory.name]=dict(json.loads(file.read_text()),execution=json.loads((directory/'execution.json').read_text()))
        for filename,target in [('worker-guard-highwater.jsonl',highwater),('project-highwater.jsonl',highwater),('worker-guard-events.jsonl',events)]:
            if (directory/filename).exists():target.extend(dict(control=directory.name,kind=filename,**json.loads(line)) for line in (directory/filename).read_text().splitlines())
    decision['resource_compliance']=dict(status='pass' if resources and all(q['aggregate_status']==q['project_status']=='pass' for q in resources.values()) else 'fail_or_missing',
        maximum_aggregate=max((q['maximum_aggregate'] for q in resources.values()),default=None),maximum_project=max((q['maximum_project'] for q in resources.values()),default=None),
        project_limit=16,aggregate_limit=72,control_helpers_included=2,historical_CT1_CT2_failures='unchanged; not overwritten or rerun')
    save('fresh-resources.json',resources);save('resource-highwater.json',highwater);save('resource-events.json',events)
    comparisons=json.loads((a.comparisons/'results.json').read_text());save('equal-time-comparisons.json',comparisons)
    decision['spatial_uncertainty']=comparisons['evolved_spatial_status']
    decision['common_physical_domain']={k:v['common_grid_domain'] for k,v in comparisons['spaces'].items()}
    for case,row in decision['cases'].items():
        if comparisons['spaces'].get(case,{}).get('common_grid_domain',{}).get('status')=='fail':
            row['pilot_completion']='failed_additional_common_grid_domain'
    decision['diagnostic_execution']='pass' if all(r['initialization']=='pass' and '0.008' in r['milestones'] and r['milestones']['0.008']['status'] in ('pass','fail') for r in decision['cases'].values()) else 'incomplete_required_case_diagnostics'
    controls=json.loads((a.root/'controls-v1/results.json').read_text());save('degree-pressure-and-coverage-controls.json',controls)
    decision['boundary_pressure_implementation']=controls['status'];decision['admission_negatives']=controls['coverage']['status']
    for f in a.root.glob('profile*/**/results.json'):
        if f.name=='results.json':
            name='profile-'+hashlib.sha256(str(f).encode()).hexdigest()[:8]+'.json';save(name,json.loads(f.read_text()));exports[name]=dict(source=str(f),sha256=sha(f))
    # Pin selected heavy reports/configurations and histories, not full matrices
    # or repeated per-step fields in Git. Stage JSON already carries array hashes.
    for f in sorted(a.root.glob('**/case.evolution.json')):
        raw[str(f.relative_to(a.root))]=sha(f)
        if (f.parent/'case.xml').exists():raw[str((f.parent/'case.xml').relative_to(a.root))]=sha(f.parent/'case.xml')
    save('native-artifact-hashes.json',raw);save('export-provenance.json',exports)
    save('decision.json',decision)
    source=a.root/'source';git=['git','-c','index.threads=1','-c','core.preloadIndex=false']
    save('source-scope.json',dict(starting_revision='db37c6b11882e7a06413a29eb6623ab212c4d4a5',
        export_revision=subprocess.check_output(git+['rev-parse','HEAD'],cwd=source,text=True).strip(),
        export_status=subprocess.check_output(git+['status','--porcelain'],cwd=source,text=True),
        governing_runtime=runtime_identity(runtime/'IncNavierStokesSolver'),observer_runtime=runtime_identity(runtime/'RANSPreservationProbe'),
        governing_runtime_build_source='a897896524b8761d3a954d491f96d1b33941d94f',
        code={str(f.relative_to(source)):sha(f) for f in (source/'solvers/IncNavierStokesSolver/Tests/RANS').glob('ct_accepted3*.py')},
        specifications={str(f.relative_to(source)):sha(f) for f in (source/'docs/literature').glob('*A93*A95*tex')},
        bridge_sha256=sha(a.ct2/'relocation-runtime-bridge.json'),scope='physical runtime reused exactly; separate timing-only runtime explicitly tested; no equations, input functions, ownership or stage order changed'))
    m.seal(a.output,[f.name for f in a.output.iterdir() if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json')])
    save('package-hashes.json',{f.name:sha(f) for f in a.output.iterdir() if f.is_file()})
    print(json.dumps(decision));return 0
if __name__=='__main__':raise SystemExit(main())
