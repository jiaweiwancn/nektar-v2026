#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT8 primitives, failures and exact provenance, after family closure."""
import argparse,copy,json,os,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_feasibility import write_json,require
from ct_accepted8_gate import pair
from ct_accepted7_gate import resources
from rans_request_manifest import RequestManifest,sha,runtime_identity

def read(p):return json.loads(Path(p).read_text())
def sealed(folder):
    f=folder/'completed-artifacts.json';s=read(f);require(all(sha(folder/n)==h for n,h in s['files'].items()),'changed sealed evidence:'+str(folder));return read(folder/'results.json')
def receipt(r):
    return {k:v for k,v in r.items() if k not in ('effective_initial','live_fixed_boundary','fixed_boundary')}
def compact_metrics(row):
    r=copy.deepcopy(row);t=r.get('traction',{});j=t.pop('interface_rows',[])
    if j:t['largest_component_jump_witness']=max(j,key=lambda q:max(map(abs,q['traction_jump'])))
    return r

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','preflight','controls','assessment','pressure','regression','runtime','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--pair',type=Path,action='append',default=[]);a=p.parse_args();s=a.root/'source'
    pre=sealed(a.preflight);control=sealed(a.controls);assessment=sealed(a.assessment);ref=sealed(a.root/'reference-v3');build=sealed(a.root/'build-v3')
    pressure=sealed(a.pressure);regression=sealed(a.regression)
    require(pressure['status']==regression['status']=='pass','required read-only diagnosis or historical equivalence incomplete')
    m=RequestManifest(a.output,dict(schema='CT8CompactEvidenceV1',reference=sha(a.root/'reference-v3/results.json'),preflight=sha(a.preflight/'results.json'),compatibility=sha(a.controls/'results.json'),assessment=sha(a.assessment/'results.json'),pressure=sha(a.pressure/'results.json'),regression=sha(a.regression/'results.json'),pairs={str(f):sha(f/'results.json') for f in a.pair},source=sha(Path(__file__))))
    scopes=dict(central_reference=ref['central_status'],extraction_sensitivity=dict(failed_roots=ref['extraction_summary']['failed_roots'],scope='central admission is not robust to all extraction corners; no favorable realization selected'),
        native_operator='pass' if pre['status']=='pass' and control['status']=='pass' else 'fail',historical_flat_equivalence=regression['status'],pressure_diagnostic_execution=pressure['status'],bounded_preservation={},terrain='not_admitted',full_domain='blocked_unchanged3000slot_guard',CT6='historical failed temporal/degree gates unchanged; read-only vector assessment only')
    write_json(a.output/'fit-and-reference.json',ref)
    small=copy.deepcopy(pre)
    for name,row in small.get('initial',{}).items():row['initialization_receipt']=receipt(row['initialization_receipt']);row['metrics']=compact_metrics(row['metrics'])
    for row in small.get('small_operator',{}).values():row['initialization_receipt']=receipt(row['initialization_receipt'])
    if small.get('projection_contrast'):
        for name in ('coarse_P3','fine_P4'):small['projection_contrast'][name]=compact_metrics(small['projection_contrast'][name])
    write_json(a.output/'native-preflight.json',small);write_json(a.output/'compatibility.json',control)
    write_json(a.output/'pressure-diagnosis.json',pressure);write_json(a.output/'flat-history-regression.json',regression)
    negatives={}
    for folder in a.pair:
        r=sealed(folder);gate=pair(r);scopes['bounded_preservation'][str(r['target'])]=gate
        out=copy.deepcopy(r)
        for name,row in out['members'].items():
            row['initialization_receipt']=receipt(row['initialization_receipt']);row['accepted_metrics']=[compact_metrics(q) for q in row['accepted_metrics']]
            histories=row.pop('traction_history',[]);initial=row.pop('actual_initial_traction_profile',{})
            row['signed_traction_history']=[dict(time=q['time'],**{k:q['profile'][k] for k in ('state_sha256','mean_wall_shear','normalized_wall_traction_error','normalized_interface_traction_RMS','wall_spatial_standard_deviation','minimum_wall','maximum_wall')}) for q in histories]
        out['recomputed_qualification']=gate;write_json(a.output/('pair-T'+str(r['target'])+'.json'),out)
        timing={}
        for name,row in r['members'].items():
            f=Path(row['directory'])/'timing.jsonl';groups={}
            if f.exists():
                for line in f.read_text().splitlines():
                    q=json.loads(line);v=groups.setdefault(q['name'],dict(calls=0,exclusive_process_cpu_seconds=0.,exclusive_wall_seconds=0.))
                    v['calls']+=q.get('calls',1)
                    for key in ('exclusive_process_cpu_seconds','exclusive_wall_seconds'):v[key]+=q[key]
            timing[name]=dict(artifact=str(f),sha256=sha(f) if f.exists() else None,exclusive_scopes=groups,execution=row['execution'],scope='exclusive scopes only; no sum of nested inclusive intervals; solver execution includes initialization/output')
        write_json(a.output/('timing-T'+str(r['target'])+'.json'),timing)
        for label,mutate,reason in [('missing_member',lambda q:q['members'].pop('0.001'),'missing_timestep_member'),('missing_epsilon',lambda q:q['members']['0.002']['budgets']['budgets']['native']['equations'].pop('outer_epsilon'),'0.002:budget_primitives'),
            ('foreign_model',lambda q:q['members']['0.002'].__setitem__('signature','A93'),'0.002:foreign_identity'),('incomplete_endpoint',lambda q:q['members']['0.002'].__setitem__('time',.007),'0.002:incomplete_or_rejected'),
            ('wrong_timestep_primitive',lambda q:q['temporal']['U']['difference'].__setitem__('weighted_RMS',1.),'temporal:U'),('wrong_matrix',lambda q:q['members']['0.002']['budgets']['captured_k_matrix'].__setitem__('maximum_defect',1.),'0.002:captured_frozen_matrix')]:
            q=copy.deepcopy(r);mutate(q);bad=pair(q);require(bad['status']=='fail' and reason in bad['reasons'],'wrong negative reason:'+label);negatives[str(r['target'])+'/'+label]=dict(status='pass',expected_reason=reason,reasons=bad['reasons'])
    write_json(a.output/'classifier-negatives.json',negatives)
    for name in ('observations','temporal','domain'):write_json(a.output/(name+'.json'),read(a.assessment/(name+'.json')))
    closed={str(f):read(f/'resource-final.json') for f in sorted(a.root.glob('launch-*')) if (f/'resource-final.json').exists()}
    resource=resources(closed);scopes['numerical_resources']=resource
    require(resource['status']=='pass','resource compliance failed')
    costs={}
    for path,r in closed.items():
        f=Path(path);e=read(f/'execution.json');costs[f.name]={k:e[k] for k in ('command','exit_code','wall_seconds','allocation_wait_seconds','process_cpu_seconds','process_cpu_scope')}
        costs[f.name].update(resource=r,execution_sha256=sha(f/'execution.json'),project_highwater=[json.loads(x) for x in (f/'project-highwater.jsonl').read_text().splitlines()],aggregate_highwater=[json.loads(x) for x in (f/'worker-guard-highwater.jsonl').read_text().splitlines()])
    write_json(a.output/'resource-costs.json',costs)
    failures={}
    for root in sorted(a.root.glob('preflight-v*')):
        f=root/'failure.json'
        if f.exists():failures[root.name]=dict(path=str(f),sha256=sha(f),body=read(f),scope='retained static preparation, no admitted primary trajectory')
    for root in ('controls-v1','assessment-v1'):
        f=a.root/root/'results.json'
        if f.exists():q=read(f);failures[root]=dict(path=str(f),sha256=sha(f),status=q.get('status'),failure=q.get('failure'),scope='retained independent diagnostic preparation')
    f=a.root/'completion-jobs-v1/results.json'
    if f.exists():
        for i,row in enumerate(read(f)['jobs']):
            if row['exit_code']!=0:
                log=f.parent/(str(i)+'.log');failures['completion-job-'+str(i)]=dict(execution=row,log=str(log),log_sha256=sha(log),log_tail=log.read_text()[-3000:],scope='independent diagnostic execution failure, not a new physical trajectory or a relabelled qualification')
    write_json(a.output/'failed-preparations.json',failures)
    write_json(a.output/'decision.json',dict(schema='CT-Accepted-8-DecisionV1',scopes=scopes,
        physical_applicability='conditional central flat mean-profile input calibration; no benchmark-prediction validation, developed precursor or terrain admission',
        historical_CT4_CT7='unchanged, including narrow A93 traction, CT4 pressure interpretation and CT6/7 failed time/degree gates',collection_resource='attach this family closed execution/resource receipts before publication; no pending physical trajectory'))
    g=['git','-c','core.preloadIndex=false','-c','index.threads=1','-c','core.filemode=false'];revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=s,text=True).strip()
    changed=subprocess.check_output(g+['diff','--name-only','6f545f53299915657a18368ac18879a7fda935e5..HEAD'],cwd=s,text=True).splitlines()
    provenance=dict(schema='CT8SourceRuntimeV1',reviewed_start='6f545f53299915657a18368ac18879a7fda935e5',source_revision=revision,source_status=subprocess.check_output(g+['status','--porcelain'],cwd=s,text=True),build=build,
        runtime={n:runtime_identity(a.runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSCurvedAcceptedKernelProbe','RANSNewModelsProbe')},
        changed_source_sha256={n:sha(s/n) for n in changed if (s/n).is_file()},specification_sha256={n:sha(s/'docs/literature'/n) for n in ('accepted_awf_appendix_A93_A95_E94_E95.tex','A93_A95_E94_E95_neutral_ABL_EP_AWF_combined.tex')},
        dataset_sha256={n:sha(s/'docs/plans/cao-tamura-validation-data'/n) for n in ('profiles.csv','snapshot.json','source_manifest.json','extraction.json')},
        receipts={str(f):dict(result=sha(f/'results.json'),request=sha(f/'request-manifest.json'),seal=sha(f/'completed-artifacts.json')) for f in [a.root/'reference-v3',a.preflight,a.controls,a.assessment,a.pressure,a.regression,*a.pair]},
        storage=subprocess.check_output(['findmnt','-n','-o','SOURCE,FSTYPE,TARGET','/media/jay/新加卷1'],text=True).strip(),threads=THREADS,
        reuse='CT7 and original failures immutable. Reference-v2 consumed by fresh physical histories; reference-v3 adds boundary negatives on RT3. Assessment-v2 geometry used pinned RT2; final metadata supplement on RT3 is separately identified. CT6 retained vectors are read-only. No old trajectory qualifies the displaced route.')
    write_json(a.output/'provenance.json',provenance)
    m.seal(a.output,[f.name for f in a.output.glob('*.json') if f.name not in ('request-manifest.json','completed-artifacts.json')]);print('CT8 compact scopes written; inspect separate primitive gates',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
