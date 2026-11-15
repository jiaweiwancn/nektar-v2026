#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Two accepted target replay steps, one rejected attempt and lossless retention."""
import argparse,json,os,resource,time
from pathlib import Path
from collections import defaultdict
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted9_controls import enable
from ct_accepted11_target import bounded_budgets
from ct_accepted12_retention import retain
from ct_accepted_feasibility import write_json,require
from ct_accepted_open import run
from ct_accepted_pilot import artifact_hash
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,sha,runtime_identity
import xml.etree.ElementTree as ET


def snapshots(report,legacy):
    check={};eligible=[]
    for now,old in zip(report['pimple_attempts'],legacy,strict=True):
        for kind in ('flow_snapshot','k_snapshot'):
            f=Path(now[kind]);g=Path(old[kind]);require(sha(f)==now[kind+'_sha256'] and sha(g)==old[kind+'_sha256'],'changed stage snapshot')
            a=n.read(f);b=n.read(g);err=[]
            for key,value in a.items():
                if key in ('frozen_flow_matrix','free_velocity_mass'):
                    require(value['sha256']==b[key]['sha256'],'captured sparse matrix changed');eligible.append(Path(value['file']))
                elif key=='pressure_operator':
                    for name in ('G','D','C'):
                        require(value[name]['sha256']==b[key][name]['sha256'],'captured complete pressure block changed');eligible.append(Path(value[name]['file']))
                elif key=='matrix_file':eligible.append(Path(value))
                else:require(value==b[key],'stage primitive changed:'+kind+':'+key)
            eligible.append(f);check[str(f)]=dict(sha256=sha(f),legacy_sha256=sha(g),all_values_and_matrix_hashes_equal=True)
    return check,eligible


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct10','ct11','runtime','probe','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    require(sha(a.runtime/'IncNavierStokesSolver')=='5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7','foreign governing runtime')
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');pair=sealed(a.ct11/'target-v1');c=pair['case'];old=pair['members']['0.001'];folder=Path(old['directory']);legacy=n.read(folder/'case.evolution.json')
    require(sha(folder/'case.evolution.json')==old['native_report_sha256'],'changed CT11 native report')
    zero=a.ct11/'target-v1/zero-0.001';z=n.read(zero/'case.evolution.json')
    m=RequestManifest(a.output,dict(schema='CT12TargetProfileRetentionV1',source=sha(Path(__file__)),contract=sha(a.contract),
        case=c,dt=.001,accepted_steps=2,deliberate_rejected_attempts=1,initial_history=artifact_hash(z['checkpoint']),
        legacy_report=sha(folder/'case.evolution.json'),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.probe/'RANSPreservationProbe'),
        scope='replay only,not a new physical qualification or longer continuation; detailed timing and external byte-preserving retention'))
    result=dict(schema='CT12TargetProfileRetentionV1',status='incomplete',members={});eligible=[];required={};reports={};template=folder/'case.xml';parent=Path(z['checkpoint'])
    for step in (1,2):
        dest=a.output/('step'+str(step));n.make(dest,a.runtime,c,.001,1,parent,template);enable(dest)
        os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_TIMING_FILE']=str(dest/'timing.jsonl');os.environ['NEKTAR_RANS_DETAILED_TIMING']='1'
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE');os.environ.pop('NEKTAR_RANS_DETAILED_TIMING')
        r=n.read(dest/'case.evolution.json');require(e['exit_code']==0 and not r['failure'] and r['time']==legacy['audit_times'][step],'target profiling replay failed')
        require(exact_vector_hash(r['state'])==exact_vector_hash(legacy['accepted_states'][step]),'detailed timing changed accepted state')
        if step==2:require(same_history(legacy,r),'segmented profiling restart changed native history')
        checks,files=snapshots(r,[legacy['pimple_attempts'][step-1]]);eligible.extend(files);eligible.append(dest/'solver.log')
        costs=defaultdict(lambda:dict(calls=0,CPU=0.,wall=0.));exclusive=0.;initial_CPU=0.
        for line in (dest/'timing.jsonl').read_text().splitlines():
            t=json.loads(line)
            if t['step']<=0:initial_CPU+=t['exclusive_process_cpu_seconds'];continue
            v=costs[t['name']];v['calls']+=t.get('calls',1);v['CPU']+=t['exclusive_process_cpu_seconds'];v['wall']+=t['exclusive_wall_seconds'];exclusive+=t['exclusive_process_cpu_seconds']
        for f in (dest/'case.evolution.json',dest/'case.xml',Path(r['checkpoint']),dest/'timing.jsonl'):required[str(f.relative_to(a.output))]=sha(f)
        result['members'][str(step)]=dict(time=r['time'],execution=e,accepted_state_sha256=exact_vector_hash(r['state']),history_sha256=artifact_hash(r['checkpoint']),configuration_sha256=sha(dest/'case.xml'),
            native_report_sha256=sha(dest/'case.evolution.json'),timing_sha256=sha(dest/'timing.jsonl'),initialization_exclusive_CPU=initial_CPU,
            accepted_exclusive_CPU=exclusive,categories=dict(costs),stage_equivalence=checks,
            legacy_exclusive_CPU=pair['first_stage_resource_plan']['timing'][0]['exclusive_step_CPU'][str(step)])
        parent=Path(r['checkpoint']);reports[step]=r;write_json(a.output/'progress.json',result);print('CT12 profile accepted',step,exclusive,flush=True)
    # The deliberate failure exercises retention/rollback, not constitutive admission.
    dest=a.output/'rollback';n.make(dest,a.runtime,c,.001,1,Path(reports[1]['checkpoint']),a.output/'step1/case.xml');enable(dest)
    root=ET.parse(dest/'case.xml').getroot();set_parameter(root,'RANSWallPimpleHardAuditRejectStep',2);write(root,dest/'case.xml');os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    e=run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log');r=n.read(dest/'case.evolution.json')
    require(e['exit_code']!=0 and 'deliberate_PIMPLE_hard_audit_rejection' in r['failure'] and same_history(reports[1],r),'rejected operation failed to retain prior native history')
    require(bool(r['rejected_state']) and r['pimple_attempts'],'missing actual rejected candidate/stage')
    checks,files=snapshots(r,[legacy['pimple_attempts'][1]]);eligible.extend(files);eligible.append(dest/'solver.log')
    for f in (dest/'case.evolution.json',dest/'case.xml',Path(r['checkpoint'])):required[str(f.relative_to(a.output))]=sha(f)
    result['rollback']=dict(status='pass',execution=e,failure=r['failure'],time=r['time'],last_accepted_state_sha256=exact_vector_hash(r['state']),
        rejected_state_sha256=exact_vector_hash(r['rejected_state']),history_sha256=artifact_hash(r['checkpoint']),native_report_sha256=sha(dest/'case.evolution.json'),stage_equivalence=checks)
    # Independent physical budgets consume these actual replay stage packets
    # before any payload becomes eligible for compressed retention.
    initial=n.read(a.ct10/'target-preflight-v1/target/readiness.json')
    budget=bounded_budgets(a.probe,a.output/'step2',initial,reports[2],c)
    require(budget['gate']['status']=='pass','actual replay physical-budget audit failed');result['budgets']=budget
    audit_receipt=a.output/'independent-audit.json';write_json(audit_receipt,budget);required[str(audit_receipt.relative_to(a.output))]=sha(audit_receipt)
    controls={'native','q48f16','q64f16','q64f32'}
    for name in controls:
        directory=a.output/'step2'/name
        required[str((directory/'case.xml').relative_to(a.output))]=sha(directory/'case.xml')
        eligible.extend(f for f in directory.iterdir() if f.is_file() and f.name!='case.xml')
    accepted_categories={}
    for step in (1,2):
        directory=a.output/('step'+str(step));sizes=dict(native=0,audit=0)
        for f in directory.rglob('*'):
            if f.is_file():sizes['audit' if f.relative_to(directory).parts[0] in controls else 'native']+=f.stat().st_size
        accepted_categories[str(step)]=sizes
    binding=dict(model='A95',case=c,governing_sha256=sha(a.runtime/'IncNavierStokesSolver'),historical_report_sha256=old['native_report_sha256'],
        initial_sha256=old['actual_initial_sha256'],last_accepted_sha256=exact_vector_hash(reports[2]['state']),native_rejected_report_sha256=sha(dest/'case.evolution.json'))
    proof=dict(status='pass',binding=binding,required_files=required,
        prerequisite='every accepted state,exact stage fields/operators,history,rejection and independent final budget audited before pruning; required reports/histories remain uncompressed')
    result['retention']=retain(a.output,a.output.parent/'retention-v1',eligible,proof)
    index=n.read(a.output.parent/'retention-v1/index.json');category={}
    for step,sizes in accepted_categories.items():
        category[step]={}
        for lane in ('native','audit'):
            subset={k:v for k,v in index['files'].items() if Path(k).parts[0]=='step'+step and (Path(k).parts[1] in controls)==(lane=='audit')}
            blobs={v['blob'] for v in subset.values()}
            compressed=sum((a.output.parent/'retention-v1'/f).stat().st_size for f in blobs)
            category[step][lane]=dict(raw_bytes=sizes[lane],eligible_bytes=sum(v['bytes'] for v in subset.values()),unique_compressed_bytes=compressed,
                retained_bytes=sizes[lane]-sum(v['bytes'] for v in subset.values())+compressed)
    result['retention']['accepted_only_categories']=category
    result['retention']['forecast_scope']='Accepted-only categories exclude the intentionally duplicate rollback payload; failure dedup cannot be counted as ordinary evolution savings.'
    result.update(status='pass',scope='two accepted target steps replayed exactly with genuine restart,one zero-acceptance rejection; physical pressure gate unchanged; output policy external only')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json','independent-audit.json']);return 0
if __name__=='__main__':raise SystemExit(main())
