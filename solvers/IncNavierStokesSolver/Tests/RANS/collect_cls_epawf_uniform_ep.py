#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Join immutable source/control readbacks, preserving each case's own h0."""
import argparse
import copy
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import load,write_json,require,sha,state_hash,ROOT
from run_cls_epawf_pimple import identity
from compare_cls_epawf_pimple import physical_differences
from continue_cls_epawf_pimple import header,process_identity
from analyse_cls_epawf_constitutive_v2 import read_report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--upstream',type=Path,required=True)
    p.add_argument('--uniform-analysis',type=Path,required=True);p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();old=a.upstream.resolve();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    u=load(a.uniform_analysis/'summary.json');require(u['status']=='analysis_completed','uniform readback incomplete')
    source=load(a.sources/'summary.json');require(source['status']=='read_only_source_diagnosis_completed','source readback incomplete')
    original={}
    paths={'initial':old/'analysis16/initial','step8':old/'analysis16/step8','step16':old/'analysis16/step16',
           'step32':old/'monitor/step-000032/step32','step64':old/'monitor/step-000064/step64','step72':old/'monitor/step-000072/step72'}
    summaries=[load(old/'analysis16/summary.json'),load(old/'monitor/step-000032/summary.json'),
               load(old/'monitor/step-000064/summary.json'),load(old/'monitor/step-000072/summary.json')]
    for s in summaries:original.update(s['states'])
    upstream_rows=[json.loads(x) for x in (old/'campaign/trajectory/steps.jsonl').read_text().splitlines()]
    uniform_rows=[json.loads(x) for x in (r/'uniform-campaign/trajectory/steps.jsonl').read_text().splitlines()]
    result=dict(schema=1,base='56657a7069af7c9600277611833bd37745d604fa',branch='codex/cls-epawf-pimple-one-sweep',
        source_revision=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        source_status=subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip(),
        root=str(r),upstream_root=str(old),collected_unix=time.time(),
        policy_checks=load(r/'policy-checks-r2/summary.json'),durability_checks=load(r/'durability.json'),
        qp_check=(r/'qp-checks.txt').read_text().strip(),sources=source,
        uniform_analysis=copy.deepcopy(u),upstream_states=original,common_time=[],within_case_drift={})
    result['policy_checks']['coverage_tripwire']={k:v for k,v in result['policy_checks']['coverage_tripwire'].items() if k!='scalar_domain'}
    for name in list(result['sources']['states']):
        s=result['sources']['states'][name]
        # Full near-cutoff profiles and weak coefficient loads remain in Linux.
        s['points']=[v for v in s['points'] if not v.get('valid') or (v.get('z')==0. and (v['variant']==-1 or v['x']==1.6))]
    result['uniform_analysis'].pop('progress',None)
    for key in ('uniform_analysis','upstream_states'):
        states=result[key].get('states',{}) if key=='uniform_analysis' else result[key]
        for s in states.values():s.pop('profiles',None)
    # Keep per-step compact increments, not all repeated inventory maps/cuts.
    for b in result['uniform_analysis']['budgets']:
        b['active_constraint_count']=len(b.pop('active_constraints',[]));b.pop('inventory_counts',None)
    for name in sorted(set(paths)&set(u['states']),key=lambda n:u['states'][n]['time']):
        up=load(paths[name]/'samples.json');un=load(a.uniform_analysis/name/'samples.json')
        require(abs(up['time']-un['time'])<1e-13,'cross-policy comparison changed time')
        own_up=load(paths['initial']/'samples.json');own_un=load(a.uniform_analysis/'initial/samples.json')
        result['common_time'].append(dict(name=name,time=un['time'],
            cross_policy_physical_L2=physical_differences(up['physical_samples'],un['physical_samples']),
            upstream_drift_from_own_h0=physical_differences(own_up['physical_samples'],up['physical_samples']),
            uniform_drift_from_own_h0=physical_differences(own_un['physical_samples'],un['physical_samples']),
            upstream_integrated_k=original[name]['global_k_integral'],uniform_integrated_k=u['states'][name]['global_k_integral'],
            upstream_near_k=original[name]['near_k_integral'],uniform_near_k=u['states'][name]['near_k_integral'],
            upstream_mean_stress=original[name]['mean_wall_stress'],uniform_mean_stress=u['states'][name]['mean_wall_stress'],
            upstream_alpha0_norms=original[name]['spatial_block_norms'],uniform_alpha0_norms=u['states'][name]['spatial_block_norms'],
            upstream_trace=original[name]['epsilon_trace_defect'],uniform_trace=u['states'][name]['epsilon_trace_defect']))
    for name,root in [('upstream',old/'campaign/trajectory'),('uniform',r/'uniform-campaign/trajectory')]:
        copy_to=a.output/(name+'-checkpoint');copy_to.mkdir()
        import hashlib,struct
        for attempt in range(20):
            progress=load(root/'progress.json');payload=(Path(progress['checkpoint'])/'State.bin').read_bytes()
            if hashlib.sha256(payload).hexdigest()==progress['checkpoint_sha256']:break
            time.sleep(.1)
        else:raise ValueError('checkpoint rotated throughout readback; no checkpoint is relabelled')
        (copy_to/'State.bin').write_bytes(payload);cp=header(copy_to,True)
        require(cp['sha256']==progress['checkpoint_sha256'] and cp['steps']==progress['accepted_steps'],'checkpoint snapshot differs')
        signature_length=struct.unpack_from('<Q',payload,24)[0]
        expected=load(root/'runtime/identity.json')['frozen'];env=os.environ.get('LD_LIBRARY_PATH')
        try:
            os.environ['LD_LIBRARY_PATH']=str(root/'runtime/lib');runtime=identity(root/'runtime/IncNavierStokesSolver')
        finally:
            if env is None:os.environ.pop('LD_LIBRARY_PATH',None)
            else:os.environ['LD_LIBRARY_PATH']=env
        require(runtime==expected,'frozen executable/library identity changed')
        result[name+'_campaign']=dict(progress=progress,runtime=runtime,checkpoint_snapshot=str(copy_to),checkpoint_sha256=cp['sha256'],
            checkpoint_steps=cp['steps'],checkpoint_time=cp['time'],checkpoint_signature_sha256=hashlib.sha256(payload[32:32+signature_length]).hexdigest(),
            native_process_now=process_identity(progress['native_identity']['pid']) if progress.get('native_identity') else None,
            driver_process_now=process_identity(progress['pid']),
            supervisor_status=load(root.parent/'status.json'))
    result['uniform_launch']=load(r/'uniform-launch.json');result['uniform_monitor_launch']=load(r/'uniform-monitor-launch.json')
    result['upstream_terminal']=load(r/'upstream-terminal/manifest.json')
    attempt=result['upstream_terminal']['first_rejection']['attempt']
    result['upstream_terminal']['first_rejection']['attempt']={k:v for k,v in attempt.items() if k not in ('input_domain','k_domain')}
    for key in ('active_constraints','inventory_counts','discovered_constraint_families'):
        result['upstream_terminal']['first_rejection']['attempt']['k_bound'].pop(key,None)
    # Exact native one-sweep counts, plus separately counted post-k work.
    for name,rows in [('upstream',upstream_rows),('uniform',uniform_rows)]:
        keys=('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves',
              'additional_governing_equation_solves','bound_calls','bound_activations','bound_active_set_iterations','bound_mass_factor_applications','bound_kkt_solves')
        result[name+'_counts']={k:sum(q['attempt'].get(k,0) for q in rows) for k in keys}
        result[name+'_counts']['accepted_steps']=len(rows)
        result[name+'_counts']['implicit_callbacks_accepted']=len(rows)
    def config(xml):
        t=ET.parse(xml).getroot()
        return dict(sha256=sha(xml),path=str(xml),solver_info={v.attrib['PROPERTY']:v.attrib['VALUE'] for v in t.findall('./CONDITIONS/SOLVERINFO/I')},
            parameters=[v.text.strip() for v in t.findall('./CONDITIONS/PARAMETERS/P')],
            boundary_xml=ET.tostring(t.find('./CONDITIONS/BOUNDARYCONDITIONS'),encoding='unicode'),
            functions=[ET.tostring(v,encoding='unicode') for v in t.findall('./CONDITIONS/FUNCTION')])
    result['case_inputs']={name:config(Path(rows[0]['segment'])/'case.xml') for name,rows in [('upstream',upstream_rows),('uniform',uniform_rows)]}
    up=ET.parse(Path(upstream_rows[0]['segment'])/'case.xml').getroot();un=ET.parse(Path(uniform_rows[0]['segment'])/'case.xml').getroot()
    result['same_spatial_and_external_problem']={part:ET.tostring(up.find(part))==ET.tostring(un.find(part)) for part in ('GEOMETRY','EXPANSIONS','CONDITIONS/BOUNDARYCONDITIONS')}
    require(all(result['same_spatial_and_external_problem'].values()),'uniform control changed mesh/spaces/external XML boundary data')
    require(result['case_inputs']['upstream']['functions']==result['case_inputs']['uniform']['functions'],'physical function changed')
    result['archived_xml']=dict(path=str(ROOT/'docs/plans/cls-epawf-revised-evidence/m8-r1/archive-session.xml'),
        sha256=sha(ROOT/'docs/plans/cls-epawf-revised-evidence/m8-r1/archive-session.xml'))
    files=[]
    for base in (r/'validation-runtime',r/'trace-probe-runtime',r/'policy-checks-r2',r/'upstream-terminal',a.sources,a.uniform_analysis):
        if not base.exists():continue
        for f in sorted(base.rglob('*')):
            if f.is_file():files.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)))
    for folder in (r/'uniform-campaign/trajectory').glob('steps-*'):
        for suffix in ('case.xml','case.evolution.json.gz','receipt.json','accepted-metrics.json.gz'):
            f=folder/suffix
            if f.exists():files.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)))
    write_json(a.output/'artifact-index.json',dict(immutable_artifacts=files,mutable_progress_is_snapshot=True))
    result['artifact_index']=dict(path=str(a.output/'artifact-index.json'),sha256=sha(a.output/'artifact-index.json'),files=len(files))
    result['classification']='short-time reference preservation comparison; no steady state or physical-time accuracy claim'
    write_json(a.output/'native.json',result);print(a.output/'native.json')


if __name__=='__main__':main()
