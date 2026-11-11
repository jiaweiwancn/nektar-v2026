#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Supplement an immutable H1 pilot with zero-step identity and stage checks."""
import argparse
import os
from pathlib import Path
import platform
import subprocess
import xml.etree.ElementTree as ET

from cao_tamura_budget_h1 import compatibility_case
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import atomic, freeze
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash, set_parameter, write


def verify_label_upgrade(old,new):
    require(same_history(old,new),'diagnostic labels changed state/history/time')
    a0=load(Path(old['pimple_attempts'][0]['k_snapshot']));a1=load(Path(new['pimple_attempts'][0]['k_snapshot']))
    for name in ('committed','post_flow','pre_k','raw','corrected','history','fixed_reference'):
        require(state_hash(a0[name])==state_hash(a1[name]),'diagnostic labels changed a governing stage: '+name)
    # The historical finite report's neutral-log ratios are retained as
    # explicitly named scale diagnostics, not experimental profile errors.
    def renamed(audits):
        rows=[]
        for audit in audits:
            audit=dict(audit)
            if 'accepted_awf_reference_departure' in audit:
                audit['neutral_log_scale_ranges_diagnostic_only']=audit.pop('accepted_awf_reference_departure')
            if 'r4_outer_gradient_witness' in audit:
                audit['neutral_log_scale_gradient_witness_diagnostic_only']=audit.pop('r4_outer_gradient_witness')
            rows.append(audit)
        return rows
    require(renamed(old['state_audits'])==new['state_audits'],'label upgrade changed audit values')
    return dict(stage_history_bitwise=True,audit_values_unchanged=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('solver','campaign','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--upgraded-solver',type=Path,help='optional finite-only diagnostic-label revision; compare one independent native step')
    p.add_argument('--verify-label-only',action='store_true',help='verify already completed native one-step reports; no setup or evolution')
    a=p.parse_args();require(platform.system()=='Linux','authoritative Linux contracts required')
    if a.verify_label_only:
        old=load(a.output/'baseline-one-step/case.evolution.json');new=load(a.output/'label-upgrade-one-step/case.evolution.json')
        checks=verify_label_upgrade(old,new)
        atomic(a.output/'label-verification.json',dict(status='passed',checks=checks,governing_solves=0,physical_commits=0,
            old_report_sha256=sha(a.output/'baseline-one-step/case.evolution.json'),
            new_report_sha256=sha(a.output/'label-upgrade-one-step/case.evolution.json'),
            scope='postprocess immutable completed verification legs; no re-solve after a Python path-type error'))
        print(a.output/'label-verification.json');return
    a.output.mkdir(parents=True,exist_ok=False)
    solver,identity=freeze(a.solver.resolve(),a.output/'runtime')
    initial=load(a.campaign/'initial/case.evolution.json');whole=load(a.campaign/'whole8/case.evolution.json')
    cp=Path(initial['checkpoint']);before=sha(cp/'State.bin')
    summary=dict(status='started',runtime=identity,parent_checkpoint_sha256=before,
        governing_solves=0,physical_commits=0,restart_rejections={},stage_checks={})
    atomic(a.output/'summary.json',summary)
    def change_mesh(tree):
        count=0
        for v in tree.findall('./GEOMETRY/VERTEX/V'):
            xyz=list(map(float,v.text.split()))
            if abs(xyz[0]-.1)<1e-13:
                xyz[0]=.09;v.text=' '.join(format(x,'.17g') for x in xyz);count+=1
        require(count>1,'mesh control must move a complete tangential row, preserving affine panels')
    for label,mutate in [('affine-mesh',change_mesh),
        ('gamma-law',lambda t:set_parameter(t,'RANSWallRegularizationScale',.9))]:
        xml,_=compatibility_case(a.output/label,a.campaign/'input/finite-inflow.txt',a.campaign/'preflight.json',0,restart=cp)
        tree=ET.parse(xml).getroot();mutate(tree);write(tree,xml)
        env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(solver.parent/'lib'),TMPDIR=str(xml.parent))
        with (xml.parent/'solver.log').open('w') as log:
            code=subprocess.run([str(solver),str(xml)],cwd=xml.parent,env=env,
                stdout=log,stderr=subprocess.STDOUT,timeout=180).returncode
        q=load(xml.with_suffix('.evolution.json'))
        require(code!=0 and 'checkpoint signature mismatch' in q['failure'],
            'restart did not reject the changed interpretation/affine geometry at its signature guard')
        # Signature failure precedes integrator construction, so the minimal
        # startup-failure receipt can legitimately omit steps/history fields.
        require(q.get('steps',0)==0 and not q.get('pimple_attempts',[]),
            'incompatible restart advanced or invoked a governing solve')
        summary['restart_rejections'][label]=dict(failure=q['failure'],xml_sha256=sha(xml),
            report_sha256=sha(xml.with_suffix('.evolution.json')),exit_code=code)
    # The retained packets are immutable production output, not diagnostic re-solves.
    eps0,eps1=map(int,whole['offsets'][3:5])
    packets=sorted((a.campaign/'whole8/pimple-states').glob('*.json'))
    require(len(packets)==8,'one saved production packet per accepted step required')
    for path in packets:
        q=load(path)
        require(state_hash(q['pre_k'][eps0:eps1])==state_hash(q['corrected'][eps0:eps1]),
            'post-k path changed the solved outer epsilon coefficients')
        require(state_hash(q['raw'])==state_hash(q['corrected']),
            'inactive bound did not return the raw native vector unchanged')
        summary['stage_checks'][path.name]=dict(packet_sha256=sha(path),
            solved_epsilon_retained=True,inactive_bound_bitwise_noop=True)
    require(sha(cp/'State.bin')==before,'parent checkpoint bytes changed')
    require(initial['finite_inflow_reference']['samples']==whole['finite_inflow_reference']['samples'],
        'prescribed reference samples changed with evolving fields')
    require(state_hash(initial['fixed_reference'])==state_hash(whole['fixed_reference']) and
        state_hash(initial['fixed_defect'])==state_hash(whole['fixed_defect']),
        'fixed reference or defect changed')
    summary.update(status='passed',reference_samples_sha256=state_hash([
        v for row in initial['finite_inflow_reference']['samples'] for v in row]),
        reference_samples_bitwise_immutable=True,reference_defect_immutable=True,
        state_advanced_only_in_retained_production=True)
    if a.upgraded_solver:
        upgraded,new_identity=freeze(a.upgraded_solver.resolve(),a.output/'runtime-upgraded')
        reports=[]
        for label,image in [('baseline-one-step',solver),('label-upgrade-one-step',upgraded)]:
            xml,_=compatibility_case(a.output/label,a.campaign/'input/finite-inflow.txt',a.campaign/'preflight.json',1,restart=cp)
            env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(image.parent/'lib'),TMPDIR=str(xml.parent))
            with (xml.parent/'solver.log').open('w') as log:
                code=subprocess.run([str(image),str(xml)],cwd=xml.parent,env=env,
                    stdout=log,stderr=subprocess.STDOUT,timeout=180).returncode
            q=load(xml.with_suffix('.evolution.json'));require(code==0 and not q['failure'],'diagnostic label regression rejected')
            reports.append(q)
        verify_label_upgrade(*reports)
        summary['diagnostic_label_upgrade']=dict(status='passed',runtime=new_identity,stage_history_bitwise=True,
            audit_values_unchanged=True,physical_steps_per_independent_verification_leg=1,
            governing_solves_per_equation_across_verification_legs=2)
    atomic(a.output/'summary.json',summary);print(a.output/'summary.json')


if __name__=='__main__':main()
