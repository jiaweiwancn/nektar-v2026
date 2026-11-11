#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded native finite-profile compatibility and checkpoint qualification.

Never launches mapped production or an experimental precursor. Numerical
rejection is preserved, not retried with modified physics or positivity gates.
"""
import argparse
from datetime import datetime,timezone
import os
from pathlib import Path
import platform
import subprocess
import time
import xml.etree.ElementTree as ET

from cao_tamura_budget_h1 import make_profile,compatibility_case
from continue_cls_epawf_pimple import CampaignLock,atomic,freeze,header
from compare_cls_epawf_pimple import same_history
from run_cls_accepted_awf import check_report
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,set_parameter,write


_output=None


def main():
    global _output
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('solver','probe','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--timeout',type=float,default=1800.)
    a=p.parse_args();require(platform.system()=='Linux','authoritative Linux native execution required')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);_output=a.output;lock=CampaignLock(a.output)
    os.environ.update(THREADS)
    for key in ('VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
    solver,identity=freeze(a.solver.resolve(),a.output/'runtime-solver')
    probe,probe_identity=freeze(a.probe.resolve(),a.output/'runtime-probe')
    start=time.monotonic();summary=dict(version='CaoTamuraBudgetH1NativeV1',status='started',
        started_UTC=datetime.now(timezone.utc).isoformat(),solver=identity,probe=probe_identity,
        rank=1,threads=dict(THREADS),physical_scope='affine finite-profile compatibility, prescribed sides/top; not developed precursor',
        declared=dict(dt=1e-5,accepted_steps=8,endpoint=.00008,maximum_steps=16,
            target='integrated one-sweep/provider/restart contracts, not experimental fit or equilibrium'),tests={})
    def save(**values):summary.update(values);atomic(a.output/'summary.json',summary)
    save()
    profile=a.output/'input';manifest=make_profile(profile)
    env=dict(os.environ,LD_LIBRARY_PATH=str(probe.parent/'lib'))
    preflight=a.output/'kernels.json'
    cmd=[str(probe),'--profile',str(profile/'finite-inflow.txt'),'--queries',str(profile/'permitted-queries.txt'),'--output',str(preflight)]
    rc=subprocess.run(cmd,cwd=a.output,env=env,stdout=(a.output/'kernels.log').open('w'),stderr=subprocess.STDOUT,timeout=180).returncode
    require(rc==0 and load(preflight)['status']=='passed','native H1 geometry/reference/preflight checks failed')
    result=load(preflight)['preflight'];atomic(a.output/'preflight.json',result)
    require(next(c for c in result['candidates'] if c['Hw']==result['selected']['Hw'])['all_declared_K_curvature_budget_cases_admissible'],'selected candidate fails explicit cutoff budget gate')
    require(result['finite_profile_covariance_minimum']>0,'fixed completed inlet fails covariance preflight')
    manifest['geometry']['h_eff_m']=result['selected']['h_eff'];manifest['selected_preflight']=result['selected']
    atomic(profile/'manifest.json',manifest)
    save(status='initializing',preflight_sha256=sha(preflight),profile_sha256=sha(profile/'finite-inflow.txt'),manifest_sha256=sha(profile/'manifest.json'))

    def affine_mesh_change(tree,folder):
        # Same-size, still-affine geometry: move the complete interior x row.
        # A single moved vertex would test the unsupported-map guard instead
        # of ordinary restart's strict geometry identity.
        for v in tree.findall('./GEOMETRY/VERTEX/V'):
            xyz=list(map(float,v.text.split()))
            if abs(xyz[0]-.1)<1e-13:
                xyz[0]=.09;v.text=' '.join(format(x,'.17g') for x in xyz)

    def native(name,steps,cp=None,reject=0,mutate=None,expect=True):
        xml,definition=compatibility_case(a.output/name,profile/'finite-inflow.txt',a.output/'preflight.json',steps,restart=cp,reject=reject)
        if mutate:
            t=ET.parse(xml).getroot();mutate(t,xml.parent);write(t,xml)
        env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),TMPDIR=str(xml.parent))
        t0=time.monotonic()
        with (xml.parent/'solver.log').open('w') as log:
            child=subprocess.Popen(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(solver),str(xml)],cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,pass_fds=(lock.stream.fileno(),))
            try:code=child.wait(timeout=a.timeout)
            except subprocess.TimeoutExpired:
                save(status='infrastructure_interrupted',native_pid=child.pid,folder=str(xml.parent));raise
        report=xml.parent/'case.evolution.json';q=load(report) if report.exists() else None
        receipt=dict(exit_code=code,seconds=time.monotonic()-t0,usage=(xml.parent/'usage.txt').read_text().strip(),
            xml_sha256=sha(xml),report_sha256=sha(report) if report.exists() else None,log_sha256=sha(xml.parent/'solver.log'),
            failure=q.get('failure') if q else None)
        atomic(xml.parent/'receipt.json',receipt);summary['tests'][name]=receipt;save()
        if expect:require(code==0 and q and not q.get('failure'),'native compatibility rejected: '+name)
        else:require(code!=0 or (q and q.get('failure')),'expected rejection missing: '+name)
        if expect:
            check_report(q,True,None)
            require(q['finite_inflow_reference']['mode']=='FiniteInflowV1' and q['evaluation_reuse']['mode']=='ExactStepV1','provider/reuse mode differs')
            require(q['audit_revision']=='NativeElementSideV2','audit side provenance differs')
        return q,receipt,definition

    try:
        zero,rz,definition=native('initial',0);cp=Path(zero['checkpoint']);base_hash=sha(cp/'State.bin')
        require(len(zero['state'])==definition['native_size_estimate'],'native dimensions differ')
        whole,rw,_=native('whole8',8,cp)
        first,_,_=native('split4',4,cp);second,_,_=native('restart4',4,Path(first['checkpoint']))
        require(same_history(whole,second),'split/uninterrupted native histories differ')
        require(whole['signature']==zero['signature'],'immutable profile identity changed during evolution')
        for x in (zero,first,second,whole):
            require(x['finite_inflow_reference']['identity']==zero['finite_inflow_reference']['identity'],'finite reference modified by evolving fields')
            require(x['finite_inflow_reference']['stress_work_relative_defect']<1e-12 and x['finite_inflow_reference']['epsilon_production_relative_defect']<1e-12,'stress-work production inconsistent')
        startup,_,_=native('startup0',0,cp)
        rejected,rr,_=native('forced-rejection1',1,cp,reject=1,expect=False)
        require(rejected and 'deliberate' in rejected['failure'].lower(),'forced rejection did not reach hard gate')
        require(same_history(startup,rejected),'rejection failed to restore accepted state/time/histories')
        require(sha(cp/'State.bin')==base_hash,'parent checkpoint bytes changed')
        def change_profile(tree,folder):
            f=folder/'changed-reference.txt';lines=(profile/'finite-inflow.txt').read_text().splitlines();q=lines[1].split();q[1]='1e-12';lines[1]=' '.join(q)
            f.write_text('\n'.join(lines)+'\n');set_info(tree,'RANSWallReferenceFile',str(f));set_info(tree,'RANSWallReferenceSHA256',sha(f))
        for name,mutate in [
            ('profile-content',change_profile),
            ('profile-digest',lambda t,f:set_info(t,'RANSWallReferenceSHA256','0'*64)),
            ('frame',lambda t,f:set_parameter(t,'RANSWallReferenceOriginZ',-1e-5)),
            ('timestep',lambda t,f:set_parameter(t,'TimeStep',5e-6)),
            ('gamma-law',lambda t,f:set_parameter(t,'RANSWallRegularizationScale',.9)),
            ('mesh',affine_mesh_change)]:
            q,rc,_=native('incompatible-'+name,0,cp,mutate=mutate,expect=False)
        attempts=whole['pimple_attempts'];require(len(attempts)==8 and whole['implicit_callback_count']==8,'implicit callback count differs')
        budget=[t['k_bound'].get('integral_total',0) for t in attempts]
        activations=sum(t['k_bound']['activations'] for t in attempts)
        history=header(Path(whole['checkpoint']),True)
        save(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat(),seconds=time.monotonic()-start,
            accepted_physical_steps=whole['steps'],physical_time=whole['time'],native_dimensions=len(zero['state']),
            initialization_work=zero['initialization_lift_calls_iterations'],solve_counts={k:sum(t[k] for t in attempts) for k in
                ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves')},
            correction=dict(activations=activations,signed_integral=sum(budget),positive_integral=sum(t['k_bound']['integrals_inner_outer_near_positive_negative'][3] for t in attempts),negative_integral=sum(t['k_bound']['integrals_inner_outer_near_positive_negative'][4] for t in attempts),absolute_integral=sum(abs(t['k_bound']['integrals_inner_outer_near_positive_negative'][3])+abs(t['k_bound']['integrals_inner_outer_near_positive_negative'][4]) for t in attempts)),
            restart_bitwise=True,rollback_bitwise=True,reference_immutable=True,
            initial_checkpoint_sha256=base_hash,final_checkpoint_sha256=history['sha256'],
            initial_state_sha256=state_hash(zero['state']),final_state_sha256=state_hash(whole['state']),
            final_native_history_sha256=[state_hash(v) for v in whole['native_history']],
            final_native_history_times=whole['native_history_times'],finite_reference_checks=whole['finite_inflow_reference'],
            stage_audits=whole.get('stage_audits'),state_audits=whole.get('state_audits'),attempt_metrics=[{k:v for k,v in t.items() if k not in ('k_bound',)} for t in attempts],
            performance=dict(whole8_seconds=rw['seconds'],usage=rw['usage']))
    except Exception as error:
        if summary['status']!='infrastructure_interrupted':save(status='native_validation_stopped',failure=str(error),seconds=time.monotonic()-start)
        raise
    print(a.output/'summary.json')


if __name__=='__main__':
    try:main()
    except Exception as error:
        if _output and (_output/'summary.json').exists():
            q=load(_output/'summary.json')
            if q.get('status') in ('started','initializing'):
                q.update(status='preparation_rejected',failure=str(error),completed_UTC=datetime.now(timezone.utc).isoformat())
                atomic(_output/'summary.json',q)
        raise
