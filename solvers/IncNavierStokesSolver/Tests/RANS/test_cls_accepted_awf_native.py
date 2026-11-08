#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native accepted-AWF transactions and exact frozen-system accounting.

All trajectories here are verification-only. Production retains one solve per
equation; optional probe contact diagnostics request no k-only solves.
"""
import argparse
import os
from pathlib import Path
import platform
import xml.etree.ElementTree as ET

from compare_cls_epawf_pimple import same_history
from run_cls_accepted_awf import CASES,COMPANIONS,case_xml,check_report
from run_cls_epawf_picard import THREADS,load,require,run,write_json,state_hash,sha,set_info,set_parameter,write
from run_cls_epawf_pimple import identity
from continue_cls_epawf_pimple import atomic,header

_created_output=None


def main():
    global _created_output
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--probe',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case',choices=tuple(CASES),default='Accepted-AWF-2G')
    p.add_argument('--phase',choices=('all','constitutive'),default='all')
    p.add_argument('--companion',choices=COMPANIONS,default='faithful');p.add_argument('--companion-reason',default='')
    p.add_argument('--perturbation',type=float,default=0.);a=p.parse_args()
    require(platform.system()=='Linux','authoritative native tests run on Linux');os.environ.update(THREADS)
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    _created_output=a.output;a.solver=a.solver.resolve()
    summary=dict(status='started',case=a.case,companion=a.companion,solver=identity(a.solver),
        production_governing_solves=0,test_governing_solves_per_equation=0,ranks=1,threads=THREADS,
        steps_are='independent verification trajectories, not hidden production work')
    def save(**values):summary.update(values);atomic(a.output/'summary.json',summary)
    save()
    def native(name,steps,cp=None,reject=0,mutate=None,case=None,companion=None):
        directory=a.output/name;xml=case_xml(case or a.case,directory,steps,restart=cp,
            companion=companion or a.companion,companion_reason=a.companion_reason or 'verification-only incompatible model',
            perturbation=a.perturbation,hard_reject=reject)
        if mutate:
            tree=ET.parse(xml).getroot();mutate(tree);write(tree,xml)
        rc=run([str(a.solver),'case.xml'],directory,'solver.log',7200);write_json(directory/'receipt.json',rc)
        report=load(directory/'case.evolution.json') if (directory/'case.evolution.json').exists() else None
        return report,rc,directory
    zero,rc,initial=native('initial',0)
    if not zero or rc['exit_code'] or zero.get('failure'):
        save(status='initialization_rejected',failure=zero.get('failure') if zero else 'native exit without report',initial_receipt=rc)
        print(a.output/'summary.json');return
    definition=load(initial/'definition.json');check_report(zero,False,definition);cp0=Path(zero['checkpoint'])
    require(zero['implicit_callback_count']==0 and zero['initial_k_bound_check']['activations']==0,'startup modified/advanced native reference')
    cp0_hash=sha(cp0/'State.bin')
    mass=None
    if a.probe:
        a.probe=a.probe.resolve()
        def probe(name,operation,**options):
            folder=a.output/name;folder.mkdir()
            write_json(folder/'input.json',dict(operation=operation,state=zero['state'],history=zero['state'],
                time=0.,alpha=0.,require_history_h0=False,inner_stress_model=definition['model_choices']['inner_stress'],**options))
            receipt=run([str(a.probe),str(initial/'case.xml'),'input.json','result.json'],folder,'solver.log',1800)
            write_json(folder/'receipt.json',receipt);require(receipt['exit_code']==0,'read-only native constitutive check failed: '+name)
            return load(folder/'result.json')
        sampling=probe('endpoint-sampling','sampling_checks')
        require(sampling['governing_solves']==0 and max(sampling[k] for k in ('endpoint_sampling_error',
            'exact_global_row_error','nonuniform_profile_K_error','dependent_epsilon_H_error'))<1e-11,'native endpoint/dependent consumer check failed')
        save(endpoint_sampling=sampling)
        repaired=probe('negative-projection','k_bound',fixture='accepted_negative')
        b=repaired['k_bound'];require(b['success'] and b['activations']==1 and b['calls']==1 and b['governing_equation_solves']==0,
            'negative diagnostic candidate was not repaired by one logical algebraic operation')
        require(b['finite_pool_complete'] and b['primal_floor_defect']==0 and b['ground_floor_defect']==0,
            'negative fixture retains a physical floor/pool defect')
        require(repaired['non_k_change']==repaired['essential_k_change']==0 and repaired['shared_cg_jump']<1e-11,
            'negative repair changed non-k/essential/shared coefficients')
        offsets=[int(x) for x in repaired['offsets']];lo,hi=offsets[2],offsets[3]
        raw,corr=repaired['raw'],repaired['corrected']
        require(state_hash(raw[:lo]+raw[hi:])==state_hash(corr[:lo]+corr[hi:]),'negative repair changed non-k bits')
        prescribed=[int(x[0]) for x in zero['reference_constraints'] if lo<=x[0]<hi]
        require(prescribed and state_hash([raw[j] for j in prescribed])==state_hash([corr[j] for j in prescribed]),
            'negative repair changed prescribed-k bits')
        require(b['stationarity_mass_dual']<1e-9*max(b['mass_norm_correction'],1e-30) and
            b['dual_feasibility_defect']<1e-12 and abs(b['complementarity'])<1e-9*max(b['primal_objective'],1e-30),
            'negative-fixture scaled KKT checks failed')
        save(active_projection_diagnostic=dict(success=True,native_step_accepted=False,governing_solves=0,
            complete_hard_admissible=repaired['complete_hard_admissible'],complete_hard_reasons=repaired['complete_hard_reasons'],
            complete_audit_failure=repaired['complete_audit_failure'],raw_minimum=b['raw_minimum'],corrected_minimum=b['corrected_minimum'],
            non_k_exact=True,prescribed_k_exact=True,shared_cg_jump=repaired['shared_cg_jump'],mass_norm=b['mass_norm_correction'],
            stationarity_mass_dual=b['stationarity_mass_dual'],complementarity=b['complementarity'],
            interpretation='QP positivity feasibility is separate from linear-stress covariance/hard-audit qualification'))
        mass=probe('mass-accounting','contact',h=[],export_mass=True)
        require(mass['diagnostic_k_solves']==0 and mass['production_governing_solves']==0,
            'consistent-mass readback performed a transport solve')
        save(same_state_native_frozen_load_agreement=mass['native_frozen_load_agreement'],diagnostic_k_solves=0)
    if a.phase=='constitutive':
        require(a.probe is not None,'constitutive phase needs the native probe')
        save(status='passed',accepted_checkpoint=header(cp0));print(a.output/'summary.json');return
    whole,rc,folder=native('whole2',2,cp0);require(rc['exit_code']==0 and not whole['failure'],'two-step reference trajectory failed')
    check_report(whole,True,definition)
    first,rc,folder=native('split1',1,cp0);require(rc['exit_code']==0 and not first['failure'],'split first step failed')
    check_report(first,True,definition);cp1=Path(first['checkpoint']);cp1_hash=sha(cp1/'State.bin')
    second,rc,folder=native('restart1',1,cp1);require(rc['exit_code']==0 and not second['failure'],'compatible restart failed')
    check_report(second,True,definition);save(restart_exact=same_history(whole,second),test_governing_solves_per_equation=4)
    startup,rc,folder=native('startup1',0,cp1);require(rc['exit_code']==0 and not startup['failure'],'accepted checkpoint startup failed')
    rejected,rc,folder=native('forced-rejection2',1,cp1,reject=2)
    require(rejected and 'deliberate' in rejected['failure'].lower(),'forced hard rejection was not exercised')
    check_report(rejected,True,definition)
    require(rejected['pimple_attempts'][-1]['k_transport_solves']==1 and rejected['pimple_attempts'][-1]['bound_calls']==1,
        'hard rejection did not follow the scheduled k solve and logical correction')
    save(rollback_exact=same_history(startup,rejected),test_governing_solves_per_equation=5)
    mismatches={}
    other=next(c for c,(H,partition) in CASES.items() if H==CASES[a.case][0] and partition!=CASES[a.case][1])
    checks=[('coefficient',lambda t:set_parameter(t,'Kinvis',1.51e-5),None,None),
            ('sampler',lambda t:set_info(t,'RANSWallKSampling','ColumnAverage'),None,None),
            ('tail',lambda t:set_parameter(t,'RANSWallTailStart',.9*CASES[a.case][0]),None,None),
            ('mesh',None,other,None),
            ('outer_model',None,None,'faithful' if a.companion!='faithful' else 'regularized')]
    for name,mutate,case,companion in checks:
        d,rc,directory=native('incompatible-'+name,0,cp1,mutate=mutate,case=case,companion=companion)
        require(rc['exit_code']!=0 or (d and d.get('failure')),'incompatible checkpoint accepted: '+name)
        mismatches[name]=dict(exit_code=rc['exit_code'],failure=d.get('failure') if d else None,
            log=str(directory/'solver.log'),log_sha256=sha(directory/'solver.log'),
            checkpoint_unchanged=sha(cp1/'State.bin')==cp1_hash,
            interpretation='explicit model/configuration guard' if name in ('sampler','tail') else 'native physical/signature compatibility')
    save(incompatible_cases=mismatches)
    require(sha(cp0/'State.bin')==cp0_hash and sha(cp1/'State.bin')==cp1_hash,'native checkpoint bytes changed')
    invariants=[]
    for t in whole['pimple_attempts']:
        q=load(Path(t['k_snapshot']));offset=[int(x) for x in whole['offsets']];lo,hi=offset[2],offset[3]
        raw,corrected=q['raw'],q['corrected'];non_k=state_hash(raw[:lo]+raw[hi:])==state_hash(corrected[:lo]+corrected[hi:])
        fixed=[int(r[0]) for r in whole['reference_constraints'] if lo<=r[0]<hi]
        require(fixed,'prescribed-k invariance fixture has no prescribed coefficients')
        essential=state_hash([raw[j] for j in fixed])==state_hash([corrected[j] for j in fixed])
        require(non_k and essential,'post-k repair changed non-k/prescribed coefficients')
        invariants.append(dict(time=t['time'],non_k_exact=non_k,prescribed_k_exact=essential,prescribed_k_count=len(fixed),
            bound_activation=t['k_bound']['activations']))
    save(corrector_invariants=invariants)
    if mass is not None:
        import numpy as np
        rows=np.array(mass['free_k_rows'],dtype=int);M=np.array(mass['free_k_mass']);accounting=[]
        for t in whole['pimple_attempts']:
            q=load(Path(t['k_snapshot']));require(np.array_equal(rows,np.array(q['transport_rows'],dtype=int)),'frozen/free mass rows differ')
            require(sha(Path(q['matrix_file']))==q['matrix_sha256'],'saved transport matrix changed')
            A=np.fromfile(q['matrix_file'],dtype='<f8').reshape(M.shape);r=np.array(q['transport_residual']);d=np.array(q['raw_increment'])
            delta=(np.array(q['corrected'])-np.array(q['raw']))[rows];linear=r+A@d;added=A@delta;storage=q['alpha']*(M@delta)
            corrected=r+A@(d+delta);identity_defect=float(max(abs(corrected-added-linear)))
            require(identity_defect<1e-11 and max(abs(linear))<1e-9,'raw/corrected frozen residual identity failed')
            require(abs(float(max(abs(corrected)))-t['k_corrected_frozen_linear_residual'])<1e-10,'corrected residual accounting differs')
            require(abs(float(max(abs(storage)))-t['k_added_storage_residual'])<1e-10,'physical consistent storage accounting differs')
            accounting.append(dict(time=t['time'],raw_residual=float(max(abs(linear))),corrected_residual=float(max(abs(corrected))),
                added_A_delta=float(max(abs(added))),storage_alpha_M_delta=float(max(abs(storage))),identity_defect=identity_defect))
        save(same_state_native_frozen_load_agreement=mass['native_frozen_load_agreement'],frozen_storage_accounting=accounting,
            diagnostic_k_solves=0,projection_mass_definition='native consistent PicardKMass Mff; transport A is distinct')
    save(status='passed',accepted_checkpoint=header(Path(second['checkpoint'])))
    print(a.output/'summary.json')


if __name__=='__main__':
    try:main()
    except Exception as error:
        # A focused-check failure must not leave its progress labelled started.
        # Preserve every candidate/result; only publish the test's terminal status.
        if _created_output is not None:
            output=_created_output
            summary=output/'summary.json'
            previous=load(summary) if summary.exists() else {}
            previous.update(status='test_failed',failure=str(error))
            atomic(summary,previous)
        raise
