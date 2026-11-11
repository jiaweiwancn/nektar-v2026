#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One admitted A95 same-vertical-space experiment, capped at T=.016."""
import argparse,copy,os,resource,subprocess,sys,time
from pathlib import Path
from qualify_rans_preservation import load,replay_request,Q
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import THREADS,write_json,exact_vector_hash,compare_samples
from rans_duration_contract import require,MEMBERS,identity
from rans_preservation_oracle import Reference,weak_balance
from diagnose_rans_balanced_spatial import run_probe
from run_rans_neutral_pilot import scaled,FIELDS
from compare_rans_balanced_evolution import compare


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--old-root',type=Path,required=True)
    p.add_argument('--contract',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    r=a.root.resolve();old=a.old_root.resolve();s=Path(__file__).parent;runtime=r/'runtime';cfg=dict(model='A95',policy='Natural',wall='Accepted')
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'),PYTHONDONTWRITEBYTECODE='1',NEKTAR_RANS_PROBE_TIMEOUT='14400')
    contract=load(a.contract);require(contract['schema']=='DirectionalResolutionExperimentV1','directional_contract_version')
    prerequisites=['directional-layout-v2/results.json','isotropic-compatibility-v1/results.json','directional-history-v3/results.json',
        'directional-geometry-v1/results.json','frozen-scalar-reconciliation-v1/results.json']
    proofs={}
    for name in prerequisites:
        data=load(r/name);require(data['status']=='pass','directional_prerequisite_failed:'+name);proofs[name]=sha(r/name)
    regression=load(r/'compatibility-v1-execution.json')
    require(len(regression['events'])==4 and all(regression['events'][i]['exit_code']==0 for i in (0,2,3)),
            'directional_shared_regressions_failed')
    # Event 1's unequal-point MatrixFree startup failure remains archived. The
    # new separately sealed supplement recomputes target-layout primitives and
    # completes signature/CLS tests; it does not relabel that failed job.
    proofs['compatibility-v1-execution.json']=sha(r/'compatibility-v1-execution.json')
    require(load(r/'directional-history-v3/results.json')['actual_unknowns']==2769,'directional_native_count_unverified')
    for name in ('runtime-resolution-v3','runtime-budget-v5'):
        if not (r/name).exists():(r/name).symlink_to(runtime,target_is_directory=True)
        require((r/name).resolve()==runtime,'directional_runtime_alias_foreign')
    binding=old/'physical-budget-bindings/A95-N-Accepted-response-2e-05-dt-0.002-0.016/results.json'
    previous=load(binding);require(previous['reconstruction_status']=='pass' and previous['accuracy_status']=='failed_or_unresolved','historical_budget_status_changed')
    neutral_proof=old/'duration-audits/A95-N-Accepted-duration-0.002-0.5/results.json'
    require(load(neutral_proof)['status']=='pass','historical_A95_neutral_not_qualified')
    oldad=load(old/'admissions/A95-N-Accepted.json');common=Path(oldad['common_directory'])
    require(sha(common/'packet.json.gz')==oldad['common_packet_sha256'],'directional_common_packet_changed')
    request=dict(schema='DirectionalResolutionCohortV1',configuration=cfg,operator='GeneratedVolumeOwnerV1',
        contract_sha256=sha(a.contract),proofs=proofs,runtime={n:runtime_identity(runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},
        members=dict(B5='same-runtime isotropic B5 (5,5)/(3,3)',B11='same B5 vertical mesh (6,5)/(4,3)',
                     **{'B5-dt2':'independently initialized isotropic half-step','B11-dt2':'independently initialized directional half-step'}),
        times=[.016],dt=[.002,.001],amplitudes=[0.,2e-5],history='fresh native IMEX1 for each member and amplitude',
        historical_neutral_proof_sha256=sha(neutral_proof),historical_failed_budget_binding_sha256=sha(binding),
        common_packet_sha256=sha(common/'packet.json.gz'),common_reuse='same GeneratedVolumeOwnerV1 isotropic p3 tests; exact old/new frozen compatibility proof',
        gates=dict(neutral_absolute_fields=1e-3,neutral_vertical=1e-4,traction=1e-3,neutral_timestep_floor=1e-7,neutral_timestep_fraction=.1,
            response_spatial_floor=5e-7,response_spatial_fraction=.1,response_temporal_floor=1e-7,response_temporal_fraction=.05,
            budget_spatial_ratio=.85,budget_floor=5e-9,budget_quadrature_multiplier=5,budget_timestep_floor=1e-9,budget_timestep_fraction=.1),
        phase_contract='same vertical space: neutral absolute accuracy, actual residuals, weak/domain/conservation and timestep bounds; no new neutral analytical-error ratio claim. The .85 total-budget gate and historical balanced-space decisions remain unchanged.',
        extension_authorized=False,source_sha256=sha(Path(__file__)))
    out=r/'directional-resolution-v1';seal=RequestManifest(out,request,a.resume)
    if seal.completed(out):return 0 if load(out/'results.json')['status']=='pass' else 1
    events=[];failures=[];receipt_files=[]
    def run(name,cmd,required=True):
        log=out/(name+'.log');cost=out/(name+'-execution.json')
        if cost.exists():
            d=load(cost)
            require([x for x in d['command'] if x!='--resume']==[str(x) for x in cmd if str(x)!='--resume'],
                    'directional_resume_command_mismatch:'+name);code=d['exit_code']
        else:
            before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
            with log.open('x') as stream:code=subprocess.run(list(map(str,cmd)),stdout=stream,stderr=subprocess.STDOUT).returncode
            after=resource.getrusage(resource.RUSAGE_CHILDREN)
            write_json(cost,dict(command=list(map(str,cmd)),exit_code=code,wall_seconds=time.monotonic()-start,
                process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
        events.append(dict(name=name,exit_code=code,cost_sha256=sha(cost)))
        write_json(out/'progress.json',dict(events=events,failures=failures))
        if code:failures.append(name)
        if required:require(code==0,'directional_execution_or_required_gate:'+name)
        return code
    try:
        meshbase=load(r/'meshes/A95-N-Accepted/Balanced5V1.json');meshes={};hypotheses={}
        for control,orders,n in [('StreamwiseBaseV1',[5,5],2355),('StreamwiseP6V1',[6,5],2769)]:
            spec=dict(meshbase,control=control,field_degrees=orders,pressure_degrees=[v-2 for v in orders])
            mp=out/(control+'-mesh.json');hp=out/(control+'-hypothesis.json')
            if not mp.exists():write_json(mp,spec)
            require(load(mp)==spec,'directional_declared_mesh_changed')
            source=old/'physical-budgets/A95-N-Accepted-response-2e-05-dt-0.002-perturbed-B5-0.016-Q2-F32/results.json'
            hypothesis=dict(schema='PhysicalBudgetAxisHypothesisV1',configuration=cfg,control=control,mesh_sha256=sha(mp),
                base_mesh_sha256=sha(r/'meshes/A95-N-Accepted/Balanced5V1.json'),source_primary_budget=str(source),source_primary_budget_sha256=sha(source),
                expected_unknowns=n,hypothesis='change only streamwise field/compatible pressure order at fixed pz and physical parameters; preserve total baseline and perturbed budget classifications')
            if not hp.exists():write_json(hp,hypothesis)
            require(load(hp)==hypothesis,'directional_hypothesis_changed');meshes[control]=mp;hypotheses[control]=hp
        tag='A95-N-Accepted-directional-dt-0.002';base=[tag+'-neutral-'+m for m in MEMBERS];pert=[tag+'-response-'+m for m in MEMBERS]
        def command(label,i,amp,initial=False):
            control='StreamwiseBaseV1' if i in (0,2) else 'StreamwiseP6V1'
            cmd=[sys.executable,'-B',s/'run_rans_budget_axis_control.py','--root',r,'--binding',binding,'--hypothesis',hypotheses[control],
                '--mesh-file',meshes[control],'--control',control,'--label',label,'--dt',.001 if i>=2 else .002,'--amplitude',amp]
            if (r/'campaigns'/label).exists():cmd+=['--resume']
            if initial:cmd+=['--initialize-only']
            return cmd
        for amp,labels in ((0.,base),(2e-5,pert)):
            for i,label in enumerate(labels):run('initialize-'+label,command(label,i,amp,True))
            q=[load(r/'campaigns'/label/'request-manifest.json')['request'] for label in labels]
            for i,j in ((0,2),(1,3)):
                for key in ('raw_neutral_sha256','reference_sha256','boundary','spatial_control'):
                    require(q[i][key]==q[j][key],'directional_same_space_initial_identity:'+key)
                if amp:require(q[i]['actual_initial']['state_sha256']==q[j]['actual_initial']['state_sha256'],'directional_actual_perturbation_differs')
        ad=dict(configuration=dict(cfg,descriptor=load(r/'campaigns'/base[0]/'request-manifest.json')['request']['descriptor']),
            common_directory=str(common),common_packet_sha256=sha(common/'packet.json.gz'),duration_diagnostic_runtime=runtime_identity(runtime/'RANSPreservationProbe'),
            meshes={m:dict(path=str(meshes[c]),sha256=sha(meshes[c]),specification=load(meshes[c])) for m,c in [('B5','StreamwiseBaseV1'),('B11','StreamwiseP6V1')]},
            space_names=['B5-IsotropicP5','B5-UniformX6Z5'])
        admission=out/'admission.json'
        if not admission.exists():write_json(admission,ad)
        require(load(admission)==ad,'directional_admission_changed')
        # Both actual initial states are audited before any trajectory advances.
        initial_checks=[]
        for i in (0,1):
            for amp,labels in ((0.,base),(2e-5,pert)):
                label=labels[i];q=load(r/'campaigns'/label/'request-manifest.json')['request'];dest=out/('initial-screen-'+label)
                if not (dest/'screen.json').exists():
                    dest.mkdir()
                    from run_cls_epawf_picard import write
                    # The immutable initialization request binds a normalized XML;
                    # use the matching actual initializer where supplied, otherwise
                    # regenerate the same declared zero-step configuration.
                    if amp:__import__('shutil').copy2(r/'campaigns'/label/'actual-initial/case.xml',dest/'case.xml')
                    else:
                        import run_rans_neutral_pilot as pilot
                        from types import SimpleNamespace
                        args=SimpleNamespace(root=r,runtime_root=runtime,model='A95',policy='Natural',wall='Accepted',dt=q['dt'],degree=5,
                            mesh_file=meshes['StreamwiseBaseV1' if i==0 else 'StreamwiseP6V1'],quadrature_index=1,cache='ExactStepV1',columns=1,velocity_amplitude=0.,level=0)
                        target=dest/'generated';pilot.configure(target,args,Reference('A95','Natural','Accepted'),0)
                        import xml.etree.ElementTree as ET
                        xml=ET.parse(target/'case.xml').getroot()
                        if i==1:
                            for e in xml.findall('./EXPANSIONS/E'):e.set('NUMMODES','5,4' if e.get('FIELDS')=='p' else '7,6')
                        write(xml,dest/'case.xml')
                    seed=load(r/'campaigns'/label/'initial-packet.json');state=q['actual_initial']['state_coefficients'] if amp else seed['reference']
                    req=replay_request(seed,[],state=state);req.update(export_basis=True,export_resolution=True)
                    packet,_=run_probe(runtime,dest,req);weak=weak_balance(Reference('A95','Natural','Accepted'),packet['observation'])
                    require(packet['physical_admissibility']=='pass' and weak['status']=='pass','directional_initial_domain_or_weak_failure')
                    write_json(dest/'screen.json',dict(status='pass',weak=weak,actual_state_sha256=exact_vector_hash(state),unknowns=sum(packet['observation']['sizes']),
                        initialization=packet['initialization_receipt'],packet_sha256=sha(dest/'packet.json.gz')))
                initial_checks.append(load(dest/'screen.json'))
        for i,label in enumerate(base):run('evolve-'+label,command(label,i,0.))
        bundle=compare(r,base,.016,inputs_only=True);ref=Reference('A95','Natural','Accepted');w=bundle['samples'][0]['grid_weights']
        differences=[scaled(compare_samples(bundle['samples'][x],bundle['samples'][y],w),ref) for x,y in ((0,1),(0,2),(1,3))]
        time_checks={f:dict(spatial_difference=differences[0][f]['weighted_RMS'],timestep_difference=max(d[f]['weighted_RMS'] for d in differences[1:]),
            limit=max(1e-7,.1*differences[0][f]['weighted_RMS'])) for f in FIELDS}
        require(all(v['timestep_difference']<=v['limit'] for v in time_checks.values()),'directional_neutral_temporal_uncertainty')
        write_json(out/'neutral-comparison.json',dict(status='pass',fields=time_checks,provenance=bundle['provenance'],
            analytical_errors=[m['reference_error'] for m in bundle['metrics']],actual_initial_drifts=[m['change_from_initial'] for m in bundle['metrics']]))
        run('neutral-audit',[sys.executable,'-B',s/'audit_rans_response.py','--root',r,'--admission',admission,'--cohort-prefix',tag+'-neutral','--time','.016'])
        for i,label in enumerate(pert):run('evolve-'+label,command(label,i,2e-5))
        run('response-comparison',[sys.executable,'-B',s/'compare_rans_nonlinear_response.py','--root',r,'--perturbed',*pert,'--baseline',*base,'--time','.016',
            '--output-directory',out/'response-comparison'],False)
        run('response-audit',[sys.executable,'-B',s/'audit_rans_response.py','--root',r,'--admission',admission,'--cohort-prefix',tag+'-response','--time','.016'],False)
        # Complete independent budgets after a failed response check as well.
        for kind,labels in (('baseline',base),('perturbed',pert)):
            for i,(m,label) in enumerate(zip(MEMBERS,labels,strict=True)):
                for q,face in (((1,16),(2,16),(2,32)) if i<2 else ((2,32),)):
                    output=tag+'-'+kind+'-'+m+'-0.016-Q'+str(q)+('-F32' if face==32 else '')
                    cmd=[sys.executable,'-B',s/'audit_rans_physical_scalar_budget.py','--root',r,'--label',label,'--time','.016','--output-label',output,
                        '--quadrature-index',q,'--face-order',face,'--runtime-name','runtime']
                    if (r/'physical-budgets'/output).exists():cmd+=['--resume']
                    run('budget-'+kind+'-'+m+'-'+str(q)+'-'+str(face),cmd,False)
        analysis=out/'budget-analysis'
        run('budget-analysis',[sys.executable,'-B',s/'analyze_rans_physical_budgets.py','--root',r,'--cohort',tag,'--time','.016','--output',analysis],False)
        run('budget-accuracy',[sys.executable,'-B',s/'classify_rans_physical_budget_accuracy.py','--source',analysis/'results.json','--output',out/'budget-accuracy'],False)
        run('budget-binding',[sys.executable,'-B',s/'bind_rans_ownership_budget_accuracy.py','--root',r,'--admission',admission,
            '--cohort',tag+'-response','--baseline-cohort',tag+'-neutral','--dt','.002','--time','.016',
            '--analysis',analysis/'results.json','--accuracy',out/'budget-accuracy/results.json','--output',out/'budget-binding'],False)
        result=dict(status='pass' if not failures else 'failed_or_unresolved',configuration=cfg,time=.016,events=events,failures=failures,
            response_comparison=load(out/'response-comparison/results.json'),physical_budget_analysis=load(analysis/'results.json'),
            physical_budget_accuracy=load(out/'budget-accuracy/results.json'),initial_checks=initial_checks,
            physical_budget_binding=load(out/'budget-binding/results.json'),
            historical_accuracy_status_unchanged='failed_or_unresolved',response_extension_eligible=False,scope=request['phase_contract'])
    except Exception as e:result=dict(status='failed_or_unresolved',failure=repr(e),configuration=cfg,time=.016,events=events,failures=failures,
        response_extension_eligible=False,scope='stopped at the concrete required gate; independent completed evidence preserved')
    write_json(out/'results.json',result)
    seal.seal(out,[p for p in out.rglob('*') if p.is_file() and p.name not in ('request-manifest.json','completed-artifacts.json')])
    print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
