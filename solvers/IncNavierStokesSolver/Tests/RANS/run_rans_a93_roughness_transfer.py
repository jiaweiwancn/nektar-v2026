#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One predeclared A93 mesh transfer; pinned runtime and sealed coarse reuse."""
import argparse,copy,gzip,json,os,shutil,subprocess,sys,time
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from qualify_rans_preservation import load,replay_request,Q
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import write_json
from rans_duration_contract import require,residual_reduction
from rans_preservation_oracle import Reference,weak_balance
from rans_spatial_resolution import analyze,common_test,test_parent_geometries,test_effective_coordinates
from diagnose_rans_balanced_spatial import run_probe
from run_a93_a95_e94e95_matrix import case_xml
from compare_rans_balanced_evolution import compare
from finish_rans_configuration_duration import finish
from sample_rans_common_endpoints import sample


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--ownership-root',type=Path,required=True)
    p.add_argument('--contract',type=Path,required=True);p.add_argument('--phase-contract',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();old=a.ownership_root.resolve();contract=load(a.contract)
    require(contract['schema']=='DirectionalResolutionExperimentV1','transfer_contract_version')
    phase=load(a.phase_contract)
    require(phase['schema']=='DirectionalInitializationPhaseV1' and phase['original_contract_sha256']==sha(a.contract) and
            phase['required_evolved_milestones']==[.016,.1,.5] and phase['evolved_maximum_ratio']==.85,
            'transfer_initialization_phase_or_evolved_gate_changed')
    require(contract['A93']['candidate_count_limit']==1 and contract['A93']['donor_node']==2.75 and
            contract['A93']['target_node']==.1375,'transfer_hypothesis_changed')
    source=Path(__file__).parent;tag='A93-N-Accepted-duration-0.002'
    base=load(old/'meshes/A93-N-Accepted/Balanced5V1.json')
    original=load(old/'meshes/A93-N-Accepted/Balanced11V1.json');mesh=copy.deepcopy(original)
    require(2.75 not in base['knots'] and .1375 not in original['knots'],'transfer_not_a_single_new_candidate')
    mesh['knots']=sorted((set(mesh['knots'])-{2.75})|{.1375})
    require(set(base['knots'])<=set(mesh['knots']) and set(mesh['retained_active_breaks'])<=set(mesh['knots']),
            'transfer_removed_required_node')
    mesh.update(schema='RoughnessTransferMeshV1',candidate='RoughnessTransferV1',donor_node=2.75,target_node=.1375,
                source_B11_sha256=sha(old/'meshes/A93-N-Accepted/Balanced11V1.json'),
                rationale='signed common rows 75/81 and varying-flux supports; one outer donor split to .1.. .175')
    for name in ('tmp','campaigns','comparisons','meshes','admissions'):(r/name).mkdir(exist_ok=True)
    for name in ('runtime','runtime-resolution-v3','runtime-budget-v5'):
        require(not (r/name).exists(),'transfer_destination_already_used')
        (r/name).symlink_to(old/name,target_is_directory=True)
    runtime=old/'runtime';analysis=(r/'runtime-resolution-v3').resolve()
    os.environ.update(LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
    newmesh=r/'meshes/RoughnessTransferV1.json';write_json(newmesh,mesh)
    execution=r/'execution';intent=dict(schema='A93RoughnessTransferExecutionV1',contract_sha256=sha(a.contract),phase_contract_sha256=sha(a.phase_contract),
        donor=2.75,target=.1375,mesh_sha256=sha(newmesh),expected_unknowns=2967,
        old_admission_sha256=sha(old/'admissions/A93-N-Accepted.json'),runtime=runtime_identity(runtime/'IncNavierStokesSolver'),
        times=[.016,.1,.5],timesteps=[.002,.001],one_candidate_only=True)
    seal=RequestManifest(execution,intent);events=[]
    def run(name,cmd):
        t=time.monotonic()
        with (execution/(name+'.log')).open('x') as f:c=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
        events.append(dict(name=name,command=list(map(str,cmd)),exit_code=c.returncode,wall_seconds=time.monotonic()-t))
        write_json(execution/'progress.json',dict(status='running',events=events))
        require(c.returncode==0,'transfer_subprocess_failed:'+name)
    def pilot(member,T,initialize=False):
        cmd=[sys.executable,'-B',str(source/'run_rans_neutral_pilot.py'),'--root',str(r),'--runtime-root',str(runtime),
            '--label',tag+'-'+member,'--model','A93','--policy','Natural','--wall','Accepted','--mesh-file',str(newmesh),
            '--dt','.001' if member.endswith('dt2') else '.002','--milestones','.016','.1','.5','--max-time',str(T),'--timeout','14400']
        if initialize:cmd+=['--initialize-only']
        if (r/'campaigns'/(tag+'-'+member)).exists():cmd+=['--resume']
        return cmd
    try:
        for m in ('B5','B5-dt2'):
            folder=old/'campaigns'/(tag+'-'+m);req=load(folder/'request-manifest.json')['request']
            for T in (.016,.1,.5):
                require(RequestManifest(folder,req,True).completed(folder/('T-'+format(T,'.4f'))),'transfer_coarse_unsealed')
                require((folder/('scalar-observations-'+str(T))/'results.json').exists(),'transfer_coarse_observer_missing')
            (r/'campaigns'/(tag+'-'+m)).symlink_to(folder,target_is_directory=True)
        run('initialize-new-fine',pilot('B11',.016,True))
        os.environ['LD_LIBRARY_PATH']=str(analysis/'lib')
        initial=load(r/'campaigns'/(tag+'-B11')/'initial-metrics.json')
        require(sum(initial['dof'])==2967,'transfer_actual_unknown_count')
        require(all(v['weighted_RMS']<=(1e-4 if f=='v' else 1e-3) for f,v in initial['accuracy']['reference_error'].items()),
                'transfer_initial_reference_accuracy')
        ad=load(old/'admissions/A93-N-Accepted.json');common=Path(ad['common_directory'])
        with gzip.open(common/'packet.json.gz','rt') as stream:cp=json.load(stream)
        ref=Reference('A93','Natural','Accepted');screens={}
        for m,knots in [('B5',base['knots']),('B11',mesh['knots'])]:
            folder=r/'initial-screen'/m;seed=load(r/'campaigns'/(tag+'-'+m)/'initial-packet.json')
            case_xml(folder,'A93','Natural','Accepted',steps=0,degree=5,quadrature=Q[1],dt=.002,
                descriptor_probe=runtime/'RANSNewModelsProbe',radial_knots=knots,boundary_policy='FrozenReferenceV1')
            packet,_=run_probe(analysis,folder,replay_request(seed,[],export_basis=True)|dict(export_resolution=True))
            cross=r/'initial-screen'/(m+'-common');cross.mkdir();shutil.copy2(common/'case.xml',cross/'case.xml')
            request=replay_request(cp,test_effective_coordinates(packet));request.update(export_test_basis=True,grid_test_geometries=test_parent_geometries(cp,packet))
            xp,_=run_probe(analysis,cross,request)
            screens[m]=dict(spatial=analyze(packet,ref),common_test=common_test(cp,packet,xp,ref),
                            analytical_weak_balance=weak_balance(ref,packet['observation']))
            require(screens[m]['analytical_weak_balance']['status']=='pass','transfer_initial_analytical_weak')
        gate=residual_reduction(screens['B5'],screens['B11'],'A93')
        write_json(r/'initial-screen/results.json',dict(status='pass',unknowns=2967,checks=screens,
            initial_residual_reduction_diagnostic=gate,phase_contract_sha256=sha(a.phase_contract),
            meaning='construction, analytical weak balance and common-test inclusion verified; evolved .85 gate applies at declared milestones'))
        ad['meshes']['B11']=dict(path=str(newmesh),sha256=sha(newmesh),specification=mesh)
        ad['duration_diagnostic_runtime']=runtime_identity(analysis/'RANSPreservationProbe')
        ad['mesh_transfer']=dict(contract_sha256=sha(a.contract),initial_screen_sha256=sha(r/'initial-screen/results.json'),
            original_failed_duration_sha256=sha(old/'duration-audits'/(tag+'-0.5')/'results.json'),
            coarse_reuse='exact sealed original requests and histories; no coarse advancement or receipt mutation')
        admission=r/'admissions/A93-N-Accepted.json';write_json(admission,ad)
        milestones=[]
        for T in (.016,.1,.5):
            for m in ('B11','B11-dt2'):run(m+'-'+str(T),pilot(m,T))
            labels=[tag+'-'+m for m in ('B5','B11','B5-dt2','B11-dt2')]
            samples=sample(r,labels,T,r/'common-sampling'/str(T),analysis)
            comparison=compare(r,labels,T,common_sampling=samples);write_json(r/'comparisons'/(tag+'-'+str(T)+'.json'),comparison)
            audit=finish(r,admission,.002,T)
            row=dict(time=T,status='pass' if comparison['status']=='pass' and audit['status']=='pass' else 'failed_or_unresolved',
                comparison=comparison,audit=audit)
            milestones.append(row);write_json(r/'results.json',dict(status=row['status'],milestones=milestones,initial_screen=gate,unknowns=2967))
            require(row['status']=='pass','transfer_milestone_gate:'+str(T))
        seal.seal(execution,['progress.json',r/'results.json']);return 0
    except Exception as e:
        write_json(execution/'failure.json',dict(status='failed_or_unresolved',reason=repr(e),events=events))
        if not (r/'results.json').exists():write_json(r/'results.json',dict(status='failed_or_unresolved',failure=repr(e),events=events))
        return 1


if __name__=='__main__':raise SystemExit(main())
