#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal the completed H27 feasibility account and publish honest unmet gates.

Run only after every admitted operation has a receipt. The explicit 300-second
control/publication provision covers this script and non-numerical transport
and audit work outside individually timed operations; it is not production
cost. No further science may use the closed account.
"""
import argparse,datetime,importlib.metadata,json,os,platform,shutil,subprocess,sys,tarfile,time
from pathlib import Path
from resources import atomic,identity,update,workers,size

A_COMMIT='529eb9337bde9eb1be3606160f0ee1617ca4b54a'
HORIZON=505.02466521278404
TERMINAL=[447.8483832973273,HORIZON]

def records(value):
    if isinstance(value,dict):
        if 'path'in value and'sha256'in value:yield value
        else:
            for v in value.values():yield from records(v)
    elif isinstance(value,list):
        for v in value:yield from records(v)

def main():
    started=time.monotonic();p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True);p.add_argument('--mac-binding',type=Path,required=True)
    a=p.parse_args();root=a.root;out=root/'publication-feasibility-v1'
    assert not out.exists();book=root/'accounting.json'
    account=json.loads(book.read_text());assert not account['closed']
    assert all(v['status']!='RESERVED'for v in account['entries'].values())
    for label in account['entries']:
        if label!='bootstrap-and-control':assert(root/'operations'/label/'receipt.json').is_file(),label
    takeover=json.loads((root/'takeover.json').read_text());old=[]
    for pin in records(takeover['old_closed_accounts']):
        measured=identity(pin['path']);assert measured['sha256']==pin['sha256']
        assert json.loads(Path(pin['path']).read_text())['closed'];old.append(measured)
    assert len(old)==3
    previous=Path('/home/jay/build/wj-h26-native-flow-20261117/source-v8.json')
    previous_manifest=json.loads(previous.read_text())
    for pin in previous_manifest['files']:
        assert identity(pin['path'])['sha256']==pin['sha256'],pin['path']
    mac=json.loads(a.mac_binding.read_text());binding=[]
    assert mac['git_source']==A_COMMIT
    for pin in mac['files']:
        measured=identity(root/pin['executed_relative_path'])
        assert measured['sha256']==pin['sha256'],pin['repo_relative_path']
        binding.append(dict(repo_relative_path=pin['repo_relative_path'],**measured))
    native=json.loads((root/'science/native-comparison-v2.json').read_text())
    assert native['diagnostic_verification']=='PASS'and native['model_gate']=='NOT_MET'
    assert native['maximum_independent_error']<=1e-10
    bound=json.loads((root/'science/planar-bound-v1.json').read_text())
    assert bound['uniform_diffusion_bound']['diffusion_change_cannot_repair_retained_inlet']
    outer=json.loads((root/'science/outer-reference-v1.json').read_text())
    assert all(v['status']=='SOLVED'for v in outer['runs'])
    proposal=json.loads((root/'science/proposal-screen-v1.json').read_text())
    assert proposal['proposal_only']and not proposal['active_model_modified']and not proposal['task3_ready']
    assert all(v['status']=='POSITIVE_FINITE_COVERAGE_ONLY'for v in proposal['cases'])
    live=workers();own=[v for v in live['processes']if str(root)in v['args']and v['pid']!=os.getpid()]
    assert not own,own
    assert live['workers']<=72
    receipts=[json.loads(p.read_text())for p in sorted((root/'operations').glob('*/receipt.json'))]
    assert max(v['peak_host_workers']for v in receipts)<=72
    source_build=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','build','runtime')))
    assert source_build<=8*1024**3 and size(root/'science')<=8*1024**3
    free=os.statvfs(root).f_bavail*os.statvfs(root).f_frsize;assert free>=15*1024**3
    general=Path('/media/jay/新加卷1/nektar-v2026')
    checkout=dict(path=str(general),revision=subprocess.check_output(['git','-C',str(general),'rev-parse','HEAD'],text=True).strip(),
        status=subprocess.check_output(['git','-C',str(general),'status','--short'],text=True))
    assert checkout['revision']=='f4b70baec10016dd098005bc5a5507b87e4a57bf'
    assert checkout['status']==' M docs/developer-guide/src/library\n'
    def close(x):
        assert not x['closed']and all(v['status']!='RESERVED'for v in x['entries'].values())
        x['entries']['control-and-final-publication']=dict(status='CONSERVATIVELY_PREPAID',phase='evidence',seconds=300.,
            scope='H27 non-numerical control/audit/transport setup outside individually timed operations, seal/publication and transport verification; idle parent waits are not independent science operations')
        used=sum(v['seconds']for v in x['entries'].values());assert used<=43200
        x.update(closed=True,closed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            closure=dict(reason='Independently confirmed scientific model/startup barrier; positive H27 closure not achieved',
                         used_seconds=used,unused_seconds=43200-used,task3_ready=False))
    account=update(book,close)
    out.mkdir();(out/'raw').mkdir();(out/'operations').mkdir()
    for f in sorted((root/'science').glob('*')):
        if f.is_file():shutil.copy2(f,out/'raw'/f.name)
    for folder in sorted((root/'operations').iterdir()):
        dest=out/'operations'/folder.name;dest.mkdir()
        for name in('request.json','receipt.json','command.log'):
            f=folder/name
            if f.exists():shutil.copy2(f,dest/name)
    shutil.copy2(book,out/'accounting-final.json');shutil.copy2(root/'takeover.json',out/'takeover.json')
    gates=dict(H27_MODEL_AND_REFERENCE='NOT_MET',H27_COUPLED_NUMERICS='NOT_MET',
        H27_PHYSICAL_PILOTS_AND_ACCURACY='NOT_MET',H27_RESTART_MONITORING_STORAGE='NOT_MET',
        H27_COMPLETE_COST='TARGET_NOT_MET',TASK3_READY='NO')
    operation_rows=[dict(label=label,**row)for label,row in account['entries'].items()]
    def root_rows(name):
        d=json.loads((root/'science'/name).read_text());rows=[]
        for c in d['cases']:
            for q in c.get('roots',[]):
                if'audit'not in q:continue
                v=q['audit'];rows.append(dict(K=c['parameters']['K'],c_over_H=v['c_over_H'],
                    budget_error_over_A=v['budget_error_over_A'],ground_flux_over_A=v['ground_flux_over_A'],
                    maximum_weak_moment=max(map(abs,v['weak_moments_over_A'])),
                    analytic_ground_flux_derivative=v['analytic_ground_flux_derivative'],
                    min_sampled_principal=min(s['min_sampled_normalized']for s in v['principal_samples'])))
        return rows
    audit=dict(status='PASS_FOR_EVIDENCE_INTEGRITY_ONLY',git_source=A_COMMIT,
        host=platform.node(),platform=platform.platform(),python=sys.version,
        python_packages={name:importlib.metadata.version(name)for name in('numpy','scipy')},
        native_source_binding=binding,mac_binding=identity(a.mac_binding),
        unchanged_old_closed_accounts=old,unchanged_H26_source_manifest=identity(previous),
        unchanged_H26_source_files_checked=len(previous_manifest['files']),
        native_inherited_inputs_checked=len(json.loads((root/'source-v2.json').read_text())['verified_native_inputs']),
        native_source_manifest=identity(root/'source-v3.json'),executable=native['executable'],
        maximum_independent_native_error=native['maximum_independent_error'],
        source_build_runtime_bytes=source_build,science_bytes=size(root/'science'),free_bytes=free,
        peak_actual_host_workers=max(v['peak_host_workers']for v in receipts),
        maximum_admission_worker_bound=max(v['allocation']['workers']+v['reservation']for v in receipts),
        peak_build_tree_rss=max(v['peak_tree_rss_bytes']for v in receipts if v['phase']=='build'),
        peak_analysis_tree_rss=max(v['peak_tree_rss_bytes']for v in receipts if v['phase']!='build'),
        live_workers=live['workers'],compiler_jobs=4,numerical_library_threads=1,
        general_checkout=checkout,unresolved_operation_outcomes=0,
        filesystem=json.loads(subprocess.check_output(['findmnt','-J','-T',str(root),'-o','TARGET,SOURCE,FSTYPE'],text=True)),
        publisher=identity(__file__))
    atomic(out/'integrity-resource-audit.json',audit)
    summary=dict(schema='WJH27FeasibilityClosureV1',positive_closure=False,gates=gates,
        checkpoint_A=A_COMMIT,reviewed_base=takeover['reviewed'],
        native_source=identity(root/'source-v3.json'),native_executable=native['executable'],
        reference_original=root_rows('reference-kin-v1.json'),reference_refined=root_rows('reference-range-v1.json'),
        reference_auxiliary_control=root_rows('reference-aux-v1.json'),
        endpoint_minus_control=root_rows('endpoint-minus-v1.json'),endpoint_plus_control=root_rows('endpoint-plus-v1.json'),
        planar_bound=bound['uniform_diffusion_bound'],native_comparisons=native['comparisons'],
        outer_reference=outer,proposal_only_screen=dict(cases=[{k:v for k,v in c.items()if k!='states'}for c in proposal['cases']],
            active_model_modified=False,qualification='Finite frozen-state coverage only; adoption and coupled qualification pending'),
        physical_reference_and_pilots=dict(performed=False,reason='Specified model fails necessary principal admissibility on required startup',
            developed_2d_reference=None,physical_adjustment_to_inlet=None,accepted_physical_steps=0),
        coupled_numerics=dict(registered_solver_implemented=False,current_coupled_iterations=None,accuracy_qualified_dt=None),
        restart_storage=dict(coupled_cold_restart_performed=False,recovery_writes=0,rotations=0,
            scientific_terminal_retention_qualified=False,feasibility_artifacts_only=True),
        complete_cost=dict(status='UNMEASURED_PREREQUISITE_NOT_MET',per_trajectory=None,total_campaign=None,
            hard_gate_seconds=86400,margin=1.25,reason='No admissible coupled accuracy-qualified recipe; H26 timing is not extended'),
        horizon=HORIZON,terminal_window=TERMINAL,operations=operation_rows,
        account=account['closure'],resource_audit=identity(out/'integrity-resource-audit.json'))
    atomic(out/'summary.json',summary)
    manifest=dict(schema='WJH27Task3LaunchHoldV1',status='NOT_READY_DO_NOT_LAUNCH',gates=gates,
        horizon=HORIZON,terminal_window=TERMINAL,model='Specified H27 tensor/reference candidate1 inadmissible; amendment is proposal only',
        source=A_COMMIT,diagnostic_executable=native['executable'],coupled_solver_executable=None,
        initial_state=None,restart_state=None,qualified_mesh_p_q_dt=None,
        planned_case_classes=['developed_control','retained_inlet_startup','finite_perturbation','second_sign_or_shape'],
        per_trajectory_complete_forecasts=None,total_campaign_forecast=None,
        launch_commands=None,abort_before_launch=['scientific tensor decision required','coupled numerics not implemented/qualified',
            'same-inlet steady reference and physical pilots missing','cold restart/rotation/retention missing','complete accuracy-qualified cost unmeasured'],
        next_review='docs/plans/wj-h27-model-amendment-proposal.md')
    atomic(out/'task3-launch-hold.json',manifest)
    atomic(out/'publication-manifest.json',dict(files=[dict(relative_path=str(f.relative_to(out)),**identity(f))for f in sorted(out.rglob('*'))if f.is_file()],
        closed_account=identity(book),publisher_seconds_before_manifest=time.monotonic()-started))
    archive=root/'h27-feasibility-evidence-v1.tar.gz'
    assert not archive.exists()
    with tarfile.open(archive,'w:gz')as t:t.add(out,arcname=out.name)
    atomic(root/'h27-feasibility-evidence-v1.json',dict(archive=identity(archive),
        manifest=identity(out/'publication-manifest.json'),closed_account=identity(book),
        actual_publication_seconds=time.monotonic()-started,provision_seconds=300))
    print(json.dumps(dict(output=str(out),archive=identity(archive),closed_account=account['closure'])))

if __name__=='__main__':main()
