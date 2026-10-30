#!/usr/bin/env python3
"""Collect compact closure evidence from the completed authoritative Linux runs."""
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def require(condition,message):
    if not condition:raise RuntimeError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--build',type=Path);args=parser.parse_args()
    require(sys.platform.startswith('linux'),'Collect authoritative evidence on Linux')
    run=args.run.resolve();build=args.build.resolve() if args.build else run/'final-build';source=Path(__file__).resolve().parents[4]
    phases={p:(run/f'final-{p}.status').read_text().strip() for p in 'ABCDEF'}
    require(all(v=='0' for v in phases.values()),'Final validation phases are incomplete or failed: '+str(phases))
    paths={
        'study':run/'final-study/results-study.json',
        'pilots':run/'final-verification/results-pilots.json',
        'restart':run/'final-verification/results-restart.json',
        'columns':run/'final-columns/results-columns.json',
        'robust':run/'final-robust/results-robust.json',
        'outflow':run/'final-outflow-equivalence/results.json',
        'standard':run/'final-standard/results-standard.json'}
    paths.update({f'm{m}':run/f'final-regressions/m{m}/validation.json' for m in (3,4,5)})
    data={key:json.loads(path.read_text()) for key,path in paths.items()}
    binaries={p.name:sha(p) for p in (build/'IncNavierStokesSolver',build/'RANSFieldAudit')}
    for key,report in data.items():
        require(report.get('binaries',report.get('binary_sha256'))==binaries,'Mixed executable versions: '+key)
        if key in ('m3','m4','m5','outflow'):require(report['status']=='passed','Regression did not pass: '+key)
    baseline={}
    for dim in (2,3):
        for model in ('standard','cls_fixed','cls_switched'):
            name=f'five_d{dim}_{model}';case=data['study']['cases'][name]
            require(case.get('engineering_target_met') and case['spec']['flow_times']>=5,'Preservation target failed: '+name)
            baseline[name]=dict(spec=case['spec'],peak_l2=case['peak_l2'],peak_linf=case['peak_linf'],final=case['final_errors'])
    refined=data['robust']['cases']['p5_3d']
    # The accepted resolution is degree four. The higher-degree experiment
    # must be completed or explicitly recorded as a measured robustness limit.
    if refined.get('status') in ('numerical_failure','profile_target_failure'):
        require(refined.get('completed_flow_times',0)>0 and refined.get('artifacts'), 'Missing failed-refinement evidence')
        require(refined['fingerprint']['binaries']==binaries,'Failed refinement used another executable')
        if refined['status']=='profile_target_failure':
            require(refined['engineering_stop']['value']>=.001,'Recorded profile target was not exceeded')
    else:
        require(refined['spec']['flow_times']>=5,'Incomplete 3D degree-five experiment')
    for field in ('u','eps_tilde','rans_nut'):
        for names in (('study_h6','five_d2_cls_switched','study_h24'),('study_p3','five_d2_cls_switched','study_p5')):
            errors=[data['study']['cases'][name]['final_errors'][field]['normalized_l2'] for name in names]
            require(errors[0]>errors[1]>errors[2],'Spatial errors did not refine: '+field+' '+str(names))
    perturbations={}
    for name in ('study_perturbed','three_perturbed'):
        case=data['study']['cases'][name]
        require(case.get('spec',{}).get('flow_times',0)>=5 and 'final_errors' in case,'Perturbation run incomplete: '+name)
        require(case['final_errors']['k']['normalized_l2']<.1*case['peak_l2']['k'],'Perturbation did not decay: '+name)
        perturbations[name]=dict(spec=case['spec'],final=case['final_errors'],peak_l2=case['peak_l2'])
    for dim in (2,3):
        name='standard_perturbed_d'+str(dim);case=data['standard']['cases'][name]
        require(case['spec']['flow_times']>=5 and case['final_errors']['k']['normalized_l2']<.1*case['peak_l2']['k'],
            'Standard perturbation did not complete/decay: '+name)
        perturbations[name]=dict(spec=case['spec'],final=case['final_errors'],peak_l2=case['peak_l2'])
    copies={}
    for dim in (2,3):
        name=f'copied_column_d{dim}';case=data['columns']['cases'][name]
        for field in ('u','k','eps_tilde'):
            require(case['peak_l2'][field]<.001,'Copied numerical-column error exceeds 0.1%: '+name+' '+field)
        copies[name]=dict(spec=case['spec'],peak_l2=case['peak_l2'],final=case['final_errors'])
    restart_differences={}
    for path in sorted((run/'final-verification').glob('compare*.log')):
        raw=path.read_text();value=json.loads(raw[raw.index('{"time"'):])
        restart_differences[path.stem]=max(f['difference_linf'] for f in value['fields'].values())
    temporal=data['restart']['temporal_rates']
    require(all(min(rates)>1.65 for fields in temporal.values() for rates in fields.values()),'Temporal convergence gate failed')
    case_summaries={}
    for name,case in data['study']['cases'].items():
        case_summaries[name]={key:case[key] for key in ('status','message','spec','engineering_target_met','peak_l2','peak_linf') if key in case}
    files=set()
    for pattern in ('solvers/IncNavierStokesSolver/**/*.cpp','solvers/IncNavierStokesSolver/**/*.h',
                    'library/SolverUtils/Turbulence/*','library/UnitTests/SolverUtils/RANS*.cpp',
                    'library/LibUtilities/TimeIntegration/*.h','solvers/IncNavierStokesSolver/Tests/RANS/*.py',
                    'solvers/IncNavierStokesSolver/Tests/RANS/*.xml','solvers/IncNavierStokesSolver/Tests/RANS/*.sh'):
        files.update(p for p in source.glob(pattern) if p.is_file())
    files.update(source/p for p in ('solvers/IncNavierStokesSolver/CMakeLists.txt','solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt'))
    ctest=(run/'final-ctest.log').read_text();require('100% tests passed' in ctest,'CTest did not pass')
    refresh=run/'audit-refresh.json'
    audit_refresh=json.loads(refresh.read_text()) if refresh.exists() else None
    if audit_refresh:require(audit_refresh['status']=='passed' and audit_refresh['binaries']==binaries,'Audit refresh is incomplete')
    trial=run/'final-robust/p5_3d-linear-1e-13'
    trial_rows=list(csv.DictReader((trial/'p5_3d.rans.csv').open())) if trial.exists() else []
    robustness_trial=dict(status='numerical_failure',reason='coupled implicit stage failed to converge',
        last_completed_step=trial_rows[-1],
        fixture_sha256=sha(trial/'p5_3d.xml'),solver_log_sha256=sha(trial/'solver.log'),
        diagnostics_sha256=sha(trial/'p5_3d.rans.csv')) if trial_rows else None
    versions={}
    for name in ('numpy','scipy','sympy','matplotlib'):
        module=__import__(name);versions[name]=module.__version__
    report=dict(schema=1,status='passed',milestone='neutral_abl_rans_awf_m6',utc=datetime.now(timezone.utc).isoformat(),
        closure_meaning='Degree-four acceptance gates passed; measured original-wall, inconsistent-boundary and higher-degree outcomes are retained separately',
        source_parent=subprocess.check_output(['git','-C',source,'rev-parse','HEAD'],text=True).strip(),
        source_sha256={str(p.relative_to(source)):sha(p) for p in sorted(files)},binary_sha256=binaries,
        platform=platform.platform(),python=platform.python_version(),packages=versions,
        compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
        cmake=subprocess.check_output(['cmake','--version'],text=True).splitlines()[0],
        mpi=subprocess.check_output(['mpirun','--version'],text=True).splitlines()[0],
        maximum_project_workers=16,maximum_mpi_ranks_per_case=4,threads={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')},
        concurrency='A<=8 plus D<=4 plus B/E/F<=4 sequentially; C waits for A',phases=phases,
        audit_refresh=dict(path=str(refresh),sha256=sha(refresh),complete_field_audits=len(audit_refresh['commands'])) if audit_refresh else None,
        retained_robustness_trial=robustness_trial,
        higher_degree_refinement={key:refined[key] for key in ('status','message','spec','planned_flow_times','completed_time',
            'completed_flow_times','engineering_target_met','peak_l2','peak_linf','engineering_stop') if key in refined},
        ctest=dict(passed=8,log_sha256=sha(run/'final-ctest.log')),
        report_sha256={key:sha(path) for key,path in paths.items()},
        report_paths={key:str(path) for key,path in paths.items()},
        baseline=baseline,perturbations=perturbations,study=case_summaries,numerical_column=copies,temporal_rates=temporal,
        restart_maximum_linf_differences=restart_differences,
        outflow_backend_maximum_linf={key:max(f['difference_linf'] for f in value['fields'].values()) for key,value in data['outflow']['comparisons'].items()},
        regression_case_counts={key:len(data[key]['cases']) for key in ('m3','m4','m5')},
        acceptance_case_counts={key:len(data[key]['cases']) for key in ('study','pilots','restart','columns','robust','standard')})
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ('status','acceptance_case_counts','regression_case_counts','binary_sha256')},indent=2))

if __name__=='__main__':main()
