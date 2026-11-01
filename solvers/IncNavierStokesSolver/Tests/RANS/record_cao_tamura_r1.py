#!/usr/bin/env python3
"""Read-only CV3-R1 evidence export; native runs remain outside the repository."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import shutil
import sys

sys.dont_write_bytecode=True
from cao_tamura_experiment import ROOT,records,sha,write_csv,write_json
from cao_tamura_routes import coverage,route_metadata
from cao_tamura_cv3_budgets import flat_momentum_budget
from terrain_continuity import audit as strip_continuity


def bound(p):return dict(path=str(p),sha256=sha(p))


def field_bound(path):
    if not path.exists():raise ValueError('Missing native field artifact: '+str(path))
    return bound(path) if path.is_file() else dict(path=str(path),files={str(f.relative_to(path)):sha(f)
        for f in sorted(path.rglob('*')) if f.is_file()})


def step_timers(path):
    """Native steady-clock integration timers, excluding later checkpoint I/O."""
    pattern=re.compile(r'^Steps:\s*(\d+)\s+Time:\s*([\d.eE+-]+)\s+CPU Time:\s*([\d.eE+-]+)s\s*$',re.M)
    rows=[dict(step=int(n),time=float(t),seconds=float(s)) for n,t,s in pattern.findall(path.read_text(errors='replace'))]
    if any(not math.isfinite(r['time']+r['seconds']) or r['seconds']<0 for r in rows):
        raise ValueError('Invalid native integration timer: '+str(path))
    first=rows[:3];consecutive=[r['step'] for r in first]==[1,2,3]
    return dict(source=bound(path),reported_intervals=len(rows),first_three=first,
        first_three_individual_steps=consecutive,
        first_three_seconds=sum(r['seconds'] for r in first) if consecutive else None,
        basis='steady-clock TimeIntegrate/pre-integration intervals; excludes subsequent checkpoint I/O',
        limitation='Compare only identical physical/mesh/dt cases and reporting intervals; allocation/load effects remain')


def reference(folder):
    result=json.loads((folder/'evidence.json').read_text())
    if result['status']!='PASS':raise ValueError('Reference control has not passed: '+str(folder))
    rows=records(folder/'reference.awf.csv')
    if not rows or any(float(r[k])!=1. for r in rows for k in ('chi_min','chi_max')):
        raise ValueError('Reference lost full active EP correction')
    for n,s in result['artifacts'].items():
        if sha(folder/n)!=s:raise ValueError('Reference artifact changed: '+n)
    result.update(evidence=bound(folder/'evidence.json'),active_chi_range=[1.,1.])
    return result


def pilot(folder):
    operation=json.loads((folder/'operation.json').read_text())
    result=dict(directory=str(folder),operation=bound(folder/'operation.json'),status=operation['status'],
        surface=operation['config']['surface'],model=operation['config']['model'],
        route=route_metadata(folder/'case.xml'),config=operation['config'],numerics=operation['numerics'],
        calibration_eligible=False,binaries=operation['binaries'])
    path=folder/'run.json'
    if path.exists():
        r=json.loads(path.read_text());result.update(run=bound(path),exit_code=r['exit_code'],
            timeout=r['timeout'],wall_time_seconds=r.get('wall_time_seconds'))
    stopped=folder/'termination.json'
    aborted=stopped.exists() and result.get('exit_code')!=0
    if stopped.exists():result['controller_stop']=dict(receipt=bound(stopped),**json.loads(stopped.read_text()))
    if aborted:result['status']='ABORTED_NOT_SELECTED'
    log=folder/'run.log'
    if log.exists():result['native_step_timers']=step_timers(log)
    if result.get('wall_time_seconds') is None and log.exists():
        timer=re.findall(r'^Execute\s+([\d.eE+-]+)',log.read_text(errors='replace'),re.M)
        if timer:
            result.update(wall_time_seconds=float(timer[-1]),runtime_basis='native Execute region')
    errors=re.findall(r'RANS: [^\n]+',log.read_text(errors='replace')) if log.exists() else []
    result['error']=errors[-1] if errors else operation.get('error','')
    native=folder/'case.rans.csv';raw=records(native) if native.exists() else []
    keys=('time','min_k','min_eps','divergence_l2','weak_divergence_l2')
    rows=[]
    for row in raw:
        try:valid=all(math.isfinite(float(row[k])) for k in keys)
        except (ValueError,TypeError,KeyError):valid=False
        if valid:rows.append(row)
        elif not aborted:raise ValueError('Incomplete/nonfinite completed native record: '+str(native))
    result['incomplete_native_rows']=len(raw)-len(rows)
    result.update(accepted_timesteps=len(rows),last_time=float(rows[-1]['time']) if rows else 0.)
    if rows:
        result.update(native=bound(native),minimum_k=min(float(r['min_k']) for r in rows),
            minimum_epsilon=min(float(r['min_eps']) for r in rows),
            maximum_strong_divergence=max(float(r['divergence_l2']) for r in rows),
            maximum_weak_divergence=max(float(r['weak_divergence_l2']) for r in rows))
        try:result['coverage']=coverage(folder/'case.xml')
        except Exception as e:result['coverage']=dict(status='UNAVAILABLE' if aborted else 'FAIL',error=str(e))
        residual=folder/'case.rans-residual.csv'
        free=[float(r['free_coefficient_l2']) for r in records(residual) if r['phase']=='1'] if residual.exists() else []
        result['maximum_free_residual']=max(free) if free else None
        try:
            result['independent_strip_continuity']=strip_continuity(folder,'case')
            if len(rows)>=3:result['independent_strip_momentum']=flat_momentum_budget(folder/'case.xml')
        except Exception as e:
            if not aborted:raise
            result['prefix_budget_limitation']=str(e)+'; incomplete exploratory prefixes do not qualify for development'
    return result


def developed(folder,solver_sha256):
    state=json.loads((folder/'development.json').read_text())
    exit_checks=[]
    unforced=[s['elapsed_flow_times']-state['config']['startup_flow_times'] for s in state['stages']]
    if (not state['complete'] or state['purpose']!='feasibility_only' or len(state['stages'])<3
            or min(unforced[-3:])<3 or unforced[-1]<5
            or any(b-a<1.-1e-9 for a,b in zip(unforced[-3:-1],unforced[-2:]))
            or not all(s['drift']['accepted'] for s in state['stages'][-2:])):
        raise ValueError('Undeveloped composite case: '+str(folder))
    for s in state['stages']:
        case=folder/f"stage-{s['index']:02d}"/'case.xml'
        r=s['run'];m=r['stages']
        manifest=json.loads((case.parent/'manifest.json').read_text())
        tf=(manifest['config']['xmax']-manifest['config']['xmin'])/manifest['reference']['Uinf']
        if abs(m['final_time']/tf-s['elapsed_flow_times'])>1e-9:
            raise ValueError('Development duration differs from the native trajectory')
        if (r['case_sha256']!=sha(case) or r['solver_sha256']!=solver_sha256 or r['exit_code'] or not r['binaries_unchanged_during_run']
                or m['maximum_guard_activations'] or m['final_startup_weight']
                or min(m['minimum_k'],m['minimum_epsilon'])<=0
                or max(m[k] for k in ('maximum_weak_divergence','maximum_projected_normal',
                    'maximum_normal_integral_defect','maximum_resolved_mass_defect'))>1e-7
                or not math.isfinite(s['maximum_free_residual']) or s['maximum_free_residual']>1e-7):
            raise ValueError('Invalid native stage in '+str(folder))
        if route_metadata(case)['kind']!='EPTransition':raise ValueError('A diagnostic cannot certify the composite route')
        exit_checks.append(dict(stage=s['index'],coverage=coverage(case),
            strip_continuity=strip_continuity(case.parent,'case'),strip_momentum=flat_momentum_budget(case)))
    return state,exit_checks


def numerical_validation(runs,build,history,ctest_log,checks_root=None):
    """Bind the new numerical map to its declared, independent checks."""
    solver=sha(build/'IncNavierStokesSolver');checks={}
    recipe=json.loads((history/'whole-recipe.json').read_text())
    pointwise=recipe.get('numerics',{}).get('RANSNewtonScalarScaling')=='Pointwise'
    for surface in ('S','R'):
        for model in ('standard_fixed','cls_fixed'):
            path=(checks_root or runs/'newton-map-checks-v2')/f'{surface}-{model}'/'evidence.json'
            value=json.loads(path.read_text())
            if (value['status']!='PASS'
                    or value['cases']['candidate']['run']['solver_sha256']!=solver
                    or max(value['scaled_field_errors'].values())>1e-6
                    or (not pointwise and value['runtime_ratio_baseline_over_candidate']<=1)):
                raise ValueError('Numerical-map equivalence/speed check failed: '+str(path))
            if pointwise and (not value.get('producer_sha256')
                    or sha(path.parent/'producer.py')!=value['producer_sha256']):
                raise ValueError('Pointwise matched-field producer is not bound')
            checks[surface+'/'+model]=dict(receipt=bound(path),**value)
    path=history/'evidence.json';value=json.loads(path.read_text())
    if (value['status']!='PASS' or value['run']['solver_sha256']!=solver
            or not value['map_mismatch_rejected'] or max(value['scaled_field_errors'].values())>1e-6
            or (pointwise and not value.get('scalar_scaling_mismatch_rejected'))):
        raise ValueError('Numerical-map history/MPI check failed')
    if sha(history/'producer.py')!=value['producer_sha256']:
        raise ValueError('History producer differs from its recorded source')
    log=ctest_log.read_text()
    if log.count('Test Passed.')!=16 or 'Test Failed.' in log:
        raise ValueError('The complete 16-test Release suite has not passed')
    return dict(status='PASS',map_equivalence=checks,history=dict(receipt=bound(path),**value),
        history_recipe=bound(history/'whole-recipe.json'),pointwise_scalar_scaling=pointwise,
        ctest=dict(log=bound(ctest_log),passed=16,failed=0))


def development_attempt(folder):
    """Retain failed, retired and unselected precursors alongside selections."""
    operation=json.loads((folder/'operation.json').read_text()) if (folder/'operation.json').exists() else {}
    state=json.loads((folder/'development.json').read_text())
    result=dict(directory=str(folder),status=operation.get('status','RUNNING'),
        development=bound(folder/'development.json'),complete=state['complete'],
        config=state['config'],numerics=state['numerics'],segments=[])
    start=Path(str(folder)+'.start.json')
    if start.exists():result['start']=dict(receipt=bound(start),**json.loads(start.read_text()))
    if operation:result['operation']=dict(receipt=bound(folder/'operation.json'),**operation)
    if (folder/'termination.json').exists():result['controller_stop']=dict(receipt=bound(folder/'termination.json'),
        **json.loads((folder/'termination.json').read_text()))
    if result.get('controller_stop') and not state['complete']:
        result['recorded_operation_status']=result['status']
        result['status']='RETIRED_NOT_DEVELOPED'
    for path in sorted(folder.glob('stage-*/run.json')):
        r=json.loads(path.read_text());segment=dict(run=bound(path),receipt=r)
        if (path.parent/'run.log').exists():segment['native_step_timers']=step_timers(path.parent/'run.log')
        native=path.parent/'case.rans.csv'
        accepted=set();incomplete=0
        for v in records(native) if native.exists() else []:
            try:valid=all(math.isfinite(float(v[k])) for k in ('time','min_k','min_eps','weak_divergence_l2','divergence_l2'))
            except (KeyError,TypeError,ValueError):valid=False
            if valid:accepted.add(float(v['time']))
            else:incomplete+=1
        if incomplete and r['complete']:raise ValueError('Invalid completed native diagnostic: '+str(native))
        segment['incomplete_native_rows']=incomplete
        if accepted:
            try:segment['coverage']=coverage(path.parent/'case.xml')
            except Exception as error:
                segment['coverage']=dict(status='UNAVAILABLE' if result.get('controller_stop') else 'FAIL',
                    error=str(error))
        if not r['complete'] and accepted:
            extrema={};minimum_kp={}
            for points in path.parent.glob('case.terrain-points-r*.csv'):
                for v in records(points):
                    if v.get('phase')!='1' or not v.get('normal_error') or float(v['time']) not in accepted:continue
                    time=float(v['time'])
                    if time not in extrema or abs(float(v['normal_error']))>abs(float(extrema[time]['normal_error'])):
                        extrema[time]={k:v[k] for k in ('time','x','mask','normal_error','normal_target','chi','reason','uH','uN','kP')}
                    if time not in minimum_kp or float(v['kP'])<float(minimum_kp[time]['kP']):
                        minimum_kp[time]={k:v[k] for k in ('time','x','kP','chi','reason','wall_re','cutoff')}
            segment['pointwise_normal_extrema']=[extrema[t] for t in sorted(extrema)]
            segment['eliminated_k_minima']=[minimum_kp[t] for t in sorted(minimum_kp)]
            try:
                segment['independent_strip_continuity']=strip_continuity(path.parent,'case')
                if len(accepted)>=3:
                    segment['independent_strip_momentum']=flat_momentum_budget(path.parent/'case.xml')
            except Exception as error:
                segment['prefix_budget_limitation']=str(error)
        result['segments'].append(segment)
    return result


def checkpoint_checks(runs,build):
    results={}
    for path in sorted(runs.glob('*checkpoint-v*/evidence.json')):
        value=json.loads(path.read_text());item=dict(receipt=bound(path),result=value)
        if value['status']=='PASS':
            if (max(value['scaled_field_errors'].values())>1e-6
                    or sha(path.parent/'producer.py')!=value['producer_sha256']):
                raise ValueError('Invalid checkpoint comparison: '+str(path))
            item['uses_selected_binary']=value['binaries']['IncNavierStokesSolver']==sha(build/'IncNavierStokesSolver')
            item['ranks']={}
            for name,case in value['cases'].items():
                start=Path(case['field']['path']).parent/'run-start.json'
                if (sha(start)!=case['run_start_sha256'] or field_bound(Path(case['field']['path']))!=case['field']):
                    raise ValueError('Producing checkpoint or receipt changed')
                producing=json.loads(start.read_text())
                if producing['solver_sha256']!=value['binaries']['IncNavierStokesSolver']:
                    raise ValueError('Checkpoint comparison has a different producing binary')
                item['ranks'][name]=producing['ranks']
        results[path.parent.name]=item
    return results


def development_figures(attempts,out):
    """Plot accepted native diagnostics only; no experimental scores."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=[];bindings={}
    for item in attempts:
        folder=Path(item['directory']);c=item['config']
        for segment in item['segments']:
            directory=Path(segment['run']['path']).parent
            path=directory/'case.rans.csv'
            if not path.exists():continue
            manifest=json.loads((directory/'manifest.json').read_text())
            tf=(c['xmax']-c['xmin'])/manifest['reference']['Uinf']
            bindings[str(path)]=sha(path)
            for r in records(path):
                try:
                    values={k:float(r[k]) for k in ('time','divergence_l2','weak_divergence_l2','min_k','min_eps')}
                    if not all(math.isfinite(v) for v in values.values()):continue
                except (KeyError,ValueError,TypeError):continue
                rows.append(dict(control=folder.parent.name,surface=c['surface'],model=c['model'],
                    stage=directory.name,unforced_flow_times=values['time']/tf-c['startup_flow_times'],
                    startup_weight=float(r.get('startup_weight',0)),**values))
    if not rows:return None
    write_csv(out/'development-native.csv',rows,list(rows[0]))
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    names=sorted({r['control'] for r in rows if r['model']=='cls_fixed'})
    colours={name:plt.get_cmap('tab20')(i%20) for i,name in enumerate(names)}
    for column,surface in enumerate(('S','R')):
        controls=sorted({r['control'] for r in rows if r['surface']==surface and r['model']=='cls_fixed'})
        for name in controls:
            values=sorted((r for r in rows if r['surface']==surface and r['model']=='cls_fixed'
                and r['control']==name and r['startup_weight']==0),key=lambda r:r['time'])
            if not values:continue
            label=name.removeprefix('development-')
            for ax,key in zip(axes[:,column],('divergence_l2','weak_divergence_l2')):
                ax.semilogy([r['unforced_flow_times'] for r in values],
                    [max(1e-16,r[key]) for r in values],label=label,color=colours[name])
        axes[0,column].set_title(surface+' CLS: accepted unforced steps')
        axes[0,column].set_ylabel('Strong divergence L2 (native units)')
        axes[1,column].set_ylabel('Weak divergence L2 (native units)')
        axes[1,column].axhline(1e-7,color='black',linestyle='--',linewidth=.8,label='Weak acceptance gate')
        axes[1,column].set_xlabel('Unforced domain flow-throughs')
        axes[0,column].legend(fontsize=7,ncol=2)
        for ax in axes[:,column]:ax.grid(alpha=.2)
    fig.suptitle('CLS feasibility diagnostics; rejected stages and forced startup times excluded')
    fig.savefig(out/'cls-development.png',dpi=180);plt.close(fig)
    return dict(csv=bound(out/'development-native.csv'),figure=bound(out/'cls-development.png'),
        native_sources=bindings,rule='Every complete finite accepted diagnostic row; forced rows retained in CSV but excluded from plot')


def newton_action_diagnostics(runs,out):
    """Retain recomputed actions; an Arnoldi estimate alone is not this check."""
    rows=[];sources={}
    paths=sorted(runs.glob('pilot-*/*/run.log'))+sorted(runs.glob('conditioning-*/run.log'))
    pattern=re.compile(r'^RANS Newton (\d+) GMRES (\d+) relative residual ([\d.eE+-]+)$')
    for path in paths:
        selected=[]
        for line,text in enumerate(path.read_text().splitlines(),1):
            match=pattern.match(text)
            if match:
                value=float(match[3])
                if not math.isfinite(value):raise ValueError('Nonfinite recomputed Newton action')
                selected.append(dict(control=str(path.parent.relative_to(runs)),line=line,
                    newton_iteration=int(match[1]),gmres_iterations=int(match[2]),relative_residual=value))
        if selected:sources[str(path)]=sha(path);rows.extend(selected)
    if not rows:return None
    write_csv(out/'newton-actions.csv',rows,list(rows[0]))
    return dict(csv=bound(out/'newton-actions.csv'),rows=len(rows),sources=sources,
        definition='Native norm of J(step)-rhs divided by rhs norm, recomputed after the GMRES solve',
        limitation='Diagnostic linear action; physical stage acceptance uses separate equation/continuity/positivity gates')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reference-root',type=Path,help='Current-binary reference controls, if stored in a separate epoch')
    p.add_argument('--native-guards',type=Path,help='Current-binary native guard evidence.json')
    p.add_argument('--map-history',type=Path,help='Passed new-map history/rank-change directory')
    p.add_argument('--map-checks',type=Path,help='Current-binary S/R standard/CLS matched-field checks')
    p.add_argument('--ctest-log',type=Path,help='Complete current-binary Release LastTest.log')
    p.add_argument('--developments',type=Path,help='S/R -> model -> completed development directory')
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    build_receipt=json.loads((a.build/'build-receipt.json').read_text())
    if (any(sha(ROOT/name)!=value for name,value in build_receipt['sources'].items())
            or any(sha(a.build/name)!=value for name,value in build_receipt['binaries'].items())
            or sha(a.build/'geometry-lib/libSpatialDomains.so')!=build_receipt['geometry_sha256']):
        raise ValueError('Build source or immutable binary set differs from its receipt')
    reference_root=a.reference_root or a.runs
    refs={k:reference(reference_root/k) for k in ('reference','perturbation','molecular')}
    solver_sha256=sha(a.build/'IncNavierStokesSolver')
    if any(r['binaries']['IncNavierStokesSolver']!=solver_sha256 for r in refs.values()):
        raise ValueError('Reference controls do not use the selected numerical binary')
    attempts=[pilot(p.parent) for p in sorted(a.runs.glob('pilot-*/*/operation.json'))]
    if any(r['status']=='RUNNING' for r in attempts):raise ValueError('Wait for live pilots before recording a final decision')
    cases={};complete=False
    if a.developments:
        paths=json.loads(a.developments.read_text())
        for s in ('S','R'):
            common=None
            for m in ('standard_fixed','cls_fixed'):
                folder=Path(paths[s][m]);state,exit_checks=developed(folder,solver_sha256);c=state['config']
                if c['surface']!=s or c['model']!=m or c['diffusion']!='Molecular':raise ValueError('Wrong primary case')
                physical={k:v for k,v in c.items() if k not in ('model','dt','steps','checkpoint')}
                if common is not None and common!=physical:raise ValueError('Closure inputs/mesh differ')
                artifacts={}
                for stage in state['stages']:
                    directory=folder/f"stage-{stage['index']:02d}"
                    for name in ('case.xml','manifest.json','run.json','incident.csv','incident.csv.json','full.csv','wall.json'):
                        artifacts[str(directory.relative_to(folder)/name)]=bound(directory/name)
                    artifacts[str(directory.relative_to(folder)/'case.fld')]=field_bound(directory/'case.fld')
                common=physical;cases[s+'/'+m]=dict(development=bound(folder/'development.json'),state=state,
                    exit_checks=exit_checks,artifacts=artifacts,statistics=bound(folder/'statistics.csv'))
        complete=True
    native=a.native_guards or a.runs/'native-guards/evidence.json'
    guards=json.loads(native.read_text()) if native.exists() else None
    if complete and (guards is None or guards['status']!='PASS'
            or guards['binary_sha256']['IncNavierStokesSolver']!=solver_sha256):
        raise ValueError('Current-binary native provenance/history checks are missing')
    solver_checks=None
    if a.map_history is not None or a.ctest_log is not None:
        if a.map_history is None or a.ctest_log is None:raise ValueError('Supply both numerical history and CTest evidence')
        solver_checks=numerical_validation(a.runs,a.build,a.map_history,a.ctest_log,a.map_checks)
    if complete and any(c['state']['numerics'].get('RANSNewtonMap')=='PreconditionedResidual' for c in cases.values()) and solver_checks is None:
        raise ValueError('Supply new-map numerical validation evidence')
    checkpoint_validation=checkpoint_checks(a.runs,a.build)
    if (complete and any(s['run']['ranks']==8 for c in cases.values() for s in c['state']['stages'])
            and not any(v.get('uses_selected_binary') and v.get('ranks')==dict(baseline=4,candidate=8)
                        for v in checkpoint_validation.values())):
        raise ValueError('Missing validated four-to-eight-rank checkpoint comparison')
    development_attempts=[development_attempt(Path(str(p)[:-len('.start.json')]))
        for p in sorted(a.runs.glob('development-*/*.start.json'))]
    if complete and any(v['status']=='RUNNING' for v in development_attempts):
        raise ValueError('Wait for or explicitly retire live unselected development attempts')
    continuation_controls=[]
    for path in sorted(a.runs.glob('conditioning-*/evidence.json')):
        value=json.loads(path.read_text())
        if value.get('producer_sha256') and sha(path.parent/'producer.py')!=value['producer_sha256']:
            raise ValueError('Continuation producer differs from its recorded source')
        item=dict(receipt=bound(path),result=value,status=value['status'],
            artifacts={n:bound(path.parent/n) for n in ('producer.py','case.xml','manifest.json',
                'run.json','run.log','case.rans.csv','case.rans-residual.csv') if (path.parent/n).exists()})
        stop=path.parent/'termination.json'
        if stop.exists():
            item['controller_stop']=dict(receipt=bound(stop),**json.loads(stop.read_text()))
            item['status']='RETIRED_DIAGNOSTIC_NOT_DEVELOPED'
        native=path.parent/'run.json'
        if native.exists():
            run_value=json.loads(native.read_text())
            item['native_solver_completed']=bool(run_value.get('complete') and run_value.get('exit_code')==0)
            if (item['native_solver_completed'] and value['status']=='FAIL'
                    and run_value.get('stages',{}).get('completed_stages',3)<3):
                item['audit_failure_kind']='minimum_three_native_stages_required_for_BDF2_inventory'
        continuation_controls.append(item)
    if complete and any(v['status']=='RUNNING' for v in continuation_controls):
        raise ValueError('Wait for or explicitly retire live conditioning diagnostics')
    result=dict(version='cao-tamura-cv3-r1-evidence-v1',status='CLOSED' if complete else 'OPEN',
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        producer_sha256=sha(__file__),command=[sys.executable,__file__,*sys.argv[1:]],
        exit_audit_sources={n:sha(Path(__file__).parent/n) for n in
            ('cao_tamura_routes.py','terrain_continuity.py','cao_tamura_cv3_budgets.py')},
        complete=complete,calibration_complete=False,references=refs,pilots=attempts,developed_cases=cases,
        native_guards=dict(receipt=bound(native),**guards) if guards else None,
        numerical_validation=solver_checks,
        checkpoint_validation=checkpoint_validation,
        attempted_developments=development_attempts,
        continuation_controls=continuation_controls,
        build=build_receipt,
        build_binding=dict(recorded_sources_match=True,recorded_binaries_match=True,geometry_matches=True),
        limitation='R1 feasibility only; no fit selection, freeze, reserved flat or hill scores; CV4 convergence pending')
    columns=['surface','model','status','accepted_timesteps','last_time','exit_code','timeout',
             'wall_time_seconds','maximum_weak_divergence','maximum_strong_divergence','error','directory']
    write_csv(out/'pilots.csv',[{k:r.get(k,'') for k in columns} for r in attempts],columns)
    for name,r in refs.items():
        write_json(out/(name+'.json'),r)
        producer=reference_root/name/'producer.py'
        if sha(producer)!=r['script_sha256']:raise ValueError('Reference producer differs from its recorded source')
        archived=out/'producers';archived.mkdir(exist_ok=True)
        shutil.copy2(producer,archived/(r['script_sha256']+'.py'))
    if solver_checks:
        producer=a.map_history/'producer.py'
        shutil.copy2(producer,out/'producers'/(solver_checks['history']['producer_sha256']+'.py'))
        for value in solver_checks['map_equivalence'].values():
            if value.get('producer_sha256'):
                producer=Path(value['receipt']['path']).parent/'producer.py'
                shutil.copy2(producer,out/'producers'/(value['producer_sha256']+'.py'))
    shutil.copy2(__file__,out/'producers'/(result['producer_sha256']+'.py'))
    for item in continuation_controls:
        producer=item['artifacts'].get('producer.py')
        if producer:
            shutil.copy2(Path(producer['path']),out/'producers'/(producer['sha256']+'.py'))
    if checkpoint_validation:
        checks=out/'checkpoint-checks';checks.mkdir()
        for name,item in checkpoint_validation.items():
            shutil.copy2(Path(item['receipt']['path']),checks/(name+'.json'))
            if item['result']['status']=='PASS':
                producer=Path(item['receipt']['path']).parent/'producer.py'
                shutil.copy2(producer,out/'producers'/(item['result']['producer_sha256']+'.py'))
    result['producer_epochs']=[]
    archived_sources=list(a.runs.glob('producer-epochs/*/*.py'))
    archived_sources+=list(a.runs.glob('development-*/*/producer-sources/*.py'))
    for producer in sorted(archived_sources):
        digest=sha(producer);shutil.copy2(producer,out/'producers'/(digest+'.py'))
        result['producer_epochs'].append(dict(source=bound(producer),archived=digest+'.py'))
    result['compiled_source_epochs']=[]
    for manifest in sorted(a.runs.glob('compiled-source-epochs/*/manifest.json')):
        epoch=manifest.parent;value=json.loads(manifest.read_text())
        target=out/'compiled-source-epochs'/epoch.name;target.mkdir(parents=True)
        for name,digest in value['sources'].items():
            source=epoch/name
            if sha(source)!=digest:raise ValueError('Archived compiled source changed: '+str(source))
            destination=target/name;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source,destination)
        shutil.copy2(manifest,target/'manifest.json')
        result['compiled_source_epochs'].append(dict(manifest=bound(manifest),
            archived=str(target.relative_to(out)),sources_verified=True))
    result['read_only_diagnostics']=[]
    for path in sorted(a.runs.glob('*diagnostic-v*/evidence.json')):
        value=json.loads(path.read_text());producer=path.parent/'producer.py'
        if sha(producer)!=value['producer_sha256']:raise ValueError('Diagnostic producer changed')
        if complete and value['status']=='RUNNING':raise ValueError('Wait for live read-only diagnostics')
        for source,digest in value.get('extra_sources',{}).items():
            if sha(Path(source))!=digest:raise ValueError('Diagnostic source changed: '+source)
            shutil.copy2(source,out/'producers'/(digest+Path(source).suffix))
        for case in value.get('cases',{}).values():
            for source,digest in case.get('source_sha256',{}).items():
                if sha(Path(source))!=digest:raise ValueError('Diagnostic native source changed')
            folder=Path(case['folder'])
            if (sha(folder/'case.rans.csv')!=case['stage_csv_sha256']
                    or sha(folder/'manifest.json')!=case['manifest_sha256']):
                raise ValueError('Diagnostic stage or manifest changed')
            if case.get('field') and field_bound(Path(case['field']['path']))!=case['field']:
                raise ValueError('Diagnostic native field changed')
        target=out/'diagnostics';target.mkdir(exist_ok=True)
        shutil.copy2(path,target/(path.parent.name+'.json'))
        shutil.copy2(producer,out/'producers'/(value['producer_sha256']+'.py'))
        for name,digest in value.get('compact_artifacts',{}).items():
            if Path(name).name!=name or sha(path.parent/name)!=digest:raise ValueError('Diagnostic compact artifact changed')
            directory=target/path.parent.name;directory.mkdir(exist_ok=True)
            shutil.copy2(path.parent/name,directory/name)
        result['read_only_diagnostics'].append(dict(receipt=bound(path),result=value))
    recipes=out/'recipes';recipes.mkdir()
    for pattern in ('pilot-*.json','development-*.json'):
        for file in a.runs.glob(pattern):shutil.copy2(file,recipes/file.name)
    shutil.copy2(a.build/'build-receipt.json',out/'build.json')
    result['reported_developed_controls']={}
    for item in development_attempts:
        if not item['complete']:continue
        folder=Path(item['directory']);state=json.loads((folder/'development.json').read_text())
        if sha(folder/'statistics.csv')!=state['statistics_sha256']:
            raise ValueError('Developed-control statistics changed: '+str(folder))
        final=Path(state['final_stage']);target=out/'developed-unfitted'/folder.parent.name/folder.name
        target.mkdir(parents=True)
        files={}
        for name,source in (('development.json',folder/'development.json'),
                ('statistics.csv',folder/'statistics.csv'),('statistics.csv.json',folder/'statistics.csv.json'),
                ('manifest.json',final/'manifest.json'),('incident.csv',final/'incident.csv'),
                ('incident.csv.json',final/'incident.csv.json'),('wall.json',final/'wall.json')):
            shutil.copy2(source,target/name);files[name]=bound(source)
        result['reported_developed_controls'][str(folder.relative_to(a.runs))]=dict(artifacts=files,
            scope='Unfitted producer-reported developed control; no common-pair selection or calibration freeze')
        control=result['reported_developed_controls'][str(folder.relative_to(a.runs))]
        if item.get('operation',{}).get('binaries',{}).get('IncNavierStokesSolver')==solver_sha256:
            _,checks=developed(folder,solver_sha256)
            control.update(independent_current_binary_exit_checks=checks,
                scope='Unfitted developed control with independent current-binary exit checks; no common-pair selection or calibration freeze')
    result['development_diagnostics']=development_figures(development_attempts,out)
    result['newton_action_diagnostics']=newton_action_diagnostics(a.runs,out)
    # Publish the exit decision only after every bound artifact/producer copy
    # has succeeded. A failed export must never leave a CLOSED receipt.
    write_json(out/'evidence.json',result)
    print('CV3-R1',result['status'],flush=True)


if __name__=='__main__':main()
