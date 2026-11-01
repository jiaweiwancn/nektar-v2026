#!/usr/bin/env python3
"""CV3 bounded flat-floor calibration; preparation is never a developed flow.

The frozen calibration/crosscheck split and physical-time gates are enforced
here. Large native fields and logs belong in the separate Linux run root.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
from cao_tamura_experiment import (ROOT, DATA, THREADS, allocation, compare, config,
    environment, execute, frozen_data, generate, observe, prepare, records, run,
    sha, statistics, write_csv, write_json)
from generate_m4 import info, parameter
from generate_m3 import element
from cao_tamura_cv3_budgets import flat_momentum_budget
from terrain_continuity import audit as strip_continuity

VERSION = 'cao-tamura-cv3-v1'
CONTROLS = ('inlet_delta', 'epsilon_scale', 'k_scale', 'ks')
CENTRE = dict(inlet_delta=.21, epsilon_scale=1., k_scale=1., ks=.004)
BRACKETS = dict(inlet_delta=(.17, .25), epsilon_scale=(.5, 2.),
                k_scale=(.75, 1.5), ks=(.002, .008))
DRIFT_LIMITS = dict(velocity=2e-4, shear=2e-3, k=1e-3,
                    epsilon=1e-3, friction=.002, delta99=.002)


def design(surface):
    """Centre plus one-axis probes, declared before viewing fit residuals."""
    if surface not in ('S', 'R'):
        raise ValueError('Expected S or R')
    centre = {**CENTRE, 'ks': 0. if surface == 'S' else CENTRE['ks']}
    result = [dict(id='centre', parameters=centre)]
    for key in CONTROLS:
        if surface == 'S' and key == 'ks':
            continue
        for name, value in zip(('low', 'high'), BRACKETS[key]):
            result.append(dict(id=key+'-'+name, parameters={**centre, key: value}))
    return result


def probes(surface, path, crosscheck=False):
    """Native heights; d scenarios move observations, never physical walls."""
    frozen_data()
    rows = []
    for p in records(DATA/'profiles.csv'):
        if p['series_id'] not in [f'f3{a}_{surface}' for a in ('a', 'b', 'c', *(['d'] if crosscheck else []))]:
            continue
        z = .04*float(p['y_over_H'])
        ident = p['series_id']+':'+p['point_id']
        rows.append(dict(id=ident, x=0., span=0., z=z))
        if p['series_id'] == 'f3c_R':
            for shift in (-.001, .001):
                rows.append(dict(id=ident+f':d{.0035+shift:.4f}', x=0., span=0., z=z+shift))
    for x in (-.5, 0., .5):
        for i in range(391):
            rows.append(dict(id=f'edge:{x:g}:{i}', x=x, span=0., z=.022+i*(.8-.022)/390))
    rows.append(dict(id='incident-UH', x=0., span=0., z=.04))
    write_csv(path, rows, ['id', 'x', 'span', 'z'])
    write_json(str(path)+'.json', dict(version=VERSION, scope='reserved-flat' if crosscheck else 'calibration-only',
        cv1_snapshot_sha256=sha(DATA/'snapshot.json'), csv_sha256=sha(path),
        displacement='R 2.5/3.5/4.5mm observational scenarios; physical wall/Yap unchanged'))
    return rows


def thickness(rows, uinf, fraction=.99, x=0.):
    points = sorted((float(r['z']), float(r['u'])) for r in rows
                    if r['id'].startswith(f'edge:{x:g}:') and r['u'] != '')
    for i, (z, u) in enumerate(points):
        if u >= fraction*uinf and all(v >= fraction*uinf for _, v in points[i:]):
            if i == 0:
                return None  # below supported range: no extrapolation
            a, b = points[i-1], points[i]
            return a[0]+(fraction*uinf-a[1])*(b[0]-a[0])/(b[1]-a[1])
    return None


def wall_summary(case, time_end=None):
    """Independent native traction, both one-sided x=0 traces retained."""
    selected = []; all_points=[]
    for path in Path(case).parent.glob(Path(case).stem+'.terrain-points-r*.csv'):
        for r in records(path):
            if int(r['phase']) != 1 or (time_end is not None and float(r['time'])>time_end+1e-10):continue
            all_points.append(r)
            if abs(float(r['x'])) < 1e-10 and abs(float(r['span'])) < 1e-10:selected.append(r)
    if not selected:
        raise ValueError('Missing native incident wall traction')
    times = sorted({float(r['time']) for r in selected})
    samples = []
    for t in times:
        rs = [r for r in selected if float(r['time']) == t]
        # Global sample element gives a reproducible one-sided trace, independent of MPI ownership.
        rs.sort(key=lambda r: (int(r['sample_element']), int(r['point'])))
        values = [float(r['traction_x']) for r in rs]
        samples.append(dict(time=t, traction=values[0], trace_min=min(values), trace_max=max(values),
                            chi=float(rs[0]['chi']), reason=int(rs[0]['reason'])))
    last = samples[-1]
    if last['traction'] <= 0:
        raise ValueError('A reversed/nonpositive incident traction cannot define positive u_tau')
    final_points=[r for r in all_points if abs(float(r['time'])-times[-1])<1e-10]
    active=[r for r in final_points if float(r['chi'])>0]
    requested=[r for r in final_points if float(r['requested_chi'])>0]
    failed=[r for r in requested if int(r['reason'])>=2]
    return dict(ustar=math.sqrt(last['traction']), final=last,
                maximum_trace_jump=max(s['trace_max']-s['trace_min'] for s in samples),
                chi_range=[min(s['chi'] for s in samples), max(s['chi'] for s in samples)],
                fallback_reasons=sorted({s['reason'] for s in samples}),
                wall_correction=dict(final_time=times[-1],total_native_traces=len(final_points),
                    active_native_trace_fraction=len(active)/len(final_points),
                    requested_native_traces=len(requested),
                    fallback_fraction_of_requested=len(failed)/len(requested) if requested else None,
                    active_x_range=[min(float(r['x']) for r in active),max(float(r['x']) for r in active)] if active else None,
                    chi_range=[min(float(r['chi']) for r in final_points),max(float(r['chi']) for r in final_points)],
                    reason_codes=sorted({int(r['reason']) for r in final_points}),
                    fallback_x=sorted({float(r['x']) for r in failed}),
                    weighting='native trace counts, including one-sided shared-face copies; not an area fraction',
                    incident='x=0 is downstream Original AWF for the prescribed upstream-only EP envelope'),
                definition='positive native tangential wall traction; deterministic smallest global N-sample element',
                samples=samples)


def objective(stats, manifest, wall, displacement=None):
    """Five equal families on fixed weights; no k/epsilon/vertical-rms fitting."""
    m = json.loads(Path(manifest).read_text())
    if m['config']['geometry'] != 'floor':
        raise ValueError('Only empty-floor data may select calibration')
    ref = m['reference']; surface = m['config']['surface']
    if displacement is None: displacement=.0035 if surface=='R' else 0.
    if surface=='S' and displacement!=0.:raise ValueError('Smooth displacement is exactly zero')
    receipt = json.loads(Path(str(stats)+'.json').read_text())
    if receipt.get('output_sha256') != sha(stats) or receipt.get('experimental_fingerprint') != m['physics_fingerprint']:
        raise ValueError('Unbound native checkpoint statistics')
    rows = records(stats); predicted = {r['id']:r for r in rows}
    out = {}; full = []
    for panel in ('a', 'b', 'c'):
        ident = f'f3{panel}_{surface}'
        points = [p for p in records(DATA/'profiles.csv') if p['series_id'] == ident and float(p['weight']) > 0]
        target_rms = math.sqrt(sum(float(p['weight'])*float(p['value'])**2 for p in points)/sum(float(p['weight']) for p in points))
        scale = .03 if panel == 'a' else .1*(ref['ustar']/ref['Uinf'])**2 if panel == 'b' else .15*target_rms
        fit = []; missing = []; excluded = []; values = []
        for p in points:
            key = ident+':'+p['point_id']
            if panel == 'c' and surface == 'R' and displacement != .0035:
                key += f':d{displacement:.4f}'
            r = predicted.get(key); value = None
            if r and r['origin'] != 'checkpoint_statistics':
                raise ValueError('Prepared fields cannot enter calibration')
            if r and r['u'] != '' and panel == 'a':
                value = float(r['u'])/ref['Uinf']
            elif r and r['R_xz'] != '' and panel == 'b':
                value = -float(r['R_xz'])/ref['Uinf']**2
            elif r and r['R_xx'] != '' and panel == 'c' and float(r['R_xx']) >= 0:
                value = math.sqrt(float(r['R_xx']))/ref['ustar']
            weight = float(p['weight']); z = .04*float(p['y_over_H'])
            if panel == 'c' and surface == 'R': z += displacement-.0035
            if value is not None:
                values.append(dict(id=key, z=z, value=value, target=float(p['value']), bound=float(p['value_bound']), weight=weight))
            if z < m['config']['point_n']:
                excluded.append(key)
            elif value is None:
                missing.append(key)
            else:
                fit.append((weight, value-float(p['value'])))
        total = sum(float(p['weight']) for p in points)
        score = sum(w*(e/scale)**2 for w, e in fit)/total if fit and not missing else None
        out[ident] = dict(score=score, planning_scale=scale, retained_weight=sum(w for w, _ in fit),
            total_weight=total, excluded_below_N=excluded, missing_resolved=missing, values=values,
            conditional_rmse=math.sqrt(sum(w*e*e for w,e in fit)/sum(w for w,_ in fit)) if fit else None,
            source_bounds='digitization bounds only; source experimental uncertainty unreported')
        full.append(score)
    delta = thickness(rows, ref['Uinf']); star = wall['ustar']
    out['friction'] = dict(value=star, target=ref['ustar'], planning_scale=.05*ref['ustar'],
                           score=((star/ref['ustar']-1)/.05)**2)
    out['thickness'] = dict(delta99=delta, delta95=thickness(rows, ref['Uinf'], .95), target_approximate=.25,
                            planning_scale=.0125, score=((delta-.25)/.0125)**2 if delta is not None else None,
                            criterion='first sustained U>=fraction*Uinf; no stated source edge criterion')
    full += [out['friction']['score'], out['thickness']['score']]
    uh = float(predicted['incident-UH']['u']); rxx = predicted['incident-UH']['R_xx']
    return dict(version=VERSION, objective=sum(full)/5 if all(s is not None for s in full) else None,
        families=out, incident=dict(UH=uh, actual_ReH=.04*uh/ref['nu'], source_ReH=ref['re_h'],
        sigma_u_over_UH=math.sqrt(float(rxx))/uh if rxx and float(rxx)>=0 else None,
        source_sigma_u_over_UH=.0893 if surface=='S' else .151),
        displacement_m=displacement, weighting='equal five families; frozen profile weights divided by full frozen weight, without fit-band renormalization',
        inference='k/span variance and epsilon inferred; fixed canonical closure; no hill or reserved f3d scores',
        statistics_sha256=sha(stats), manifest_sha256=sha(manifest), cv1_snapshot_sha256=sha(DATA/'snapshot.json'))


def drift(old, new, ref, full_old, full_new, wall_old, wall_new):
    a = {r['id']:r for r in records(old)}; b = {r['id']:r for r in records(new)}
    if a.keys() != b.keys(): raise ValueError('Probe set changed during development')
    fa = records(full_old); fb = records(full_new)
    if len(fa) != len(fb) or any((r['element'],r['point']) != (s['element'],s['point']) for r,s in zip(fa,fb)):
        raise ValueError('Native full-field ordering/coverage changed')
    def diff(key, scale):
        return max((abs(float(b[k][key])-float(a[k][key]))/scale for k in a if a[k][key] and b[k][key]), default=0.)
    def fdiff(keys, scale):
        return max(abs(float(s[k])-float(r[k]))/scale for r,s in zip(fa,fb) for k in keys)
    def relative(key):
        if any(float(r[key])<=0 or float(s[key])<=0 for r,s in zip(fa,fb)):
            raise ValueError('Nonpositive field in development drift check')
        return max(abs(float(s[key])/float(r[key])-1) for r,s in zip(fa,fb))
    result = dict(velocity=max(diff('u',ref['Uinf']),diff('w',ref['Uinf']),fdiff(('u','v'),ref['Uinf'])),
        shear=diff('R_xz',ref['ustar']**2),
        k=relative('k'),epsilon=relative('eps_tilde' if 'eps_tilde' in fa[0] else 'eps'),
        friction=abs(wall_new['ustar']/wall_old['ustar']-1))
    da, db = thickness(list(a.values()),ref['Uinf']), thickness(list(b.values()),ref['Uinf'])
    result['delta99'] = abs(db/da-1) if da and db else None
    return dict(values=result, limits=DRIFT_LIMITS, accepted=all(v is not None and v<=DRIFT_LIMITS[k] for k,v in result.items()))


def assess(case, fields, directory, build, cpus):
    """Assess completed continuous native evolution using explicit checkpoints.

    This does not rerun, reset, or admit a failed/partial solve. At least three
    checkpoints must cover 3TF then two one-TF intervals after any warmup.
    """
    case=Path(case).resolve();directory=Path(directory).resolve()
    m=json.loads((case.parent/'manifest.json').read_text());c=m['config']
    receipt=json.loads((case.parent/'run.json').read_text())
    if (c['geometry']!='floor' or m['case_sha256']!=sha(case) or receipt.get('exit_code')!=0
            or not receipt.get('binaries_unchanged_during_run') or receipt.get('case_sha256')!=sha(case)
            or receipt.get('solver_sha256')!=sha(Path(build)/'IncNavierStokesSolver')):
        raise ValueError('A completed, bound native empty-floor run is required')
    if len(fields)<3:raise ValueError('Supply at least three evolved checkpoints')
    gate=receipt['stages']
    if (gate['final_startup_weight'] or gate['maximum_guard_activations'] or
            min(gate['minimum_k'],gate['minimum_epsilon'])<=0 or
            max(gate[k] for k in ('maximum_weak_divergence','maximum_projected_normal',
                'maximum_normal_integral_defect','maximum_resolved_mass_defect'))>1e-7):
        raise ValueError('Native evolution gate failed')
    residual=records(case.with_suffix('.rans-residual.csv'))
    free=max(float(r['free_coefficient_l2']) for r in residual if r['phase']=='1')
    if not math.isfinite(free) or free>1e-7:raise ValueError('Native full weak stage residual failed')
    continuity=strip_continuity(case.parent,case.stem)
    momentum=flat_momentum_budget(case)
    environment(build);allocation(1,cpus);directory.mkdir(parents=True,exist_ok=False)
    probes(c['surface'],directory/'points.csv')
    tf=(c['xmax']-c['xmin'])/m['reference']['Uinf']
    params={n.text.split('=',1)[0].strip():float(n.text.split('=',1)[1])
            for n in ET.parse(case).getroot().findall('CONDITIONS/PARAMETERS/P')}
    startup=params.get('RANSStartupDuration',0.)/tf
    state=dict(version=VERSION,status='assessing_continuous_run',complete=False,config=c,
        producing_case=str(case),producing_run_sha256=sha(case.parent/'run.json'),
        cv1_snapshot_sha256=sha(DATA/'snapshot.json'),stages=[])
    previous=None
    for i,field in enumerate(fields):
        field=Path(field).resolve();out=directory/f'frame-{i:02d}.csv'
        observe(case,field,directory/'points.csv',out,build,1,cpus)
        rows=records(out);time=float(rows[0]['time'])
        if any(r['origin']!='checkpoint' for r in rows) or time/tf<startup+3-1e-10:
            raise ValueError('Checkpoint precedes three unforced domain flow times')
        if previous and time-previous['time']<tf-1e-10:
            raise ValueError('Checkpoints must be separated by at least one domain flow time')
        full=directory/f'full-{i:02d}.csv'
        execute(allocation(1,cpus)+[Path(build)/'RANSWallFieldAudit',case,field,full],directory,f'full-{i:02d}')
        wall=wall_summary(case,time);write_json(directory/f'wall-{i:02d}.json',wall)
        item=dict(index=i,time=time,elapsed_flow_times=time/tf,checkpoint=str(field),
            checkpoint_sha256=sha(field) if field.is_file() else None,run=receipt,
            wall={k:v for k,v in wall.items() if k!='samples'},maximum_free_residual=free,
            strip_continuity=continuity,strip_momentum=momentum)
        if previous:item['drift']=drift(previous['out'],out,m['reference'],previous['full'],full,previous['wall'],wall)
        state['stages'].append(item);write_json(directory/'development.json',state)
        previous=dict(time=time,out=out,full=full,wall=wall)
    stable=all(s.get('drift',{}).get('accepted') for s in state['stages'][-2:])
    if stable:
        statistics([directory/f'frame-{i:02d}.csv' for i in range(len(fields)-3,len(fields))],directory/'statistics.csv')
        for d in ([.0025,.0035,.0045] if c['surface']=='R' else [0.]):
            write_json(directory/f'objective-d{d:.4f}.json',objective(directory/'statistics.csv',case.parent/'manifest.json',wall,d))
        state.update(status='developed_provisional_2D',complete=True,final_stage=str(case.parent),
            statistics_sha256=sha(directory/'statistics.csv'),
            stationarity='steady native drift check; CV4 convergence/long-time/3D gates pending')
    else:state['status']='not_developed_at_assessed_checkpoints'
    write_json(directory/'development.json',state)
    return state


def develop(value, directory, build, ranks, cpus, timeout=7200, max_extensions=10, numerics=None,
            probe_builder=None, stage_auditor=None, score=True,initial_unforced_flow_times=3,checkpoint_interval=None):
    """Native fixed-step 3TF + 1TF extensions; never accepts short/failing runs."""
    environment(build); allocation(ranks,cpus); directory=Path(directory).resolve(); directory.mkdir(parents=True,exist_ok=False)
    value=config(value)
    if value['geometry']!='floor' or value['dimension']!=2:
        raise ValueError('CV3 provisional calibration uses explicitly declared 2D empty floors')
    if value['inlet_reference_q'] is None or value['inlet_reference_kappa'] is None:
        raise ValueError('Declare common physical inlet reference Q/kappa independently of the closure')
    if initial_unforced_flow_times not in (1,3):raise ValueError('Initial unforced leg is one or three flow times')
    if checkpoint_interval is not None and (type(checkpoint_interval) is not int or checkpoint_interval<1):
        raise ValueError('Checkpoint interval must be a positive step count')
    numerics=numerics or {}
    allowed={'GlobalSysSoln','RANSNonlinearSolver','RANSNewtonMap','RANSNewtonDifferenceStep','RANSNewtonKrylovDimension','RANSNewtonMaxIterations',
             'IterativeSolverTolerance','RANSMatrixFreeTolerance','RANSNormalAugmentation','RANSPressureRelaxation','RANSPicardRelaxation','RANSGradDiv','PressureInverseTolerance','RANSNonlinearTolerance','RANSPressureTolerance'}
    if set(numerics)-allowed: raise ValueError('Unsupported undeclared numerical controls')
    (probe_builder or probes)(value['surface'],directory/'points.csv')
    state=dict(version=VERSION,status='developing',config=value,numerics=numerics,stages=[],
               cv1_snapshot_sha256=sha(DATA/'snapshot.json'),complete=False,
               purpose='calibration' if score else 'feasibility_only',initial_unforced_flow_times=initial_unforced_flow_times,
               checkpoint_interval=checkpoint_interval)
    write_json(directory/'development.json',state)
    previous=None; stable=0
    for index in range(max_extensions+1):
        folder=directory/f'stage-{index:02d}'; folder.mkdir()
        if previous is None:
            # Reference units are resolved by generate before physical duration is selected.
            folder.rmdir(); m=generate(value,folder); tf=(value['xmax']-value['xmin'])/m['reference']['Uinf']
            count=math.ceil((value['startup_flow_times']+initial_unforced_flow_times)*tf/value['dt'])
            tree=ET.parse(folder/'case.xml')
            for key,val in numerics.items():
                if key=='PressureInverseTolerance':
                    group=tree.getroot().find('CONDITIONS/GLOBALSYSSOLNINFO')
                    if group is None:group=element(tree.getroot().find('CONDITIONS'),'GLOBALSYSSOLNINFO')
                    field=element(group,'V',VAR='p')
                    element(field,'I',PROPERTY='IterativeSolverTolerance',VALUE=str(val))
                else:(info if key in ('RANSNonlinearSolver','GlobalSysSoln','RANSNewtonMap') else parameter)(tree.getroot(),key,val)
            parameter(tree.getroot(),'RANSStageAudit',1)
            parameter(tree.getroot(),'NumSteps',count); parameter(tree.getroot(),'IO_CheckSteps',min(count,checkpoint_interval or count))
            tree.write(folder/'case.xml',encoding='utf-8',xml_declaration=True)
            m['config'].update(steps=count,checkpoint=min(count,checkpoint_interval or count))
            m.update(case_sha256=sha(folder/'case.xml'),numerics=numerics)
            write_json(folder/'manifest.json',m); prepare(folder/'case.xml',build)
        else:
            duration=2 if index==1 and initial_unforced_flow_times==1 else 1
            count=math.ceil(duration*tf/value['dt'])
            restart=ROOT/'solvers/IncNavierStokesSolver/Examples/NeutralABL/restart.py'
            histories=[p for p in previous.glob('case.rans-*') if p.is_dir()]
            if not histories: raise ValueError('Missing full IMEX history')
            history_option=['--history-path',histories[0]] if len(histories)==1 else []
            # With frequent checkpoints, native field metadata identifies the
            # exact final history. Never guess by directory name or mtime.
            execute([sys.executable,restart,previous/'case.xml',previous/'case.fld',folder/'case.xml',
                '--steps',count,*history_option],folder,'restart')
            tree=ET.parse(folder/'case.xml'); parameter(tree.getroot(),'IO_CheckSteps',min(count,checkpoint_interval or count))
            tree.write(folder/'case.xml',encoding='utf-8',xml_declaration=True)
            m=json.loads((previous/'manifest.json').read_text());m['case_sha256']=sha(folder/'case.xml')
            m['config'].update(steps=count,checkpoint=min(count,checkpoint_interval or count))
            write_json(folder/'manifest.json',m)
        case=folder/'case.xml'
        try:
            receipt=run(case,build,ranks,cpus,timeout)
        except Exception as error:
            state.update(status='failed_or_timed_out',error=str(error));write_json(directory/'development.json',state)
            raise
        stage=receipt['stages']
        if stage['final_startup_weight'] or stage['maximum_guard_activations'] or stage['minimum_k']<=0 or stage['minimum_epsilon']<=0 or max(stage[k] for k in ('maximum_weak_divergence','maximum_projected_normal','maximum_normal_integral_defect','maximum_resolved_mass_defect'))>1e-7:
            raise ValueError('Native positivity/continuity gate failed')
        residual=records(case.with_suffix('.rans-residual.csv'))
        free=max(float(r['free_coefficient_l2']) for r in residual if r['phase']=='1')
        if not math.isfinite(free) or free>1e-7:raise ValueError('Native full weak stage residual failed')
        observe(case,folder/'case.fld',directory/'points.csv',folder/'incident.csv',build,1,cpus.split(',')[0])
        execute(allocation(1,cpus.split(',')[0])+[Path(build)/'RANSWallFieldAudit',case,folder/'case.fld',folder/'full.csv'],folder,'full')
        wall=wall_summary(case);write_json(folder/'wall.json',wall)
        item=dict(index=index,elapsed_flow_times=stage['final_time']/tf,run=receipt,wall={k:v for k,v in wall.items() if k!='samples'})
        item['unforced_elapsed_flow_times']=item['elapsed_flow_times']-value['startup_flow_times']
        item.update(maximum_free_residual=free,strip_continuity=strip_continuity(folder,'case'),
                    strip_momentum=flat_momentum_budget(case))
        if stage_auditor:item['route_audit']=stage_auditor(case)
        if previous:
            prior=json.loads((previous/'wall.json').read_text())
            item['drift']=drift(previous/'incident.csv',folder/'incident.csv',m['reference'],previous/'full.csv',folder/'full.csv',prior,wall)
            # Only one-TF extensions following at least three unforced TF
            # qualify for the two successive development drift checks.
            prior_flow=state['stages'][-1]['elapsed_flow_times']-value['startup_flow_times']
            stable=stable+1 if item['drift']['accepted'] and prior_flow>=3.-1e-10 else 0
        state['stages'].append(item);write_json(directory/'development.json',state)
        if stable>=2:
            inputs=[directory/f'stage-{j:02d}'/'incident.csv' for j in range(index-2,index+1)]
            statistics(inputs,directory/'statistics.csv')
            if score:
                for d in ([.0025,.0035,.0045] if value['surface']=='R' else [0.]):
                    write_json(directory/f'objective-d{d:.4f}.json',objective(directory/'statistics.csv',folder/'manifest.json',wall,d))
            state.update(status='developed_provisional_2D',complete=True,final_stage=str(folder),
                statistics_sha256=sha(directory/'statistics.csv'),stationarity='steady drift check; CV4 full stationarity/3D/convergence pending')
            write_json(directory/'development.json',state)
            return state
        previous=folder
    state.update(status='not_developed_within_budget');write_json(directory/'development.json',state)
    return state


def select(campaign):
    """Select only developed common candidates, retaining every attempted point."""
    campaign=Path(campaign); table=[]
    for surface in ('S','R'):
        for candidate in design(surface):
            item=dict(surface=surface,**candidate,models={})
            for model in ('standard_original','cls_original'):
                folder=campaign/surface/candidate['id']/model
                state_path=folder/'development.json'
                status=json.loads(state_path.read_text()) if state_path.exists() else {'status':'unexecuted','complete':False}
                if 'config' in status and any(status['config'][k]!=v for k,v in candidate['parameters'].items()):
                    raise ValueError('Recorded candidate differs from the declared design: '+str(folder))
                result=dict(status=status['status'],complete=status['complete'],directory=str(folder))
                if status['complete']:
                    result['objective']=json.loads((folder/('objective-d0.0035.json' if surface=='R' else 'objective-d0.0000.json')).read_text())
                    result['development_sha256']=sha(state_path)
                item['models'][model]=result
            js=[r['objective']['objective'] for r in item['models'].values() if r['complete']]
            item['common_objective']=sum(js)/2 if len(js)==2 and all(j is not None and math.isfinite(j) for j in js) else None
            table.append(item)
    selections={}
    for surface in ('S','R'):
        usable=[r for r in table if r['surface']==surface and r['common_objective'] is not None]
        selections[surface]=min(usable,key=lambda r:r['common_objective'])['id'] if usable else None
    return dict(version=VERSION,candidates=table,selected=selections,
        design_executed=all(r['status'] in ('developed_provisional_2D','failed_or_timed_out','failed_input_or_audit','not_developed_within_budget') for t in table for r in t['models'].values()),
        criterion='common five-family objective averaged across canonical standard and CLS; developed cases only',
        bounded_search='centre plus declared coordinate probes; failed points are inadmissible; no global optimum claim',
        cv1_snapshot_sha256=sha(DATA/'snapshot.json'))


def freeze(campaign, transfers, output, calibration_version):
    """Seal a completed selection and all wall transfers before reserved checks.

    transfers.json maps S/R -> model name -> developed case directory. Every
    primary policy must be tested with the same surface/inlet controls. This
    deliberately cannot freeze an undeveloped or failed numerical campaign.
    """
    selection=select(campaign)
    if not selection['design_executed'] or any(v is None for v in selection['selected'].values()):
        raise ValueError('Execute the bounded design and develop a common S/R selection before freezing')
    output=Path(output)
    if output.exists(): raise ValueError('Use a new immutable calibration version directory')
    transfers=json.loads(Path(transfers).read_text()); sealed={}
    expected=('standard_original','cls_original','standard_fixed','cls_fixed','standard')
    # Validate everything before creating a reviewable freeze.
    for surface in ('S','R'):
        selected=next(r for r in selection['candidates'] if r['surface']==surface and r['id']==selection['selected'][surface])
        if set(transfers[surface])!=set(expected):
            raise ValueError('Supply all five primary wall/closure transfers for each surface')
        base=None; numerical_base=None; cases={}
        for model in expected:
            directory=Path(transfers[surface][model]);state=json.loads((directory/'development.json').read_text())
            if not state['complete'] or state['status']!='developed_provisional_2D' or len(state['stages'])<3:
                raise ValueError('Undeveloped transfer '+surface+'/'+model)
            if not all(t['drift']['accepted'] for t in state['stages'][-2:]):
                raise ValueError('Missing successive physical-time drift gates')
            startup=state['config'].get('startup_flow_times',0.)
            if state['stages'][0]['elapsed_flow_times']<startup+3 or state['stages'][-1]['elapsed_flow_times']<startup+5:
                raise ValueError('Short physical-time runs cannot be frozen as developed precursors')
            last=Path(state['final_stage']);m=json.loads((last/'manifest.json').read_text())
            producing=state['stages'][-1]['run']
            if (m['case_sha256']!=sha(last/'case.xml') or
                    producing.get('case_sha256')!=m['case_sha256'] or
                    producing.get('exit_code')!=0 or not producing.get('binaries_unchanged_during_run')):
                raise ValueError('Unbound or changed producing case')
            floor=ET.parse(last/'case.xml').getroot().find("CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF='0']")
            if any(n.attrib.get('VALUE')!='EpsilonWall' for n in floor if n.attrib.get('VAR') in ('eps','eps_tilde')):
                raise ValueError('Case predates the corrected startup wall trace')
            if m['config']['surface']!=surface or m['config']['model']!=model or m['config']['geometry']!='floor':
                raise ValueError('Wrong transfer surface/model/geometry')
            if any(m['config'][k]!=v for k,v in selected['parameters'].items()):
                raise ValueError('Transfer retuned selected physical parameters')
            common={k:v for k,v in m['config'].items() if k not in ('model','steps','checkpoint')}
            if base is None:base=common
            elif common!=base:raise ValueError('Transfer inlet, surface or numerical mesh differs')
            if numerical_base is None:numerical_base=m.get('numerics',{})
            elif numerical_base!=m.get('numerics',{}):raise ValueError('Transfer nonlinear/stabilization controls differ')
            if model.endswith('_original') and directory.resolve()!=Path(selected['models'][model]['directory']).resolve():
                raise ValueError('Original transfer must be the selected candidate, not another realization')
            if sha(directory/'statistics.csv')!=state['statistics_sha256']:
                raise ValueError('Native statistics changed since development')
            statistics_receipt=json.loads((directory/'statistics.csv.json').read_text())
            if statistics_receipt.get('experimental_fingerprint')!=m['physics_fingerprint']:
                raise ValueError('Statistics belong to another experimental input')
            m.update(calibration_complete=True,calibration_version=calibration_version+'-'+surface,
                calibration_role='shared surface/inlet; canonical model and policy-specific epsilon trace',
                development_sha256=sha(directory/'development.json'),statistics_sha256=state['statistics_sha256'],
                numerical_status='provisional 2D development; CV4 h/p/dt/domain/3D uncertainty remains pending',
                selected_parameters=selected['parameters'],native_development_directory=str(directory),
                heldout_hill_scores_used=False,reserved_f3d_used_for_selection=False)
            cases[model]=m
        sealed[surface]=dict(calibration_version=calibration_version+'-'+surface,calibration_complete=True,
            selected_candidate=selected['id'],parameters=selected['parameters'],common_config=base,
            numerics=numerical_base,cases=cases,conditional_fits=[],canonical_model_misfit='retained in native objective and crosschecks',
            uncertainty=dict(displacement_m=[.0025,.0035,.0045] if surface=='R' else [0.],
                incoming_k_multiplier=list(BRACKETS['k_scale']),epsilon_multiplier=list(BRACKETS['epsilon_scale']),
                numerical='pending CV4; bounded parameter ranges are not confidence intervals'),
            reserved_f3d_used_for_selection=False,heldout_hill_scores_used=False)
    output.mkdir(parents=True)
    write_json(output/'selection.json',selection)
    for surface,manifest in sealed.items():
        write_json(output/(surface+'.json'),manifest)
        for model,m in manifest['cases'].items():write_json(output/(surface+'-'+model+'.json'),m)
    write_json(output/'snapshot.json',dict(version=calibration_version,
        cv1_snapshot_sha256=sha(DATA/'snapshot.json'),
        files={p.name:sha(p) for p in sorted(output.glob('*.json'))},
        sealed_before_reserved_flat_crosschecks=True))
    return sealed


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    d=sub.add_parser('design');d.add_argument('surface',choices=('S','R'));d.add_argument('output',type=Path)
    q=sub.add_parser('develop');q.add_argument('configuration',type=Path);q.add_argument('output',type=Path)
    q.add_argument('--build',type=Path,required=True);q.add_argument('--ranks',type=int,default=4)
    q.add_argument('--cpus',required=True);q.add_argument('--timeout',type=int,default=7200)
    q.add_argument('--max-extensions',type=int,default=10);q.add_argument('--numerics',type=Path)
    q=sub.add_parser('select');q.add_argument('campaign',type=Path);q.add_argument('output',type=Path)
    q=sub.add_parser('assess');q.add_argument('case',type=Path);q.add_argument('output',type=Path)
    q.add_argument('fields',type=Path,nargs='+');q.add_argument('--build',type=Path,required=True);q.add_argument('--cpus',required=True)
    q=sub.add_parser('freeze');q.add_argument('campaign',type=Path);q.add_argument('transfers',type=Path)
    q.add_argument('output',type=Path);q.add_argument('--version',required=True)
    args=p.parse_args()
    if args.command=='design':
        write_json(args.output,dict(version=VERSION,surface=args.surface,candidates=design(args.surface)))
    elif args.command=='select':
        write_json(args.output,select(args.campaign))
    elif args.command=='freeze':
        freeze(args.campaign,args.transfers,args.output,args.version)
    elif args.command=='assess':
        result=assess(args.case,args.fields,args.output,args.build,args.cpus)
        print(json.dumps(dict(status=result['status'],complete=result['complete'],
                             drifts=[s.get('drift') for s in result['stages']]),indent=2))
        return 0 if result['complete'] else 1
    else:
        result=develop(json.loads(args.configuration.read_text()),args.output,args.build,args.ranks,args.cpus,
            args.timeout,args.max_extensions,json.loads(args.numerics.read_text()) if args.numerics else None)
        print(json.dumps(result,indent=2));return 0 if result['complete'] else 1
    return 0

if __name__=='__main__':
    raise SystemExit(main())
