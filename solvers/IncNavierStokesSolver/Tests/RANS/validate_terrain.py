#!/usr/bin/env python3
"""Bounded Linux verification of terrain frames, wall constraints and ramps."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
import os
import signal
from queue import Queue
from pathlib import Path
import subprocess
import sys
import threading
import tempfile
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
os.environ.update(THREADS)
from generate_terrain import plane, hill
from generate_m4 import parameter, info
from generate_m6 import original_j0
from terrain_initial import prepare as prepare_initial
from terrain_continuity import audit as audit_continuity


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def bounded_solver(command,directory,log,timeout):
    """Stop the entire owned MPI process group on timeout; retain a receipt."""
    with subprocess.Popen(command,cwd=directory,env=os.environ.copy(),stdout=log,
                          stderr=subprocess.STDOUT,start_new_session=True) as process:
        try:
            return process.wait(timeout=timeout),False
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL)
                process.wait()
            return 124,True


def rows(path):
    with Path(path).open() as f:
        return [{k:v if k == 'field' else float(v) for k,v in row.items()} for row in csv.DictReader(f)]


def fields(path):
    return {(int(r['element']),int(r['point'])):r for r in rows(path)}


def corner_compatibility(directory):
    """Expose both one-sided continuity targets at each compact C1 foot."""
    groups = {}
    for path in Path(directory).glob('neutral.terrain-points-r*.csv'):
        for r in rows(path):
            if r['phase'] == 1 and abs(abs(r['x'])-.1) < 1e-10:
                groups.setdefault((r['time'],round(r['x'],10)),[]).append(r)
    result = []
    for (time,x),points in sorted(groups.items()):
        targets = [p['normal_target'] for p in points]
        actual = [p['normal_target']+p['normal_error'] for p in points]
        jump = max(targets)-min(targets)
        maximum = max(abs(p['normal_error']) for p in points)
        if max(actual)-min(actual) > 1e-9 or maximum+1e-9 < jump/2:
            raise RuntimeError('Inconsistent C1 corner compatibility audit')
        result.append(dict(time=time,x=x,target_jump=jump,
            unavoidable_maximum_for_continuous_trace=jump/2,
            measured_maximum=maximum,actual_trace_jump=max(actual)-min(actual)))
    if len(result) < 2: raise RuntimeError('Both C1 hill feet must be audited')
    return result


def compare(a, b, angle=0., dimension=2):
    """Compare every saved physical point, rotating Cartesian vectors/tensors."""
    first, second = fields(a), fields(b)
    if first.keys() != second.keys():
        raise RuntimeError('Field grids differ')
    error, scale = {}, {}
    c, s = math.cos(angle), math.sin(angle)
    vertical = 'v' if dimension == 2 else 'w'
    for key, reference in first.items():
        actual = dict(second[key])
        # Convert the rotated result back into the first problem's frame.
        actual['u'], actual[vertical] = c*actual['u']+s*actual[vertical], -s*actual['u']+c*actual[vertical]
        R = [[actual['rans_R_xx'],actual['rans_R_xy'],actual['rans_R_xz']],
             [actual['rans_R_xy'],actual['rans_R_yy'],actual['rans_R_yz']],
             [actual['rans_R_xz'],actual['rans_R_yz'],actual['rans_R_zz']]]
        Q = [[c,0,s],[0,1,0],[-s,0,c]]
        rotated = [[sum(Q[i][k]*R[k][l]*Q[j][l] for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
        for name,i,j in [('xx',0,0),('yy',1,1),('zz',2,2),('xy',0,1),('xz',0,2),('yz',1,2)]:
            actual['rans_R_'+name] = rotated[i][j]
        for name,value in reference.items():
            if name in ('element','point','x','y','z','time'):
                continue
            error[name] = max(error.get(name,0),abs(actual[name]-value))
            scale[name] = max(scale.get(name,1),abs(value))
    return {name:value/scale[name] for name,value in error.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--suite',choices=('planes','curved','periodic-curved','transitions','pilot','pilot-2d','preservation'),default='planes')
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--publish-to',type=Path,help='Publish completed pilot-2d cases as immutable cache aliases for a concurrent full pilot')
    args = p.parse_args()
    if not sys.platform.startswith('linux') or not 1 <= args.jobs <= 2:
        raise RuntimeError('Use Linux and at most two concurrent cases (at most eight ranks each)')
    cpus = sorted(os.sched_getaffinity(0))[:16]
    os.sched_setaffinity(0,cpus)
    cpu_leases = Queue()
    group_size = len(cpus)//args.jobs
    required = 8 if args.suite == 'pilot' else 4
    if group_size < required:
        raise RuntimeError('Each case needs at least '+str(required)+' available CPUs; reduce --jobs or increase the allocation')
    for j in range(args.jobs):
        cpu_leases.put(cpus[j*group_size:(j+1)*group_size])
    build, output = args.build.resolve(), args.output.resolve()
    geometry_library = build/'geometry-lib/libSpatialDomains.so'
    if not geometry_library.exists():
        raise RuntimeError('Build the source geometry library before terrain validation')
    os.environ['LD_LIBRARY_PATH'] = str(geometry_library.parent)+':'+os.environ.get('LD_LIBRARY_PATH','')
    source = Path(__file__).resolve().parents[4]
    if output.is_relative_to(source):
        raise RuntimeError('Store runs outside the source checkout')
    output.mkdir(parents=True,exist_ok=True)
    if args.publish_to:
        args.publish_to = args.publish_to.resolve()
        if args.suite != 'pilot-2d' or args.publish_to.is_relative_to(source) or args.publish_to == output:
            raise RuntimeError('Publish pilot-2d cases to a distinct run directory outside the source checkout')
        args.publish_to.mkdir(parents=True,exist_ok=True)
    binaries = {name:digest(build/name) for name in ('IncNavierStokesSolver','RANSWallFieldAudit')}
    binaries['libSpatialDomains.so'] = digest(geometry_library)
    cases, comparisons = {}, {}
    lock = threading.Lock()

    def run_case(name, kind='plane', ranks=1, modify=None, expected_rejection=None, cpu_ids=None, **options):
        if ranks not in (1,2,4,8) or len(cpu_ids) < ranks:
            raise RuntimeError('Unsupported rank allocation')
        directory = output/name; directory.mkdir(exist_ok=True)
        session = directory/'neutral.xml'
        receipt = directory/'receipt.json'
        candidate = tempfile.TemporaryDirectory(prefix='definition-',dir=directory) if receipt.exists() else None
        definition = Path(candidate.name)/'neutral.xml' if candidate else session
        preparation = (plane if kind == 'plane' else hill)(definition,**options)
        tree = ET.parse(definition)
        parameter(tree.getroot(),"RANSDebug",1)
        if modify:
            modify(tree.getroot())
            tree.write(definition,encoding='utf-8',xml_declaration=True)
        tree.write(definition,encoding='utf-8',xml_declaration=True)
        record = dict(preparation=preparation,ranks=ranks,threads=THREADS,generator_fixture_sha256=digest(definition),binaries=binaries)
        if receipt.exists():
            old = json.loads(receipt.read_text())
            if any(old.get(k) != v for k,v in record.items()) or not old.get('passed'):
                raise RuntimeError('Use a fresh run directory; retain failed/stale evidence: '+str(directory))
            if digest(session) != old['fixture_sha256']:
                raise RuntimeError('Accepted fixture changed: '+str(session))
            if old.get('initial_audit') and digest(directory/'initial.fld') != old['initial_audit']['field_sha256']:
                raise RuntimeError('Accepted initial field changed: '+name)
            candidate.cleanup()
            with lock:
                cases[name] = old
            return old
        allocation = ['taskset','-c',','.join(map(str,cpu_ids[:ranks]))]
        command = allocation+([str(build/'IncNavierStokesSolver'),'neutral.xml'] if ranks == 1 else
                   ['mpirun','--bind-to','none','-np',str(ranks),str(build/'IncNavierStokesSolver'),'neutral.xml'])
        record['cpus'] = cpu_ids[:ranks]
        if kind != 'plane' and not expected_rejection:
            try:
                record['initial_audit'] = prepare_initial(session,build,cpu_ids[0])
            except Exception as error:
                preparation_command = ['taskset','-c',str(cpu_ids[0]),str(build/'RANSWallFieldAudit'),
                    str(session),'--terrain-initial',str(directory/'initial.fld')]
                record.update(passed=False,failure_phase='initial preparation',
                    command=list(getattr(error,'cmd',preparation_command)),
                    exit_code=getattr(error,'returncode',124 if isinstance(error,subprocess.TimeoutExpired) else 1),error=str(error))
                receipt.write_text(json.dumps(record,indent=2)+'\n')
                raise
        record['fixture_sha256'] = digest(session)
        # Full-hex hill stages are substantially more expensive, especially
        # with centered CLS Jacobian probes. Runtime limits do not alter the
        # independently checked physical convergence criteria.
        timeout = 7200 if kind == 'hill' and options.get('dimension',2) == 3 else (
            1800 if options.get('model','').startswith('cls') or options.get('dimension',2) == 3
            or (kind == 'hill' and options.get('refinement',1) >= 4) else 600)
        record['timeout_seconds'] = timeout
        started = time.monotonic()
        with (directory/'solver.log').open('w') as log:
            code,timed_out = bounded_solver(command,directory,log,timeout)
        record.update(command=command,exit_code=code,timeout_expired=timed_out,elapsed=time.monotonic()-started,passed=False)
        receipt.write_text(json.dumps(record,indent=2)+'\n')
        if expected_rejection:
            if code == 0 or timed_out or expected_rejection not in (directory/'solver.log').read_text():
                raise RuntimeError('Expected rejection did not occur: '+name)
            record['passed'] = True; record['expected_rejection'] = expected_rejection
        else:
            if code:
                raise RuntimeError('Solver failed: '+str(directory/'solver.log'))
            command = ['taskset','-c',str(cpu_ids[0]),str(build/'RANSWallFieldAudit'),str(session),str(directory/'neutral.fld'),str(directory/'field.csv')]
            with (directory/'audit.log').open('w') as log:
                subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=120)
            record['audit_command'] = command
            stages = rows(directory/'neutral.rans.csv')
            record['completed_stages'] = len(stages)
            record['maximum_weak_divergence'] = max(r['weak_divergence_l2'] for r in stages)
            if record['maximum_weak_divergence'] > 1e-7 or min(r['min_k'] for r in stages) <= 0 or min(r['min_eps'] for r in stages) <= 0:
                raise RuntimeError('Stage/positivity gate failed: '+name)
            wall = directory/'neutral.terrain-wall.csv'
            if wall.exists():
                data = rows(wall)
                record['normal_projected_max'] = max(r['normal_projected_max'] for r in data)
                record['normal_pointwise_max'] = max(r['normal_pointwise_max'] for r in data)
                record['normal_l2'] = max(r['normal_l2'] for r in data)
                record['normal_away_feet_max'] = max(r['normal_away_feet_max'] for r in data)
                record['epsilon_projection_max'] = max(r['epsilon_projection_max'] for r in data)
                record['remote_samples'] = max(r['remote_samples'] for r in data)
                record['wall_ranks'] = max(r.get('wall_ranks',0) for r in data)
                record['remote_without_wall_samples'] = max(r.get('remote_N_without_wall_samples',0) for r in data)
                record['projection_audit'] = max(r['projection_audit'] for r in data)
                record['N_inverse_residual'] = max(r['N_inverse_residual'] for r in data)
                record['tangent_flux_residual_l2'] = max(r['tangent_flux_residual_l2'] for r in data)
                record['k_flux_residual_l2'] = max(r['k_flux_residual_l2'] for r in data)
                record['normal_integral_defect'] = max(abs(r['normal_integral_defect']) for r in data)
                record['resolved_mass_defect'] = max(abs(r['resolved_mass_defect']) for r in data)
                if record['normal_projected_max'] > 1e-7:
                    raise RuntimeError('Normal constraint gate failed: '+name)
                if kind != 'plane':
                    expected_volume = (preparation.get('length',1.6)*.3-preparation['amplitude']*.1)*(.032 if preparation['dimension'] == 3 else 1)
                    record['volume_relative_error'] = max(abs(r['resolved_volume']+r['unresolved_volume']-expected_volume)/expected_volume for r in data)
                    if record['volume_relative_error'] > 1e-8 or max(record['normal_integral_defect'],record['resolved_mass_defect']) > 1e-7:
                        raise RuntimeError('Geometric volume/mass budget gate failed: '+name)
            residuals = directory/'neutral.rans-residual.csv'
            if residuals.exists():
                audit = [r for r in rows(residuals) if r['phase'] == 1]
                record['maximum_free_stage_residual'] = max(r['free_coefficient_l2'] for r in audit)
                if record['maximum_free_stage_residual'] > 1e-7:
                    raise RuntimeError('Native free-row stage audit failed: '+name)
            record['passed'] = True
        receipt.write_text(json.dumps(record,indent=2)+'\n')
        with lock:
            cases[name] = record
            (output/'progress.json').write_text(json.dumps(dict(cases=cases,comparisons=comparisons),indent=2)+'\n')
        print(name,'PASS',flush=True)
        return record

    def run(name, **options):
        if args.suite in ('curved','periodic-curved'):
            options['dt'] = 2e-4 if options.get('model') == 'standard' else 2e-5
            options['periodic'] = args.suite == 'periodic-curved'
        allocation = cpu_leases.get()
        try:
            record = run_case(name,cpu_ids=allocation,**options)
            if args.publish_to:
                try:
                    # Atomic, no replacement: never overwrite a case already
                    # started by the other runner. Original paths stay valid.
                    (args.publish_to/name).symlink_to(output/name,target_is_directory=True)
                except FileExistsError:
                    pass
            return record
        finally:
            cpu_leases.put(allocation)

    if args.suite == 'planes':
        models = ('standard','standard_original','standard_fixed','standard_switched','cls_original','cls_fixed','cls_switched')
        jobs = []
        for model in models:
            # Original CLS has a rapid unbalanced wall-dissipation transient.
            # Resolve it with the same timestep as the hill pilot; retain the
            # complete primary/derived-field covariance gate for every model.
            step = {'dt':.0002} if model == 'cls_original' else {}
            jobs += [(model+'_flat',dict(model=model,**step)),(model+'_rotated',dict(model=model,angle=.35,**step))]
        for model in ('standard','standard_fixed','cls_fixed'):
            jobs.append((model+'_legacy',dict(model=model,legacy=True)))
        for model in ('standard_fixed','cls_fixed'):
            jobs += [(model+'_3d_flat',dict(model=model,dimension=3,nx=2,layers=4)),
                     (model+'_3d_rotated',dict(model=model,dimension=3,nx=2,layers=4,angle=.35))]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run,name,**options) for name,options in jobs]
            for future in futures:
                future.result()
        for model in models:
            errors = compare(output/(model+'_flat')/'field.csv',output/(model+'_rotated')/'field.csv',angle=.35)
            comparisons[model+'_rotation'] = errors
            if max(errors.values()) > 1e-6:
                raise RuntimeError('Rotation covariance failed: '+model+': '+str(errors))
        for model in ('standard','standard_fixed','cls_fixed'):
            errors = compare(output/(model+'_legacy')/'field.csv',output/(model+'_flat')/'field.csv')
            comparisons[model+'_flat_reduction'] = errors
            if max(errors.values()) > 1e-6:
                raise RuntimeError('Flat reduction failed: '+model+': '+str(errors))
        for model in ('standard_fixed','cls_fixed'):
            errors = compare(output/(model+'_3d_flat')/'field.csv',output/(model+'_3d_rotated')/'field.csv',angle=.35,dimension=3)
            comparisons[model+'_3d_rotation'] = errors
            if max(errors.values()) > 1e-6:
                raise RuntimeError('3D rotation covariance failed: '+model+': '+str(errors))
    if args.suite in ('curved','periodic-curved'):
        def verification(root):
            parameter(root,'RANSStageAudit',1)
            parameter(root,'RANSNonlinearTolerance',1e-8)
            parameter(root,'RANSPressureTolerance',1e-8)
        jobs = []
        for model in ('standard','standard_fixed'):
            for refinement in (1,2,4):
                jobs.append((model+'_h'+str(refinement),dict(kind='hill',model=model,steps=1,refinement=refinement,ranks=4,modify=verification)))
        for model in ('standard_original','cls_original','cls_fixed'):
            jobs.append((model+'_curved',dict(kind='hill',model=model,steps=1,ranks=4,modify=verification)))
        for augmentation in (1.,4.):
            def control(root, value=augmentation):
                verification(root); parameter(root,'RANSNormalAugmentation',value)
            jobs.append(('standard_augmentation_'+str(int(augmentation)),dict(kind='hill',model='standard_fixed',steps=1,ranks=4,modify=control)))
        def iterative_projection(root):
            verification(root); parameter(root,'RANSNormalDirectLimit',0)
        def half_difference(root):
            verification(root); parameter(root,'RANSNewtonDifferenceStep',5e-7)
        jobs.append(('standard_iterative_projection',dict(kind='hill',model='standard_fixed',steps=1,ranks=4,modify=iterative_projection)))
        jobs.append(('cls_half_difference',dict(kind='hill',model='cls_fixed',steps=1,ranks=4,modify=half_difference)))
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run,name,**options) for name,options in jobs]
            for future in futures:
                future.result()
        for model in ('standard','standard_fixed'):
            first, second, third = (cases[model+'_h'+str(level)] for level in (1,2,4))
            metrics = ('normal_pointwise_max','normal_l2','tangent_flux_residual_l2','k_flux_residual_l2')
            comparison = {name:dict(coarse=first[name],middle=second[name],fine=third[name],
                coarse_to_fine=first[name]/third[name],last_order=math.log2(second[name]/third[name])) for name in metrics}
            comparisons[model+'_refinement'] = comparison
            required = metrics if args.suite == 'periodic-curved' else metrics[1:]
            if any(comparison[name]['fine'] >= min(comparison[name]['coarse'],comparison[name]['middle']) for name in required):
                raise RuntimeError('Curved normal/flux refinement did not decrease: '+model+': '+str(comparison))
            if args.suite == 'curved':
                comparisons[model+'_C1_corner_compatibility'] = {
                    'pointwise_maximum_decreases':comparison['normal_pointwise_max']['fine'] < min(
                        comparison['normal_pointwise_max']['coarse'],comparison['normal_pointwise_max']['middle']),
                    'interpretation':'C1-foot one-sided targets and a continuous velocity trace; maximum convergence is not claimed',
                    'levels':{str(level):corner_compatibility(output/(model+'_h'+str(level))) for level in (1,2,4)}}
        for augmentation in (1,4):
            errors = compare(output/'standard_fixed_h1/field.csv',output/('standard_augmentation_'+str(augmentation))/'field.csv')
            comparisons['augmentation_'+str(augmentation)] = errors
            if max(errors.values()) > 1e-6:
                raise RuntimeError('Normal augmentation changed the converged wall problem: '+str(errors))
        for baseline,control in [('standard_fixed_h1','standard_iterative_projection'),('cls_fixed_curved','cls_half_difference')]:
            errors = compare(output/baseline/'field.csv',output/control/'field.csv')
            comparisons[control] = errors
            if max(errors.values()) > 1e-6:
                raise RuntimeError('Converged iteration-control agreement failed: '+control+': '+str(errors))
    if args.suite == 'transitions':
        from terrain_campaigns import transition_campaign
        transition_campaign(run,cases,comparisons,output,args.jobs,compare)
    if args.suite in ('pilot','pilot-2d'):
        from terrain_campaigns import pilot_campaign
        pilot_campaign(run,cases,comparisons,output,args.jobs,compare,include_3d=args.suite == 'pilot')
    if args.suite == 'preservation':
        # Five flow-through clocks on a periodic flat domain test stationarity
        # of the new terrain adapter; the inlet/outlet benchmark is separate.
        for model in ('standard','standard_fixed','cls_fixed','standard_switched','cls_switched'):
            q = 3.190675081624288 if model.startswith('cls') else 10/3
            kappa = .440553905186659 if model.startswith('cls') else math.sqrt(.48*1.3/q)
            matching = .3/kappa*math.log(.008/.0002) if model == 'standard' else .09*original_j0(q*.09,.008,.002,1.5e-5)
            duration = 5*.4/(matching+.3/kappa*math.log(10))
            def settings(root):
                parameter(root,'IO_CheckSteps',40)
            run(model,model=model,nx=1,layers=12,dt=duration/40,steps=40,ranks=4,modify=settings)
            root = ET.parse(output/model/'neutral.xml').getroot()
            params = {p.text.split('=',1)[0].strip():float(p.text.split('=',1)[1]) for p in root.findall('CONDITIONS/PARAMETERS/P')}
            errors = dict(u=0.,v=0.,k=0.,epsilon=0.)
            for row in rows(output/model/'field.csv'):
                z = row['y'];star = params['Ustar'];kap = params['Kappa']
                expected = dict(u=params['MatchingVelocity']+star/kap*math.log(z/params['H']),
                    v=0.,k=params['Qstar']*star**2,epsilon=star**3/(kap*z))
                for name,value in expected.items():
                    actual = row[('eps_tilde' if model.startswith('cls') else 'eps') if name == 'epsilon' else name]
                    errors[name] = max(errors[name],abs(actual-value)/(star if name == 'v' else abs(value)))
            comparisons[model+'_five_flow_reference'] = dict(duration=duration,flow_clocks=5,
                periodic_streamwise=True,normalized_linf=errors,engineering_target=1e-3)
            if max(errors.values()) > 1e-3:
                raise RuntimeError('Terrain-plane neutral preservation exceeded 0.1 percent: '+model+': '+str(errors))
            if model != 'standard' and not all(row['chi_min'] == row['chi_max'] == 1
                    for row in rows(output/model/'neutral.terrain-wall.csv')):
                raise RuntimeError('Terrain reference left the correction plateau: '+model)
        for closure in ('standard','cls'):
            errors = compare(output/(closure+'_fixed')/'field.csv',output/(closure+'_switched')/'field.csv')
            comparisons[closure+'_five_flow_fixed_switched'] = errors
            if max(errors.values()) > 1e-6:
                raise RuntimeError('Fixed/switched terrain reference fields disagree: '+closure)
    for name,record in cases.items():
        if not record.get('expected_rejection') and (output/name/'neutral.terrain-wall.csv').exists():
            comparisons[name+'_continuity'] = audit_continuity(output/name)
    if args.suite in ('curved','periodic-curved'):
        for model in ('standard_fixed','cls_fixed'):
            if comparisons[model+'_curved_continuity' if model == 'cls_fixed' else model+'_h1_continuity']['maximum_transpiration'] < 1e-6:
                raise RuntimeError('Curved AWF did not exercise nonzero transpiration: '+model)
    if args.suite == 'transitions':
        for model in ('standard','cls'):
            if comparisons[model+'_span_continuity_continuity']['maximum_span_flux_divergence'] < 1e-6:
                raise RuntimeError('Independent audit did not exercise spanwise strip-flux divergence: '+model)
    evidence = dict(status='PASS',suite=args.suite,cases=cases,comparisons=comparisons,binaries=binaries,
                    source_revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip())
    (output/'evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print('PASS',len(cases),'cases;',len(comparisons),'field comparisons',flush=True)


if __name__ == '__main__':
    main()
