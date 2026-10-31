#!/usr/bin/env python3
"""Authoritative Linux acceptance of standard k-epsilon with analytical wall strips.

At most four jobs, each with at most four single-threaded MPI ranks. Cached
solver data require identical fixtures and binaries; all checks are reapplied.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import threading
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
from generate_m6 import generate_abl, original_j0
from generate_m4 import info, parameter
from generate_m5 import reference


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    if not path.exists():
        return []
    return [{k: v if k == 'field' else float(v) for k, v in r.items()}
            for r in csv.DictReader(path.open())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--suite', choices=('checks', 'study', 'all'), default='all')
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    require(sys.platform.startswith('linux'), 'Use the authoritative Linux environment')
    require(1 <= args.jobs <= 4, 'Use 1 to 4 jobs; at most 16 MPI workers')
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[-16:])
    source = Path(__file__).resolve().parents[4]
    build, output = args.build.resolve(), args.output.resolve()
    require(not output.is_relative_to(source), 'Keep generated results outside the source tree')
    output.mkdir(parents=True, exist_ok=True)
    solver, audit = build/'IncNavierStokesSolver', build/'RANSFieldAudit'
    binaries = {p.name: digest(p) for p in (solver, audit, build/'RANSKernelProbe')}
    commands, cases, comparisons, oracle_errors = [], {}, {}, {}
    lock = threading.Lock()
    verification = False

    def execute(command, directory, label, rejection=None, timeout=7200):
        command = list(map(str, command))
        with lock:
            commands.append(dict(command=command, cwd=str(directory), timeout=timeout))
        log = directory/(label+'.log')
        with log.open('w') as stream:
            result = subprocess.run(command, cwd=directory, env=os.environ.copy(),
                                    stdout=stream, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        raw = log.read_text()
        require((result.returncode != 0 and rejection in raw) if rejection else result.returncode == 0,
                'Inspect '+str(log))
        return raw

    def audit_json(command, directory, label):
        raw = execute(command, directory, label)
        return json.loads(raw[raw.index('{"time"'):])

    def close(actual, expected, label, tolerance=3e-9):
        error = abs(actual-expected)/max(1., abs(actual), abs(expected))
        with lock:
            oracle_errors[label] = max(oracle_errors.get(label, 0.), error)
        require(error < tolerance, 'Independent '+label+' mismatch')

    def check_points(record, root):
        params = {n.text.split('=', 1)[0].strip(): float(n.text.split('=', 1)[1])
                  for n in root.findall('CONDITIONS/PARAMETERS/P')
                  if n.text.split('=', 1)[0].strip() in (
                      'RANSStripHeight', 'RANSStripPointP', 'RANSStripPointN', 'RANSRoughnessHeight',
                      'Kinvis', 'RANSCmu', 'RANSCepsilon1', 'RANSCepsilon2', 'RANSSigmaEpsilon',
                      'RANSAWFOn', 'RANSAWFOff', 'RANSAWFFixedWeight')}
        get = lambda key, default: params.get(key, default)
        q = 1/math.sqrt(get('RANSCmu', .09))
        kappa = math.sqrt((get('RANSCepsilon2', 1.92)-get('RANSCepsilon1', 1.44))*get('RANSSigmaEpsilon', 1.3)/q)
        h, p, n = get('RANSStripHeight', 2), get('RANSStripPointP', 1), get('RANSStripPointN', 4)
        nu, rough = get('Kinvis', 1.5e-5), get('RANSRoughnessHeight', .1)
        settings = {n.attrib['PROPERTY']: n.attrib['VALUE'] for n in root.findall('CONDITIONS/SOLVERINFO/I')}
        points = record['points']
        # Include initial/final phases and independent tangential locations.
        selected = points[:12]+points[-12:]
        @lru_cache(maxsize=None)
        def local(k):
            return reference(k, h, p, rough, nu, q, kappa)
        for row in selected:
            ref = local(row['kP'])
            force = [row['forcing_x'], row['forcing_y']]
            traction = [row[key]/ref['j0']-ref['j1']/ref['j0']*f for key, f in zip(('uH', 'vH'), force)]
            production = (sum(t*t for t in traction)*ref['production_weight']
                          +2*sum(t*f for t, f in zip(traction, force))*ref['production_weight1']
                          +sum(f*f for f in force)*ref['production_weight2'])/h
            close(row['traction_x'], traction[0], 'traction_x')
            close(row['traction_y'], traction[1], 'traction_y')
            close(row['production'], production, 'production')
            if settings['RANSWallModel'] == 'AWFEquilibrium':
                for key, value in [('reference_A', ref['numerator']), ('reference_u', ref['matching_velocity']),
                                   ('reference_P0', ref['p0']), ('lambda', ref['lam']), ('cutoff', ref['cutoff'])]:
                    close(row[key], value, key)
            chi = row['chi']
            close(row['dissipation'], (1-chi)*ref['epsilon_mean']+chi*ref['p0'], 'mean_dissipation')
            close(row['epsP'], (1-chi)*ref['epsilon_p']+chi*ref['numerator']/p, 'epsilon_P')
            close(row['epsH'], (1-chi)*ref['epsilon_h']+chi*ref['numerator']/h, 'epsilon_H')
            if settings['RANSWallModel'] == 'AWFEquilibrium' and settings['RANSAWFWeight'] == 'Switched':
                us = math.hypot(ref['ustar'], nu/h)
                sensors = [math.hypot((row['uH']-ref['matching_velocity'])/us, row['vH']/us),
                           math.hypot((row['uN']-row['uP']-ref['ustar']/kappa*math.log(n/p))/us,
                                      (row['vN']-row['vP'])/us),
                           (row['kN']-row['kP'])/row['kP'],
                           (row['epsN']-ref['numerator']/n)/(ref['numerator']/n),
                           h*math.hypot(row['pressure_x'], row['pressure_y'])/us**2,
                           h*math.hypot(row['convection_x'], row['convection_y'])/us**2, row['wP']/us]
                for key, value in zip(('r_U', 'r_DeltaU', 'r_K', 'r_Epsilon', 'r_Pressure', 'r_Convection', 'r_Normal'), sensors):
                    close(row[key], value, key)
        record['independently_checked_points'] = len(selected)

    def run(name, model='standard_switched', dimension=2, ranks=1, modify=None, rejection=None,
            preservation=False, point_audit=False, prepared=None, **options):
        directory = output/name
        directory.mkdir(exist_ok=True)
        session = directory/(name+'.xml')
        if prepared:
            spec = prepared(session)
        else:
            base = dict(model=model, dimension=dimension, nx=1, ny=1, degree=4, layers=12,
                        dt=.5, steps=4, backend='matrixfree', boundary='outflow', tolerance=1e-8, monitor=True)
            base.update(options)
            spec = generate_abl(session, **base)
            tree = ET.parse(session)
            root = tree.getroot()
            parameter(root, 'IterativeSolverTolerance', 1e-15 if verification else 1e-13)
            parameter(root, 'RANSAWFPointDiagnostics', int(point_audit))
            if modify:
                modify(root)
            tree.write(session, encoding='utf-8', xml_declaration=True)
        root = ET.parse(session).getroot()
        # Modifiers may change the closure, wall anchor or duration. Report
        # the actual fixture values rather than the generator's initial ones.
        keys = ('TimeStep', 'NumSteps', 'Qstar', 'Kappa', 'Ustar', 'MatchingVelocity',
                'Uref', 'L', 'H', 'RANSRoughnessHeight')
        values = {n.text.split('=', 1)[0].strip(): float(n.text.split('=', 1)[1])
                  for n in root.findall('CONDITIONS/PARAMETERS/P')
                  if n.text.split('=', 1)[0].strip() in keys}
        spec = dict(spec, q=values['Qstar'], kappa=values['Kappa'], dt=values['TimeStep'],
                    steps=int(values['NumSteps']), friction_velocity=values['Ustar'],
                    matching_velocity=values['MatchingVelocity'], uref=values['Uref'],
                    roughness=values['RANSRoughnessHeight'],
                    physical_time=values['TimeStep']*values['NumSteps'],
                    flow_times=values['TimeStep']*values['NumSteps']*values['Uref']/values['L'],
                    log_z0=math.log(values['H'])-values['Kappa']*values['MatchingVelocity']/values['Ustar'])
        fingerprint = dict(fixture=digest(session), binaries=binaries, ranks=ranks, rejection=rejection)
        cache = directory/'case.json'
        cached = json.loads(cache.read_text()) if cache.exists() else {}
        cached_fingerprint = {k:v for k,v in cached.get('fingerprint', {}).items() if k != 'validator'}
        if cached and cached_fingerprint == fingerprint:
            record = cached
            require(all((directory/p).is_file() and digest(directory/p) == sha
                        for p,sha in record['artifacts'].items()), 'Cached artifact changed: '+name)
            print('Rechecking cached solver data '+name, flush=True)
        else:
            print('Running '+name, flush=True)
            command = (['mpirun', '--bind-to', 'none', '--map-by', 'slot', '-np', str(ranks)] if ranks > 1 else [])
            execute(command+[solver, session.name], directory, 'solver', rejection)
            record = dict(fingerprint=fingerprint, spec=spec, ranks=ranks,
                          session=str(session), field=str(directory/(name+'.fld')))
            if rejection:
                record['expected_rejection'] = rejection
            else:
                record['audit'] = audit_json([audit, session.name, name+'.fld', name+'.fld', 'profiles.csv'], directory, 'audit')
                record['stages'] = rows(directory/(name+'.rans.csv'))
                record['wall'] = rows(directory/(name+'.awf.csv'))
                record['errors'] = rows(directory/(name+'.rans-errors.csv'))
                record['points'] = [r for path in sorted(directory.glob(name+'.awf-points-r*.csv')) for r in rows(path)]
            record['artifacts'] = {str(p.relative_to(directory)): digest(p) for p in directory.rglob('*')
                                   if p.is_file() and p.suffix in ('.xml', '.log', '.csv', '.fld', '.chk')}
            cache.write_text(json.dumps(record, indent=2)+'\n')
        if not rejection:
            stages = record['stages']
            expected = next(int(float(n.text.split('=', 1)[1])) for n in root.findall('CONDITIONS/PARAMETERS/P')
                            if n.text.split('=', 1)[0].strip() == 'NumSteps')
            require(len(stages) == expected, 'Incomplete timesteps: '+name)
            for row in stages+record['wall']+record['points']:
                require(all(math.isfinite(v) for v in row.values()), 'Nonfinite diagnostic: '+name)
                require(row.get('guard_activations', 0) == 0, 'Guard activated: '+name)
            require(all(r['min_k'] > 0 and r['min_eps'] > 0 for r in stages), 'Nonpositive outer state: '+name)
            require(all(r['min_kP'] > 0 and r['min_epsP'] > 0 for r in record['wall']), 'Nonpositive wall state: '+name)
            if point_audit:
                require(record['points'], 'Missing point diagnostics')
                check_points(record, root)
            if record['errors']:
                fields = {r['field'] for r in record['errors']}
                last = max(r['time'] for r in record['errors'])
                record['peak_l2'] = {f: max(r['normalized_l2'] for r in record['errors'] if r['field'] == f) for f in fields}
                record['peak_linf'] = {f: max(r['pointwise_relative_linf'] for r in record['errors'] if r['field'] == f) for f in fields}
                record['final_errors'] = {r['field']: r for r in record['errors'] if r['time'] == last}
                required = ['u', 'k', 'eps', 'rans_nut', 'rans_R_xx', 'rans_R_yy', 'rans_R_zz', 'rans_R_xz']
                record['engineering_target_met'] = all(record['peak_l2'][f] < .001 and record['peak_linf'][f] < .001 for f in required)
                if preservation:
                    require(spec['flow_times'] >= 5 and record['engineering_target_met'], 'Five-flow preservation failed: '+name)
                    require(all(r['chi_min'] == 1 and r['chi_max'] == 1 and
                                r['cutoff_fallback_points'] == r['regime_fallback_points'] == r['range_fallback_points'] == 0
                                for r in record['wall']), 'Reference left correction plateau: '+name)
            print('Completed '+name, flush=True)
        record['spec'] = spec
        record['revalidated_with'] = digest(Path(__file__))
        with lock:
            cases[name] = record
            (output/('progress-'+args.suite+'.json')).write_text(json.dumps(cases, indent=2)+'\n')
        return record

    def compare(name, first, second, tolerance=1e-7):
        value = audit_json([audit, first['session'], first['field'], second['field']], output, name)
        maximum = max(f['difference_linf'] for f in value['fields'].values())
        require(maximum < tolerance, name+' field difference too large')
        comparisons[name] = dict(maximum_absolute_linf=maximum, tolerance=tolerance, fields=value['fields'])

    def configured(root):
        q = 1/math.sqrt(.075)
        kappa = math.sqrt(.5*1.17/q)
        matching = .16*original_j0(q*.16, 2, .1, 1.5e-5)
        for key, value in dict(RANSCmu=.075, RANSCepsilon1=1.4, RANSCepsilon2=1.9,
                               RANSSigmaK=.85, RANSSigmaEpsilon=1.17, Qstar=q, Kappa=kappa,
                               MatchingVelocity=matching, Uref=matching+.4/kappa*math.log(5)).items():
            parameter(root, key, value)
        for node in root.findall('CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF="2"]/N[@VAR="eps"]'):
            node.attrib['VALUE'] = '-Ustar^4/(1.17*Z)'

    kernel = {}
    if args.suite in ('checks', 'all'):
        verification = True
        raw = execute(['ctest', '--test-dir', build, '--parallel', '1', '--output-on-failure'], output, 'ctest', timeout=900)
        require('100% tests passed' in raw and 'out of 9' in raw, 'Expected nine CTests')
        raw = execute([build/'RANSKernelProbe', '--output', output/'kernels.json'], output, 'kernels', timeout=120)
        kernel = json.loads((output/'kernels.json').read_text())
        import validate_m2 as independent
        for case in kernel['standard_walls']:
            mu, c1, c2, _, se, _, _ = case['model_parameters']
            q = 1/math.sqrt(mu)
            kappa = math.sqrt((c2-c1)*se/q)
            independent.near(case['q'], q, 'standard_q')
            independent.near(case['kappa'], kappa, 'standard_kappa')
            independent.check_wall(case, q, kappa)
        kernel = dict(native=kernel['summary'], independent_standard_wall_cases=len(kernel['standard_walls']),
                      independent_errors=independent.errors)
        fixed = run('fixed_2d', model='standard_fixed', point_audit=True)
        switched = run('switched_2d', point_audit=True)
        compare('fixed_switched', fixed, switched)
        run('switched_3d', dimension=3, ranks=4, point_audit=True)
        run('configured', model='standard_fixed', dt=.025, steps=4, modify=configured, point_audit=True)
        original = run('original', model='standard_original', dt=.001, steps=2, point_audit=True)
        def zero(root):
            info(root, 'RANSAWFWeight', 'Fixed'); parameter(root, 'RANSAWFFixedWeight', 0)
        disabled = run('zero_weight', dt=.001, steps=2, modify=zero, point_audit=True)
        compare('original_recovery', original, disabled)
        def depart(root):
            node = root.find('CONDITIONS/FUNCTION[@NAME="InitialConditions"]/E[@VAR="eps"]')
            node.attrib['VALUE'] = '1.5*('+node.attrib['VALUE']+')'
            for suffix in ('U', 'DeltaU', 'K', 'Pressure', 'Convection', 'Normal'):
                parameter(root, 'RANSAWFWeight'+suffix, 0)
        # A reference inlet pins epsilon_N to equilibrium at its wall corner.
        # A periodic column makes this a genuinely uniform switch-off control.
        off = run('switched_off', boundary='periodic', dt=.025, steps=2, modify=depart, point_audit=True)
        off_original = run('original_departure', model='standard_original', boundary='periodic',
                           dt=.025, steps=2, modify=depart, point_audit=True)
        compare('switch_off_recovery', off, off_original)
        require(all(p['chi'] == 0 for p in off['points']), 'Off-equilibrium switch stayed active')
        for label, key, value, message in [
                ('shifted', 'RANSHeightShift', .02, 'without RANSHeightShift'),
                ('independent_z0', 'RANSZ0', .02, 'do not prescribe RANSZ0'),
                ('cutoff', 'RANSStripPointP', .001, 'cutoff/high-Re preflight')]:
            run('reject_'+label, steps=1, modify=lambda r,k=key,v=value: parameter(r,k,v), rejection=message)
        def fallback(root):
            info(root, 'RANSAWFFailurePolicy', 'UseOriginal')
            info(root, 'RANSAWFWeight', 'Fixed')
            # Reduce the reference Kappa until the capped family is inadmissible,
            # keeping P/H/N fixed so this tests fallback, not a stiff geometry.
            parameter(root, 'RANSCepsilon2', 1.46)
        fb = run('cutoff_fallback', dt=.025, steps=2, modify=fallback, point_audit=True)
        fb_original = run('original_inadmissible_reference', model='standard_original', dt=.025, steps=2,
                          modify=lambda r: parameter(r, 'RANSCepsilon2', 1.46), point_audit=True)
        compare('cutoff_original_recovery', fb, fb_original)
        require(all(p['chi'] == 0 for p in fb['points']), 'Inadmissible cutoff did not fall back')
        whole = run('restart_whole', steps=8)
        prefix = run('restart_prefix', steps=4)
        helper = source/'solvers/IncNavierStokesSolver/Examples/NeutralABL/restart.py'
        def restart(session):
            execute([sys.executable, helper, prefix['session'], prefix['field'], session, '--steps', 4], session.parent, 'prepare')
            return prefix['spec']
        resumed = run('restart_1_to_4', ranks=4, prepared=restart)
        compare('restart_complete_fields', whole, resumed)
        whole3 = run('restart_whole_3d', dimension=3, ranks=4, steps=4)
        prefix3 = run('restart_prefix_3d', dimension=3, ranks=4, steps=2)
        def restart3(session):
            execute([sys.executable, helper, prefix3['session'], prefix3['field'], session, '--steps', 2], session.parent, 'prepare')
            return prefix3['spec']
        resumed3 = run('restart_4_to_2', dimension=3, ranks=2, prepared=restart3)
        compare('restart_complete_fields_3d', whole3, resumed3, tolerance=1e-6)
        def periodic(root):
            parameter(root, 'RANSStageAudit', 1)
            # The large-volume constant-test budgets need a tighter state
            # fixed point than the 0.1% profile-preservation experiment.
            parameter(root, 'RANSNonlinearTolerance', 1e-12)
            parameter(root, 'RANSPressureTolerance', 1e-8)
        column = run('periodic_budget', boundary='periodic', dt=.1, steps=4, modify=periodic, point_audit=True)
        require(all(r['budget_available'] == 1 and abs(r['k_budget']) < 1e-6 and abs(r['momentum_budget_x']) < 1e-5
                    for r in column['wall']), 'Conservative wall/outer budget failed')
        # Original numerical paths remain covered by the existing independent suites.
        execute([sys.executable, Path(__file__).with_name('validate_m2.py'), '--build', build,
                 '--output', output/'legacy-kernels'], output, 'legacy-kernels', timeout=900)
        execute([sys.executable, Path(__file__).with_name('validate_m5.py'), '--build', build,
                 '--output', output/'legacy-cls', '--pilot'], output, 'legacy-cls', timeout=1800)

    if args.suite in ('study', 'all'):
        verification = False
        jobs = []
        for dim in (2, 3):
            for mode in ('fixed', 'switched'):
                jobs.append((f'five_d{dim}_{mode}', dict(model='standard_'+mode, dimension=dim,
                            ranks=1 if dim == 2 else 4, flow_times=5, preservation=True)))
        for label, updates in [('h6', dict(layers=6)), ('h24', dict(layers=24)),
                               ('p3', dict(degree=3)), ('p5', dict(degree=5)), ('dt025', dict(dt=.25))]:
            jobs.append((label, dict(flow_times=5, **updates)))
        jobs.append(('dt025_3d', dict(dimension=3, ranks=4, dt=.25, flow_times=5, preservation=True)))
        jobs.append(('perturbed', dict(perturbation=.01, flow_times=5)))
        def richards_norris(root):
            # Match the SAME U(10)=10, z0=.02 standard reference used by the
            # original log-wall control by solving the AWF momentum anchor.
            from scipy.optimize import brentq
            kappa = math.sqrt(.48*1.3*.3)
            star = kappa*10/math.log(10/.02)
            matching = star/kappa*math.log(2/.02)
            roughness = brentq(lambda h: star**2*original_j0((10/3)*star**2, 2, h, 1.5e-5)-matching,
                              0, 1.99, xtol=1e-13)
            for key,value in dict(Ustar=star, MatchingVelocity=matching, Uref=10,
                                  RANSRoughnessHeight=roughness).items():
                parameter(root,key,value)
        jobs.append(('richards_norris', dict(steps=300, modify=richards_norris, preservation=True)))
        def launch(job):
            name, options = job
            try:
                return run(name, **options)
            except Exception as error:
                with lock:
                    cases[name] = dict(status='failed', message=str(error), options=options)
                    (output/('progress-'+args.suite+'.json')).write_text(json.dumps(cases, indent=2)+'\n')
                print('FAILED '+name+': '+str(error), flush=True)
                return cases[name]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            list(pool.map(launch, jobs))
        require(all(c.get('status') != 'failed' for c in cases.values()), 'Study contains failed cases; inspect progress JSON')
        rates = {}
        for field in ('u', 'k', 'eps', 'rans_nut'):
            herrors = [cases[n]['final_errors'][field]['normalized_l2'] for n in ('h6', 'five_d2_switched', 'h24')]
            perrors = [cases[n]['final_errors'][field]['normalized_l2'] for n in ('p3', 'five_d2_switched', 'p5')]
            if field != 'k':
                require(herrors[0] > herrors[1] > herrors[2] and perrors[0] > perrors[1] > perrors[2],
                        'Spatial errors did not decrease: '+field)
            rates[field] = dict(h_errors=herrors, h_orders=[math.log(a/b, 2) for a,b in zip(herrors,herrors[1:])], p_errors=perrors)
        perturbation = cases['perturbed']
        require(perturbation['final_errors']['k']['normalized_l2'] < .1*perturbation['peak_l2']['k'], 'Perturbation did not decay')
    else:
        rates = {}
    # Compact committed evidence; all step/point data and fields remain in run directories.
    compact = {name: {k: v for k, v in c.items() if k not in ('errors', 'stages', 'wall', 'points')}
               | dict(maximum_weak_divergence=max((r['weak_divergence_l2'] for r in c.get('stages', [])), default=0))
               for name, c in cases.items()}
    paths = [p for p in (source/'library/SolverUtils/Turbulence').glob('*.h') if p.name != 'HITGenerator.h']
    paths += [source/p for p in ('library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp',
        'library/SolverUtils/Turbulence/RANSModels.cpp', 'library/UnitTests/SolverUtils/RANSKernelProbe.cpp',
        'solvers/IncNavierStokesSolver/EquationSystems/RANSAnalyticalWallStrip.cpp',
        'solvers/IncNavierStokesSolver/EquationSystems/RANSAnalyticalWallStrip.h',
        'solvers/IncNavierStokesSolver/EquationSystems/RANSVelocityCorrectionScheme.cpp',
        'solvers/IncNavierStokesSolver/Tests/RANS/generate_m6.py')]
    report = dict(status='passed', suite=args.suite, utc=datetime.now(timezone.utc).isoformat(),
        source_parent=subprocess.check_output(['git', '-C', source, 'rev-parse', 'HEAD'], text=True).strip(),
        source_sha256={str(p.relative_to(source)): digest(p) for p in paths+[Path(__file__)]},
        binaries=binaries, platform=platform.platform(), threads=THREADS, maximum_workers=args.jobs*4,
        commands=commands, kernel=kernel, cases=compact, comparisons=comparisons,
        oracle_maximum_scaled_errors=oracle_errors, spatial_refinement=rates)
    (output/('validation-'+args.suite+'.json')).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='passed', cases=len(cases), comparisons=comparisons, spatial_refinement=rates), indent=2), flush=True)


if __name__ == '__main__':
    main()
