#!/usr/bin/env python3
"""Finite rational interval witnesses for prescribed-gamma incompatibility.

Numerical searches only nominate candidates. An infeasibility conclusion uses
two necessary conditions, proved with outward-rounded intervals at rational
admissible states: a rank-one incompressible diffusion lower requirement and a
modeled-covariance Rayleigh upper requirement on the same fixed GammaScale.
Failure to obtain disjoint bounds is inconclusive. This does not qualify any
smaller envelope, adaptive coefficient, changed C(q), PDE, or terrain model.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                  'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.optimize import differential_evolution, minimize

F = Fraction
RADIUS, RATIOS, SEED = F(1, 4), (F(1, 2), F(3, 2)), 823147
RATIONAL_DENOMINATORS = (1024, 4096)
RANDOM_STATES, DE_ITERATIONS, DE_POPULATION, LOCAL_ITERATIONS = 384, 90, 8, 120
SEARCH_RAYLEIGH_COEFFICIENT = 1/64
G0 = np.zeros((3, 3)); G0[0, 2] = 1
SREF = G0+G0.T


def stamp():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), *arguments], text=True).strip()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def plan(target):
    return dict(schema=1, milestone='M8', frozen_before_numerical_execution=True,
        purpose='necessary fixed-GammaScale compatibility on one declared broad constitutive box',
        beta=['2/5', '9/20'], normalized_gradient_center=G0.astype(int).tolist(),
        normalized_gradient_entry_radius=str(RADIUS), trace_free=True,
        k_over_k0=list(map(str, RATIOS)), epsilon_over_epsilon0=list(map(str, RATIOS)),
        normalization='A*=A/G0; k*=Qstar*(k/K0); epsilon*=epsilon/epsilon0; nu_t0=u_star^2/G0=1 in these units',
        gamma='fixed GammaScale=gamma/nu_t0; Sref*=ex tensor ez+ez tensor ex, independent of the state',
        principal_target_over_nut0=str(target),
        necessary_conditions=['GammaScale >= target-D_beta(H)/(|a|^2*|xi|^2), H=a tensor xi, a dot xi=0',
                              'GammaScale <= v^T R_beta v / v^T(S-Sref)v when the denominator is positive'],
        candidate_search=dict(seed=SEED, random_states=RANDOM_STATES,
            propagation_grid='primitive integer vectors with components -2..2, sign duplicates removed',
            covariance_vectors='positive generalized-Rayleigh direction when base stress is SPD; finite integer controls otherwise',
            minimum_search_Rayleigh_gamma_coefficient='1/64 times |v|^2, a search conditioning choice only; proof checks the exact positive coefficient',
            differential_evolution=dict(workers=1, population_multiplier=DE_POPULATION, max_iterations=DE_ITERATIONS, polish=False),
            local_refinement=dict(method='Powell', max_iterations=LOCAL_ITERATIONS),
            rational_denominators=list(RATIONAL_DENOMINATORS), retained_candidates_per_condition=2,
            interpretation='search only; sampled extrema are not mathematical certificates'),
        proof='outward Decimal interval stress and complete active-branch tangent at finite rational witnesses; exact rational rank-one and trace tests',
        success_meaning='screen observations and proofs reproduced; an inconclusive screen is not a feasibility certificate',
        infeasibility_scope='only this declared state box, fixed nonnegative scalar GammaScale and prescribed Sref, separately for each fixed beta',
        excluded=['scan-as-proof', 'global model impossibility', 'invariance of the state box', 'adaptive gamma', 'changed C(q)', 'coupled PDE or terrain qualification'])


def gradient(parameters):
    a = G0.copy()
    a[0, 0], a[1, 1] = parameters[:2]
    total = a[0, 0]+a[1, 1]
    if total > .25:
        a[1, 1] = .25-a[0, 0]
    elif total < -.25:
        a[1, 1] = -.25-a[0, 0]
    a[2, 2] = -a[0, 0]-a[1, 1]
    for index, offset in zip(((0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)), parameters[2:8]):
        a[index] += offset
    return a


def state_parameters(a, kr, er):
    return np.array([a[0, 0], a[1, 1], a[0, 1], a[0, 2]-1,
                     a[1, 0], a[1, 2], a[2, 0], a[2, 1], kr, er], dtype=float)


def unit_vector(theta, phi):
    return np.array([math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)])


def transverse_frame(xi):
    axis = np.eye(3)[int(np.argmin(np.abs(xi)))]
    first = np.cross(xi, axis); first /= np.linalg.norm(first)
    return np.column_stack((first, np.cross(xi, first)))


def rank_one(parameters):
    xi = unit_vector(parameters[10], parameters[11]); basis = transverse_frame(xi)
    a = basis@np.array([math.cos(parameters[12]), math.sin(parameters[12])])
    return a, xi


def primitive_directions():
    output = {}
    for values in itertools.product(range(-2, 3), repeat=3):
        common = math.gcd(math.gcd(abs(values[0]), abs(values[1])), abs(values[2]))
        if not common:
            continue
        key = tuple(v//common for v in values)
        if next(v for v in key if v) < 0:
            key = tuple(-v for v in key)
        output[key] = np.asarray(key, dtype=float)/np.linalg.norm(key)
    return list(output.values())


def covariance_direction(stress, strain_change, controls):
    values, vectors = np.linalg.eigh(stress)
    candidates = list(controls)+[vectors[:, 0]]
    if values[0] > 1e-12:
        eig, vec = eigh(strain_change, stress, check_finite=True)
        candidates.append(vec[:, -1])
    result = None
    for vector in candidates:
        denominator = float(vector@strain_change@vector)
        if denominator > SEARCH_RAYLEIGH_COEFFICIENT*float(vector@vector):
            value = float(vector@stress@vector)/denominator
            if result is None or value < result[0]:
                result = value, np.asarray(vector)/np.linalg.norm(vector)
    return result


def search_candidates(beta, q, oracle):
    rng = np.random.default_rng(SEED+int(1000*beta))
    controls = primitive_directions(); candidates = []
    for shear in (.75, 1., 1.25):
        for kr, er in itertools.product((.5, 1., 1.5), repeat=2):
            a = G0.copy(); a[0, 2] = shear
            candidates.append(state_parameters(a, kr, er))
    for _ in range(RANDOM_STATES):
        values = rng.uniform(-.25, .25, 8)
        kr, er = rng.uniform(.5, 1.5, 2)
        if rng.random() < .5:
            kr, er = rng.choice((.5, 1.5), 2)
        candidates.append(np.r_[values, kr, er])
    principal, covariance = [], []
    def retain(records, record):
        records.append(record); records.sort(key=lambda item: item['value']); del records[2:]
    for parameters in candidates:
        a = gradient(parameters); k, ep = q*parameters[8], parameters[9]
        tensor = oracle.stress_tangent(a, k, ep, beta).reshape(3, 3, 3, 3)
        for xi in controls:
            basis = transverse_frame(xi)
            acoustic = np.einsum('j,ijmn,n->im', xi, tensor, xi)
            reduced = basis.T@(.5*(acoustic+acoustic.T))@basis
            values, vectors = np.linalg.eigh(reduced)
            retain(principal, dict(value=float(values[0]), parameters=parameters.tolist(),
                                   a=(basis@vectors[:, 0]).tolist(), xi=xi.tolist(), origin='seeded_state_and_integer_direction_screen'))
        stress = oracle.indexed_stress(a, k, ep, beta).real
        upper = covariance_direction(stress, a+a.T-SREF, controls)
        if upper is not None:
            retain(covariance, dict(value=upper[0], parameters=parameters.tolist(), v=upper[1].tolist(), origin='seeded_state_and_Rayleigh_screen'))
    bounds = [(-.25, .25)]*8+[(.5, 1.5)]*2
    def principal_objective(parameters):
        a = gradient(parameters); pol, xi = rank_one(parameters); h = np.outer(pol, xi)
        stress = oracle.indexed_stress(a.astype(complex)+1j*1e-30*h, q*parameters[8], parameters[9], beta)
        return float(-np.sum(h*stress.imag)/1e-30)
    def covariance_objective(parameters):
        a = gradient(parameters); vector = unit_vector(parameters[10], parameters[11])
        denominator = float(vector@(a+a.T-SREF)@vector)
        if denominator <= SEARCH_RAYLEIGH_COEFFICIENT:
            return 1e6+1e3*abs(denominator)
        stress = oracle.indexed_stress(a, q*parameters[8], parameters[9], beta).real
        return float(vector@stress@vector)/denominator
    optimization = []
    for name, objective, extra, seed in (
            ('principal', principal_objective, [(0., math.pi), (0., 2*math.pi), (0., 2*math.pi)], SEED+17),
            ('covariance', covariance_objective, [(0., math.pi), (0., 2*math.pi)], SEED+29)):
        result = differential_evolution(objective, bounds+extra, seed=seed, popsize=DE_POPULATION,
            maxiter=DE_ITERATIONS, workers=1, updating='immediate', polish=False, tol=1e-8)
        local = minimize(objective, result.x, method='Powell', bounds=bounds+extra,
                         options=dict(maxiter=LOCAL_ITERATIONS, xtol=1e-8, ftol=1e-10))
        chosen = local if local.fun < result.fun else result
        optimization.append(dict(condition=name, differential_evolution_success=bool(result.success),
            differential_evolution_message=str(result.message), de_evaluations=int(result.nfev),
            local_success=bool(local.success), local_message=str(local.message), local_evaluations=int(local.nfev),
            selected_value=float(chosen.fun), interpretation='floating candidate only'))
        record = dict(value=float(chosen.fun), parameters=chosen.x[:10].tolist(), origin='bounded_DE_and_Powell_candidate_search')
        if name == 'principal':
            pol, xi = rank_one(chosen.x); record.update(a=pol.tolist(), xi=xi.tolist()); retain(principal, record)
        else:
            record['v'] = unit_vector(chosen.x[10], chosen.x[11]).tolist(); retain(covariance, record)
    return dict(beta=beta, principal_candidates=principal, covariance_candidates=covariance,
                random_states=RANDOM_STATES, deterministic_states=27, optimization=optimization,
                interpretation='No floating extremum in this object is a mathematical bound')


def rational(value, denominator):
    return F(round(float(value)*denominator), denominator)


def rational_state(record, denominator, require_smooth=False):
    parameters = [rational(value, denominator) for value in record['parameters']]
    a = np.full((3, 3), F(0), dtype=object); a[0, 2] = F(1)
    a[0, 0], a[1, 1] = parameters[:2]
    total = a[0, 0]+a[1, 1]
    if total > RADIUS:
        a[1, 1] = RADIUS-a[0, 0]
    elif total < -RADIUS:
        a[1, 1] = -RADIUS-a[0, 0]
    a[2, 2] = -a[0, 0]-a[1, 1]
    for index, value in zip(((0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)), parameters[2:8]):
        a[index] += value
    def gap(matrix):
        s, o = matrix+matrix.T, matrix-matrix.T
        return sum((v*v for v in s.flat), F(0))-sum((v*v for v in o.flat), F(0))
    if require_smooth and gap(a) == 0:
        # A02>=3/4, so changing A20 by a nonzero rational step changes the
        # invariant gap. Choose an inward step at an A20 box boundary.
        a[2, 0] += F(-1 if a[2, 0] > 0 else 1, denominator)
    if sum(a[i, i] for i in range(3)) != 0:
        raise AssertionError('Rational state is not exactly trace-free')
    for index in np.ndindex(a.shape):
        center = F(1) if index == (0, 2) else F(0)
        if abs(a[index]-center) > RADIUS:
            raise AssertionError('Rational witness left the preregistered gradient box')
    if not all(RATIOS[0] <= ratio <= RATIOS[1] for ratio in parameters[8:10]):
        raise AssertionError('Rational witness left the preregistered scalar box')
    return a, parameters[8], parameters[9], gap(a)


def cross(a, b):
    return np.array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]], dtype=object)


def strings(array):
    array = np.asarray(array, dtype=object)
    return [str(x) for x in array] if array.ndim == 1 else [[str(x) for x in row] for row in array]


def interval_witness(record, beta, denominator, kind, certificate, oracle, target):
    Iv = certificate.Iv
    a, kr, er, gap = rational_state(record, denominator, kind == 'principal')
    q = Iv(certificate.Q_LO, certificate.Q_HI)
    state = dict(gradient=strings(a), k_over_k0=str(kr), epsilon_over_epsilon0=str(er),
                 exact_trace=str(sum(a[i, i] for i in range(3))), invariant_squared_gap=str(gap))
    g, k, ep = certificate.interval_array(a), q*kr, Iv(er)
    result = dict(kind=kind, beta=str(beta), rational_denominator=denominator, state=state,
                  nominated_by=record['origin'], floating_nomination=record['value'])
    numeric_g = np.asarray(a, dtype=float); numeric_k = float((Decimal(certificate.Q_LO)+Decimal(certificate.Q_HI))/2)*float(kr)
    if kind == 'principal':
        xi = np.array([rational(value, denominator) for value in record['xi']], dtype=object)
        xi_square = sum((v*v for v in xi), F(0))
        if xi_square == 0:
            raise AssertionError('Rational propagation vector vanished')
        numeric_eta = -np.cross(np.asarray(xi, dtype=float), np.asarray(record['a']))/float(xi_square)
        eta = np.array([rational(value, denominator) for value in numeric_eta], dtype=object)
        pol = cross(xi, eta); h = np.outer(pol, xi)
        dot = sum((x*y for x, y in zip(pol, xi)), F(0))
        norm = sum((entry*entry for entry in h.flat), F(0))
        shift = sum((x*y for x, y in zip(h.flat, (h+h.T).flat)), F(0))
        if dot != 0 or norm <= 0 or shift != norm or gap == 0:
            raise AssertionError('Rational principal witness is not a nonzero incompressible rank-one direction on a smooth branch')
        branch = 'strain' if gap > 0 else 'rotation'
        tensor, invariant = certificate.interval_tangent(g, k, ep, Iv(beta), branch)
        linear = h.ravel(); action = Iv(0)
        for i in range(9):
            for j in range(9):
                action += tensor[i, j]*(linear[i]*linear[j])
        normalized = action/norm; necessary = Iv(target)-normalized
        H = np.asarray(h, dtype=float)
        numeric = float(H.ravel()@oracle.stress_tangent(numeric_g, numeric_k, float(er), float(beta), branch=branch)@H.ravel())/float(norm)
        middle = .5*(float(normalized.lo)+float(normalized.hi))
        if abs(numeric-middle) > 1e-9*max(1., abs(middle)):
            raise AssertionError('Independent indexed directional tangent disagrees with the interval witness')
        result.update(a=strings(pol), xi=strings(xi), exact_a_dot_xi=str(dot), rank_one=strings(h),
            norm_squared=str(norm), exact_gamma_shift=str(shift), physical_branch=branch,
            branch_invariant_interval=invariant.pair(), principal_form_interval=action.pair(),
            principal_over_norm_interval=normalized.pair(), necessary_gamma_threshold_interval=necessary.pair(),
            certified_necessary_lower=str(necessary.lo), adverse_principal_certified=normalized.hi < 0,
            independent_indexed_value=numeric)
    else:
        vector = np.array([rational(value, denominator) for value in record['v']], dtype=object)
        if not any(vector):
            raise AssertionError('Rational covariance vector vanished')
        reference = np.full((3, 3), F(0), dtype=object); reference[0, 2] = reference[2, 0] = F(1)
        change = a+a.T-reference
        denominator_exact = sum((vector[i]*change[i, j]*vector[j] for i in range(3) for j in range(3)), F(0))
        result.update(v=strings(vector), exact_positive_gamma_coefficient=str(denominator_exact))
        if denominator_exact <= 0:
            result.update(usable_upper_bound=False, reason='Rationalized Rayleigh gamma coefficient is not positive')
            return result
        stress = certificate.interval_stress(g, k, ep, Iv(beta), branch='physical')
        rayleigh = sum((stress[i, j]*(vector[i]*vector[j]) for i in range(3) for j in range(3)), Iv(0))
        upper = rayleigh/denominator_exact
        v = np.asarray(vector, dtype=float)
        numeric = float(v@oracle.indexed_stress(numeric_g, numeric_k, float(er), float(beta)).real@v)/float(denominator_exact)
        middle = .5*(float(upper.lo)+float(upper.hi))
        if abs(numeric-middle) > 1e-9*max(1., abs(middle)):
            raise AssertionError('Independent indexed covariance disagrees with the interval witness')
        result.update(usable_upper_bound=True, covariance_rayleigh_interval=rayleigh.pair(),
            covariance_stress_intervals=[[entry.pair() for entry in row] for row in stress],
            necessary_gamma_threshold_interval=upper.pair(), certified_necessary_upper=str(upper.hi),
            independent_indexed_value=numeric)
    return result


def prove_candidates(search, beta, target, certificate, oracle):
    witnesses = []
    for denominator in RATIONAL_DENOMINATORS:
        for kind, key in (('principal', 'principal_candidates'), ('covariance', 'covariance_candidates')):
            for record in search[key]:
                witnesses.append(interval_witness(record, beta, denominator, kind, certificate, oracle, target))
    lowers = [v for v in witnesses if v['kind'] == 'principal' and v['adverse_principal_certified']]
    uppers = [v for v in witnesses if v['kind'] == 'covariance' and v.get('usable_upper_bound')]
    result = dict(beta=str(beta), witnesses=witnesses, status='inconclusive', fixed_gamma_nonnegative=True,
        principal_target_over_nut0=str(target), numerical_scan_is_not_proof=True,
        interpretation='Failure to find incompatible witnesses is not a feasibility or all-state certificate')
    if not lowers or not uppers:
        result['reason'] = 'No certified adverse principal witness or usable covariance upper witness'
        return result
    low = max(lowers, key=lambda v: Decimal(v['certified_necessary_lower']))
    high = min(uppers, key=lambda v: Decimal(v['certified_necessary_upper']))
    lower, upper = Decimal(low['certified_necessary_lower']), Decimal(high['certified_necessary_upper'])
    effective_lower = max(Decimal(0), lower)
    gap = certificate.Iv(effective_lower)-certificate.Iv(upper)
    result.update(selected_principal_witness=witnesses.index(low), selected_covariance_witness=witnesses.index(high),
        necessary_gamma_lower=str(effective_lower), necessary_gamma_upper=str(upper),
        certified_gap_lower=str(gap.lo), empty_intersection=effective_lower > upper)
    if effective_lower > upper:
        result.update(status='infeasible_on_declared_box',
            interpretation='The two interval-certified necessary GammaScale constraints have empty intersection; no fixed nonnegative scalar GammaScale satisfies both properties at every state in this box')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate-source', type=Path, required=True, help='Immutable source checkout containing the qualified public interval module')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--target', default='0', help='Rational principal lower target in nu_t0 units (default0 tests strict-positive compatibility)')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific feasibility searches and interval evaluation belong on Linux')
    target = F(args.target)
    if target < 0:
        raise ValueError('Principal target must be nonnegative')
    source, proof_source, output = Path(__file__).resolve().parents[4], args.certificate_source.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output/'experiment-plan.json').write_text(json.dumps(plan(target), indent=2)+'\n')
    report = dict(schema=1, milestone='M8', status='failed', success=False, started_utc=stamp(),
        experiment_plan_sha256=digest(output/'experiment-plan.json'), script_sha256=digest(__file__),
        source_commit=git(source, 'rev-parse', 'HEAD'), source_clean_start=not git(source, 'status', '--porcelain'),
        certificate_source=str(proof_source), certificate_source_commit=git(proof_source, 'rev-parse', 'HEAD'),
        certificate_source_clean_start=not git(proof_source, 'status', '--porcelain'),
        allocated_workers=1, cpu_affinity=sorted(os.sched_getaffinity(0)), thread_environment=THREADS,
        versions=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        command=[sys.executable, *sys.argv], results=[], failures=[], production_or_solver_execution=False,
        scope='fixed scalar gamma feasibility on the preregistered constitutive box only')
    def save():
        (output/'validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    start = time.monotonic(); save()
    try:
        if not report['source_clean_start'] or not report['certificate_source_clean_start']:
            raise RuntimeError('Source and immutable certificate checkout must both start clean')
        module_dir = proof_source/'solvers/IncNavierStokesSolver/Tests/RANS'
        paths = {name: module_dir/(name+'.py') for name in ('cls_epawf_m8_certificate', 'cls_epawf_m8_oracle')}
        report['immutable_modules'] = {name: dict(path=str(path), sha256=digest(path)) for name, path in paths.items()}
        certificate = load_module(paths['cls_epawf_m8_certificate'], '_m8_gamma_certificate')
        oracle = load_module(paths['cls_epawf_m8_oracle'], '_m8_gamma_indexed')
        lo, hi = certificate.Iv(certificate.Q_LO), certificate.Iv(certificate.Q_HI)
        root_lo = certificate.cmu_interval(lo)*lo.square()-1
        root_hi = certificate.cmu_interval(hi)*hi.square()-1
        if not root_lo.hi < 0 < root_hi.lo:
            raise RuntimeError('Reference Q interval does not bracket its defining root')
        report['reference_root'] = dict(bracket=[certificate.Q_LO,certificate.Q_HI],
            Cq2_minus1_at_lower=root_lo.pair(), Cq2_minus1_at_upper=root_hi.pair(), precision=certificate.PREC)
        q = .5*(float(certificate.Q_LO)+float(certificate.Q_HI))
        for beta in (F(2, 5), F(9, 20)):
            search = search_candidates(float(beta), q, oracle)
            (output/('search-beta-'+str(beta).replace('/', '-')+'.json')).write_text(json.dumps(search, indent=2)+'\n')
            result = prove_candidates(search, beta, target, certificate, oracle)
            report['results'].append(result); save()
        report.update(success=True, status='screen_verified', conclusion='infeasible_on_declared_box' if
            all(r['status'] == 'infeasible_on_declared_box' for r in report['results']) else 'inconclusive_for_at_least_one_variant',
            qualification='no fixed-gamma model or trajectory is qualified by this screen')
    except Exception as error:
        report['failures'].append(str(error))
    try:
        report.update(source_commit_end=git(source, 'rev-parse', 'HEAD'), source_clean_end=not git(source, 'status', '--porcelain'),
            certificate_source_commit_end=git(proof_source, 'rev-parse', 'HEAD'),
            certificate_source_clean_end=not git(proof_source, 'status', '--porcelain'))
        if not report['source_clean_end'] or not report['certificate_source_clean_end'] or \
                report['source_commit_end'] != report['source_commit'] or report['certificate_source_commit_end'] != report['certificate_source_commit']:
            raise RuntimeError('A pinned source checkout changed during feasibility qualification')
        for value in report.get('immutable_modules', {}).values():
            if digest(value['path']) != value['sha256']:
                raise RuntimeError('An immutable interval/indexed module changed during qualification')
    except Exception as error:
        report['success'] = False; report['status'] = 'failed'; report['failures'].append(str(error))
    report.update(finished_utc=stamp(), seconds=time.monotonic()-start); save()
    print(json.dumps({key:report.get(key) for key in ('success','status','conclusion','failures','seconds')}), flush=True)
    return 0 if report['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
