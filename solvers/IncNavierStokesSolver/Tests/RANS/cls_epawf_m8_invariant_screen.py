#!/usr/bin/env python3
"""Research-only interval-witness screen of the two approved C(q) tails.

The original response is retained for q<=33/10. Candidate A continues C with
the C1 power tail and fixes beta=beta0. Candidate B additionally uses the
specified state-dependent beta tail. Complex-step searches differentiate C
and beta as parts of the law. They only nominate candidates; a gate fails
only after a rational state/direction has a directed-interval counterexample.
No production model is adopted, and absence of a witness is inconclusive.
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
from scipy.optimize import differential_evolution, minimize
import cls_epawf_m8_gamma_feasibility as common

F, QC = Fraction, Fraction(33, 10)
G0, SEED = common.G0, 830219
RATIONAL_DENOMINATORS = (1024, 4096)
RANDOM_STATES, DE_ITERATIONS, DE_POPULATION, LOCAL_ITERATIONS = 384, 90, 8, 120


def stamp():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), *arguments], text=True).strip()


def load_module(path):
    spec = importlib.util.spec_from_file_location('_m8_invariant_intervals', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def plan():
    return dict(schema=1, milestone='M8', research_screen_only=True, frozen_before_numerical_execution=True,
        production_adoption=False, candidates=['A', 'B'], beta0=['2/5', '9/20'], gamma=0,
        qc=str(QC), C_original='0.3/(1+0.35*q^(3/2))*(1-exp(-0.36*exp(0.75*q)))',
        alpha='-qc*C_original_prime(qc)/C_original(qc)',
        candidate_A='C=C_original below qc, Cc*(q/qc)^(-alpha) above; beta=beta0',
        candidate_B='same C; beta=beta0 below qc, beta0*(qc*Cc/(q*Ctail))^2 above',
        derivative_policy='differentiate the entire C(q) and beta(q); no frozen adaptive coefficient',
        join_policy='C continuation is C1; B beta is continuous with potentially distinct one-sided derivatives; certify strictly active sides only',
        normalized_gradient_center=G0.astype(int).tolist(), gradient_entry_radius='1/4', exact_trace_free=True,
        k_over_k0=['1/2', '3/2'], epsilon_over_epsilon0=['1/2', '3/2'],
        normalization='A*=A/G0; k*=Qstar*(k/K0); epsilon*=epsilon/epsilon0; stress in u_star^2, principal form in nu_t0 units',
        unchanged_components=['inner capped k repair', 'outer transport constants', 'invariant definition max(strain,rotation)', 'original C below qc'],
        search=dict(seed=SEED, random_states=RANDOM_STATES,
            propagation_grid='primitive integer directions -2..2; incompressible acoustic minimum over the perpendicular plane',
            differential_evolution=dict(population_multiplier=DE_POPULATION, max_iterations=DE_ITERATIONS, workers=1, polish=False),
            local_refinement=dict(method='Powell', max_iterations=LOCAL_ITERATIONS),
            retained_nominees_per_gate=2, rational_denominators=list(RATIONAL_DENOMINATORS),
            rational_side_preservation='at most 8 inward A20 grid steps and 8 scalar-ratio grid steps preserve the nominated smooth invariant/tail sides',
            explicit_controls=['A02=1,A20=1/4096,k/K0=5/4,epsilon/epsilon0=1/2,H=ey tensor ex',
                               'A02=5/4,k/K0=3/2,epsilon/epsilon0=1/2,v=(1,0,5)']),
        proof='directed interval full tangent or stress at exact rational trace-free admissible states and rational incompressible rank-one directions',
        counterexample_gate='principal/norm has negative upper endpoint, or covariance Rayleigh quotient has negative upper endpoint',
        absence_interpretation='inconclusive, never a qualification of this box or a production model',
        excluded=['sample scans as proof', 'production adoption', 'global model impossibility', 'state-box invariance', 'coupled PDE stability', 'terrain accuracy'])


def original_coefficient(q):
    with np.errstate(over='ignore', under='ignore'):
        return .3/(1+.35*q**1.5)*(-np.expm1(-.36*np.exp(.75*q)))


def tail_constants():
    qc = float(QC); a = .35*qc**1.5; da = .525*math.sqrt(qc)
    u = .36*math.exp(.75*qc); e = math.exp(-u)
    Cc = .3/(1+a)*(1-e)
    prime = -.3*da/(1+a)**2*(1-e)+.225*u*e/(1+a)
    return Cc, -qc*prime/Cc


def indexed_stress(gradient, k, epsilon, beta0, candidate, branch='physical'):
    """Independent dimensional indexed tensors, including adaptive beta'."""
    g = np.asarray(gradient, dtype=complex)
    s, o, eye = g+g.T, g-g.T, np.eye(3)
    ss, oo = np.sum(s*s), np.sum(o*o)
    qs, qo = k/epsilon*np.sqrt(ss/2), k/epsilon*np.sqrt(oo/2)
    q = qs if branch == 'strain' else qo if branch == 'rotation' else (qs if qs.real >= qo.real else qo)
    beta = beta0+0*q
    if q.real <= float(QC):
        C = original_coefficient(q)
    else:
        Cc, alpha = tail_constants(); logarithm = np.log(q/float(QC))
        C = Cc*np.exp(-alpha*logarithm)
        if candidate == 'B':
            beta = beta0*np.exp(-2*(1-alpha)*logarithm)
    nt = C*k*k/epsilon
    t1 = np.einsum('ik,jk->ij', s, s)-ss*eye/3
    t2 = np.einsum('ik,jk->ij', o, s)+np.einsum('jk,ik->ij', o, s)
    t3 = np.einsum('ik,jk->ij', o, o)-oo*eye/3
    t4 = np.einsum('ki,lj,kl->ij', s, o, s)+np.einsum('kj,li,kl->ij', s, o, s)
    quadratic = k*nt/epsilon*(-.1*t1+.1*t2+.26*t3)
    cubic = nt*k*k/epsilon**2*(-10*C*C*t4+5*C*C*(oo-ss)*s)
    return 2*k*eye/3-nt*s+quadratic+beta*cubic


def tangent(gradient, k, epsilon, beta0, candidate, branch='physical'):
    gradient = np.asarray(gradient, dtype=float); result = np.empty((9, 9))
    for column in range(9):
        perturbed = gradient.astype(complex); perturbed.flat[column] += 1j*1e-30
        result[:, column] = -indexed_stress(perturbed, k, epsilon, beta0, candidate, branch).imag.ravel()/1e-30
    return result


def physical_branches(a, k, epsilon):
    s, o = a+a.T, a-a.T
    qs, qo = k/epsilon*math.sqrt(float(np.sum(s*s))/2), k/epsilon*math.sqrt(float(np.sum(o*o))/2)
    branch = 'strain' if qs >= qo else 'rotation'
    q = max(qs, qo)
    return branch, 'original' if q <= float(QC) else 'tail', q


def controls():
    a = G0.copy(); a[2, 0] = 1/4096
    principal = dict(value=None, parameters=common.state_parameters(a, 1.25, .5).tolist(),
        a=[0.,1.,0.], xi=[1.,0.,0.], origin='preregistered_transverse_control')
    a = G0.copy(); a[0, 2] = 1.25
    covariance = dict(value=None, parameters=common.state_parameters(a, 1.5, .5).tolist(),
        v=[1.,0.,5.], origin='preregistered_parallel_shear_covariance_control')
    return principal, covariance


def nominate(beta, candidate, qstar):
    rng = np.random.default_rng(SEED+int(1000*beta)+(0 if candidate == 'A' else 10000))
    directions = common.primitive_directions(); states = []
    for shear in (.75, 1., 1.25):
        for kr, er in itertools.product((.5, 1., 1.5), repeat=2):
            a = G0.copy(); a[0, 2] = shear; states.append(common.state_parameters(a, kr, er))
    for _ in range(RANDOM_STATES):
        kr, er = rng.uniform(.5, 1.5, 2)
        if rng.random() < .5:
            kr, er = rng.choice((.5, 1.5), 2)
        states.append(np.r_[rng.uniform(-.25, .25, 8), kr, er])
    principal, covariance = [], []
    def retain(records, record):
        a = common.gradient(record['parameters'])
        branch, tail, invariant = physical_branches(a, qstar*record['parameters'][8], record['parameters'][9])
        record.update(invariant_branch=branch, tail_branch=tail, q=invariant)
        records.append(record); records.sort(key=lambda item:item['value']); del records[2:]
    for parameters in states:
        a = common.gradient(parameters); k, epsilon = qstar*parameters[8], parameters[9]
        linear = tangent(a, k, epsilon, beta, candidate).reshape(3, 3, 3, 3)
        for xi in directions:
            basis = common.transverse_frame(xi)
            acoustic = np.einsum('j,ijmn,n->im', xi, linear, xi)
            reduced = basis.T@(.5*(acoustic+acoustic.T))@basis
            values, vectors = np.linalg.eigh(reduced)
            retain(principal,dict(value=float(values[0]),parameters=parameters.tolist(),
                a=(basis@vectors[:,0]).tolist(),xi=xi.tolist(),origin='seeded_state_and_acoustic_direction_screen'))
        stress = indexed_stress(a,k,epsilon,beta,candidate).real
        values, vectors = np.linalg.eigh(stress)
        retain(covariance,dict(value=float(values[0]),parameters=parameters.tolist(),v=vectors[:,0].tolist(),origin='seeded_state_covariance_screen'))
    bounds = [(-.25,.25)]*8+[(.5,1.5)]*2
    def principal_objective(parameters):
        a = common.gradient(parameters); pol, xi = common.rank_one(parameters); h = np.outer(pol,xi)
        stress = indexed_stress(a.astype(complex)+1j*1e-30*h,qstar*parameters[8],parameters[9],beta,candidate)
        return float(-np.sum(h*stress.imag)/1e-30)
    def covariance_objective(parameters):
        a = common.gradient(parameters)
        return float(np.linalg.eigvalsh(indexed_stress(a,qstar*parameters[8],parameters[9],beta,candidate).real)[0])
    optimizations = []
    for name, objective, extra, seed in (
        ('principal',principal_objective,[(0.,math.pi),(0.,2*math.pi),(0.,2*math.pi)],SEED+37),
        ('covariance',covariance_objective,[],SEED+53)):
        result = differential_evolution(objective,bounds+extra,seed=seed,popsize=DE_POPULATION,
            maxiter=DE_ITERATIONS,workers=1,updating='immediate',polish=False,tol=1e-8)
        local = minimize(objective,result.x,method='Powell',bounds=bounds+extra,
                         options=dict(maxiter=LOCAL_ITERATIONS,xtol=1e-8,ftol=1e-10))
        chosen = local if local.fun < result.fun else result
        optimizations.append(dict(gate=name,de_success=bool(result.success),de_message=str(result.message),
            de_evaluations=int(result.nfev),local_success=bool(local.success),local_message=str(local.message),
            local_evaluations=int(local.nfev),selected_value=float(chosen.fun)))
        record = dict(value=float(chosen.fun),parameters=chosen.x[:10].tolist(),origin='bounded_DE_and_Powell_nomination')
        if name == 'principal':
            pol, xi = common.rank_one(chosen.x); record.update(a=pol.tolist(),xi=xi.tolist()); retain(principal,record)
        else:
            a = common.gradient(chosen.x)
            _, vectors = np.linalg.eigh(indexed_stress(a,qstar*chosen.x[8],chosen.x[9],beta,candidate).real)
            record['v'] = vectors[:,0].tolist(); retain(covariance,record)
    explicit_principal, explicit_covariance = controls()
    for record, kind in ((explicit_principal,'principal'),(explicit_covariance,'covariance')):
        a = common.gradient(record['parameters']); k, ep = qstar*record['parameters'][8],record['parameters'][9]
        branch, tail, q = physical_branches(a,k,ep); record.update(invariant_branch=branch,tail_branch=tail,q=q)
        if kind == 'principal':
            h = np.outer(record['a'],record['xi']); record['value'] = float(h.ravel()@tangent(a,k,ep,beta,candidate)@h.ravel())/float(np.sum(h*h))
        else:
            v = np.asarray(record['v']); record['value'] = float(v@indexed_stress(a,k,ep,beta,candidate).real@v)/float(v@v)
    return dict(candidate=candidate,beta0=beta,principal=principal+[explicit_principal],covariance=covariance+[explicit_covariance],
        optimizations=optimizations,interpretation='Floating nominations only; no sample is an all-state certificate')


def invariant_gap(a):
    s, o = a+a.T,a-a.T
    return sum((v*v for v in s.flat),F(0))-sum((v*v for v in o.flat),F(0))


def strict_state(record, denominator, principal, certificate):
    a, kr, er, _ = common.rational_state(record,denominator,False)
    adjustments = []
    if principal:
        sign = 1 if record['invariant_branch'] == 'strain' else -1
        for _ in range(8):
            if sign*invariant_gap(a) > 0:
                break
            value = a[2,0]+F(sign,denominator)
            if abs(value) > F(1,4):
                return None,'Cannot preserve the smooth nominated invariant side inside the box'
            a[2,0] = value; adjustments.append('A20:'+str(F(sign,denominator)))
        if sign*invariant_gap(a) <= 0:
            return None,'Nominated invariant side remains unresolved after bounded rational steps'
    Iv = certificate.Iv
    def q_interval():
        qs, qo = certificate.invariants_interval(certificate.interval_array(a),Iv(certificate.Q_LO,certificate.Q_HI)*kr,Iv(er))
        return qs if invariant_gap(a) > 0 else qo if invariant_gap(a) < 0 else Iv(max(qs.lo,qo.lo),max(qs.hi,qo.hi))
    q = q_interval(); qc = Iv(QC)
    if principal:
        above = record['tail_branch'] == 'tail'
        for _ in range(8):
            if (q.lo > qc.hi if above else q.hi < qc.lo):
                break
            epsilon_step = F(-1 if above else 1,denominator)
            k_step = -epsilon_step
            if F(1,2) <= er+epsilon_step <= F(3,2):
                er += epsilon_step; adjustments.append('epsilon_ratio:'+str(epsilon_step))
            elif F(1,2) <= kr+k_step <= F(3,2):
                kr += k_step; adjustments.append('k_ratio:'+str(k_step))
            else:
                return None,'Cannot preserve nominated tail side inside the scalar box'
            q = q_interval()
        if not (q.lo > qc.hi if above else q.hi < qc.lo):
            return None,'Tail join remains unresolved after bounded rational steps'
    return (a,kr,er,q,adjustments),None


def pairs(values):
    return {key:value.pair() if hasattr(value,'pair') else str(value) for key,value in values.items()}


def prove(record, beta, candidate, kind, denominator, certificate):
    result = dict(candidate=candidate,beta0=str(beta),kind=kind,rational_denominator=denominator,
                  nomination=record,usable=False,counterexample_certified=False)
    state, reason = strict_state(record,denominator,kind == 'principal',certificate)
    if state is None:
        result['reason'] = reason; return result
    a,kr,er,q,adjustments = state; Iv = certificate.Iv
    g,k,ep = certificate.interval_array(a),Iv(certificate.Q_LO,certificate.Q_HI)*kr,Iv(er)
    gap = invariant_gap(a); branch = 'strain' if gap > 0 else 'rotation' if gap < 0 else 'physical'
    tail = 'tail' if q.lo > Iv(QC).hi else 'original' if q.hi < Iv(QC).lo else 'physical'
    result.update(usable=True,state=dict(gradient=common.strings(a),k_over_k0=str(kr),epsilon_over_epsilon0=str(er),
        exact_trace=str(sum(a[i,i] for i in range(3))),invariant_squared_gap=str(gap)),
        q_interval=q.pair(),invariant_branch=branch,tail_branch=tail,rational_side_adjustments=adjustments)
    numeric_g = np.asarray(a,dtype=float)
    numeric_k = .5*(float(certificate.Q_LO)+float(certificate.Q_HI))*float(kr)
    if kind == 'principal':
        xi = np.array([common.rational(value,denominator) for value in record['xi']],dtype=object)
        xi_norm = sum((v*v for v in xi),F(0))
        if xi_norm == 0:
            raise AssertionError('Rational propagation vector vanished')
        numeric_eta = -np.cross(np.asarray(xi,dtype=float),np.asarray(record['a']))/float(xi_norm)
        eta = np.array([common.rational(value,denominator) for value in numeric_eta],dtype=object)
        pol = common.cross(xi,eta); h = np.outer(pol,xi)
        dot = sum((v*w for v,w in zip(pol,xi)),F(0)); norm = sum((v*v for v in h.flat),F(0))
        if dot != 0 or norm <= 0 or branch == 'physical' or tail == 'physical':
            raise AssertionError('Unresolved or inadmissible rational principal witness')
        linear,invariant = certificate.interval_tangent(g,k,ep,beta,candidate,branch=branch,tail_branch=tail)
        flat = h.ravel()
        action = sum((linear[i,j]*(flat[i]*flat[j]) for i in range(9) for j in range(9)),Iv(0))/norm
        H = np.asarray(h,dtype=float)
        numeric = float(H.ravel()@tangent(numeric_g,numeric_k,float(er),float(beta),candidate,branch)@H.ravel())/float(norm)
        result.update(a=common.strings(pol),xi=common.strings(xi),rank_one=common.strings(h),exact_a_dot_xi=str(dot),
            norm_squared=str(norm),principal_over_norm_interval=action.pair(),
            coefficients=pairs(certificate.coefficient_intervals(invariant,beta,candidate,tail_branch=tail)))
    else:
        vector = np.array([common.rational(value,denominator) for value in record['v']],dtype=object)
        norm = sum((v*v for v in vector),F(0))
        if norm <= 0:
            raise AssertionError('Rational covariance direction vanished')
        stress = certificate.interval_stress(g,k,ep,beta,candidate,branch='physical',tail_branch='physical')
        action = sum((stress[i,j]*(vector[i]*vector[j]) for i in range(3) for j in range(3)),Iv(0))/norm
        v = np.asarray(vector,dtype=float)
        numeric = float(v@indexed_stress(numeric_g,numeric_k,float(er),float(beta),candidate).real@v)/float(norm)
        result.update(v=common.strings(vector),norm_squared=str(norm),covariance_rayleigh_interval=action.pair(),
            stress_intervals=[[entry.pair() for entry in row] for row in stress])
    midpoint = .5*(float(action.lo)+float(action.hi))
    if abs(numeric-midpoint) > 1e-9*max(1.,abs(midpoint)):
        raise AssertionError('Independent fully differentiated indexed oracle disagrees with interval witness')
    result.update(independent_indexed_value=numeric,counterexample_certified=action.hi < 0)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate-source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Numerical invariant-response search and interval proofs belong on Linux')
    source,proof_source,output = Path(__file__).resolve().parents[4],args.certificate_source.resolve(),args.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    (output/'experiment-plan.json').write_text(json.dumps(plan(),indent=2)+'\n')
    report = dict(schema=1,milestone='M8',success=False,status='failed',started_utc=stamp(),failures=[],results=[],
        source_commit=git(source,'rev-parse','HEAD'),source_clean_start=not git(source,'status','--porcelain'),
        certificate_source=str(proof_source),certificate_source_commit=git(proof_source,'rev-parse','HEAD'),
        certificate_source_clean_start=not git(proof_source,'status','--porcelain'),
        command=[sys.executable,*sys.argv],script_sha256=digest(__file__),
        helper_sha256=digest(Path(common.__file__)),experiment_plan_sha256=digest(output/'experiment-plan.json'),
        allocated_workers=1,cpu_affinity=sorted(os.sched_getaffinity(0)),thread_environment=THREADS,
        versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        production_model_adopted=False,qualification='research constitutive screen only; no PDE or robust-model qualification')
    def save():
        (output/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    start = time.monotonic(); save()
    try:
        if not report['source_clean_start'] or not report['certificate_source_clean_start']:
            raise RuntimeError('Both source checkouts must be pinned and clean')
        module_dir = proof_source/'solvers/IncNavierStokesSolver/Tests/RANS'
        paths = {name:module_dir/(name+'.py') for name in ('cls_epawf_m8_invariant_certificate','cls_epawf_m8_certificate')}
        report['immutable_modules'] = {name:dict(path=str(path),sha256=digest(path)) for name,path in paths.items()}
        sys.path.insert(0,str(module_dir)); certificate = load_module(paths['cls_epawf_m8_invariant_certificate'])
        report['tail_parameters'] = pairs(certificate.tail_parameters())
        lo,hi = certificate.Iv(certificate.Q_LO),certificate.Iv(certificate.Q_HI)
        import cls_epawf_m8_certificate as base_certificate
        flo,fhi = base_certificate.cmu_interval(lo)*lo.square()-1,base_certificate.cmu_interval(hi)*hi.square()-1
        if not flo.hi < 0 < fhi.lo:
            raise RuntimeError('Original reference Q bracket is not certified')
        report['reference_root'] = dict(bracket=[certificate.Q_LO,certificate.Q_HI],lower=flo.pair(),upper=fhi.pair())
        qstar = .5*(float(certificate.Q_LO)+float(certificate.Q_HI))
        for candidate in ('A','B'):
            for beta in (F(2,5),F(9,20)):
                search = nominate(float(beta),candidate,qstar)
                (output/f'search-{candidate}-beta-{str(beta).replace("/","-")}.json').write_text(json.dumps(search,indent=2)+'\n')
                witnesses = [prove(record,beta,candidate,kind,denominator,certificate)
                    for denominator in RATIONAL_DENOMINATORS for kind in ('principal','covariance') for record in search[kind]]
                principal = [i for i,w in enumerate(witnesses) if w['kind']=='principal' and w['counterexample_certified']]
                covariance = [i for i,w in enumerate(witnesses) if w['kind']=='covariance' and w['counterexample_certified']]
                report['results'].append(dict(candidate=candidate,beta0=str(beta),witnesses=witnesses,
                    status='counterexample_certified' if principal or covariance else 'inconclusive',
                    principal_gate='failed_by_explicit_interval_witness' if principal else 'inconclusive',
                    covariance_gate='failed_by_explicit_interval_witness' if covariance else 'inconclusive',
                    negative_principal_witnesses=principal,negative_covariance_witnesses=covariance,
                    absence_of_witness_is_not_qualification=True)); save()
        report.update(success=True,status='screen_verified',conclusion='declared_box_gates_failed_by_explicit_witnesses' if
            all(r['status']=='counterexample_certified' for r in report['results']) else 'inconclusive_for_at_least_one_candidate')
    except Exception as error:
        report['failures'].append(str(error))
    try:
        report.update(source_commit_end=git(source,'rev-parse','HEAD'),source_clean_end=not git(source,'status','--porcelain'),
            certificate_source_commit_end=git(proof_source,'rev-parse','HEAD'),
            certificate_source_clean_end=not git(proof_source,'status','--porcelain'))
        if not report['source_clean_end'] or not report['certificate_source_clean_end'] or \
                report['source_commit_end']!=report['source_commit'] or report['certificate_source_commit_end']!=report['certificate_source_commit']:
            raise RuntimeError('A source checkout changed during the invariant screen')
        for value in report.get('immutable_modules',{}).values():
            if digest(value['path']) != value['sha256']:
                raise RuntimeError('An immutable proof module changed during screening')
        if digest(Path(common.__file__)) != report['helper_sha256']:
            raise RuntimeError('The rational/search helper changed during screening')
    except Exception as error:
        report.update(success=False,status='failed'); report['failures'].append(str(error))
    report.update(finished_utc=stamp(),seconds=time.monotonic()-start); save()
    print(json.dumps({key:report.get(key) for key in ('success','status','conclusion','failures','seconds')}),flush=True)
    return 0 if report['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
