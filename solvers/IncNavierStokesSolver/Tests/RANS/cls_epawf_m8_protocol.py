"""Immutable M8 candidate/observation contract, separate from legacy CV3 fits.

Only metadata are read here.  The original-first CV3 selector/freezer is never
called for a repaired CLS model.  Physical scores remain conditional on a
complete qualified native trajectory and a versioned flat-only calibration.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


VERSION = 'cls-epawf-m8-robust-v1'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def contract(source, mode, parameters):
    source = Path(source)
    data = source/'docs/plans/cao-tamura-validation-data'
    snapshot = json.loads((data/'snapshot.json').read_text())
    for name, expected in snapshot['files'].items():
        if sha(data/name) != expected:
            raise ValueError('Frozen Cao--Tamura observation changed: '+name)
    scores = json.loads((data/'scores.json').read_text())
    legacy = {}
    for name in ('cao_tamura_campaign.py', 'cao_tamura_calibration.py',
                 'cao_tamura_experiment.py', 'cao_tamura_observables.py'):
        path = source/'solvers/IncNavierStokesSolver/Tests/RANS'/name
        legacy[name] = sha(path)
    return dict(
        schema=1, version=VERSION, milestone='M8', mode=mode,
        recorded_utc=datetime.now(timezone.utc).isoformat(),
        prerequisite='User approved explicit common-cubic and reference-preserving regularized variants; original M0 contract remains frozen',
        selected_parameters=dict(parameters),
        variants=dict(
            Original='Unchanged original CLS/AWF failure control; no repair or model qualification inferred',
            CubicScaled='beta*(c4,c6,c7), original C(q), paired bounded inner destruction repair; reference-local candidate',
            Regularized='CubicScaled minus prescribed gamma(x)*(S-Sref), same source repair and adapted baseline'),
        constitutive_contract=dict(
            reference='Incoming fixed Cartesian flow/normal, analytic Gref=us/(kappa*r) flow tensor normal; no current-state or hill-normal feedback',
            gamma='GammaUniform+GammaScale*kappa*us*physical wall distance; both coefficients fixed independently of evolving fields',
            field_coupling='Same selected Reynolds stress in momentum, turbulence production and total interface traction',
            required_evidence='All active/limiting branch tangents bounded over declared state envelope; a finite scan alone cannot certify that envelope',
            local_certificate_only='Common beta certificate at neutral reference does not qualify unrestricted hill states'),
        source_contract=dict(
            theta=parameters['theta'], delta=parameters['delta'],
            limiter='C2 quintic' if parameters['smooth'] else 'piecewise linear min(1,k/(theta*K))',
            normalization='b=int(J*Dpre*limiter)/int(J*Dpre); m=max(b,delta); lambda=b/m',
            definitions='Dpos=Dpre*limiter/m; Epos=lambda*Epre; Cpos=Epos-Dpos; net k source=P_A-Dpos',
            plateau='Original source and all its derivatives retained where every ray value has k/K>theta',
            depletion='lambda<1 explicitly changes total dissipation and dependent epsilon trace; no double-counted redistribution',
            baseline='Same repaired source at chi=0 and in EP transition; Original remains separate',
            discrete_limit='Continuum quasi-positivity is not a high-order CG/IMEX positivity theorem'),
        verification_gates=dict(
            isolated_identity=1e-10, scaled_coupled_residual=1e-9,
            mpi_and_restart_physical_difference=1e-8,
            positivity='k>=0 at every declared source/mass/trace/transfer/native/dense finite evaluation set; outer k,epsilon and physical K strictly positive',
            covariance='Every finite evaluated covariance minimum >=-1e-10*K0 for the state range claimed as robust',
            source_mechanism='Depletion counterexample is recorded separately; passing its source balance cannot override a covariance failure',
            long_duration='Only complete predeclared time intervals in named flow-through and turnover units count; a solver stop leaves this gate failed',
            refinement='dt, p and mesh varied independently at fixed physical wall inputs; separate H/transition-width model studies'),
        physical_validation=dict(
            status='conditional; blocked until candidate and native long-time/refinement gates pass',
            backend='Current revised backend is native replicated dense CG; old CV strip trajectories do not qualify this wall representation',
            geometry=dict(hill_height=.04, hill_half_length=.1, roof=.8,
                          x=[-1., 1.2], matching_height=.011,
                          physical_k_sampler=.0055,
                          sensor_distances=[.0165, .022],
                          source='Cao--Tamura frozen geometry; sampler/sensor values are predeclared M8 wall-model choices'),
            routes=['spatial EP/Upstream', 'all-baseline same repaired operator'],
            paired_required='Floor/hill pairs with identical physical inlet, surface, matching/sampling parameters, variant and gamma prescription',
            first_surfaces=['S', 'R'],
            calibration_roles=scores['data_roles']['calibration'],
            reserved_flat=scores['data_roles']['reserved_flat_crosscheck'],
            reserved_hill='All hill observations and dimensional holdouts; never select beta,gamma,H,width or inlet parameters by their scores',
            development='At least three unforced domain flow-throughs and two successive one-flow-through extensions; steady drift or preregistered stationary sampling required',
            steady_maximum_drift=dict(velocity_over_Uinf=2e-4,
                                     shear_over_ustar2=2e-3,
                                     relative_k=1e-3, relative_epsilon=1e-3,
                                     friction_and_delta99=.002),
            physical_targets=scores['primary'],
            missing_observations='Report pressure and wall stress as paired CFD observables where frozen dataset has no direct experimental counterpart; no invented data or pooled replacement scores',
            coverage='Preserve unsupported heights/weights and data uncertainties; incomplete ranges remain conditional'),
        frozen_data=dict(version=snapshot['version'], snapshot_sha256=sha(data/'snapshot.json'),
                         source_pdf_sha256=snapshot['source_sha256'], files=snapshot['files']),
        legacy_tool_hashes=legacy,
        new_selection_policy='A passed M8 qualification receipt plus explicit candidate model identity is required; old original-first select/freeze commands are not invoked',
        failure_policy='Preserve every native failure report, accepted prefix, command and source hash. No rejected/unfinished trajectory enters developed averages or physical scores.',
        resources=dict(maximum_linux_workers=72, numerical_threads_per_rank=1,
                       default_ranks=1, sequential_cases=True))
