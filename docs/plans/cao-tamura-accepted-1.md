# CT-Accepted-1: benchmark feasibility and conservative open floor

This package starts clean at `b127d307edb423f93a93aa9d861ce4b6a486573c`,
matching the fetched reviewed remote head. It is a separate Accepted-only
capability and feasibility package. The directional package's failed A93,
response and budget decisions are unchanged and are not admission criteria for
this new boundary problem.

The pre-run [contract](cao-tamura-accepted-1-evidence/contract.json) pins both
specifications, the unchanged Cao–Tamura v1.1 dataset and the source paper hash.
It declares Reynolds-matched scaling, four wall heights, three K factors,
common and once-fitted model-conditioned roughness scenarios, independent flat
and metric-weighted column checks, the new conservative open policy and bounded
native verification/pilot gates. No hill data enter a fit. Neither k nor epsilon
is uniquely measured. Curved columns do not establish native hill capability.

Native development and validation use a new isolated Linux source/build/runtime
and run area under `/media/jay/新加卷/nektar-v2026-runs/ct-accepted-1-20261111`.
The project allocation is at most 16 workers including helpers, at most one
native family at a time, four compiler jobs, numerical-library threads one,
16-worker startup headroom and a fresh aggregate census below 72. Existing CLS
and WJ installations/campaigns are preserved.

The conservative boundary equations pass their independent native checks.
The intended A95 rough-floor pilot fails its actual constitutive-domain gate
at attempted nondimensional time **0.004**, before the planned **0.2** endpoint.
A93 and E94/E95 complete their two-step smokes. Experimental prediction
validation is unattempted. The aggregate resource audit also records a failed
72-worker gate, discussed below. These outcomes remain separate in the
[required-scope decision](cao-tamura-accepted-1-evidence/decision.json).

## Benchmark and feasibility contract

Lengths are scaled by the paper's hill height, 0.04 m; velocity is scaled by
the incident velocity at that height. Thus `Hh=U(Hh)=1` and
`nu=1/Re_H`. The central frozen incident profiles give:

| Surface | Re_H | Uinf/U(Hh) | u_tau/U(Hh) |
|---|---:|---:|---:|
| Smooth | 12000 | 1.260222109 | 0.041587330 |
| Rough | 10500 | 1.456234077 | 0.074267938 |

No dimensional velocity or air viscosity is treated as measured. The exact
compact cosine-squared hill retains half-length 0.10 m, incident depth about
0.25 m and eventual roof 0.8 m. Paper vertical height and AWF normal distance
are separate quantities. The physical ground stays at zero. Cube height
5 mm, reported aerodynamic z0 (0.004/0.2 mm), assumed rough displacement
3.5 mm and AWF equivalent input h_r retain their distinct meanings and source
statuses. SS/RR are the future primary cases; mixed surfaces remain deferred.

The frozen v1.1 auditor passes: 98 series, 1151 points and 39 scalar facts.
Observations, extraction bounds, weights and calibration/holdout roles are
unchanged. Only visible flat Fig. 3a data enter this package's fit. The rough
16 mm anchor is interpolated in log(y-d); the rough 8 mm value is extrapolated
below the first visible point at 10.2409 mm. Its U/Uinf estimate is 0.3979593,
with extraction envelope 0.0105288 and two-versus-three-point method difference
0.0043913. Displacement and experimental uncertainties remain unquantified.
The 8 mm extrapolation is displayed, not used as a measured fitting anchor.

[ct_accepted_feasibility.py](../../solvers/IncNavierStokesSolver/Tests/RANS/ct_accepted_feasibility.py)
implements the appendix independently of the production kernel. The pre-run
contract fixes H_A=2/4/8/16 mm, K/K_ref=0.8/1/1.2, common h_r=1/5 mm controls,
and two separately labelled model-conditioned fits at 16 mm. One fits visible
velocity; the other fits the reported z0 as an **unshifted** intercept scenario.
Each fit is held fixed in every subsequent H/K sensitivity. Model-conditioned
h_r is an effective input, not a new measurement of the surface.

All q and Natural kappa values come from defining roots. Wall constants are
`c_mu=q^-2`, `c_l=kappa*q^(3/2)`, `alpha=kappa/sqrt(q)`. Neither inversion nor
cutoff uses conventional kappa=0.41. The selected velocity-fit references give:

| Model | q | kappa | h_r, mm | Lambda at 16 mm | cutoff/H_A | unshifted z0, mm |
|---|---:|---:|---:|---:|---:|---:|
| A93 | 3.401265417 | 0.428323629 | 3.377941 | 3.844739 | 0.0581494 | 0.175766 |
| A95 | 3.822000426 | 0.404060995 | 4.235842 | 3.712025 | 0.0664022 | 0.226942 |
| E94E95 | 3.324506456 | 0.412087472 | 4.032219 | 3.758775 | 0.0633693 | 0.208546 |

The compatible matching velocity is 0.7822045101 and total traction
0.00551572660 in these three references. Minimum reference k/K is 1 at the
endpoint. Inner covariance minima/K are 0.374835/0.407070/0.368180; outer
reference minima/K are 0.381443/0.289908/0.365870. These are flat-reference
statements, not global realizability claims.

The [complete flat table](cao-tamura-accepted-1-evidence/flat-feasibility.csv)
has 180 candidates: 102 admissible references, 42 geometric exclusions and
36 negative-k references. It includes J0, J2, finite-nu production, cutoff
margin, plateau, traction, intercept, Q0 endpoint, k minimum/location,
covariance and mechanical-work defects. Quadrature is split at the physical
breakpoints and viscosity-scale subpanels; orders 48/96 assess uncertainty.
All k stationary points are bracketed through Q0 roots on monotone intervals.
Native flat-kernel comparisons have maximum scaled defect 1.11e-14.
Molecular-only cutoff failures and h_r>=H rejections have their intended
native reasons. Production/dissipation and finite-nu work identities pass.

The required smooth counterexample is reproduced at H_A=8 mm, h_r=0 and
K=q*u_tau^2 using source-derived u_tau/nu:

| Model | Lambda | ground k/K |
|---|---:|---:|
| A93 | 2.744266120 | -0.323602373 |
| A95 | 2.692138199 | -0.057156647 |
| E94E95 | 2.706834352 | -0.388142326 |

Cutoff existence does not make these references admissible. This does not
exclude every forced, finite-depth or nonequilibrium smooth solution. No sink,
cutoff or k value is clipped, and no EP fallback or redistribution is added.

The [936 curved columns](cao-tamura-accepted-1-evidence/curved-columns.csv)
use 39 stations per model/surface/height: crest, inflections, intermediate
points, both foot limits and outside-flat controls. Exact-foot entries select
the outside branch; inside limits retain the curvature jump. They solve the
monotone J-weighted cutoff equation with counterfactual p0(K), then use
`Q_J=integral(J*(E-p0))` and `k_col=K-integral(Q_J/(J*D))`. Zero curvature
recovers the independent flat cutoff and k profile. No flat covariance formula
is used to certify curved flow.

Local stretch, sufficient global reach and thin-layer adequacy are distinct.
The sufficient normal-map reach is `1/sup|h''|=1.266515 Hh`; all declared
heights satisfy it and the local stretch bound. Whole-hill thin adequacy
requires max(|curvature|*H_A)<=0.1. With the once-fitted rough inputs, A93 at
4 mm passes all sampled column/geometry checks. The A95/E94E95 4 mm layers
fail h_r<H_A; their 8/16 mm layers fail the whole-hill thin criterion.
Positive frozen columns elsewhere do not constitute terrain equilibria or
curved covariance/stability certificates. Native mapped evolution stays guarded.

## Conservative open-floor implementation

`RANSWallBoundaryPolicy=AcceptedOpenFlatV1` is opt-in, Accepted-only, canonical
affine serial flat 2D, with uniform isotropic Modified_A spaces and the retained
3000-unknown limit. `RANSOpenTop=DirichletCompatibleV1` and
`RANSOpenBackflow=RejectV1` are explicit. Evaluation caching is Off; selecting
ExactStepV1 is rejected. Existing closed behavior and signatures are preserved.

Inlet velocity/k/outer epsilon and compatible top values are prescribed;
ground velocity is zero with natural k diffusion. The outlet has free
velocity/scalar rows, prescribed **full total traction** and outward scalar
diffusive data. With `un=U.n`, the added weak loads are

```
momentum: + integral_outlet v . (U*un - t)
scalar:   + integral_outlet phi * (f*un - q_f)
```

Thus zero q_f retains advective transport. These terms occur in the original
residual, frozen one-sweep matrices/residuals, resolution components and scalar
reconciliation. Volume/interface ownership remains GeneratedVolumeOwnerV1.
Participating epsilon rows use the outlet face's native generating geometry.

Absolute total traction fixes the pressure level. All physical pressure test
functions are free. The legacy mean equation/multiplier is inactive; its unused
zero slot is excluded from the solve. No outlet pressure is imposed.
Changing traction by `-c*n` shifts Pi by c. The pressure correction still uses
the consistent free-velocity mass inverse, including the inherited complete
pressure action. Flow -> epsilon -> k ordering and frozen prescribed sink are
unchanged. Dependent epsilon takes precedence at its compatible inlet/H
corner; incompatible data are rejected. There is no inner epsilon equation or
extra two-sided epsilon-flux condition.

Detected outlet inflow below `-1e-12*U(Hh)` is rejected with time, point, side
and geometry. Signed flux is retained in assembly. The check is a declared
Gaussian face inventory, not a continuum sign certificate. Internal reversal
and negative wall shear are not outlet-policy failures.

The boundary identity binds actual inlet/top, traction, diffusive and source
functions and parameters into native restart and replay identity. Restart also
binds timestep, method and quadrature. Read-only preservation identity retains
its deliberate quadrature independence. Explicit/implicit live boundary maps,
actual/effective initialization and native histories are recorded.

Source mapping:

| Responsibility | Source |
|---|---|
| Policy/data identity, outlet loads/status/exports | RANSWallEvolutionOpenBoundary.cpp |
| Constraints, pressure level, original residual, guards | RANSWallEvolutionPatch.cpp |
| Frozen matrices and actual mixed-system export | RANSWallEvolutionPicard.cpp |
| Stage outflow checks, unchanged sweep order | RANSWallEvolutionPimple.cpp |
| Full pressure action, physical term reconstruction | RANSWallEvolutionActions.cpp / Resolution.cpp |
| Initialization/restart and physics identity | RANSEquilibriumWallEvolution.cpp / RANSWallPhysicsSignature.h |
| Read-only sampling and physical scalar reconstruction | RANSPreservationProbe.cpp |
| Raw divergence/isotropic-work witnesses | RANSWallEvolutionFeasible.cpp |

## Native verification

[Boundary receipts](cao-tamura-accepted-1-evidence/boundary_verification.json)
compare independent tensor/source weak loads and matrices with production
assembly, not residual subtraction. The 188-slot two-column p3/p1 fixture has
187 physical unknowns and a 72-row free mixed system. Mixed rank is 72;
pressure action, divergence and actual consistent-mass PIMPLE operator each
have rank 12. No pressure-level null mode remains (relative cutoff 1e-11).
This is finite-case qualification, not an inf-sup theorem.

| Check | Maximum observed defect | Limit |
|---|---:|---:|
| Independent raw weak load, scaled | 4.90e-16 | 2e-10 |
| Frozen matrix entries | 1.16e-14 | 2e-10 |
| Global Stokes / volume | 1.73e-15 / 3.56e-15 | 2e-9 |
| Physical term sum versus original residual | 1.89e-15 | 2e-10 |
| Independent manufactured free weak equations | 2.44e-15 | 2e-9 |
| Evolved c=0.125 pressure offset, coefficients | 4.86e-15 | 2e-10 |

The manufactured test prescribes its correct, nonzero full traction and uses
declared verification forcing. It is absent from the physical pilot. Its
independent actual storage/source/physical-gradient-flux defects are
8.65e-14 (inner k), 1.84e-14 (outer k), 4.35e-15 (outer epsilon).
Prescribed-row reactions are not substituted for face fluxes.

[History controls](cao-tamura-accepted-1-evidence/history_controls.json) pass
bitwise same-method restart and complete rejected-step rollback; changed
traction histories are rejected. Ordinary closed CLS states/histories are
bitwise equal under old/new runtimes and cache on/off. Existing native
diagnostic negative controls pass. Outlet backflow, unsupported cache use and
incompatible inlet/dependent-epsilon corners fail for the intended reasons.
The registered `CTAcceptedOpenConservativeNative` CTest passes in 54.69 s.

## Unforced rough-floor pilot and first failure

The declared A95 pilot uses two x columns `[0,.125,.25]`, H_A=.4, roof=20,
incident depth=6.25, p3/p1 and **1085 assembled slots (1084 physical)**.
It has 492 free velocity and 48 physical pressure rows. Physical nu is
1/10500. Initial/inlet mean velocity is a C1 Hermite finite-depth inference
matching the wall profile, U(Hh)=1 and Uinf at the incident depth. Inner k is
the nonconstant independent wall reference's Hermite representation and is
transported thereafter. Outer k tapers to `1.5*(.003*Uinf)^2`; epsilon uses
model-resolved length-scale samples with a finite-depth Hermite cap. Neither
is claimed uniquely measured. Continuous cubic-extremum audits and native
inventories confirm positive initial k/epsilon, including between samples.

Initial Pi is derived from complete normal stresses. Initial analytical total
outlet traction is held fixed, with zero scalar diffusive data and active
advection. No body/scalar forcing is present. Eta and corrections are zero;
there is no reference reset or epsilon-trace repair. Independent zero-step
native histories bind each dt; same-space initial coefficients are bitwise
identical. Wall/other orders are 32/8/10/10/10. All input and consumed-boundary
identities are in [pilot.json](cao-tamura-accepted-1-evidence/pilot.json).
The [consumed-input/stage supplement](cao-tamura-accepted-1-evidence/source-stage-supplement.json)
asserts equality of zero-step, replay and consumed initial coefficients, both
live boundary maps and outlet loads. It retains every small native stage and
accepted-state audit, including continuity, source/storage and refreshed-trace
defects. Closed old/new signature strings are also bitwise equal.

| Model / dt | Requested T | Last accepted T | Result |
|---|---:|---:|---|
| A95 / .002 | .2 | .002, 1 step | A95 radicand failure at attempted .004 |
| A95 / .001 | .2 | .003, 3 steps | A95 radicand failure at attempted .004 |
| A93 / .001 | .002 | .002, 2 steps | Smoke passes |
| E94E95 / .001 | .002 | .002, 2 steps | Smoke passes |

Both A95 failures follow momentum prediction and pressure correction, before
epsilon transport. Rejected states, accepted histories, stage/domain receipts
and [independent witnesses](cao-tamura-accepted-1-evidence/a95-failure-witnesses.json)
are retained. On unchanged actual outer-volume rows, the minima are
**-0.04632764 / -0.14882060**, at `(x,z)=(.248369158,6.298925259)`, geometry
25, strictly inside the outlet-adjacent cell above the incident-layer depth.
There is no endpoint roundoff or owner mismatch. For dt=.002, k=2.86359e-5,
epsilon=8.20472e-9 and raw divergence=-.00133918. Full A95 production/stress is
undefined there; no clipped tensor is reported. The formal isotropic-work
component is recorded separately.

The original unforced pilot remains failed. A fresh two-attempt timed
reproducer has bitwise identical accepted state/history and rejected
coefficients. It uses no old evolved checkpoint. The intended T=.2 would
exceed the Uinf bulk transit .1717; accepted prefixes reach only about
1.2%/1.7% of that transit. A full transit or developed precursor is not achieved.

## Budgets, uncertainty and limitations

The [postmortem](cao-tamura-accepted-1-evidence/postmortem.json) compares both
A95 trajectories at their common accepted T=.002. Scaled physical weighted
RMS dt differences are U=2.38e-7, w=5.74e-8, k=4.92e-7,
epsilon=2.26e-7 and Pi=2.46e-6 (limit 1e-3). Pressure retains one absolute
traction-determined level, with **no mean subtraction**. This is a common-prefix
uncertainty check, not a passing planned endpoint or measured temporal order.

Independent physical budgets use consecutive accepted states and actual
frozen flow/epsilon/pre-k stages, retaining dependent boundary storage and
current/frozen coefficient splitting. Whole-k external totals exclude both
internal-interface faces; the two one-sided interface fluxes remain visible.
The outer epsilon budget includes its own H boundary. With outward total flux
`F_out=advective_outward-diffusive_outward`, the defect is
`storage-source+F_out`. The native packet's net boundary term has the opposite
sign. No scalar*div(U) correction is introduced.

The physical-rate normalization is `scalar_scale * patch_volume * U(Hh)/Hh`,
with K0 or epsilon(H_A) as scalar scale. The frozen-coefficient controls use
wall orders 48/64, other orders 12/16 and physical-face orders 16/32, on the
same coefficients and boundaries. Reconstruction and finite-defect accuracy
are separately classified. At the common prefix or smoke endpoints:

| Configuration | Scaled mass defect | Pointwise div RMS / max | Largest scaled scalar defect across controls |
|---|---:|---:|---:|
| A95 dt=.002, T=.002 | 3.66e-16 | 1.34e-4 / 8.25e-3 | 8.96e-5 |
| A95 dt=.001, T=.002 | 1.39e-17 | 1.50e-4 / 9.13e-3 | 8.88e-5 |
| A93 T=.002 | 4.22e-16 | 3.70e-4 / 1.25e-2 | 8.37e-5 |
| E94E95 T=.002 | 3.52e-16 | 5.63e-5 / 1.01e-2 | 9.79e-5 |

Mass limit is 1e-8; scalar target is .02. These finite defects and quadrature
sensitivity pass the declared targets at the stated states. They do not prove
spatial convergence. Small weak continuity does not imply small pointwise
divergence. Native ranges and signed isotropic production-work witnesses are
retained, with independent full 3x3 R/Q/production and model-specific epsilon
source comparisons. Outer nu/nu_t and sqrt(nu*epsilon)/k indicators are retained
for a later finite-Re assessment. Outer diffusion remains high-Re, and the
E-model destruction denominator remains k.

No coupled A95 endpoint at .2 exists. No longer response, A93 perturbation,
EP/transition, common-kappa, mixed-surface, native curved/3D/MPI, open hill,
developed precursor or experimental validation is claimed. Old directional
failures and all eight historical unresolved epsilon-order observations are
unchanged. Preparation failures (coordinate declaration, diagnostic corner,
manufactured traction, expression syntax, endpoint extension and directory
hashing) remain separately archived; none is relabelled a physical success.

## Provenance, resource incident and reproduction

The governing/diagnostic executable sources were built at
`a909093955f2c29fbeab500a4d0693f2f008de70`; the
freeze records clean source `8f90f99e7b90e931aacf81ab1f96dd8ade9a729a` (the
additional change is test-only). Later commits add diagnostics/scripts and
evidence, with no subsequent governing C++ change. GNU 11.4 Release uses
`-O3 -DNDEBUG -std=gnu++17`, MPI enabled, runtime rank one. The isolated
ldd closure is frozen; actual dynamically loaded MPI plugins are separately
hashed. Full source/library/configuration/history identities are in
[source-scope](cao-tamura-accepted-1-evidence/source-scope.json),
[runtime identity](cao-tamura-accepted-1-evidence/runtime-identity.json),
[loaded maps](cao-tamura-accepted-1-evidence/actual-loaded-library-maps.json)
and [raw artifact hashes](cao-tamura-accepted-1-evidence/raw-artifact-hashes.json).

```
IncNavierStokesSolver SHA256
216cce2edae987dd137b10393e7ab3aa3abf8b7daeb39cf40f505d8e1667dd39
RANSPreservationProbe SHA256
1b04fcf83f09170e061cdff698677554163936dd6e56bd4105bb805120b8b464
```

The minimal A95 reproducer takes 9.82 s wall / 8.76 s process CPU and peaks
at 197148 KiB RSS. Exclusive instrumented assembly totals .134 CPU s;
dense factorization/triangular solves total .021 CPU s; domain computation
2.986 CPU s and object initialization 3.416 CPU s dominate. Inclusive and
exclusive timings are not added. Smokes cost about 11.26 process CPU s each;
their peak RSS is about 198-201 thousand KiB. All heavy fields, histories,
matrices and logs remain under the Linux experiment root.

**Aggregate resource compliance failed.** Booked launches used at most one
native family, four compiler workers, explicit single-thread libraries and
16-worker startup headroom; their recorded external counts were <=29.
However the monitor recorded **79 active threads**, above the hard 72 limit.
An ad hoc read-only NumPy inspection with an unset thread environment is the
likely cause. The old monitor did not retain that peak's PID/time, so exact
attribution cannot be established. No pause events were recorded. The failed
high-water receipt is retained in
[resource-control-incident.json](cao-tamura-accepted-1-evidence/resource-control-incident.json).
Subsequent numerical helpers explicitly set all library threads to one; the
corrected NumPy helper observes exactly one thread. The owned command wrapper
also propagates a deliberately failing exit code 7 unchanged, with a fresh
observed census of 18 threads. The guard now logs peak
PID/thread/ownership inventories even without pause events; its self-test
passes. The historical maximum is not reset. Native numerical results are
retained with this procedural failure, not presented as full resource compliance.

For exact retained-runtime reproductions on Linux (use fresh destinations):

```bash
CT_ROOT='/media/jay/新加卷/nektar-v2026-runs/ct-accepted-1-20261111'
CT_TESTS="$CT_ROOT/source/solvers/IncNavierStokesSolver/Tests/RANS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$CT_ROOT/tmp" LD_LIBRARY_PATH="$CT_ROOT/runtime-v2/lib"

# Each control destination must be new; the wrapper books the whole family,
# starts its owner-only monitor and preserves the numerical exit code.
ct_run() {
  CT_RUN_TAG="$1"
  shift
  python3 -B "$CT_TESTS/ct_accepted_command.py" \
    --control "$CT_ROOT/reproduce-control-$CT_RUN_TAG" -- "$@"
}
ct_run data python3 -B "$CT_TESTS/cao_tamura_data.py" \
  --data "$CT_ROOT/source/docs/plans/cao-tamura-validation-data" \
  --report "$CT_ROOT/reproduce-data-audit.json"
ct_run boundary python3 -B "$CT_TESTS/ct_accepted_open.py" --runtime "$CT_ROOT/runtime-v2" \
  --output "$CT_ROOT/reproduce-boundary"
ct_run feasibility python3 -B "$CT_TESTS/ct_accepted_feasibility.py" \
  --contract "$CT_ROOT/source/docs/plans/cao-tamura-accepted-1-evidence/contract.json" \
  --data "$CT_ROOT/source/docs/plans/cao-tamura-validation-data" \
  --native-probe '/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-directional-20261111/runtime/RANSNewModelsProbe' \
  --output "$CT_ROOT/reproduce-feasibility"
ct_run regression ctest --test-dir "$CT_ROOT/build" -R '^CTAcceptedOpenConservativeNative$' \
  -j 1 --output-on-failure
```

The exact gated native launch commands and measured execution receipts are
committed beside this report. To reproduce the physical failure, copy
[minimal-a95-radicand.xml](cao-tamura-accepted-1-evidence/minimal-a95-radicand.xml)
as `case.xml` into a fresh Linux run directory, then invoke the frozen solver
through `ct_run witness "$CT_ROOT/runtime-v2/IncNavierStokesSolver" case.xml`
with the environment above. It initializes a fresh history, accepts one step
and exits nonzero at the second attempt. Its only portable-export change is
the checkpoint output path. The original tested configuration hash is retained.
Python >=3.10, NumPy/SciPy and the source-compatible Nektar SDK are required;
the frozen-data auditor itself uses only the standard library. Exact isolated
build commands/flags are in the build-execution and environment receipts.

## Next bounded gate

Do not extend this A95 precursor or start hill evolution. The smallest useful
next experiment is one declared refinement of the outlet/incident-depth
supports, keeping these physical inlet functions and equations fixed: split
the last x interval at .1875 and add z=6.5 above delta. Its predicted stored
count is 1669 (must be confirmed natively), below 3000. Fresh .002/.001
histories through the first failure interval can test whether local normal
strain/divergence and pressure/normal-stress representation move the domain
loss. A low inferred free-stream epsilon and the original A95 domain remain
independent physical/model questions; no regularization is authorized.

For a later scalable-backend/curved-coupling package, use this 1085-slot case
and the 188-slot matrix fixture as exact small references. Port conservative
volume/interface/external-face operators, full pressure action and consistent
free-mass PIMPLE correction; preserve absolute traction pressure level,
dependent scalar storage and strict method/history identity. Sparse assembly
will allow useful multicolumn/terrain spaces, but cannot by itself fix the
observed constitutive failure or the dominant diagnostic cost. Before native
curvature, independently test one shared 5 mm rough-layer height with the
already fixed h_r fits, then qualify the full J-weighted Cartesian coupling
and vertical-height observations. A passing frozen column is insufficient.
