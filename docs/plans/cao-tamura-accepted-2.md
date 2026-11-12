# CT-Accepted-2: A95 domain diagnosis and rough-wall operating height

This package starts clean at `eedb3738017f46f9c71b857bf814213ad070d7f7`,
identical to the fetched GitHub SSH head. HTTPS fetch failed before that
successful check. The [pre-run contract](cao-tamura-accepted-2-evidence/contract.json)
fixes the physical functions, five bounded mesh controls, histories, diagnostics,
gates and resource allocation. CT-Accepted-1, the directional results and the
historical eight unresolved epsilon-order observations are unchanged.

The governing equations and one-sweep numerical operator are unchanged.
The C++ change preserves computed flow-stage diagnostics before an A95
constitutive exception. A separate diagnostic correction replaces the hardcoded
streamwise length in scalar-budget normalization with the actual geometry.
It gives the identical scale on the original patch.

The four primary layouts remain physically inadmissible. The single predeclared
below-delta companion passes the .008/.02 admissions and both unchanged-method
histories through **T=.2**, including the final physical-domain, mass, scalar-budget
and timestep checks. This is a bounded small-patch result with unresolved full
evolved spatial convergence. The independent wall study supports a 5 mm frozen
reference under the retained thin criterion, with physical roughness-layer
qualifications below. **Resource compliance fails**: private helper allocation
was undercounted, and a later foreign debugger caused an aggregate excursion.
Neither is erased by the corrective controls. No developed precursor or
experimental prediction is admitted.

## Source, storage and physics

Mac source/Git work uses the existing implementation worktree. Authoritative
Linux work is isolated under
`/media/jay/新加卷/nektar-v2026-runs/ct-accepted-2-20261112`.
The server rebooted after the eight primary trajectories and their comparison
collector error. The disk was temporarily unmounted, then mounted at
`/media/jay/新加卷1`. The user explicitly selected that actual mount.

The [storage bridge](cao-tamura-accepted-2-evidence/storage-relocation.json)
does not edit an old request, signature, report, coefficient or history.
Read-only access maps only that disk prefix and verifies artifact hashes.
Executable and complete linked-library identities agree after the same prefix
mapping, with all bytes unchanged. Remaining work uses a new detached
`source-after-reboot` checkout and fresh histories on the actual mount.
The old worktrees, runtime, CLS/WJ projects and shared installations are preserved.

The retained physics is Natural A95/A93/E94E95, original high-Re outer
equations, zero outer molecular diffusion, finite inner nu, Accepted budget
cutoff, transported nonconstant inner k, total epsilon, full three-dimensional
stresses, signed raw-gradient production and DeviatoricStrainV1. The inner
prescribed sink has zero local reaction. There is no inner epsilon equation,
reference restoration, forcing, clipping, regularization, bound projection or
post-k epsilon-trace repair.

AcceptedOpenFlatV1 retains its compatible prescribed top, inlet data,
full-total-traction outlet, scalar advective fluxes and diffusive data,
backflow rejection, absolute pressure level and consistent-free-mass PIMPLE.
GeneratedVolumeOwnerV1, the unsupported-route guards and 3000-slot limit remain.
The physical initial/inlet/top/outlet FUNCTION expressions are copied unchanged.
In particular, epsilon is never regenerated using a refined mesh knot list.

## Computed-stage diagnosis

`RANSWallEvolutionPimple.cpp` now publishes the already computed linear
residuals and algebraic defects before calling the complete corrected-state
residual. With the existing opt-in snapshot flag it writes:

- The actual input, history, predictor, corrected velocity/pressure and fixed reference.
- Frozen flow matrix, free-velocity mass, pressure operator, native free-row IDs and hashes.
- Predictor increment, pressure increment, subtracted mass correction and full pressure action.
- Predictor/pressure linear residuals, continuity and momentum-splitting identities.

This adds no governing solve. The diagnostic does not evaluate an inadmissible
A95 tensor to obtain those quantities. Native generating geometry, coordinates,
field rows and side metadata are exported separately at the valid input and
applied to the saved coefficients. The outer volume rows are fixed for each
mesh; no rounded-coordinate owner substitution is introduced.

The baseline dt=.002/.001 accepted states, rejected states, failure and native
histories are bitwise identical to CT1. Full stress and production are marked
undefined at negative-radicand states. Raw G, divergence, k/epsilon, strain
eigenvalues, strain/rotation norms, Cmu and radicand remain calculable.

A separate admissible nonzero-roughness control uses h_r=.1 and an independently
projected pressure direction p=x. It includes both pressure stress work and
the below-roughness volume load sampled from the outer H-side tangential
pressure gradient. Maximum action defect is 5.56e-16. Omitting the latter term
gives 2.25625e-3 on free velocity rows. The old zero-roughness fixture cannot
establish this coverage. No pressure-action source defect was found.

The saved binary arrays also support matrix-vector replay of predictor,
mass correction, pressure, continuity and splitting identities, without another
linear solve. The timed snapshot-Off reproducer and high-order normal-load
controls are separate read-only/verification receipts. Binary matrix-vector
replay defects are at most 5.43e-19 (predictor), 1.02e-20 (mass correction),
5.39e-21 (pressure), 6.53e-21 (continuity), and 1.51e-19 (splitting).
Snapshot-Off plus timing gives bitwise identical accepted/rejected states,
failure and native history. A separate reporting-only rebuild corrects the
open Functions projection description from the radial Gauss-32 rule to its
actual Gauss-max(8,p+3) mass-correction rule. Numerical observations and the
minimal trajectory/history are bitwise identical under both binaries.
The sealed governing trajectories keep their original runtime identity.

## Mesh matrix and pressure representation

All cases are isotropic p3/p1, serial flat 2D with the unchanged physical patch,
two or three streamwise columns, H_A=.4, delta=6.25, roof=20 and nu=1/10500.
Each timestep starts from its own fresh zero-step native history; same-space
initial coefficients are identical. Actual native counts match the predictions.

| Layout | Added nodes | Stored slots | Initial dPi/dz at original witness |
|---|---|---:|---:|
| B | none | 1085 | 2.877213723554e-5 |
| X | x=.1875 | 1545 | 2.877213723554e-5 |
| Z+ | z=6.5 | 1172 | 9.468164691212e-4 |
| XZ+ | both | 1669 | 9.468164691212e-4 |
| Z− | z=5.125 | 1172 | 5.009747526694e-6 |

The analytical supplied Pi=-Qzz is zero above delta. Its continuous P1
projection is not. An independent one-dimensional assembled P1 mass/load
calculation reproduces the native gradients, including both supplied hypotheses.
It uses native projection order 8 and the actual supplied pressure function;
it does not call the production constitutive kernel.

The common weak test is
`phi=x*(.25-x)*exp(-((z-6.25)/1.5)^2)`, with unchanged normalization.
Its native outer pressure and nonlinear-stress contributions are:

| Layout | Pressure load | Q load | Signed sum |
|---|---:|---:|---:|
| B / X | -2.099334432842e-6 | 1.739497758381e-6 | -3.59836674461e-7 |
| Z+ / XZ+ | -1.853215064860e-6 | 1.739497758381e-6 | -1.13717306479e-7 |
| Z− | -1.765873528010e-6 | 1.739498933861e-6 | -2.63745941489e-8 |

These are signed physical term loads, not an equilibrium-restoring source.
The initial receipts' V1 keys called raw coefficient-load Euclidean norms
`free_w`; those norms include constrained rows. The final source/export labels
them `unmasked_w`. They are not mass-dual norms and are not compared as nested
native-space norms. The fixed physical test avoids changing test spaces.
A smaller weak load can coexist with a worse pointwise pressure gradient.
Independent 48/96 integration of the supplied analytical Qzz and represented
P1 pressure gives a normal-load quadrature difference <=4.56e-15. The native
versus tight-load defect is <=1.26e-12 in the primary meshes and 7.46e-14 in Z−,
well below the interpreted load changes. This diagnoses the initialized
representation, not an assertion that an evolved field is an equilibrium.

| Layout / dt | Last accepted T | First actual rejection |
|---|---:|---|
| B / .002 | .002 | attempted .004, A95 radicand |
| B / .001 | .003 | attempted .004, A95 radicand |
| X / .002 | 0 | attempted .002, A95 radicand |
| X / .001 | 0 | attempted .001, constitutive inventory |
| Z+ / .002 | 0 | attempted .002, A95 radicand |
| Z+ / .001 | 0 | attempted .001, constitutive inventory |
| XZ+ / .002 | 0 | attempted .002, A95 radicand |
| XZ+ / .001 | 0 | attempted .001, A95 radicand |
| Z− / .002, .001 | .2 | none through .008, .02 and .2 |

The [Z− activation](cao-tamura-accepted-2-evidence/zminus-activation.json)
was recorded before evolution. It is the sole extra below-delta candidate;
there is no iterative search for a survivor.

A fixed physical observation grid resolves the union of partitions and active
breaks. A denser grid covers x near the artificial outlet and z around delta,
including the original (.2483691580323,6.2989252590303) witness. Interfaces keep
their native side/field provenance, and pressure is never regauged by subtracting
a mean. Native and independent inventories are both retained.

In all failed primary attempts, sampled loss is already present in the
momentum predictor. For B/.001 at .003, that predictor is outside the domain
but pressure correction returns the flow to it; the state is accepted.
The method deliberately uses coefficients frozen from its admissible input
during those linear solves. That intermediate excursion is distinct from
a corrected state where the actual constitutive law is required and undefined.
It is not evidence that pressure correction alone caused the rejection.

At the original geometry-25 point in B/.002's failed .004 attempt:

| Stage at that same physical point | div(U) | A95 radicand |
|---|---:|---:|
| Input | -7.83903e-4 | .474751 |
| Momentum predictor | -1.36457e-3 | -.0693987 |
| Pressure-corrected flow | -1.33918e-3 | -.0463276 |

The corrected raw physical tensor is

~~~
G = [[-1.337583966672e-3, 0, 8.272750987418e-6],
     [ 0,                0, 0],
     [ 1.732964491329e-5, 0,-1.597262199523e-6]]
~~~

Here k=2.863594597571e-5, epsilon=8.204721671112e-9,
tau=3490.17884, Cmu=.08949297, strain norm=.00109163137,
rotation norm=6.40419e-6, and strain eigenvalues are
(-.00089131287,.00044491913,.00044639374). The pressure gradient is
(-3.93925e-5,0,2.76553e-5). The predictor adds -5.79686e-4 to Gxx;
the mass correction adds +2.53805e-5, partly reducing the compression.
Planar trace removal gives radicand .462620 here, strictly as a diagnostic.
Full actual stress, production and traction at this inadmissible state are
undefined; no formal limiting tensor is substituted. The initial and accepted
outlet traction reconstructions, including their nonzero prescribed data, are
retained separately in the audit bodies.

## Mechanism and limitations

The projection tail, represented-versus-prescribed outlet traction and local
gradient response are coupled. X changes the response despite an unchanged
initial vertical pressure profile. Z+ improves the chosen weak sum while
amplifying the pressure gradient and failing earlier. Z− refines the interval
that generates the tail, reducing both the gradient and the weak sum without
changing any physical input. These support a representation/coupling mechanism;
they do not establish a unique cause or a general spatial-convergence theorem.

The retained local k and epsilon remain positive at the primary failure.
The inferred low free-stream epsilon makes k/epsilon large, so a small
bulk-scaled gradient error matters to the constitutive invariants. No epsilon
increase is used. The short outlet distance and one-sweep sensitivity remain
independent questions; this matrix does not vary the domain length or method.

Planar trace removal is diagnostic only: subtract div(G)/2 from Gxx and Gzz,
then recompute invariants. It is never applied to a trajectory or production.
Exact incompressibility itself cannot guarantee A95 admissibility:
G=diag(-a,0,a), k=epsilon=1 gives
`rad=1-18*a^2/(6.5+3*a)^2`, negative for
a>6.5/(sqrt(18)-3). The a=6 trace-free negative control has no physical tensor.
Thus neither weak continuity nor zero divergence is a general domain certificate.

## Conditional pilot and physical budgets

The .008 admission checks actual accepted endpoint coefficients on native and
fixed grids, the original covariance tolerance, mass, independent physical
scalar budgets under native/48/64 volume and 16/32 face controls, the .001
scaled weighted-RMS timestep limit, and supported pressure representation
improvement. Each continuation uses genuine unchanged-method restart history.

Both Z− members pass through .2 (100/200 total accepted steps), with matched
histories at each milestone. The independent final native outer-volume
radicand minima are .0670346/.0667780; the fixed-grid minima are
.179950/.178927. The broader native audit inventory gives the smaller
.0391853/.0389589 at (x,z)=(0,6.109375), illustrating why one inventory is
insufficient. All remain admissible. Fixed-grid covariance minima are
8.89794e-6/8.90186e-6. Corrections and additional governing solves are zero.

| Final scaled weighted-RMS dt difference | Value | Limit |
|---|---:|---:|
| U | 2.45436e-7 | 1e-3 |
| w | 1.70890e-7 | 1e-3 |
| k | 2.22369e-6 | 1e-3 |
| outer total epsilon | 2.43947e-6 | 1e-3 |
| augmented Pi, absolute level | 3.06579e-7 | 1e-3 |

These are equal-time physical weighted RMS values, with scales
U(Hh), K0, epsilon(H_A), and U(Hh)^2. They bound timestep differences and do
not establish temporal order. The .2 scaled change from the actual initialized
field is approximately 6.23e-5 for U, 3.97e-5 for w, 9.37e-4 for k,
9.53e-4 for epsilon and 1.012e-4 for Pi. These are numerical initial-state
drifts, not claims that the supplied finite-depth inference is an exact equilibrium.

Final scaled mass defects are 2.36e-16/2.44e-16. Sampled global divergence
weighted RMS is 8.88e-4/8.93e-4, and maximum is .04540/.04566; the broader
native inventory reaches .05625/.05657. Native maximum absolute outer
isotropic-production work is about 2.857e-4, at the outer H trace. Weak
continuity therefore does not justify dropping divergence contributions.
Maximum refreshed epsilon-trace lags are 1.44145e-7/7.39990e-8, retained as
splitting defects. The frozen imposed trace is not repaired after k transport.

Across all final volume/face controls, the largest scaled scalar-rate defect
is 3.262e-5/3.282e-5, below .02. Quad/face sensitivity and full region/source/
storage/flux primitives are retained in the complete audit bodies. Passing
reconstruction is separate from this finite-defect accuracy check.
Full evolved spatial uncertainty remains unresolved because the original
coarser and alternative primary layouts cannot supply admissible .2 endpoints.

T=.2 exceeds the nominal Uinf bulk transit .1717 (about 1.165 transits), but
not the interface-speed transit .3196. It is not a full-layer transit or a
developed precursor. The .02-to-.2 continuation alone costs 317.8/630.0 process
CPU seconds and 318.5/630.8 elapsed seconds for 90/180 additional steps, with
solver peak RSS about 362/368 MiB. The dense limit remains 3000.

Scalar budgets use actual consecutive accepted states, frozen flow/epsilon/pre-k
stages, dependent epsilon-boundary storage, and independently reconstructed
physical gradients on external and interface faces. Whole-k totals exclude
internal-interface fluxes; region totals retain their two contributions.
Storage-source+outward advective-minus-diffusive flux is reported with
current/frozen splitting separately. There is no scalar*div(U) correction and
no replacement of physical gradients by nominal zero Neumann data or reactions.

The rate scale uses actual patch volume times K0 or epsilon(H_A), with
Hh/U(Hh)=1. On this patch the volume is 5. A synthetic doubled-length unit check
must double the scale and halve the normalized defect; it does not alter a
physical result. Mass and scalar targets remain 1e-8 and .02 respectively.

## Rough-wall operating heights and extraction sensitivity

The independent study has 17 predeclared deterministic scenarios: central,
four height-only joint anchor extremes, four velocity-only extremes, four
two-axis joint extremes and four coherent shared-axis shifts. Every interpolation
rebrackets the observations; extrema include cases where an observation height
crosses the 16 or 40 mm anchor. The disjoint source-height envelopes preserve
point order, allowing enumeration of all adjacent candidate pairs and boundary
vertices. Normalization at 40 mm and calibration at 16 mm are propagated together.

The central values are pinned exactly to CT1. Its y/H and published-axis
coordinate columns were rounded separately; their tiny numerical difference
is retained explicitly instead of changing the central roughness fit.
Other scenarios use the frozen axis-coordinate and value extraction envelopes.

Each model/scenario is fitted once at 16 mm and that h_r is held fixed at
H_A=5/8/11 mm and every hill column. There is no hill fit or per-height refit.
These are effective model-conditioned roughness inputs, not new measured surface
properties or confidence intervals. The assumed 3.5 mm displacement, approximate
delta/Re, calibration and experimental measurement uncertainty are separately
unquantified. Coherent shifts cover a shared extraction component; no
statistical independence or root-sum-square model is assumed.

The study uses the unchanged independent finite-nu flat and metric machinery:
counterfactual p0(K), monotone weighted cutoff, all Q/Q_J stationary roots,
J-weighted normal diffusion, one-sided foot limits, zero-curvature controls,
and tighter quadrature. Flat covariance is reported only for its declared
flat velocity/strain field. Frozen curved k positivity is not a terrain
equilibrium or covariance/stability certificate.

The historical thin ratio <=.1 remains a separate flag. A <=.25 sensitivity
flag is also shown; it does not replace the old gate or supply a model-error
bound. Normal-map reach and positive local Jacobians are distinct.
Five millimetres is at the nominal cube tops; 8/11 mm are above them.
An eligible kernel at those heights does not demonstrate roughness-sublayer
homogenization. The first visible rough mean point is about 10.24 mm.
The 5/8 mm matching values therefore have no direct visible velocity anchor.

The study completes 153 flat candidates and 5967 frozen hill columns.
All 17 scenarios are reference-admissible for every model/height. Minimum
column k/K is 1 over the tested stations; actual field realizability remains
separate. The extraction-only model-conditioned ranges are:

| Model | h_r range, mm | Smallest clearance below H_A=5 mm |
|---|---:|---:|
| A93 | 2.89538–3.97258 | 1.02742 mm |
| A95 | 3.67941–4.93517 | .06483 mm |
| E94E95 | 3.50322–4.69375 | .30625 mm |

The narrow A95 clearance is not robust to unspecified physical/measurement
uncertainties merely because the extraction envelope fits below 5 mm.
All local Jacobian and sufficient global reach checks pass. Maximum thin
ratios are .0986960, .157914 and .217131 for 5/8/11 mm. Thus 5 mm passes .1;
8/11 mm retain failure of that historical threshold while passing the separate
.25 flag. Neither flag supplies a curved-flow error bound. Matching velocity,
traction, budgets, cutoff, flat stress margins and all interior k extrema are
in the adjacent tables; native central flat-kernel controls pass.

Five millimetres is a useful **mathematical frozen-reference candidate**, with
a physical interface decision still needed at the cube tops. It has not been
substituted into the H_A=16 mm native pilot. The kernel matching velocity and
traction in curved-column tables are its counterfactual flat-compatible scales,
not actual hill traction or a curved momentum equilibrium. The old smooth 8 mm negative-k references remain an unresolved
SS preparation issue, despite Lambda>1. RR eligibility does not qualify SS.


## Exact identities and reproduction

The governing/diagnostic runtime was built from clean
`a897896524b8761d3a954d491f96d1b33941d94f` (the computed-flow C++ change is
`3c6b826e`). The reporting-only runtime was built at
`8ba59b9d2016eb289006261474d80f1959963ba1`. No governing cohort is relabelled
as a result of that reporting binary.

| Executable | Governing runtime SHA256 |
|---|---|
| IncNavierStokesSolver | 858c9cdb8f343d0f68f109d4992123021f4d87d6e1a32d9cf4a7cb902f11984b |
| RANSPreservationProbe | 43ecbaed68ecc2e6424e05df1b8090646f58cf9914ca82abd2758f6862d17475 |
| RANSNewModelsProbe | 34dab79530e4e2f51b510d7d5709dbd727bd77b86da5c7a3f728cb6da0c6a669 |

The reporting solver SHA256 is
`c8b508dd879ff5d4b93a55f764910d63c769c97d28821644209e04a1211fa35f`.
Its probe SHA256 is
`cbbe767ceeef0a58e72ec348d9e643887915acbb3e73758521d15eb9d1af3155`.
Full loaded-library bytes/paths, descriptor, mesh, basis, actual quadrature,
live explicit/implicit boundary hashes, consumed initialization, checkpoint and
native-history identities are retained in the runtime, initialization,
coupling, request/seal and raw-artifact receipts.

The independent scalar/constitutive kernel sources and both specifications are
unchanged from CT1. The data auditor again passes 98 series, 1151 points and
39 reported facts. The source paper SHA remains
`adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df`.
All fieldwise comparisons include terminal times, physical coordinates, sides,
weights and one absolute traction pressure gauge.

Use a NEW destination for each command. The existing directories below are
immutable inputs/evidence, not destinations to overwrite:

~~~bash
CT=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-2-20261112
OLD=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-1-20261111
S=$CT/source-after-reboot/solvers/IncNavierStokesSolver/Tests/RANS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR=$CT/tmp LD_LIBRARY_PATH=$CT/runtime/lib
ct_run() {
  tag=$1; shift
  python3 -B "$S/ct_accepted_command.py" --control "$CT/NEW-control-$tag" -- "$@"
}
ct_run initial python3 -B "$S/ct_accepted2_diagnosis.py" \
  --runtime "$CT/runtime" --previous "$OLD" --mode initial \
  --output "$CT/NEW-initial"
ct_run matrix python3 -B "$S/ct_accepted2_diagnosis.py" \
  --runtime "$CT/runtime" --previous "$OLD" --mode advance \
  --preparation "$CT/NEW-initial" --output "$CT/NEW-matrix"
ct_run zminus-initial python3 -B "$S/ct_accepted2_diagnosis.py" \
  --runtime "$CT/runtime" --previous "$OLD" --mode initial --layouts Zminus \
  --output "$CT/NEW-zminus-initial"
ct_run zminus python3 -B "$S/ct_accepted2_diagnosis.py" \
  --runtime "$CT/runtime" --previous "$OLD" --mode advance --layouts Zminus \
  --preparation "$CT/NEW-zminus-initial" --output "$CT/NEW-zminus"
ct_run wall python3 -B "$S/ct_accepted2_wall.py" \
  --runtime "$CT/runtime" --previous "$OLD" \
  --data "$CT/source-after-reboot/docs/plans/cao-tamura-validation-data" \
  --output "$CT/NEW-wall"
~~~

Before each continuation, audit its actual parent cohort. For example:

~~~bash
ct_run audit python3 -B "$S/ct_accepted2_audit.py" \
  --runtime "$CT/runtime" --preparation "$CT/NEW-zminus-initial" \
  --trajectories "$CT/NEW-zminus" --output "$CT/NEW-zminus-audits"
ct_run continue02 python3 -B "$S/ct_accepted2_continue.py" --root "$CT" \
  --preparation "$CT/NEW-zminus-initial" --baseline-initial "$CT/NEW-initial" \
  --parent "$CT/NEW-zminus" --audits "$CT/NEW-zminus-audits" \
  --target .02 --output "$CT/NEW-continue02"
~~~

Audit that .02 cohort before using it as the .2 parent. The optional explicit
preparation arguments preserve the tested defaults and bind the fresh
preparation/baseline hashes in new requests. Never substitute a timestep in
a checkpoint. The complete executed argv and original defaults are retained
in the evidence; interface-only reproduction additions do not relabel old runs.

The original failed comparison collector and its read-only recovery remain
separate requests. The recovered comparisons use the original eight report
hashes; no native trajectory was rerun merely to fix that collector.
The two-axis study's first precision-preflight failure and the first Git
resource-control failure also remain archived. Exact executed argv, CPU/wall
costs and current source hashes are in the complete receipts. Original build
argv is preserved; a rebuild at another path has a new binary identity and
requires the comparison control rather than an edited old manifest.

The timed minimal reproducer consumes 8.82 process-CPU seconds and 9.33 elapsed
seconds. Its dense triangular solves total .0084 exclusive CPU seconds,
momentum factorization .000873 and consistent-mass factorization .000336.
These startup/two-attempt timings are not a scalable performance model; they
reinforce the need to retain and profile full physical checks.

## Verification, provenance and next gate

The open weak/matrix, pressure-offset, history/rollback, ordinary closed CLS,
cache and native diagnostic controls pass after the C++ instrumentation change.
Unaffected closure/wall evidence is reused through exact source identities.
No old trajectory is recast as a result of this runtime.

The isolated build uses GNU 11.4 Release, MPI enabled/rank one, four compiler
workers, and one-thread numerical libraries. Each numerical family, including
read-only NumPy/SciPy work, runs through the owned launcher. The historical
CT1 79-thread incident remains failed. The final aggregate audit records **147**, with an unowned 129-thread gdb
(PID 22924) contributing 129 threads and external total 131 during the later
reporting build. At 1794451796.1175 the guard paused only our 16-process family;
it resumed it at 1794451798.9317, a 2.814 s family suspension. The physical
cohort was already complete. No foreign process was controlled. The hard
72-worker gate remains failed; earlier partial statements based on the
pre-build receipt are superseded explicitly.

A separate private-allocation failure is retained: leased-family counts up to
16 plus the launcher/monitor exceed the declared 16. V1 omitted those two
helpers and falsely passed that subgate. V2 counts them, preserves the original
high-water evidence, and reports aggregate and project failures separately.
The first Git correction also failed because index preloading is distinct from
index reader threads. The completed corrective test uses index.threads=1,
core.preloadIndex=false and pack.threads=1, observes one Git thread, and forces
a real two-job native-test build. Its observed project peak is 14, including
control helpers. This repairs subsequent allocation but does not erase either
failed gate. Process CPU and elapsed time are separately recorded; the 2.814 s
suspension is not reported as compute CPU.

Preparation failures are retained: missing element metadata on test rows,
a list/array comparison collector error, central coordinate-column precision,
and the reboot/mount interruption. Primary native jobs are reused only through
a new read-only recovery receipt binding runtime and native artifact hashes.
Neither collection nor a passing process exit advances a physical cohort.

Heavy state/matrix/checkpoint data stay on Linux. The adjacent compact evidence
contains full comparison/classification bodies, raw hashes, requests, tables,
native counts, source/runtime/library/consumed-input identities and exact commands.
The numerical small-patch result is ready to serve as an equivalence reference;
its procedural resource failures remain explicit. The smallest next numerical
control is one predeclared combination of x=.1875 and z=5.125 (predicted 1669
slots, to confirm natively), with unchanged inputs and .002/.001 histories to
the same horizon. This tests streamwise/outlet coupling after repairing the
pressure representation; it is not another unbounded mesh search. No epsilon
retuning is justified by this package.

Before scalable flat-flow development, qualify that evolved representation/
pointwise-divergence control and preserve this dense reference under any sparse
backend or immutable-row cache. Do not reuse A95 validity with incomplete
velocity-gradient dependencies. Measured solve cost alone is not the dominant
cost in the small prototype; full checks must remain.

Before native curved coupling, implement and verify the complete Cartesian
metric momentum/scalar/pressure/interface and roughness-force equations,
normal-map geometry and vertical-height observations. The current mapped guard
remains. Before S/R approach-flow qualification, settle the physical rough-wall
interface and inferred turbulence profiles, and independently resolve the
smooth reference's negative-k issue. Before hill validation, admit those
approach flows, a scalable multicolumn domain and native curved coupling, then
use the untouched frozen validation roles and scores. No crest, reattachment,
EP/transition, terrain, MPI/3D or experimental prediction study was attempted.\n