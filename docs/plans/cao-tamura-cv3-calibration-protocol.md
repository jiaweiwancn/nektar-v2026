# CV3 protocol v2: EP-AWF transition calibration first

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: computational experiment plan amendment
- Origin Date: 2026-11-01
- Verification Status: R1 reference/native and standard checks pass; CLS developed-floor gates unmet; R1 OPEN and CV3 IN PROGRESS
- Version Label: cao-tamura-cv3-protocol-v2
- Supersedes: [original-first protocol v1](cao-tamura-cv3-calibration-protocol-v1.md), preserved byte-for-byte from commit `476ee5b`

## Decision and retained evidence

The primary calibration route is now the actual proposed wall-treatment
arrangement on an empty floor:

```text
EP-AWF upstream -> upstream transition -> original AWF downstream
```

Successful development of an all-original-AWF precursor is no longer a
prerequisite for this route. Uniform original AWF remains a separate baseline
and diagnostic track. Compatible standard log-law and switched-AWF comparisons
also have their own validity/status records; primary EP calibration does not
certify them automatically.

The [v1 attempt](cao-tamura-cv3-calibration-attempt.md), its 32-case design,
13 developed standard cases, unsuccessful CLS cases, later controls, source
hashes and raw results remain unchanged. They describe the original-first
route, not an execution of this revision. No primary composite precursor is
yet calibrated or frozen. This amendment does not close CV3 or reclassify a
failed baseline as a successful simulation.

## Revised sequence

| Task | Deliverable / gate | Dependency | Status |
|---|---|---|---|
| CV3-R1 — Reference recovery and composite-floor feasibility | Reproduce the documented CLS EP-AWF reference on the current solver; establish valid developed composite-floor candidates for S/R and both closures, with audited active correction, transition and original section | CV1–CV2 and delivered CV3 fixes | OPEN: reference/standard checks pass; CLS development gates unmet; [execution record](cao-tamura-cv3-r1-reference-feasibility.md) |
| CV3-R2 — Primary composite calibration and freeze | Route-aware bounded common-parameter search, uncertainty/residual assessment and versioned primary S/R manifests | CV3-R1 | PLANNED |
| CV3-R3 — Reserved flat crosschecks | Score reserved Fig.3d only for the frozen, validated primary policies; record unsupported observations | CV3-R2 freeze | PLANNED |
| Baseline / transfer track | Original-AWF and compatible-log comparisons with transferred physical inputs; successful, failed and unexecuted cases distinguished | Primary parameters available for transfers; diagnostic originals can run independently | SEPARATE; NOT A PRIMARY EP GATE |

CV3 closes after R1–R3 pass their gates, including published common/conditional
fits, sensitivity and incident residuals. Every requested baseline retains an
explicit status. A baseline failure can remain unresolved while the primary EP
study advances; its accuracy or stability cannot then be claimed. CV4–CV8 keep
their numerical, sampling and observable-level requirements for each case used.

## Physical route and comparison reference

Retain the canonical standard k–epsilon and high-Re CLS equations and
coefficients, physical wall/Yap distance, integrated AWF strip continuity,
0.8 m roof and primary resolved `Molecular` diffusion. Hh=0.040 m is the hill
height; hm/P/N=11/5.5/22 mm are physical wall-strip heights. Smooth ks is
exactly zero. Rough equivalent sand-grain ks, observational displacement and
reported z0 remain separate quantities. This revision changes the study order,
not the closure or wall physics.

Use this provisional calibration anchor before the CV5 width study:

- Domain x=-1.0..1.2 m, empty Plane floor.
- Fixed correction chi_fixed=1 where admissible; the C2 spatial envelope is
  one upstream of xa and zero downstream of xb.
- Width 2Hh=0.080 m; xa=-0.280 m, xb=-0.200 m. The full ramp ends before
  the hill's upstream foot at x=-0.100 m.
- Observe the incident profile at x=0, inside the downstream original-AWF
  section. Keep the full original-AWF section in the native solve.

The anchor is a declared starting geometry, not a selected width or a
width-independence result. Every later hill case must have a floor counterpart
with the identical inlet/surface parameters, mask width/placement and physical
boundary conditions. Changing width retains the frozen physical inputs and
requires a matching floor reference; its incident changes are measured rather
than removed by retuning. Material changes after numerical refinement or width
selection reopen the affected calibration and require a new version.

## Next task: CV3-R1

This is a verification and feasibility study, not another 32-case fitting
campaign and not a hill-validation run. Store its cases under a fresh Linux
directory such as `calibration/cv3/ep-transition-v2/r1/` outside the source tree.
Bind the current revision/source, binary, geometry library, configuration and
native observation hashes. Preserve the v1 run directories and evidence.

### R1.1 — Prepare route-aware cases and receipts

Reuse the native generator, preparer, development/restart and audit routines.
Make the primary route explicit in case/manifest/observation metadata. Distinguish
the composite envelope, a uniform-EP diagnostic and an all-original baseline;
the last two cannot masquerade as composite results. Check the exact chi(x),
endpoint ordering, requested/active correction and original recovery at zero
weight. Inlet/reference initialization and policy epsilon collars must agree.

The delivered v1 `cao_tamura_campaign.py` and the `select`/`freeze` commands in
`cao_tamura_calibration.py` hard-code original-first selection and all-policy
freezing. They remain legacy tools. Adapt/version the primary workflow before
using it for R2; do not run the old defaults as a v2 campaign or bypass their
guards by manually setting calibration-complete flags. R1 can reuse the generic
native case/develop operations with explicit composite configurations.

### R1.2 — Recover the known CLS EP-AWF reference

Reproduce the documented accepted degree-four, five-flow-through reference
using the current solver and its recorded physical/outer boundary data. Use the
CLS-compatible Q/kappa and the original benchmark's resolved-diffusion policy;
experimental common-inlet reference constants do not replace them. Keep chi=1
uniformly for this reference test. Require the prior 0.1% maximum-over-time
profile gates, positive fields, converged physical residuals and no fallback.
Repeat the selected small-perturbation recovery control with its own recorded
decay criterion; its initial departure deliberately exceeds the unperturbed gate.

Use actual unforced URANS evolution, without field resets, frozen turbulence
coefficients or reference residual forcing. Record a regression if the former
accepted case no longer passes. A Molecular bridge is assessed against its
finite-error result; exact finite-nu log equilibrium is not assumed. A passing
reference is a prerequisite for attributing experimental failures to the new
inlet, transition or constitutive behaviour.

### R1.3 — Four unfitted experimental composite-floor cases

| Surface | Standard k–epsilon | High-Re CLS |
|---|---|---|
| S | `standard_fixed`, ks=0 | `cls_fixed`, ks=0 |
| R | `standard_fixed`, ks=0.004 m | `cls_fixed`, ks=0.004 m |

Use the declared provisional anchor and the v1 centre physical inlet controls
(delta_in=0.21 m, k/epsilon multipliers one, MeasuredRms/ShearLength inference)
as unfitted starting points. Keep common physical inputs above N across the
closures, including the declared common inlet Q/kappa. Each closure keeps its
canonical coefficients and required H..N policy epsilon collar. Source Re_H
scaling, frozen observations, rough d scenarios and positive outer turbulence
remain as in CV1–CV2/v1. None of these four definitions is a calibrated case.
The centre is a starting probe, not a required successful parameter point.
If it is inadmissible, retain the failure and declare a small feasibility
bracket before further runs, with common physical inputs across the closures.
R1 needs one valid common starting pair per surface; it does not require the
arbitrary centre to succeed or a full parameter search to be repeated first.

Use a ramp-resolved mesh rather than the v1 all-original two-column mesh.
Align elements with xa/xb and begin with at least four elements across the
ramp, degree-four velocity/scalars, degree-two pressure and eight GLL points.
Declare all normal/outer layers and quadrature; short unscored throughput/CFL
checks select a fixed-step recipe before the long runs. Verify its sensitivity
where needed to distinguish a numerical rejection. Full h/p/dt/domain and
observable uncertainty remain CV4 work. A short solve cannot establish them.

Advance in stages: native initialization/first steps, one flow-through, then
at least three unforced domain flow-throughs and two further one-flow-through
extensions. Steady candidates use the existing pointwise/full-field drift
criteria. A valid bounded unsteady solution needs a stationarity/sampling
protocol using the plan's sampling/uncertainty requirements declared before
its assessment. R1 remains open until that assessment passes. Failed timesteps
cannot be counted as physical unsteadiness or averaged into acceptance.

### R1.4 — Audit the whole route and its failure modes

Record profiles and native wall diagnostics in the full-EP region, through the
ramp, just downstream, at the upstream foot and at x=0. Retain both one-sided
traces where applicable. Audit positivity, full free-row stage residuals, weak
divergence, projected/strong normal errors, resolved mass, integrated strip
normal-plus-side flux and independent strip momentum identities.

Apply the existing 1e-7 native free/weak/projected/mass gates and independent
oracle tolerances. For steady development retain the v1 drift limits: velocity
2e-4 U_inf, shear 2e-3 u_tau^2, relative k/epsilon 1e-3 at every native point,
and friction/delta99 0.2%, for two successive one-flow-through extensions.
Also publish strong divergence and pointwise wall mismatch; weak convergence
does not certify them. Carry out the existing representative MPI/history
checks after any implementation change affecting evolution or provenance.

For the reference require active chi=1 with no fallback. For the experimental
route require the prescribed envelope where reference admissibility holds and
report every admissibility/fallback reason and location. Prescribed chi=0
downstream is not fallback. A fully reverted upstream region cannot qualify
as an EP-AWF result. Set any additional quantitative coverage acceptance rule
before calibration; do not lower it after examining fit scores.

If a composite case fails, retain its first growing-error locations and last
valid full-history checkpoints. A matched uniform-EP run on the same mesh,
inlet and numerical recipe can isolate the transition/original section. Its
successful evolution does not certify the composite route. Likewise, failure
of an all-original baseline alone does not reject the EP composite route.
Do not introduce a closure regularization or replace integrated continuity
without a separately declared derivation and validation track.

### R1.5 — Exit decision

Pass R1 only when the reference controls and all four composite cases meet
their declared validity/development/coverage gates, with native receipts and
incident profiles. Experimental agreement is not an R1 feasibility gate.
Publish the outcome, failures, numerical recipe, runtime and limits; commit
and push the task evidence. A timeout, unbounded growth, inadmissible correction
coverage or an incomplete case leaves R1 open. No all-original success is
required to pass R1.

### R1 execution preregistration — 2026-11-01

Before experimental evolution, use the recipe in `cao_tamura_r1.py`: p4/P2,
eight GLL points, twelve logarithmic outer layers with explicit H/N and outer
blend endpoints plus three blend cells. The initial x nodes are -1, -.6,
-.32, -.28, -.26, -.24, -.22, -.20, -.16, -.10, 0, .4, 1.2 m. Thus the
four ramp elements have equal 20 mm widths. Start with six unscored dt=.01 s
implicit-advection steps, Newton–Krylov, grad-div=.01, current-operator
tolerance 1e-11, native inverse tolerance 1e-8 and pressure inverse tolerance
1e-13. Pilot timings and convergence select the long-run step before scoring.
All numerical revisions get fresh directories and immutable receipts.

Use a strict coverage gate: at every accepted post-startup native wall trace,
actual chi must equal the fixed weight times the independently evaluated C2
envelope within 5e-9. No requested-positive trace may fall back for cutoff,
regime or representability. Every accepted time must sample full-EP, ramp and
downstream regions. Downstream prescribed chi=0 with SwitchInactive reason is
valid. Counts include one-sided duplicates and are not area fractions.

The preparer's scalar expressions now evaluate at the physical wall-foot x.
In H..N the epsilon collar uses the actual wall policy: the linear spatial
blend of Original and EP wall epsilon, rather than full EP epsilon everywhere.
The initial floor Dirichlet trace uses the same expression. Velocity retains
the existing solenoidal streamfunction, and all physical inlet profiles at
and above N remain common across closures. This is an initialization correction;
it changes neither the native time-evolving AWF nor the closure equations.

Further unscored numerical controls are declared in hashed JSON files in the
R1 run root **before** launching the corresponding cases. They test: dt=.001;
uniform EP at that timestep; a current-operator tolerance of 5e-14 with a
physical Newton difference step of 1e-4; compact x nodes -1, -.28, -.26,
-.24, -.22, -.2, 0, 1.2 with the same four ramp elements, p4/P2/eight points
and twelve normal layers; matched serial-direct/four-rank-iterative solves
at dt=.05 and a dt=.01 MPI control. Unit Picard-map and pressure relaxation
is a separate Newton-conditioning control; converged stage equations and
all physical residual gates remain unchanged. These are numerical checks,
not a fit search or a claim of mesh/timestep convergence.

A startup-continuation pilot uses the already verified fixed initial weak
load, fading to exactly zero after 0.2 domain flow times, dt=.05 and twenty
steps. Its forced times cannot count as development or enter statistics.
Any use in development retains at least three subsequent unforced flow
times and two further one-flow-time drift passes. The R1 staged developer
supports an initial one-unforced-flow-time check, continuation to three,
and then one-flow-time extensions. The accepted analytical reference controls
remain unforced throughout. A separate five-flow molecular reference uses
the same 0.1% engineering error gate, reporting its finite errors explicitly.

Before further physical-input probes, declare the two-point feasibility
bracket epsilon_scale={1,2}, with delta_in=.21, k_scale=1, ks_S=0 and
ks_R=.004 unchanged. The second point is the existing v1 upper epsilon bound;
it tests whether a smaller turbulent-diffusion scale permits positive,
well-converged evolution on the ramp-resolved mesh. Apply it to both closures
and surfaces without inspecting objective scores. It is an unfitted starting
point, not a calibration selection or evidence that the centre has no solution.

One additional feasibility bracket tests outer_length={.025,.010} m, retaining
the measured mean/k inlet family and all other physical controls. The
ShearLength epsilon inference has interior valleys where fitted mean shear
changes sign; changing this unmeasured length tests their role without
changing closure coefficients. Declare the bracket before its native probes,
use the same value across closures for any accepted surface pair, and retain
the original .025 m controls. This does not establish an experimentally
preferred length or numerical convergence. R2 must retain this inference
uncertainty if the alternative is used.

Fixed timesteps may be selected separately for each surface/closure from the
unscored numerical pilots; no timestep changes during a trajectory are allowed.
Common physical inputs, mesh, mask and native validity/development gates remain
mandatory. A smaller or larger timestep is not a retuned inlet parameter;
steady temporal terms vanish in the zero-drift limit. This provision does not
certify temporal stability or accuracy, which still requires CV4 sensitivity.

The optional [preconditioned-residual Newton map](cao-tamura-cv3-r1-newton-map.md)
has its own fixed-point derivation and validation contract. After its matched
S/R short-run comparisons, `pilot-residual-map-dt05.json` preregisters a
three-step dt=.05 check on the compact mesh for all four cases, using the same
strict physical gates. The dt=.01 development trajectories remain recorded.
Any faster trajectory starts afresh at a fixed timestep; histories cannot
change Newton maps or splice failed/retired segments into acceptance.

After the centre rough-CLS rejection and its exact-history repetition, declare
`development-graddiv01.json`: increase the existing consistent grad-div
coefficient from .01 to .1, retaining centre physical inputs, the compact
mesh, dt=.01, integrated continuity and every acceptance gate. This tests
the role of strong-divergence growth left by the P4/P2 weak constraint; it
does not prove its cause. It changes a consistent discretization term, not
the CLS stress, wall law or positivity guards. An accepted surface pair
retains common physical inputs; report the numerical coefficient explicitly
and retain its sensitivity as CV4 work. Use a fresh trajectory and preserve
the failed .01 control. No fit or reserved observation selects this control.

`development-full-mesh.json` retains dt=.01, grad-div=.01 and centre physical
inputs, restoring the original preregistered x nodes in `cao_tamura_r1.BASE`.
The compact mesh directly joins a .72 m upstream element to the .02 m ramp
elements, and a .20 m downstream element follows the ramp. The original mesh
adds approach and wake nodes. Its developed control tests whether this mesh
gradation contributes to the growing strong wall error; a short throughput
pass on either mesh is insufficient. It is a targeted rejection diagnostic,
not the later three-level convergence study. Any accepted closure pair must
share the selected mesh.

After the rough centre, epsilon-only, grad-div and original-mesh rejections,
declare two separate controls on the compact mesh at dt=.01, grad-div=.01.
`development-startup1.json` uses the already verified fixed initial weak load
with its native C2 fade over one domain flow time. The first leg then includes
one **unforced** flow; later legs retain at least three unforced flows and two
one-flow drift passes. No forced time or forced observation can qualify as
development. All analytical reference controls remain unforced throughout.
`development-corner.json` is an unforced feasibility point with k_scale=1.5
and epsilon_scale=2, the upper values of the existing v1 input bounds; delta,
roughness and canonical model coefficients remain fixed. It tests sensitivity
to unmeasured turbulence inference, without using fit scores. Any accepted
surface pair must use the same physical corner in both closures. These are
finite declared controls, not a new calibration search or a relaxation of
positivity, coverage, residual or drift limits.

`development-pressure3.json` is a separate pressure-resolution control: retain
the compact p4 velocity/scalar mesh, dt=.01, grad-div=.01 and centre physical
inputs, increasing pressure from p2 to p3. It tests whether the lower-order
weak continuity space contributes to rejection. Native free/weak/projected
and mass gates, pointwise-error reporting, positivity and correction coverage
are unchanged. It is not a pressure-convergence result; any selected pair
must share its mesh and pressure space, with refinement left to CV4.

`development-krylov512.json` targets the recorded centre-CLS rejection, whose
late 240-vector Newton GMRES action residual remains above its .01 request.
Use a fresh compact p4/p2 trajectory with centre inputs, dt=.01 and grad-div=.01,
increasing only `RANSNewtonKrylovDimension` to its supported maximum of 512.
This tests the linearized-solve capacity without relaxing nonlinear, accepted
positivity or physical-equation gates. The original rejection and exact-history
diagnostics remain separate evidence.

`development-ramp8.json` halves only the streamwise element widths across
the fixed physical ramp: eight 10 mm elements instead of four 20 mm elements.
Keep the compact outer x columns, p4/p2, original normal mesh, dt=.01,
grad-div=.01 and centre physical inputs. The original centre's largest
pointwise normal mismatch lay inside that ramp; the full-x-mesh control
left its four ramp elements unchanged. This is a targeted resolution check,
not an eight-element acceptance shortcut or a completed convergence study.
Any selected closure pair must use the same resulting mesh.

After these rejections, `pilot-uniform-ep-matched.json` uses the centre rough
inputs, compact p4/p2 mesh, dt=.01, grad-div=.01 and the verified residual map,
changing only the route to `UniformEP`. Request one domain flow-through (152
steps) without fitting. This matched diagnostic removes both the transition
and downstream Original section; it cannot qualify as a composite precursor.
`pilot-ramp8-debug.json` repeats only the refined-ramp first stage with native
debug logging, to distinguish its numerical plateau from an accepted flow.
No reference forcing or weakened acceptance criterion is introduced.

The refined-ramp debug control shows an unresolved linear-action residual
despite positive scalar trials. `pilot-ramp8-capacity.json` therefore tests
four fresh dt=.01 stages on that same mesh with 512 Krylov vectors and up to
200 Newton iterations. Every tolerance and independent physical gate remains
unchanged. This is a numerical iteration-budget check, not a physical-model
revision or development evidence; only successful staged continuation can
qualify a later precursor.

`development-direct-inverse.json` isolates the accuracy of the cached native
inverse in the residual map: use serial `DirectMultiLevelStaticCond` on the
original compact p4/p2 centre trajectory. Retain dt=.01, grad-div=.01, 240
Krylov vectors, 60 Newton iterations and all physical gates. The earlier
development trajectories used an iterative native inverse with tolerance
1e-8; short A/B tests used a serial direct inverse. The measured discrepancy
between GMRES's target and its recomputed action residual motivates this
matched long-run check. Direct inversion changes the numerical map accuracy,
not its fixed-point equations or the closure/wall physics.

Retain a matched legacy-map check beyond the rough centre's stopping point:
`conditioning-legacy-current-solve.json` continues the preserved, valid t=.45
full history with its original immutable CurrentSolve binary for ten steps.
The earlier CurrentSolve development was deliberately retired at t=.48; it
had not rejected there. This distinguishes a later numerical-map rejection
from the behaviour of the original fully solved map. It is a diagnostic, not
a current-binary developed precursor or permission to splice trajectories.

After four accepted refined-ramp stages with the larger iteration budget,
`pilot-ramp8-capacity-dt05.json` repeats three fresh stages at dt=.05. This
uses the previously declared larger-step bracket on the refined mesh, before
any long trajectory is selected. Retain the dt=.01 receipts; passing this
pilot supplies no timestep-convergence or development claim.

The larger-step refined-ramp pilot rejects before its first accepted step.
Use `development-ramp8-capacity.json` for fresh dt=.01 development of S/R and
both closures, with the successful 512-vector/200-iteration recipe. Keep the
same physical centre inputs and eight-element ramp in each closure pair.
Request the existing one-/three-/four-/five-flow staged gates, extending for
two qualified drift passes; permit 14400 seconds per stage because this
measured recipe is slower. The rough CLS allocation is eight MPI ranks, smooth
CLS four, each standard counterpart one, and the control allocation one CPU;
the separate serial smooth diagnostic can share the remaining single CPU.
Never exceed the 16-worker project cap. Frequent full histories bind any
early rank-comparison diagnostics; no pilot/retired prefix is spliced into
these fresh trajectories. This is a feasibility assessment at the stated
resolution, not a convergence result or a calibrated selection.

`pilot-ramp8-inverse12.json` tests the cached inverse's finite-difference
accuracy: repeat three refined-ramp rough-CLS steps in one rank with native
inverse tolerance 1e-12, 240 Krylov vectors, the original 60-iteration budget
and native debug logging. The one-rank 1e-8 baseline needed 55 Newton
iterations for its first stage, and its recomputed GMRES action residual
missed the linear target. This control changes inverse accuracy only; require
the same physical gates and matched-checkpoint fields within 1e-6. It supplies
no development result. Retire the unselected serial smooth direct-inverse
assessment explicitly after its valid one-flow audit to release that CPU;
retain its accepted history and distinguish retirement from rejection.

The inverse-accuracy debug control still exhausts its 240-vector linear
space. The [pointwise scalar-scaling control](cao-tamura-cv3-r1-newton-map.md)
therefore has a separate declared derivation and validation contract. It
changes Newton weighting and finite-difference normalization, retaining the
same current weak equations and positivity domain. Default field-maximum
scaling remains available. No R1 selection may use its new binary before
default regressions, current references, matched-field and MPI/history/scaling
mismatch checks pass; any selected development starts afresh with bound mode
metadata and the existing physical/coverage/drift gates.

After its four matched-field checks and all three current-binary analytical
controls pass, declare `pilot-compact-pointwise55.json`: the common compact
four-ramp-element centre, fresh unforced dt=.01 evolution for 55 steps,
pointwise scalar scaling and the already declared 512-vector/200-iteration
capacity. This crosses the recorded t=.49 rough-CLS rejection without
splicing a history or changing physical data. It is a conditioning pilot,
not a developed precursor. A selected long trajectory must start afresh
after current-binary native guards and history/mode-mismatch checks pass.

Those native and 1-to-4-rank history checks pass, including explicit scaling
mismatch rejection. `development-pointwise-scalars.json` therefore declares
fresh centre composite assessments using the same compact four-ramp-element
mesh, fixed dt=.01, pointwise scaling and 512-vector/200-iteration capacity.
All stages retain the prescribed unforced one/three/five-flow progression
and unchanged physical/coverage/drift gates. Smooth CLS and the standard
counterparts can proceed while the bounded rough-CLS conditioning pilot is
assessed; no failed or pilot segment can be spliced into a selected trajectory.

The rough compact pointwise pilot also rejects its t=.49 stage. Before any
new search or closure change, declare one final finite physical-input corner:
`k_scale=.75`, `epsilon_scale=2`, the existing v1 lower-k and upper-epsilon
bounds, with delta, roughness, canonical coefficients, mesh, mask and all
validity/development gates retained. This reduces the prescribed incoming
kinematic stress scale; it is a sensitivity of unmeasured inference, not a
fit or a proof of stability. Use identical inputs in both closures for each
surface. Start with the failed rough surface; the same declared corner is
available for smooth only if its centre assessment fails. Record the two
starting pairs separately and retain all outcomes without objective scores.
Any selected trajectory is fresh. A five-step four-rank lower-corner pilot
precedes its eight-rank rough-CLS assessment, allowing an independently
time-bound four-to-eight-rank field comparison at t=.05. Success at this
checkpoint is not a developed-flow or convergence result.

The pointwise compact smooth centre and rough lower corner also reject stages.
The remaining declared numerical comparison uses the eight-element ramp with
pointwise scaling for all four fresh centre trajectories, in
`development-ramp8-pointwise-scalars.json`. Its rough three-step pointwise
pilot already passes; the earlier long eight-element assessments used
field-maximum scaling. Keep all centre physical inputs, canonical equations,
dt=.01, 512-vector/200-iteration capacity and all validity/development gates.
Use four ranks for each CLS case and three for each standard case (14 total),
with one numerical thread per rank. This is one matched resolution/scaling
comparison; it is not a spatial convergence study or a parameter-fit search.

Both refined pointwise standard dt=.01 assessments reject a line search.
`development-ramp8-pointwise-standard-dt005.json` declares a fresh fixed
dt=.005 standard-only control on the identical eight-element mesh, with all
physical inputs, model/wall coefficients, solver capacities and acceptance
gates retained. The ongoing CLS trajectories retain dt=.01. This applies the
previous per-case fixed-timestep provision; it does not change a timestep
inside an accepted trajectory or establish temporal convergence.

Both smaller-step standard controls also reject. Before changing another
control, `pilot-ramp8-standard-debug.json` repeats only the first smooth
dt=.005 stage, with the same three-rank allocation and RANSDebug=1. Retain
recomputed GMRES action residuals and the physical state/continuity evolution
to distinguish linear capacity/accuracy from the later line-search rejection.
This logged replay cannot count as developed flow or calibration.

The replay stops GMRES before its capacity but reports recomputed action
residuals up to 3.48, including 3.02 before line-search rejection; weak
divergence remains 5.62e-6, above the physical gate. Standard uses forward
differences, whereas CLS uses centered differences. Before adding another
solver formulation, `pilot-ramp8-standard-fd-accuracy.json` repeats the same
first stage with the standard default difference step 1e-7 and native inverse
tolerance 1e-13. Retain dt=.005, three ranks and debug output. This controls
forward-difference truncation and inverse noise without changing an equation
or the physical acceptance thresholds. Both settings are bound explicitly.

That accuracy control accepts its first stage, with recomputed action residuals
.00854--.01614 instead of order-one discrepancies. Declare
`development-ramp8-standard-fd-accuracy.json`: fresh standard S/R dt=.01
trajectories on the same eight-element mesh, with difference step 1e-7 and
native inverse tolerance 1e-13. Retain pointwise scaling, current residual map,
canonical physics and all independent gates. CLS keeps its centered 1e-4
probe and existing trajectory. This is a model-specific numerical accuracy
setting, not a temporal convergence result or proof of long-time feasibility.

The refined pointwise rough-CLS trajectory rejects its t=.69 line search
after 68 accepted stages. Before attributing this rejection, replay exactly
one step from its native t=.68 field and complete IMEX history, with the same
four ranks, centered difference step, inverse tolerances and all physical
controls. Enable only RANSDebug, archive the recomputed action residuals and
use a separate `conditioning-ramp8-pointwise-cls-debug-v1` directory. This
is a conditioning diagnostic, not a developed or calibrated trajectory.

The one-step native replay accepts t=.69, but its wrapper audit fails because
the shared map audit requires three native stages for the independent BDF2
momentum inventory. Retain that wrapper failure and native exit-zero receipt
distinctly. Add an up-front three-step minimum to the continuation helper,
then repeat from the same t=.68 checkpoint for three steps in
`conditioning-ramp8-pointwise-cls-debug-v2`, retaining the same four ranks,
settings and 1800 s timeout. This tests whether refreshing numerical caches
at restart permits evolution beyond the rejected step; it is not a developed
candidate or evidence that a unique root exists.

The three-step replay accepts t=.69 and .70. At the following step its
recomputed action residual rises above .05 while GMRES stops before 512
vectors. Declare a matched three-step accuracy control from the same native
t=.68 checkpoint, `conditioning-ramp8-pointwise-cls-inverse-v1`: change only
the native inverse tolerance from 1e-8 to 1e-13, retaining centered 1e-4
differences, four MPI ranks on CPUs 69--72, debug logging and a 3600 s timeout.
Run it after the existing replay releases those CPUs. The continuation helper
records the explicit numerical override in its receipt and actual XML/manifest;
the complete physical IMEX history and configuration signature are retained.
This diagnostic does not supply a development or calibration flag.

At t=.71 the tighter-inverse replay still reports recomputed action errors
above .03 before the 512-vector limit. Before another long assessment,
declare `development-ramp8-cls-fd-accuracy.json` for a fresh rough-CLS
trajectory: use the CLS default centered difference step 1e-6 and native
inverse tolerance 1e-13, retaining the same eight-ramp-element mesh, physical
centre, fixed dt=.01, pointwise metric and all physical gates. This jointly
controls central-difference truncation and inverse noise; no causal attribution
to either setting alone is made. Do not change the difference step inside an
existing full-history trajectory. Use four ranks on CPUs 69--72 after the
conditioning control finishes, checkpoint every step, twelve extensions and
a 43200 s per-leg timeout so an expensive valid two-flow leg is not stopped
by the earlier four-hour operational limit. No other physical input probe is
introduced. The ongoing smooth CLS and standard trajectories retain their
declared settings. If the conditioning replay stalls with an action residual
above one at capacity, retire it explicitly rather than treating its prefix
as a physical failure or a developed candidate.

After the fresh accuracy trajectory accepts t=.03, use the existing matched
checkpoint auditor to compare its complete native field with the earlier
pointwise refined-ramp trajectory at t=.03 (`case_3.chk` in each producing
stage). Use one audit worker on CPU 79 and `fd-accuracy-checkpoint-v1`, require
native time .03, identical physical/mesh/dt data, unchanged producing hashes,
positive scalars and the existing 1e-6 complete-field comparison gate. The
baseline's later rejection remains a failure; this early accepted-checkpoint
comparison supplies neither a developed-flow nor a convergence claim.

Observed current rough-CLS step costs are about 70--100 s on four iterative
inverse ranks. The archived compact direct-inverse rough control used a
median 19.0 s on one rank, but has a different mesh/scaling and later rejects;
it is not a matched speed or feasibility result. Declare the unscored
`pilot-ramp8-cls-direct-accuracy.json`: clone the fresh refined accuracy recipe,
change only GlobalSysSoln to DirectMultiLevelStaticCond, run three dt=.01
steps with one MPI process/one thread on spare CPU 79 and timeout 1800 s.
Retain pointwise scaling, centered 1e-6 probes, canonical physics and all
independent equation/strip/coverage checks. Compare its complete native t=.03
field with the ongoing four-rank accuracy checkpoint using the existing
1e-6 gate and time/producer binding. Current trajectories continue untouched.
Any observed timing ratio remains an operational, allocation-dependent
throughput result, rather than an isolated performance benchmark.

Local source inspection confirms that MPI builds compile and link the
existing XxtMultiLevelStaticCond backend (ThirdPartyMPI and MultiRegions
CMake definitions). Declare `pilot-ramp8-cls-xxt-accuracy.json`, identical to
the direct-accuracy pilot except for that backend. Use one process on CPU 64,
three steps and timeout 1800 s, after the periodic Linux monitor exits. The
four trajectories use fourteen workers and the two single-process pilots
use the remaining two; no independent Linux audit or compiler is started
concurrently. Each pilot's sequential audits reuse its own CPU. Monitoring
is from the Mac. Require the same native physical/coverage/strip gates and
time-bound t=.03 complete-field agreement before considering any long-run
recipe change. This does not introduce a new closure or wall treatment.

Both single-process pilots pass their independent physical audits and native
field comparisons. Declare `pilot-ramp8-cls-xxt-mpi2-accuracy.json`: retain
the XXT recipe and run three steps on two ranks restricted to CPUs 64 and 79,
after the CPU-79 periodic monitor exits. The four ongoing trajectories plus
this pilot use exactly sixteen compute workers, with one numerical thread
per rank and no independent audit/compiler allocation. Sequential post-run
audits use one of the released pilot CPUs. Require the same full-field and
physical gates. Compare the sum of the first three reported step timers with
the existing four-rank iterative trajectory; these timers use steady-clock
elapsed time and exclude subsequent checkpoint I/O. Record native Execute
and whole-run times separately. A two-rank timing is measured directly; it
does not establish four-rank XXT scaling or long-time feasibility.

Declare a read-only constitutive diagnostic using the canonical kernel linked
by the current solver. Use immutable initial/last-accepted fields from the
failed refined rough trajectory and initial/first-flow fields from the refined
smooth trajectory. Reconstruct both velocity gradients on each affine Plane
quadrilateral from its native GLL data; verify polynomial differentiation and
independent degree-four Legendre reconstruction. Sample the four interior GLL
points nearest each element centre plus H/N and inlet traces, retaining
one-sided copies. Analyze actual gradients and a separately labelled frozen
trace-free projection; no evolved field is changed.

Evaluate T=nu_eff(G+G^T)-R_NL. At fixed k/epsilon minimize
`a_i n_j [T(G+s*h*a*n)-T(G)]_ij/(s*h)` over 361 in-plane directions,
`a.n=0`, s=+/-1, and h=1e-5/1e-6 times max(1,max|G|). Include the actual
molecular viscosity. Require standard-model agreement with nu+0.09*k^2/epsilon
within 1e-6 relative error, and agreement with the existing parallel-gradient
oracle. Report robust negatives only when both widths are below -1e-12 and
agree within 2% of their magnitude. Bind complete native fields, sources,
gradients, compiler and output hashes. A negative local tangent is a
constitutive warning, not proof of the global instability mechanism or a
development pass. Use one compiler/audit worker on CPU 79 within the existing
sixteen-worker allocation; no CFD equation, closure, wall condition or gate
is changed by this diagnostic.

The refined accuracy rough centre rejects t=.69 after accepting through .68,
reproducing the earlier last state despite smaller probes and a tighter
inverse. The refined smooth 1e-4/1e-8 trajectory rejects t=3.15 after accepting
through 3.14. Rough standard develops after two passing drift extensions;
smooth standard continues. These CLS rejections leave R1 open.

Before further runs, declare two remaining comparisons with the accurate
recipe. Run a fresh smooth CLS centre with centered 1e-6 probes and native
inverse 1e-13, keeping the same mesh and physical inputs. Repeat the already
declared rough k_scale=.75/epsilon_scale=2 corner on the eight-element ramp
with accurate CLS settings and a matched standard counterpart (forward
1e-7/native inverse 1e-13). This revisits the same physical corner previously
tested with the compact mesh and coarser numerical accuracy; it introduces
no new physical parameter point or fit score. All three use fresh dt=.01
trajectories, checkpoint interval one, twelve extensions and 43200 s per-leg
timeout. Canonical coefficients, physical wall/Yap distance, integrated strip
continuity, mask and every validity/development gate are retained.

Use CPUs 65--68/four ranks for the fresh smooth CLS, 69--72/four ranks for
the rough CLS corner and 76--78/three ranks for its standard counterpart.
The ongoing smooth standard retains 73--75/three ranks. Thus fourteen native
workers plus the diagnostic/compiler worker on 79 and supervision on 64
remain within the sixteen-worker limit. No failed or diagnostic prefix is
promoted or spliced into these fresh assessments. The diagnostic producer's
first CSV column-selection error is retained in its failed v1 directory;
the corrected producer runs in a fresh v2 directory.

The corrected diagnostic finds robust negative CLS transverse tangents in
initial and accepted sampled fields, including after the labelled trace-free
projection; the standard controls remain positive. Before interpreting the
contrast with the passing benchmark, extend the parallel control to the
canonical compatible Q returned by RANSNeutralReference::CLS at zero resolved
nu. Also sample immutable initial/last-accepted fields of the previously
declared compact rough corner and the last accepted refined smooth/accurate
rough fields. Select checkpoints from native output logs and verify their
native times. Retain v2 and write this extension to fresh v3 output. This is a
constitutive comparison, not a new physical search point, a global stability
proof or a change to any running trajectory.

At the first completed accurate smooth-CLS flow-through, compare its immutable
t=1.75 field with the earlier refined smooth-CLS first-flow field using the
existing time-bound complete-field comparator and 1e-6 primary/derived gate.
Write a fresh `smooth-first-flow-accuracy-checkpoint-v1` receipt before using
the result, with one postprocessing worker on CPU 64 within the sixteen-worker
budget. This checks accumulated numerical agreement at the same physical
point and time; it introduces no new trajectory, fit score or development pass.

The accurate refined rough CLS corner rejects the t=1.07 line search after
accepting through t=1.06, before its first flow-through. The declared rough
bracket therefore supplies no qualifying common pair. Retain the rejection
and all accepted-prefix audits. Its unselected standard counterpart may be
retired after the current three-flow-through stage passes its independent
audits; retain any earlier native rejection under its own status. Record
process ownership, private output path, binary, CPU allocation, last accepted
state and explicit retirement reason before stopping the owned MPI launcher.
That retirement is not a numerical failure or a developed result. Continue
the accurate smooth CLS assessment to its declared native outcome. R1 remains
open, with canonical equations, strip continuity and all gates retained.

If the accurate smooth trajectory reaches the earlier coarse-accuracy last
accepted time t=3.14, compare those immutable native checkpoints with the same
time/source-bound comparator and 1e-6 gate, in fresh
`smooth-last-accuracy-checkpoint-v1` output. The comparison time is fixed by
the already recorded baseline endpoint, before examining field differences.
Use one worker on CPU 64. This tests late numerical agreement, not developed
feasibility or a new physical parameter point; retain either comparison outcome.

The endpoint comparison passes. The accurate smooth trajectory subsequently
exhausts its 200-iteration nonlinear limit at t=3.15 after t=3.14, matching the
coarse-accuracy endpoint. Both new rough/smooth CLS failures have native exit
one and no timeout. The unselected rough standard counterpart is retired after
its passed physical/coverage/strip audit at 3.011757 flows, with a partial
extension retained through t=4.81; its one-to-three-flow drift fails and it is
not developed. All R1 runs have now ended or been explicitly retired. R1 remains
open because the four required developed cases are not available; no selection,
freeze or reserved score is unlocked.

## R2/R3 calibration and freeze contract

After R1, preregister the v2 bounded design before parameter selection. Keep
the five equally weighted flat families and frozen height weights on z>=N;
retain excluded/missing observations and common/conditional fits. Parameter
bounds may be revised from permitted flat-data evidence with recorded rationale
before new selection; the v1 upper-bound minima do not establish an optimum.
Publish k/epsilon/span-variance assumptions, d=2.5/3.5/4.5 mm scenarios, Re_H,
friction, thickness-criterion sensitivity, intensity and profile residuals.

Select common physical inlet/surface parameters across the two primary
composite closures separately for S and R. Freeze four validated primary
policy manifests with mask, geometry/reference convention, source/data hashes,
native development/statistics provenance and numerical limitations. Surface
summary manifests must enumerate validated policies and pending/failed
comparisons. Only a validated policy-specific manifest may unlock its reserved
flat or hill scores. Original and compatible-log transfers require their own
native evidence; they never inherit an EP calibration-complete flag.

Run reserved Fig.3d crosschecks only after this freeze. Hill scores remain
outside CV3 selection. All-original and standard-log transfers use the frozen
physical inputs; retain any incident mismatch and explicit failed/unexecuted
status instead of retuning or blocking the primary EP route. Switched and
selected resolved-viscosity sensitivities follow the fixed-mask study.

## Resources and retained requirements

Use the private Linux source/build and large-volume run root specified by
the project. Preserve the unrelated main checkout/build/installation. Never
exceed 16 project CPU workers, accounting for compilation, MPI ranks, test
concurrency and active numerical threads together. Set OMP/OPENBLAS/MKL/BLIS
threads to one for MPI runs, and record actual concurrent allocations, commands,
timeouts, exit status, source/binary hashes and runtime. Do not infer runtime
from the short WG4 pilots.

CV1 data/scoring definitions, the canonical equations, physical wall distance,
integrated strip continuity, strict accepted-state positivity and residual
gates remain in force. The revised route still has to demonstrate successful
evolution through the original-AWF downstream section. The sequence amendment
alone supplies no new numerical or experimental validation.
