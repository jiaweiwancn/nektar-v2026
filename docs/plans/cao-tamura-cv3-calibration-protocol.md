# CV3 protocol v2: EP-AWF transition calibration first

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: computational experiment plan amendment
- Origin Date: 2026-11-01
- Verification Status: sequence change authorized by the user; revised-route execution pending; CV3 remains IN PROGRESS
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
| CV3-R1 — Reference recovery and composite-floor feasibility | Reproduce the documented CLS EP-AWF reference on the current solver; establish valid developed composite-floor candidates for S/R and both closures, with audited active correction, transition and original section | CV1–CV2 and delivered CV3 fixes | NEXT; NOT EXECUTED |
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
