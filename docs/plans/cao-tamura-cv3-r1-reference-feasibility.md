# CV3-R1: reference recovery and composite-floor feasibility

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: native computational verification and feasibility
- Origin Date: 2026-11-01
- Verification Status: reference/native checks pass; both standard centres developed; CLS feasibility assessment IN PROGRESS
- Version Label: cao-tamura-cv3-r1-v1

**CV3-R1 remains open while the four developed composite precursors are being
assessed.** Passing short runs does not close this task. CV3-R2 selection,
freezing and reserved flat/hill scores remain locked. The acceptance contract
is [protocol v2](cao-tamura-cv3-calibration-protocol.md).

The current immutable `c9c6` solver has the following unfitted results on the
eight-element ramp (187 quadrilaterals, 17 actual normal elements):

| Case | Current result | R1 implication |
|---|---|---|
| S standard centre | Developed at t=8.75 s, 5.01225 flow-throughs; two successive drift passes | Valid individual feasibility control |
| R standard centre | Developed at t=9.11 s, about 6.03 flow-throughs; two successive drift passes | Valid individual feasibility control |
| S CLS centre, centered 1e-4/native inverse 1e-8 | Accepts through t=3.14; next stage exhausts its nonlinear limit | Not developed |
| R CLS centre, centered 1e-6/native inverse 1e-13 | Accepts through t=.68; next line search rejects, reproducing the coarser-accuracy rejection | Not developed |
| S CLS centre, centered 1e-6/native inverse 1e-13 | Fresh assessment running | No development claim yet |
| R lower-k/upper-epsilon corner, accurate standard/CLS | Fresh matched assessments running, k_scale=.75/epsilon_scale=2 | No common developed pair yet |

The standard cases use forward 1e-7 probes and native inverse target 1e-13.
All listed cases retain dt=.01, common physical inputs per surface, canonical
coefficients and integrated strip continuity. The tighter CLS settings change
numerical accuracy only. No failed prefix is spliced into a fresh trajectory.
The completed standard centres are independently rechecked against every
declared stage, coverage and strip-budget gate; they do not certify CLS or
close R1 by themselves.

## Reference controls

The current numerical solver reproduces the documented degree-four CLS
inlet/outlet benchmark for 415 unforced steps, dt=0.5 s, t=207.5 s (5.00873
flow-throughs). It uses the benchmark's canonical CLS Q/kappa, H/P/N=2/1/4 m,
compatible outer fluxes and prescribed outlet traction. These physical data
are distinct from the experimental measured-profile inlet and Plane-wall
composite route. No field resets or reference-residual forcing are used.

| Control | Outcome | Quantitative check |
|---|---|---|
| Unperturbed HighRe, fixed EP | PASS | Largest maximum-over-time L2/Linf preservation error over U, k, epsilon, viscosity and modeled normal/shear stresses: 0.00487844%, below 0.1% |
| Selected smooth perturbation, switched EP | PASS | k L2 error decays by 1531.16; final relevant-field L2 errors below 0.1% |
| Unperturbed Molecular, fixed EP | PASS | Largest maximum-over-time L2/Linf error: 0.0204559%, below 0.1%; this is a finite-error result |

All three controls retain positive k/epsilon, zero physical guard activations
and no admissibility fallback. Native chi is one throughout, including the
switched small-perturbation control. Maximum weak divergence is below 1e-8.
The molecular case does not establish an exact finite-viscosity logarithmic
equilibrium. This is a reproduction of the accepted resolution, rather than
a new reference convergence study.

## Implementation and native checks

Experimental cases now explicitly distinguish `EPTransition`, `UniformEP`,
`Original` and `CompatibleLog`. The producing XML, manifest, observation and
statistics receipts carry the route, and incompatible declarations or mixed
routes are rejected. R1 accepts only fixed weight one with the prescribed C2
mask, active correction at every requested-positive accepted trace and no
admissibility fallback. Prescribed downstream zero is recorded separately.
Every accepted time must sample full EP, ramp and downstream regions; native
trace counts retain one-sided duplicates and are not area fractions.
The independent coverage checker also compares wall-frame times with every
complete finite accepted solver row. A wholly missing wall frame fails;
extra wall output from an aborted, unaccepted stage is reported separately.
The missing-frame guard and eight completed native control prefixes pass,
and all 16 registered tests pass with the stricter checker. Exact producer
sources loaded before that audit change are archived separately.

The initial epsilon collar previously used full EP epsilon even downstream
of the ramp. It now uses the same spatial linear mixture of Original and EP
wall epsilon as the native AWF. The high-precision preparer evaluates scalar
collars at physical wall-foot x; the initial floor Dirichlet trace uses the
same expression. The solenoidal velocity initialization and physical inlet
profiles at/above N are retained. The change affects initialization, not the
closure equations, time-evolving wall kernel or integrated strip continuity.

The isolated Release Linux build passes all 16 registered CTests. Native
checks also pass ten surface/policy inlet preparations, independent expression
and positive-projection checks, Q1/Q2 initial residual checks, startup fade,
post-startup 1-to-4-rank history equivalence and rejection of mid-startup
restart. The first validated residual-map solver binary is
`24bdf3f7f1ba26c0639d89bf7b31ce86276b9d738d4af80ad8ebe3fbded8735f`.
The initializer/auditor is
`d71d5658efb45ee17d783ef1497c8c438ae5d39f81569b9d4776c20ede52b50d`.
The reference controls and native checks were repeated on this binary.
An optional pointwise scalar-scaling candidate has subsequently compiled in
the same isolated Release tree and passed all 16 registered CTests. Its binary
is `c9c640b86a5c4c414762e27c2cce57802f4b2be2e7263b81f68e3e17c575dfbb`;
all three analytical controls and native inlet/startup guards also pass on
this binary. The four S/R standard/CLS matched-field checks pass, with maximum
primary-field difference 2.64e-9 and maximum derived-field difference 1.87e-8
against field-maximum scaling. The new scaling is slower in those short
controls (baseline/candidate Execute ratios .528--.725); no speedup is claimed.
Its purpose is to assess scalar conditioning over large heightwise ranges.
Its full-history 1-to-4-rank check passes (maximum derived-field difference
2.23e-8), including explicit map and scalar-scaling mismatch rejection. The
complete suite also passes after adding pointwise scaling to the registered
CLS smoke test. Timing comparisons include allocation/load differences;
development robustness is assessed separately. The immutable first binary
and its compiled-source snapshot are retained.

The delivered history helper also repeats its full-history 1-to-4-rank and
map-mismatch check successfully in `residual-map-history-v4`. A separate
refined-ramp checkpoint at t=.04 matches four and eight MPI ranks: maximum
primary-field scaled difference is 4.30e-10 and maximum derived-field
difference 1.19e-8, below 1e-6. That checkpoint comparison certifies neither
long-time development nor spatial/temporal convergence.
The current pointwise binary also passes a separately time-bound lower-corner
four-to-eight-rank checkpoint at t=.05: maximum primary difference 3.17e-10
and derived difference 7.01e-9. Older checkpoint evidence is retained under
its producing binary and cannot certify a different selected binary.

## Experimental development assessment

The four centre physical inputs are unfitted: delta_in=.21 m, k/epsilon
multipliers one, smooth ks=0 and rough ks=.004 m. Both closures retain common
physical inputs above N with MeasuredRms/ShearLength inference, while their
canonical coefficients and required H..N policy collars differ. The mask
runs from -.28 to -.20 m and ends before the proposed hill foot at -.10 m.
The actual downstream Original section and x=0 incident station are retained.

The first numerical recipe failed Newton line searches. An independent action
check showed that the default finite-difference/inner-solve combination could
produce a substantially larger action residual than the linear solver's
estimate. A 5e-14 current-operator target and physical difference step 1e-4
allowed all four centre cases to complete three dt=.01 timesteps on the full
ramp-resolved pilot mesh. Each retained positive fields and the complete
prescribed correction envelope. These short runs establish neither developed
flow nor experimental agreement.

The optional [preconditioned-residual Newton map](cao-tamura-cv3-r1-newton-map.md)
then reproduced all four accepted solutions within 1e-6 scaled error, including
full primary fields and derived stresses. The largest primary-field difference
was 1.74e-8 and the largest derived-field difference 1.97e-7. All independent
physical residual, positivity, correction-coverage and strip-budget checks
passed. The identical native Execute timer measured speedups of 5.74/5.88 for
S standard/CLS and 7.87/7.09 for R standard/CLS. A five-step compact-mesh
comparison against the retained CurrentSolve trajectory also passed. This
changes the nonlinear preconditioner, retaining the current assembled equations
and default CurrentSolve path. It supplies no mesh/timestep convergence claim.
Both sides used one MPI rank and one numerical thread. Baselines were pinned
to CPUs 1--4 on socket 0 and new-map cases to CPU 64 on socket 1 while an
unrelated 48-rank study used socket 0. These observed runtime ratios therefore
include allocation/load effects and are not an isolated performance benchmark
of the map. The field and independent-equation comparisons are the numerical
equivalence evidence.

Larger-step, MPI, relaxation, augmentation and startup controls retain their
own receipts and outcomes. Exploratory runs stopped to release resources are
marked as stopped, and cannot count as developed solutions. The numerical and
unmeasured-input feasibility brackets were recorded before their native
probes; none used fit scores or reserved observations.

The staged primary recipe uses p4 velocity/scalars, p2 pressure, eight GLL
points, twelve normal layers, three outer blend cells, four 20 mm ramp
elements and compact streamwise nodes [-1,-.28,-.26,-.24,-.22,-.2,0,1.2] m.
The twelve-layer generator input plus inserted strip/blend endpoints yields
17 actual normal elements (119 quads on the seven-column compact mesh).
It uses Molecular diffusion, implicit current convection, dt=.01, fixed EP,
Newton--Krylov with PreconditionedResidual, grad-div=.01, normal augmentation 2 and Picard-map/pressure
relaxation .1. Native inverse and pressure inverse targets are 1e-8 and 1e-13.
It starts unforced and retains full history every five steps. Stages cover
one flow-through, continuation to three, and one-flow-through extensions;
two drift passes after three unforced flows are mandatory. All free/weak,
projected-normal and mass gates remain 1e-7, with independent strip inventories.

The initial CurrentSolve rough-CLS assessment was deliberately stopped after
48 valid steps, t=.48 s, to adopt the faster verified map. Its full histories,
native residuals, producer sources and explicit stop receipt remain archived;
it is neither developed nor evidence of physical infeasibility. New-map
development uses fresh trajectories with one numerical thread per MPI rank.
A maximum of 16 project CPUs is retained, including the control/diagnostic
allocation. No hill solve or calibration freeze is implied by this
work-in-progress record.

The new-map rough-CLS centre completed 48 valid steps before the t=.49
stage exhausted 60 nonlinear iterations. An exact full-history continuation
from t=.45 reproduced the rejection at .49 after accepting .46--.48; refreshing
the native caches through restart did not remove it. Before rejection the
fields stayed positive, guards/fallback were zero and maximum weak divergence
was 9.18e-10. Strong divergence grew from .00186 to .0545 and maximum pointwise
normal mismatch reached .1546 U_H at x=-.24 m in the ramp. Small projected
residuals do not certify these strong errors. The failure is retained, and
does not establish global nonexistence or a causal instability mechanism.
At the rejected .49 stage, a diagnostic run reaches a far-above-tolerance
plateau: epsilon map change .00623, weak divergence .000330 and projected wall
error 4.21e-5. Its late 240-vector GMRES action residual is about .065--.084,
above the requested .01, and a linear epsilon image is negative. This is a
nonlinear-stage rejection, not a bounded physical oscillation that can enter
statistics. It does not prove that no other positive solution exists.
The declared epsilon_scale=2 feasibility point uses the same fixed .01 step,
mesh, mask, canonical coefficients and physical gates. Its standard counterpart
develops; CLS rejects after t=1.03. Any accepted rough starting pair must use
the same physical input in both closures.
Experimental times and velocities here use the declared native U_H=1 scaling
at the source Re_H; they are not a determination of the tunnel's dimensional
velocity or physical duration. The analytical benchmark has its own recorded
dimensional reference data.

Current development outcomes are:

| Surface/closure/control | Status |
|---|---|
| S standard, centre | Developed at 8.0196 flow-throughs; final two one-flow drift gates pass |
| R standard, centre | Developed at 10.0546 flow-throughs; final two one-flow drift gates pass |
| S CLS, centre | First one-flow audit passes; subsequent evolution rejects after t=3.13 |
| R CLS, epsilon_scale=2 | Rejects after valid evolution through t=1.03 |
| S CLS, grad-div=.1 | Rejects after t=3.67 |
| R CLS, grad-div=.1 | Rejects after t=.93; last strong divergence .00945 |
| R CLS, original full x mesh | Rejects after t=.54; last strong divergence .1294 |
| S CLS, original full x mesh | Rejects after t=2.70 |
| R CLS, dt=.05 | Rejects after eight valid steps, t=.40 |
| R standard, epsilon_scale=2 | Developed at 9.049 flow-throughs; cannot validate the failed CLS counterpart |
| R CLS, k_scale=1.5 and epsilon_scale=2 | Rejects after t=.36 |
| R CLS, one-flow startup fade | Rejects after t=1.47, before the fade ends; no unforced development evidence |
| R CLS, p3 pressure | Rejects after t=.95 |
| S CLS, p3 pressure | Rejects after t=3.15 |
| S standard, p3 pressure | Rejects after t=.01 |
| R CLS, 512-vector GMRES | Rejects at the same next stage after t=.48 |
| R CLS, eight ramp elements, 240-vector/60-iteration capacity | First long-run stage exhausts the nonlinear limit; a separate serial debug control accepts its first step after 55 Newton iterations |
| R CLS, eight ramp elements, 512-vector/200-iteration capacity | Four-step pilot passes; unselected long assessment deliberately retired after t=.68 |
| S CLS, eight ramp elements, 512-vector/200-iteration capacity | Unselected long assessment deliberately retired after t=.88 |
| R/S standard, eight ramp elements, 512-vector/200-iteration capacity | Line searches reject after t=.27/.59 |
| R CLS, refined ramp, dt=.05 | Rejects before the first accepted step |
| R CLS, uniform EP, matched compact inputs/dt=.01 | Rejects after t=.96; removing the ramp/downstream Original section does not by itself remove rejection |
| R CLS, serial direct cached inverse | Rejects at the same next stage after t=.48 |
| S CLS, serial direct cached inverse | Deliberately retired at t=2.50 after a valid first one-flow audit; not a numerical failure or developed case |
| R CLS, pointwise scaling, compact centre, 512-vector/200-iteration capacity | Rejects the t=.49 line search after 48 accepted steps; changing scalar weighting does not rescue this centre |
| S CLS, pointwise scaling, compact centre, 512-vector/200-iteration capacity | Rejects the t=3.14 line search after valid evolution through t=3.13 |
| R CLS, pointwise scaling, refined ramp | Three-step serial pilot passes in 598.99 s; not developed |
| R CLS, pointwise scaling, refined ramp, centered 1e-4/native inverse 1e-8 | Fresh long assessment rejects t=.69 after 68 accepted steps through t=.68; not a timeout |
| R CLS, pointwise scaling, lower-k/upper-epsilon corner | Rejects the t=1.09 line search after 108 accepted steps, t=1.08; no first one-flow audit |

After the repeated rough-centre rejection, the final declared input corner
uses k_scale=.75 and epsilon_scale=2, within the original uncertainty bounds.
Both closures use this same physical corner on the rough surface, retaining
delta=.21, ks=.004, the compact mesh, mask and canonical equations. A five-step
four-rank CLS pilot passes, but its fresh eight-rank trajectory rejects at
1.09 (861.04 s, no timeout). Up to t=1.08 its maximum weak divergence is
9.96e-10, projected normal error 6.75e-10, and minimum k/epsilon remain
2.68e-5/9.14e-7 with no physical guard or fallback. Strong divergence reaches
.07614 and maximum strong normal mismatch .2162 U_H. These are separate
resolution diagnostics, not removed by the small weak errors. The unselected
standard counterpart is retired at t=3.38 after its valid first one-flow audit.
The unselected current-binary smooth standard centre is retired at t=5.19
after its first audit; neither retirement is a numerical failure or a developed
case. The two current-binary compact CLS rejections do not supply a common
developed pair. No fit score is used to choose these controls.

The remaining matched numerical comparison is fresh S/R standard/CLS centre
development with the eight-element ramp and pointwise scaling. It retains the
same physical data, canonical coefficients, dt=.01 and 512-vector/200-iteration
capacity. Four ranks per CLS and three per standard case use 14 workers total.
The prior long eight-element runs used field-maximum scaling, so their outcomes
do not execute this combination. All four are assessed by the same unforced
duration, independent equation/strip, exact correction-coverage and drift gates.

With the original common 1e-4 probe, refined pointwise standard rejects after
t=.05 (S) and .35 (R). Halving dt to .005 rejects before the first S step and
after t=.015 for R. A logged replay of the first S stage finds an order-one
linear-action inconsistency: recomputed GMRES relative residuals reach 3.48,
with 3.02 after 160 vectors immediately before rejection, while weak divergence
is still 5.62e-6. This is above the physical gate, rather than a reason to
relax that gate. It does not identify all prior CLS failures.

The standard model uses forward differences; CLS uses centered differences.
Restoring the standard default 1e-7 probe and tightening the cached native
inverse to 1e-13 accepts that first stage, with recomputed action residuals
.00854--.01614. Fresh standard dt=.01 assessments use those accuracy settings;
the two CLS trajectories retain their centered 1e-4 probe and 1e-8 inverse.
The standard/cubic equations and all physical acceptance gates are retained.
Both refined standard centre trajectories subsequently meet the declared
development gates, as recorded above. This does not establish spatial/temporal
convergence or rescue the matched CLS centres.

The eight-element smooth CLS trajectory completes its first full-flow audit
at t=1.75 (1.00245 flow-throughs), with maximum free-row residual 1.27e-10;
the continuation later exhausts its nonlinear limit at t=3.15 after accepting
through t=3.14. The rough counterpart
rejects its t=.69 line search after 68 accepted steps (3663.68 s measured run
time, no timeout). Its accepted prefix has maximum weak divergence 9.57e-10,
zero physical guards, and positive k/epsilon. Last accepted strong divergence
is .15379 and maximum pointwise normal mismatch is .645704 U_H. This does not
establish developed rough flow, spatial accuracy or physical nonexistence.

An exact-history diagnostic from t=.68 accepts t=.69 and .70 after restart.
The first, one-step diagnostic has a native exit-zero acceptance but a wrapper
audit failure: the shared independent BDF2 inventory requires three native
frames. That wrapper error is retained separately, and the helper now rejects
such a request before launching a solver. In the three-step replay, action
residuals grow above four at the 512-vector capacity while t=.71 stalls.
A matched 1e-13 native-inverse replay accepts the same first two steps but its
t=.71 action residual also grows above one at capacity. Both diagnostics are
explicitly retired with verified MPI ownership and native-prefix receipts;
their termination is not recorded as a numerical or physical failure.

A further rough assessment starts fresh, with the same physical centre and
mesh but the CLS default centered probe 1e-6 and native inverse target 1e-13.
No difference step is changed inside an existing history. This jointly tests
finite-difference truncation and inverse accuracy, retaining fixed dt=.01 and
every physical/coverage/development gate. It later rejects the same next t=.69
step after t=.68. No part of a diagnostic is promoted to development.
At the first common three-step checkpoint, this fresh rough trajectory agrees
with the earlier accepted solution: maximum primary-field scaled difference
1.99e-9 and derived-field difference 9.56e-8, below the existing 1e-6 gate.
Both native field times are .03; producing case, binary, geometry, manifest,
complete-field and native observation hashes are bound. This comparison does
not reclassify the baseline's later rejection or establish developed flow.

Two serial backend pilots on that same refined accuracy case also pass their
independent free/continuity/positivity, strip-budget and coverage audits. Their
t=.03 field comparisons have maximum scaled differences 5.86e-8 (direct) and
6.44e-8 (XXT), below 1e-6. The native integration timer uses a steady clock:

| Native inverse backend | MPI ranks | First three integration steps (s) |
|---|---:|---:|
| IterativeMultiLevelStaticCond | 4 | 396.1263 |
| DirectMultiLevelStaticCond | 1 | 1084.512 |
| XxtMultiLevelStaticCond | 1 | 949.879 |
| XxtMultiLevelStaticCond | 2 | 594.529 |

These measured timers exclude subsequent checkpoint I/O. The allocations and
concurrent loads differ; neither serial pilot demonstrates a wall-time gain
over the existing four-rank case. They do not establish parallel XXT scaling
or developed feasibility. The separate two-rank XXT comparison also passes
the independent audits and field gate (maximum scaled difference 5.80e-8).
It uses the two CPUs released by the periodic monitor while the four
trajectories continue, retaining the sixteen-worker project limit. No tested
alternative allocation is faster than the existing four-rank iterative case;
the current long trajectories are retained.

A separate replay of the original CurrentSolve full history at t=.45 accepts
.46--.48 and also exhausts its nonlinear limit at .49 (988.21 s, no timeout).
Thus the corresponding compact rough rejection is reproduced with both
Newton maps, rather than established as a new-map-only regression. Tightening
the cached inverse from 1e-8 to 1e-12 on the refined ramp accepts three steps
in 1339.15 s; its 240-vector linear solves still reach their iteration capacity.
The 512-vector/200-iteration recipe accepts four steps in 390.75 s on its
recorded allocation. Neither short control qualifies as developed flow.

The registered numerical controls retain integrated strip continuity and the
canonical stresses. Their status does not change the requirement for four
valid developed cases. No rejected timestep enters an average or a drift pass.

An independent last-accepted-frame trace decomposition locates the largest
normal errors at element nodes: x=-.24 for the original rough centre, and
x=-.20 for the pointwise rough lower corner and smooth centre. Their one-sided
normal-target jumps are .3037/.2349/.5762 U_H, while the normal-velocity copies
agree within 7e-18. The target differentiates a continuous, piecewise-polynomial
eliminated U_P, whose derivative need not be continuous across elements. The
constraint uses the full continuous velocity trace, rather than the lower-order
pressure trace. Maximum interior errors remain .0709/.0814/.2373 U_H, so the
discrepancy cannot be dismissed as a node-only effect. These accepted-frame
observations expose a resolution concern; they do not establish the cause of
Newton rejection or authorize an impermeable-interface replacement.

For the earlier `24bd` compact-mesh developed standard centres, maximum free-row residuals are
1.56e-10 (S) and 1.33e-10 (R), and maximum weak divergence is below 1e-9.
The final pointwise relative k/epsilon drifts are 1.85e-4/9.47e-5 for S and
1.60e-4/3.25e-4 for R, below the declared 1e-3 limit; their preceding one-flow
extensions also pass every drift gate. Maximum strong normal mismatch over
the trajectories is .00434/.00750 U_H for S/R. These nonzero strong errors
remain a CV4 convergence concern; the small weak residuals do not remove them.

The current `c9c6` refined standard centres pass independent exit rechecks at
5.012247 (S) and 6.030133 (R) flow-throughs. The R four-flow extension fails
its k/epsilon drift test (.00161017/.00168529); its subsequent five- and
six-flow extensions both pass. This failed drift is retained, rather than
accepting the case after an arbitrary five-flow duration. Final pointwise
relative k/epsilon drifts are .00010602/.00018513 (S) and
.00010767/.00009971 (R). No experimental fit or common CLS pair is selected.

## Read-only constitutive comparison

`validate_cao_tamura_r1_tangent.py` compiles the standalone
`r1_general_stress_tangent.cpp` against the canonical RANS kernel and samples
immutable complete native fields. Tensor-product differentiation is checked
against analytic polynomials and an independent P4 Legendre reconstruction
(maximum relative gradient discrepancy 2.19e-11). The general-gradient probe
matches the existing parallel oracle, and standard-model tangents agree with
nu+0.09*k^2/epsilon within 7.25e-10 relative error. The first CSV producer
error is retained in v1; corrected v2 and extended v3 diagnostics complete.
Neither diagnostic changes the solver or evolved fields.

At frozen positive k/epsilon and gradient G, the probe evaluates the
directional transverse coefficient of
`T(G)=nu_eff*(G+G^T)-R_NL(G)` in directions `a tensor n`, with a.n=0.
It includes resolved molecular nu, 361 in-plane angles and both directional
signs. A reported robust negative requires both relative widths 1e-5/1e-6
to give values below -1e-12 and to agree within 2%. Actual gradients and a
separately labelled frozen trace-free projection are evaluated; the latter
does not project or modify the native solution.

| Native CLS sample | Robust negatives, actual / trace-free | Minimum actual / trace-free tangent (m2/s) |
|---|---:|---:|
| R centre initial, refined ramp | 470/470 of 1145 each | -4.1093e-5 / -4.1093e-5 |
| R centre t=.68, refined ramp, accurate recipe | 492/485 of 1145 each | -6.2652e-4 / -5.4586e-4 |
| S centre initial, refined ramp | 499/499 of 1145 each | -2.0182e-5 / -2.0182e-5 |
| S centre t=3.14, refined ramp | 249/239 of 1145 each | -4.4443e-4 / -4.4421e-4 |
| R corner initial, compact ramp | 58/58 of 777 each | -2.2732e-5 / -2.2732e-5 |
| R corner t=1.08, compact ramp | 136/128 of 777 each | -6.0468e-4 / -3.4497e-4 |

Samples comprise four interior GLL points per element plus H/N/inlet traces,
including one-sided duplicates; these counts are not volume fractions.
Every standard counterfactual evaluated at the same native k/epsilon/G has
positive tangents. The molecular viscosities are 3.33333e-6 (S) and
3.80952e-6 m2/s (R), already included in the table.

The compatible canonical parallel control also has a robust negative
transverse tangent: Q=3.19067508, kappa=.440553905, k=.01, shear=1,
epsilon=k/Q and zero resolved nu give -2.11983e-4 m2/s. This algebraic
control is not the complete benchmark trajectory. In particular, negative
local tangents are not unique to the experimental inlet, and their presence
does not contradict the passing preservation or selected perturbation test.
Those tests do not establish stability to every spatial mode. These frozen
local results identify a constitutive concern; they do not prove the cause
of a global nonlinear rejection, physical nonexistence or developed feasibility.
The accurate smooth and matched rough-corner trajectories therefore retain
their original development and independent physical gates.

## Further numerical hypothesis, not an established cause

The P4 velocity/P2 pressure pair enforces continuity weakly. In a parallel
CLS flow the nonlinear normal stress R_NL,nn is constant for the analytical
constant-k log reference, but generally varies with height for the measured
finite-depth inlet family. Normal momentum then requires
`pi + R_NL,nn = constant`, where `pi = p_mechanical/rho + 2k/3`.
A lower-order pressure space need not represent that nonpolynomial balance.
This is a possible source of parasitic normal motion, rather than proof of
the cause of the recorded rejections. The canonical CLS transverse-tangent
findings in the v1 record remain a separate concern.

A candidate consistent pressure representation is `Pi=pi+phi`, with
`phi=n.R_NL.n`, and `R_NL*=R_NL-phi I`. Then
`nu_eff*(grad U+grad U^T)-R_NL*-Pi I` equals the original total stress
pointwise. Mechanical wall pressure must be reconstructed as
`Pi-phi-2k/3`; integrated strip continuity would be retained. This changes
the finite-dimensional pressure approximation, and is not an equivalence
claim for the existing P2 discrete solution. No such pressure extension has
been implemented or validated in the present R1 solver.

Before using that representation, test genuinely evolving, streamwise-periodic
parallel profiles with nonconstant k/epsilon: v=0 should persist without field
resets, frozen coefficients or reference residual forcing. Bind the pressure
convention to history/observations, verify physical traction and mechanical
wall-pressure reconstruction, and repeat reference, MPI/history and all four
composite development gates. Passing that diagnostic alone would not certify
the experimental route or remove the need for CV4 convergence evidence.

## Reproduction and scope

Source parent is `f04f7f751ef5ee5b1c1309bb9f0e8bf82f30c879` plus the
recorded R1 changes. Use the private Linux source/build under
`/media/jay/新加卷/nektar-v2026-runs/cao-tamura-calibrated-validation`,
leaving the unrelated main checkout/build/installation intact. Native runs
and large artifacts are under `calibration/cv3/ep-transition-v2/r1/`.
The reference helper, `cao_tamura_r1.py` and `record_cao_tamura_r1.py` bind
source/binary/geometry/configuration hashes and native command receipts.

The native Release suite is reproduced with `ctest --test-dir <private-build>
--output-on-failure --parallel 1`, pinned to the declared control CPU. The
reference helper accepts `--case reference|perturbation|molecular --build
<immutable-build> --output <fresh-folder> --cpus <one-CPU>`. Experimental
operations use `cao_tamura_r1.py develop --surface S|R --model
standard_fixed|cls_fixed --build <immutable-build> --output <fresh-folder>
--cpus <explicit-list> --ranks <count> --override <preregistered-recipe>`.
Activate the installed MPI environment first and set OMP/OPENBLAS/MKL/BLIS
threads to one and TMPDIR to the large-volume experiment directory.
Actual commands, CPU lists, native runtimes and exit codes are retained in
the run receipts; no runtime estimate is substituted for a measured solve.

Closure of R1 requires passing reference controls **and all four** developed,
audited composite cases. Numerical convergence, full stationarity/3D readiness,
calibration quality and experiment validation remain later tasks. A failed or
timed-out segment remains a failed or incomplete assessment; it cannot be
averaged into an unsteady-flow acceptance claim.
