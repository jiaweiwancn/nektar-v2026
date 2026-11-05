# Segregated CLS/EP-AWF Picard R3: sequential reduced directions

**Development status: the original M8 problem remains unresolved. No qualified
physical stage has been accepted; accepted physical steps/time remain 0/0.**
The designated volume was restored and native Linux development resumed.
A later concurrent task switched the shared Mac checkout; subsequent work uses
an isolated Mac worktree and preserves that task's PIMPLE branch. Only the
Picard commits are delivered to the requested remote branch. The
long original-start validation and an extended diagnostic solve are still
running while this report is being assembled; their final results will replace
the interim entries below before publication of the final evidence.

R3 now constructs a constrained k direction from the response of converged
velocity, pressure/gauge, and epsilon/trace solves. It no longer chooses a
conditional Jkk direction and only optionally settles the other equations.
The R2 conditional method, old backends, verification repair, restricted k
analytic assembly, original residual, and historical evidence remain intact.
The R2 projected map could be stationary in its frozen A_k^T Mff^-1 A_k
metric at ground contact while its original k residual stayed nonzero. R3
therefore scores the original reduced residual, not that projected-map defect.
No full all-field Jacobian or hidden coupled direction solver was added.

The current producing solver is `9a5b58d` (radius-constrained reduced models).
The measured-response enrichment is `bdcdc48`.
The inward restoration search target is `ef83196`;
`3cd9dca` supplies the checked half-response warm starts.
The pure-restoration control is `8687dc3f6d7c9cf7dc6e1a1cc3e3f7fd069f84da`.
Its mixed-search and relative-restoration controls are `1c23b84` and `ad1dbc3`;
the completed earlier radius/filter controls use `8dc862f`.
Its initial implementation was `819c135177a5381d9332013609378eb943c8e0b9`;
native-driven repairs include QR handling of degenerate active rows (`380bad3`),
mixed restoration (`7848ad4`), bounded ground arithmetic corrections (`645d05b`),
physical-history-bound trial packets (`e7313f7`), and a filter-consistent direction
constraint (`8dc862f`). The continuation and spatial diagnostic harnesses are
`315cbee` and `0bd00ae`. R1/R2 reports and evidence are preserved.

## Target and acceptance authority

The primary target is unchanged: original outer6 M8 mesh and physical history,
2,165 unknowns, p3/p1 spaces, native GLL8 rule and original source quadrature,
SaturatedR2, capped destruction with Linear smoothing, original transition,
closures, normalization, diffusion, boundary data and zero prescribed forcing.
The physical implicit coefficient is 400, stage time .0025, and original
complete-residual tolerance 2e-12. Loading a numerical trial changes the
nonlinear guess only. It never replaces h0, history slots, stage time or alpha.

The original qualification requires maximum normalized outer-gradient departure
<=.001, local outer k and sampled K ratios in [.99,1.01], actual outer epsilon
and projected matching trace in [.98,1.02], covariance and normalization-lambda
checks, the native scalar evaluation domain, and continuous-ground admissibility.
Finite scalar-inventory validity is not a full-element positivity theorem.
Ground admissibility uses the actual continuous cubic polynomial and its extrema,
not signs of Modified_A coefficients.

The initial represented h0 remains bitwise SHA256
`b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71`.
The native callback's numerically equivalent history (canonicalized zero signs)
is `f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca`.
The harness asserts target size, alpha, time, tolerance, h0 and history binding.
A converged affine model, vanishing projected update or eta-positive root is
never an acceptance certificate.

## Controlled non-k operation

`EquilibrateNonK` fixes every independent k coefficient. It preserves the const
physical history, stage coefficients/time and external data, and uses separate
velocity prediction, pressure Schur correction, and epsilon transport/nonlinear
trace solves. H contains original momentum, both regional continuity equations,
common gauge, outer epsilon and matching trace. Its merit excludes k.

Tolerance and update budget are arguments. The return includes convergence status,
attained block norms, work and failure reason. An unsuccessful operation returns
the caller's incoming state, rather than exposing a partial iterate as Y(k).
Expected factorization, invalid-domain and exhausted-iteration failures permit
retry. Shape and fixed-k violations remain programming errors. No history advances.
Every assessed reduced trial fraction requires H convergence before scoring.
The tighter half-response solve now starts from the midpoint of the converged
non-k response, while retaining the exact prescribed half-step k coefficients.
It must attain the same tighter H tolerance. Failure retries the canonical warm
start. In verification mode, the first used response of every model is also
solved from the canonical warm start; inconsistent local responses are rejected.
This reduces repeated nonlinear equilibration work without substituting a partial
H solve for Y(k). At the nearly restored checkpoint, the first half response
takes one update versus nine canonically, with relative action difference
3.64e-6 and maximum normalized observation-action difference4.94e-6.
The bounded original-h0 control passes the operator/history checks. Its first
14-column update uses162 flow/epsilon pairs versus239 canonically, including
the extra independent warm-start check. The finite-accuracy optimization path
is not bitwise identical: maximum k-coefficient difference is.00310, residuals
64.185405 versus64.208867 and qualification excesses.5230437 versus.5230460.
Both trials are scored using their actual rebuilt residual and inventory.

Verification repeats the same k, changes the y warm start, tightens H, and forces
a zero-budget failure. Representative measured reduced-residual differences are:

| State | Same warm start | Disturbed/tighter warm start |
| --- | ---: | ---: |
| Original h0, after H | 0 | 5.315e-12 |
| R2 endpoint, after H | 0 | 9.770e-8 |
| Ground diagnostic endpoint | 0 | 5.180e-8 |
| Same ground endpoint with H tolerance 2e-11 | 0 | 1.860e-10 |

All k coefficients stayed fixed; the deliberate insufficient-accuracy operation
failed explicitly and rolled back. These differences and two-scale response
checks are empirical local accuracy evidence, not global uniqueness of Y(k).
For example, the six-update mixed original run used 177 converged H evaluations,
1,352 velocity/pressure pairs and 1,352 epsilon solves. Maximum attained H error
was 1.986e-9 for declared 2e-9 accuracy. Tighter half-response solves use 1e-10;
the expanded diagnostic uses H=2e-11 and half-response H=2e-12.

## Actual reduced action and direction construction

The true conditional analytic free-k matrix, including exact alpha mass, supplies
one seed. The rate-Picard matrix is another seed. Neither is substituted for
Jred. Each used column of B is a one-sided difference of original k residuals
after converged H solves at fixed candidate k. All dependent wall profiles,
moving quadrature/basis evaluations, complete-ray normalization, cutoffs,
mechanical forcing, stress, production/destruction, diffusivity, sensors, chi,
and nonlinear trace quantities are rebuilt by native residual/evaluation calls.

Each column is checked at h and h/2 with tighter H. Unusable columns are rejected
(3% residual-action and 5% maximum individual-observation action limits). Typical
h is 2e-5 in dimensionless coefficient norm, and 5e-6 in the 32-column diagnostic.
Probes at ground contact are one-sided and inward; no undefined negative-k
constitutive state is used. Logs retain branch-change flags, cutoff ranges,
normalization margins and moving-inventory counts at both perturbation sizes.

Two-pass orthogonalization tests independent k seeds and response rank. Seeds
cover ground value/tangent modes, conditional/rate corrections, sampled K and
normalization on both sides, epsilon lower-bound witnesses, spatially distinct
signed qualification witnesses, active cuts, kinetic defect modes and previous
displacements. The model expands toward 32 columns when needed. Thirty-two
columns are still a tested subspace, not the entire 504-dimensional free-k space.

The optional response enrichment (`bdcdc48`, `RANSWallReducedKrylov`, zero
retains the previous space) adds up to four preconditioned missing-response
seeds in the tested configuration. With conditional k preconditioner P and a
verified reduced action Bv, the next seed is -P^-1(Bv-K0 P v)/K0. It extends
the preconditioned Krylov space using the measured non-k response; the rate
Picard block is a fallback preconditioner. A seed is generated only when the
missing response exceeds ten times its empirical two-scale action uncertainty.
Every generated seed is independently subjected to feasible H solves and the
same two-size action/observation checks before use. These seeds occupy the
declared response budget; they do not assemble an all-field matrix. Matched
18-column plain/enriched controls and an enriched original-h0 run are recorded
separately. At the same qualified saved state, the enriched model accepts
fraction .5 (displacement .2) and reduces residual20.232366 ->20.179175 with
zero qualification excess; the plain model accepts1/32 (displacement .0125),
reaching20.226733. After three enriched updates, residual is20.127916 and the
original qualification still passes. This is useful local progress, not stage
convergence or a nonlinear stationarity certificate. A small exact H/F model verifies that this enrichment recovers a
known reduced root missed by the single conditional seed; it is not M8 evidence.

Feedback changes signed directions materially. At the original-run eight-update
checkpoint, the sampled-K direction on the negative-x side has conditional merit
slope +1382.779 but reduced slope -503.291. Actual finite differences give
-502.366 at h and -503.158 at h/2. The reduced action norm is 304.223 and its
feedback component norm 196.358. Near the ground diagnostic endpoint, the
inward-feasible conditional-Newton seed has conditional merit slope +.320858
but reduced slope -.594301, with measured slopes -.589468 and -.592603.
Thus feedback enters the actual model before direction selection and can reverse
the local predicted merit sign.

At the original h0, feedback is much smaller: conditional-Newton action norm
6.822192 versus reduced 6.822456, with feedback norm .00047165; the frozen Picard
action norm is 14.912641. Fixed-point consistency does not make these actions equal.
An independent component experiment at the original-run checkpoint constructs
its sampled-K direction from the native field map at(-2.604649609,1), validates
that point functional against the native sampler, and independently assembles
the free tensor mass. Its native norm difference is4.26e-14. Controlled H solves
use2e-11 and2e-12 accuracy at perturbations1e-5 and5e-6; the ground polynomial
is unchanged by this direction. All conditional partial readbacks verify that
independent coefficients are not silently lifted.

At the smaller perturbation, the conditional-k action has merit projection
+1380.249, while the reduced action has projection-486.135. Direct chain-rule
contributions along the complete non-k response are:

| Independent y field | k-action norm | Original merit projection |
| --- | ---: | ---: |
| Velocity | 5.43352 | +16.7322 |
| Regional pressures/common gauge | 196.56889 | -1883.1109 |
| Outer epsilon/trace | .021309 | +.00574 |

The pressure contribution matches the inner-destruction contribution to about
2.2e-9 in action norm; other k source components have zero direct pressure
response in this test. Source inspection routes pressure through mechanical
wall forcing and the reconstructed destruction profile/complete-ray treatment.
This measurement does not isolate the normalization part from the profile part.
The sum of conditional y contributions leaves remainder .005528 at the larger
perturbation and .002763 at the smaller, against total feedback194.488.
The reduced two-scale action difference is7.05e-5. These are local directional
measurements, not full-space stationarity. Epsilon and velocity also influence
the solved pressure response, so the direct component sizes are not a causal
argument for dropping their equations.

Small-step merit differences include their quadratic term: at the R2 endpoint,
a nearly tangent sampled-K direction has negative first-order slope but positive
finite-h merit changes. Its large action norm explains that curvature; an action
norm check alone would not establish usable finite-step descent.

The small quadratic model uses the original k mass-dual norm, displacement scaling
S=I/K0, regularization 1e-3||S d||² plus a 1e-10 coefficient diagonal. It includes
individual signed observations, both signs, for every applicable OuterBoundsR5
row, including local outer k, sampled K, gradient components, actual epsilon and
projected trace. Domain cuts cover all current scalar evaluation points. Actual
candidates rebuild the full moving inventory, so newly active observations can
be discovered away from the previously worst point. Covariance and lambda remain
complete final-audit gates.

Outside qualification, the objective combines normalized original-residual merit
and a minimax qualification slack. Restoration can raise residual while reducing
qualification violation. The filter tracks these quantities separately. The
model now includes the filter's allowed qualification ceiling, rather than
rejecting every proposal before an actual line search when its full-step
violation exceeds that ceiling. Inside qualification, the model uses all original signed
bounds while reducing the original residual. The actual search filter permits
a normalized numerical excess below1e-8; every positive excess triggers
restoration on the next model. This search allowance never qualifies a physical
state: the complete original final audit has its unchanged exact bounds.

Two accepted steps with less than 1% relative qualification improvement reduce
the residual multiplier by ten and expand the response budget. The default mixed
control retains a minimum multiplier1e-4. A separate explicit
`RANSWallReducedRestorationWeight=0` control performs pure qualification
restoration outside the box, then restores the full original-residual objective
inside it. This tests whether a finite mixed penalty is itself holding a state
outside a feasible box. A known scalar model demonstrates such a soft-penalty
fixed point; it does not establish an M8 infeasibility result.

Qualification progress and filter comparisons now use relative violation
reduction, including close to the original bound. The previous absolute
comparisons could reject a 50% improvement when the normalized violation was
1e-5. Every positive violation remains outside the original envelope; the H
accuracy is tightened with the remaining violation, down to the original stage
tolerance. These are solver controls. No physical bound or stopping tolerance
is widened.

A native mixed-restoration example reduces violation .2720621 -> .2121202 while
raising residual14.06670 ->15.17449; the next update reaches violation .1601001
at residual16.30274. This verifies actual restoration with a residual increase.
The original-start mixed run has also reached violation .01280786 at residual
19.60144 (still an interim state, not a qualified root).

Active cuts use two-pass QR/null-space solves rather than the ill-conditioned
normal Gram system. The final combined correction is bounded after cut handling,
and predictions are recomputed. This is a radius-bounded model search; radial
scaling is not claimed to preserve the unscaled affine model's optimum. Logged
small affine KKT defects are not nonlinear stationarity certificates. Original
merit gradients and empirical two-scale action uncertainty are logged separately.

A small inward coefficient correction handles arithmetic contact only when the
computed continuous minimum is within 256 double-epsilon times the field scale.
The correction lies in an already measured ground-lift column, is re-bounded by
the radius and rescanned as a continuous polynomial. Observed coefficients were
1.776e-15: minima near -1.6e-17 became about 1.18e-16, above their arithmetic
scale 3.04e-17. This changes the numerical proposed polynomial; it is not clipping,
a k floor, a source change, or relaxation of residual/qualification acceptance.

A near-envelope restoration model adds an inward search margin of
min(1e-6,1e-4 times the original k residual norm) in normalized observation
units when the positive excess is below1e-4. This addresses the measured
unbuffered progression1.28e-5 ->1.78e-9 ->4.00e-12, whose last predicted
excess was2.44e-13: aiming at exact contact reaches the scale of inner error.
The margin shifts only the restoration model, shrinks with original residual
progress and is zero for an inside-envelope model. Its first native saved-state
update passes the original complete qualification at residual20.240190578:
gradient .0009999990002, epsilon minimum ratio .9846519511, sampled-K maximum
ratio1.0099999900, continuous-ground minimum .0005205281. Predicted residual
20.240190575 agrees with actual20.240190578 at displacement4.35512e-7
within radius .4. This is a qualified constitutive state, not an algebraic root
or an accepted physical stage. Actual trial scoring and
complete physical acceptance retain the original bounds. A standalone tiny-
excess test verifies crossing the unchanged bound instead of asymptotic
outside contact; it passes Linux ASan/UBSan and bounds assertions.

The enlarged original-h0 response space exposed a radius-handling weakness:
its raw affine model predicted qualification excess .411883, but post-solve
radial scaling to radius .4 raised the prediction to1.700588. The smaller
space's in-radius direction had predicted .456267. Thus adding useful response
directions could degrade the bounded proposal. The radius now participates in
the small constrained model through a nonnegative multiplier on V^T V. A bounded
search keeps the converged inside-radius bracket and reports the multiplier,
raw/final field norms, solve count and radius complementarity. The actual final
correction is still re-bounded after ground arithmetic handling and rescanned.
This multiplier changes direction construction, not the physical mass or source.
Anisotropic known-answer tests show the advantage over radial scaling, preserve
an independent slack, enforce an active cut and use the actual field norm for
an oblique response basis. They pass ASan/UBSan. Native controls confirm the improvement. At original h0 the raw correction
norm3.77048 becomes an in-model correction .399804 within radius .4, with
multiplier .014856 and predicted/actual qualification excess .450799/.494633.
The radius model uses16 QP solves (4,486 affine iterations); its radius
complementarity is2.33e-6, separately from affine stationarity1.48e-15.
At the matched qualified endpoint, predicted/actual residuals are
19.995316/19.995620 at a half step, compared with20.179175 for post-solve radial
scaling. After three radius-aware updates the residual is19.887304, gradient
.0009999545 and ground minimum .001078865; qualification passes, full algebraic
convergence does not. The original-start complete run remains in progress.

The buffered model's ground-cut restart and saved replay initial point also now
include the inward buffer. The earlier artifact could encode an infeasible
starting slack despite a correctly solved first native model. Replay evidence
uses the Hessian including the radius multiplier and the corrected feasible
initial point. ASan/UBSan replay passes both the radius-active and positive-buffer native
packets. The latter has affine primal/stationarity defects4.30e-23/2.03e-22.
An affine replay certificate remains a model check.

## Native progress and nonlinear trial agreement

The original unmodified h0 first has gradient departure .000500195, inside the
bound. At the first fully equilibrated non-k state with k still equal to h0 it
is .002857296, witness (x,z)=(1.625,2), component (0,2), geometry34/outer side.
Thus the violation is already present before the first reduced k update. The
current logs locate this first complete H state, not its individual internal
velocity/epsilon sub-iteration.

The following are numerical stage iterations only; none is an accepted physical
step. The final primary-run entry is pending completion.

| Experiment | Complete residual | Gradient maximum | Role |
| --- | ---: | ---: | --- |
| Historical best R2 endpoint | 2.908833889 | .003826405 | retained control |
| R3 mixed original h0, six updates (`82ca255`) | 32.866814956 | .001660092 | bounded native control, fresh readback |
| R3 pre-filter-repair original h0, eight updates (`e7313f7`) | 24.415622070 | .001273039 | filter mismatch; no convergence |
| R3 repaired filter, three more updates from that saved guess | 18.566628956 | .001264727 | additional warm-start diagnostic |
| R3 adaptive mixed h0, 17 numerical updates (`1c23b84`) | 19.494577545 | .001011874 | historical progress snapshot |
| Pure restoration from mixed checkpoint, first update (`8687dc3`) | 20.239805178 | .001000012818 | recorded first-update diagnostic |
| Inward saved-state restoration (`ef83196`), first update | 20.240190578 | .000999999000 | passes qualification only; not an algebraic root |
| Radius-aware saved-state search (`9a5b58d`), three updates | 19.887303979 | .000999954475 | passes qualification, not an algebraic root |
| Current original-h0 runs (`9a5b58d`) | pending | pending | primary final tests |

The saved R2 endpoint is used only for diagnostics. Its first reduced restoration
step changed residual 2.908834 -> 3.271267 while improving gradient .003826405 ->
.002337314 and epsilon lower ratio .974454851 -> .980310747. The gradient witness
migrated from x=1.666667 to x=-.871740. The complete inventory detected that
migration; differentiating only the old maximum would miss it. The ground
arithmetic repair then permitted further updates to residual 3.384774 and
gradient .002150084, still outside the original qualification.

Representative predicted/actual reduced trial norms and normalized maximum
qualification excesses are:

| Trial | Predicted residual | Actual residual | Predicted violation | Actual violation |
| --- | ---: | ---: | ---: | ---: |
| Original h0 first k update, radius .4 | 60.0947648 | 60.0946233 | 1.7674851 | 1.7671654 |
| Original h0 third k update, radius1 | 49.3828514 | 49.3919176 | 1.4092333 | 1.4092871 |
| Saved R2 endpoint first restoration | 3.2702170 | 3.2712669 | 1.3373250 | 1.3373141 |
| Repaired checkpoint first update, radius1 | 23.2713918 | 22.9391808 | .2674437 | .2675366 |
| Expanded ground response, 32 columns/rank32, radius .4 | 3.3114674 | 3.3149908 | 1.1442756 | 1.1444590 |

Individual signed observations also agree. At a migrated ground endpoint, the
negative epsilon lower-bound observation changed by predicted -.000268710415
versus measured -.000268710153; the sampled-K lower observation changed by
-.000364407914 in both checks. The negative gradient witness at x=-.871740
changed by predicted -5.111307e-5 versus measured -5.111284e-5 on that response.
The full inventories, both signed bounds, branch summaries and two-scale actions
are in the compact evidence, not only the maxima in this table.

A later useful full step exposed finite-step source curvature: predicted residual
16.80015 versus actual21.72800, full vector model error10.82331. Native component
reconstruction attributed10.82317 to inner destruction, with remaining error
.013539. The half step passed the actual filter and reduced residual to20.42408.
This supports safeguarding actual trials, not importing the coupled R8/R9
optimizer or treating a better source prediction as evidence of a useful direction.
The source model, normalization and destruction remain fully nonlinear in trials.

The supplementary coupled R9 iterate is a separate diagnostic starting state.
A fresh original readback matches all2,165 saved residual coefficients exactly.
It is inside the original box but has non-k residual as large as .220237.
Holding every k coefficient fixed and converging H raises k residual2.910254 ->
11.879869 and gradient .000999985 -> .003322874. Ground remains5.93348e-6.
It therefore supplies neither a controlled Y(k) nor a qualified root before
that equilibration. Its k state is not substituted for the original-h0 test,
and no coupled direction matrix or optimizer is imported.

## Algebraic continuation and the qualification obstruction

The separate opt-in continuation evaluates
F_eta=F_target+eta M_D(Q-h0), with original physical alpha400, h0, t=.0025 and
external data fixed. The retained segregated conditional solver is used as a
control. The wrapper checks the exact mass identity (relative defects about
1e-16) and always rejects physical advancement, even at eta0. Positive-eta roots
are diagnostic only. The eta-zero/original acceptance tolerance stays2e-12.

A strict eta5000 attempt stopped at about9.34e-12 and did not meet the
2e-12 algebraic tolerance. A separate explicitly declared eta-positive tolerance2e-11 allowed
following the nearby branch. Original residuals on these states remain tens,
not zero; each continued and original residual was reconstructed independently.

| Eta | Gradient departure | Ground minimum | Original qualification |
| ---: | ---: | ---: | --- |
| 5000 | .000500195 | .40511246 | passes |
| 2500 | .000582927 | .33036459 | passes |
| 1250 | .000967567 | .22662976 | passes |
| 1209.375 | .000988819 | .22106282 | passes |
| 1168.75 | .001011025 | .21525707 | fails gradient |
| 1087.5 | .001058584 | .20286116 | fails gradient |
| 925 | .001168619 | .17435135 | fails gradient |
| 600 | .001476320 | .09561097 | fails gradient |
| 500 | .001606738 | .07065555 | fails gradient |

The followed branch exits the original gradient envelope between eta1168.75
and1209.375, at witness (1.604649609,2), component(0,2), geometry34/outer side.
Ground is still well positive there. Other scalar ratios, covariance and lambda
pass on these bracket roots. The eta1209.375 root repeated from h0 and from a
nearby branch root differs by at most9.326e-14 in coefficients, with identical
history. This supports the same nearby branch in that comparison, not global
uniqueness or global absence of another qualifying target root.

The retained control failed at eta400 (a separate larger ladder also failed at
eta200), before certifying an eta-zero root. At eta400 its last accepted numerical
state still had positive ground minimum .05153 and substantial residual; this
failure is not a ground-contact bracket. No claim of nonexistence or a floating-
point lower bound follows from the failed line searches.

A separate spatial diagnostic splits only [0,.01] at .005, increasing unknowns
2165 ->2315 and regenerating that mesh's own represented h0/history. Physical
parameters, spaces/degrees, quadrature, functions and timestep stay unchanged.
It is explicitly a different discrete problem and never an original M8 success.
The auxiliary eta600 root has ground .07856770 and gradient .001480973; eta400
converges there with ground .03993422. Eta0 still failed to certify a root.
This one refinement changes the near-ground response and the retained control's
convergence behavior, but does not remove the gradient failure at eta600 or prove
a general spatial positivity obstruction. A mesh-convergence theorem is not claimed.

## Verification, provenance, and outstanding work

The fixed-state readback now optionally requests the original complete
qualification and continuous-ground audit in a fresh process. A negative physical
ground vertex is rejected before evaluating the candidate scalar model. The
focused guard check passes; the successful fixed-state audit reproduces the
previous residual norm vector bitwise. Optional audit flags leave historical
probe behavior intact. This is an observation, not native physical acceptance.


The `_GLIBCXX_ASSERTIONS` Linux solver exercises the repaired predictor with
verification On. Verified and unverified two-sweep states and norms are bitwise
identical. Predictor continuity/gauge is [8.54316e-4,2.53250e-4,1.02830e-16];
after correction it is [6.87347e-16,8.68001e-16,6.21722e-17]. Direct saddle
coefficient difference is3.60822e-16. Block sizes remain velocity980, pressure85,
free-k504 and epsilon/trace254. This verification defect is repaired independently
of the unresolved k/qualification problem.

Six focused Linux CTests pass: RANSReducedModelChecks, RANSEquilibriumTimeNative,
RANSEquilibriumEvolutionNative, RANSEquilibriumPatchNative,
RANSEquilibriumWallEvolutionSmoke and RANSEquilibriumWallStageSmoke.
The standalone model additionally passes Linux ASan/UBSan and bounds assertions,
including known feedback correction, duplicate/near-dependent cuts, constraint
discovery, epsilon lower observation, restoration, actual final radius, tiny
corrections, filter ceiling, failure and programming-shape cases. Model checks
are not labelled M8 convergence. Native response-checker results explicitly use
`passed_operator_checks`, retaining failed stage results separately.

Heavy artifacts are under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard-r3`.
The shared Linux checkout and installation were preserved; the shared revision
is `f4b70baec10016dd098005bc5a5507b87e4a57bf` and its pre-existing symlink content
status is recorded in the build manifest. Source/build trees and immutable copied
executables are isolated. Eight explicit compiler workers, one test at a time,
one rank and one thread per BLAS/OpenMP runtime were used, within the72-worker
limit including the simultaneous bounded experiments. No support claim extends
beyond the existing serial/flat tensor-cubic ground implementation.

The current original-start command can be reproduced in the isolated Linux tree
(the retained run directories must not be overwritten):

```bash
R=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard-r3
S="$R/continuation-source"
B="$R/continuation-build"
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$R/tmp"
export OMPI_MCA_orte_tmpdir_base=/tmp PMIX_MCA_gds=hash
cmake -S "$S/solvers/IncNavierStokesSolver/Tests/RANS" -B "$B" \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/home/jay/.local/opt/nektar-v2026-mpi \
  '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake' \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$B" --parallel 8 --target IncNavierStokesSolver RANSWallResolutionProbe RANSReducedModelChecks
python3 "$S/solvers/IncNavierStokesSolver/Tests/RANS/run_cls_epawf_picard.py" \
  --solver "$R/binaries-9a5b58d/IncNavierStokesSolver" \
  --resolution-probe "$R/binaries-9a5b58d/RANSWallResolutionProbe" \
  --reduced --reduced-verify --reduced-basis 18 --reduced-krylov 4 --reduced-radius .4 \
  --reduced-restoration-weight 0 --sweeps 40 --timeout 18000 \
  --output "$R/radius-original-reproduction"
```

The harness deliberately returns nonzero unless the original stage completes;
run the operator checker separately on bounded failed-stage controls. Each
native numerical command is one rank. Logs and manifests distinguish solver,
readback, expected rejected-stage and startup return codes.

[Build and sanitized regression evidence](cls-epawf-segregated-picard-r3-evidence/native/build-validation.json),
[pressure verification](cls-epawf-segregated-picard-r3-evidence/native/verification.json),
[development controls](cls-epawf-segregated-picard-r3-evidence/native/development-controls.json),
[updated native checkpoints](cls-epawf-segregated-picard-r3-evidence/native/development-updated.json),
[continued roots and independent readbacks](cls-epawf-segregated-picard-r3-evidence/native/continuation-aux.json),
[current build/control evidence](cls-epawf-segregated-picard-r3-evidence/native/build-validation-current.json),
[8687dc3 regressions and pressure verification](cls-epawf-segregated-picard-r3-evidence/native/build-validation-8687dc3.json),
[inward-search build and sanitizers](cls-epawf-segregated-picard-r3-evidence/native/build-validation-ef83196.json),
[enriched-source regressions and verification](cls-epawf-segregated-picard-r3-evidence/native/build-validation-bdcdc48.json),
[radius-source build and replays](cls-epawf-segregated-picard-r3-evidence/native/build-validation-9a5b58d.json),
[matched response-space controls](cls-epawf-segregated-picard-r3-evidence/native/response-enrichment-completed.json),
[completed radius controls](cls-epawf-segregated-picard-r3-evidence/native/radius-completed.json),
[sanitized native model replays](cls-epawf-segregated-picard-r3-evidence/native/radius-replay-checks.json),
[independent restored-state audit](cls-epawf-segregated-picard-r3-evidence/native/audit-inward-qualified-1.json),
[independent feedback components](cls-epawf-segregated-picard-r3-evidence/native/feedback-components.json),
[completed native controls](cls-epawf-segregated-picard-r3-evidence/native/development-completed.json),
[fresh residual/qualification audits](cls-epawf-segregated-picard-r3-evidence/native/native-audits.json),
[continued branch below eta600](cls-epawf-segregated-picard-r3-evidence/native/continuation-near-ground.json),
[repeated auxiliary root](cls-epawf-segregated-picard-r3-evidence/native/continuation-repeat.json),
and [changed-mesh diagnostic](cls-epawf-segregated-picard-r3-evidence/native/ground-mesh.json)
record producing revisions, executable/library hashes, exact numerical commands,
XML hashes, work, elapsed time, block histories, qualification witnesses and
independent residual reconstructions. Large raw inventories and binary replay
packets remain on Linux. Earlier pre-mount compiler evidence is preserved.

The unchanged first stage is still unresolved. Three qualified native steps to
.0075, small admissible nonlinear-guess perturbations of that physical target,
and final qualified-trajectory readback remain outstanding. Deliberate rejection
and completed failed reduced stages preserve accepted h0, native time and history;
killed processes without native receipts are recorded as unknown, not as verified
rejections. No eta-positive or changed-mesh state was used to advance the original
physical trajectory. The current evidence does not establish global infeasibility,
full-space reduced stationarity, or existence of a different qualifying root.
