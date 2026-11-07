# Segregated CLS/EP-AWF Picard R3: sequential reduced directions

**Development status: the original M8 problem remains unresolved. No qualified
physical stage has been accepted; accepted physical steps/time remain 0/0.**
The designated volume was restored and native Linux development resumed.
A later concurrent task switched the shared Mac checkout; subsequent work uses
an isolated Mac worktree and preserves that task's PIMPLE branch. Only the
Picard commits are delivered to the requested remote branch. The
completed original-start validations have not converged. Larger
response-space and saved-endpoint follow-ups continue. Completed controls and explicitly dated progress snapshots are
distinguished below; no running experiment is counted as a successful stage.

R3 now constructs a constrained k direction from the response of converged
velocity, pressure/gauge, and epsilon/trace solves. It no longer chooses a
conditional Jkk direction and only optionally settles the other equations.
The R2 conditional method, old backends, verification repair, restricted k
analytic assembly, original residual, and historical evidence remain intact.
The R2 projected map could be stationary in its frozen A_k^T Mff^-1 A_k
metric at ground contact while its original k residual stayed nonzero. R3
therefore scores the original reduced residual, not that projected-map defect.
No full all-field Jacobian or hidden coupled direction solver was added.

## Current-base action preservation (ff1c4fc)

The completed original-h0 control `current-original160-3step` now reaches
**13.4621750005** after twelve numerical updates. All 3,936 controlled non-k
evaluations converge, using 15,400 sequential flow and epsilon solves each.
A fresh process reproduces every original residual norm bit. Ground minimum
is **3.91948e-5**, gradient departure **0.0009999990412**, and actual epsilon
lower ratio **0.9800122691**; complete qualification passes. The stage still
fails the original residual tolerance, accepts **zero physical steps**, and
preserves accepted state and every physical-history slot. This supersedes the
17.7529 eight-update result as the best completed original-start result in this
collection; numerical updates are not physical steps.
[Completed original-start readback and rejection evidence](cls-epawf-segregated-picard-r3-evidence/native/followup-completed-review.json).

The first e9bbda5 reference-bridge attempt exposes a concrete direction-model
defect at original-start numerical point two (residual **49.2816311106**).
Independent same-k H solves converge near 1e-12 and agree, and two-size
direction errors stay between **4.3e-7 and 1.56e-5**. Yet the current-direction
model error grows **0.04018, 0.04346, 0.09588, 0.20879, 0.44902** across four
rank-one refreshes. Each newest action is interpolated near machine precision,
but subsequent correlated directions lose earlier measurements. The guard
correctly discards that history and rebuilds fresh responses; no accepted
bridged direction or M8 success is inferred from this failed cache attempt.

A separate one-thread reconstruction reads only the small matrix/vector
prefixes of the native affine packets, solves their recorded active faces,
and reproduces all five displacement norms within **1.6e-15**. Successive new
direction components outside the earlier span are only **1.34%, 1.09%, 0.0943%
and 1.07%**. This quantitatively supports interference between correlated
rank-one calibrations; it is not a reduced stationarity or infeasibility claim.
[Independent native direction reconstruction](cls-epawf-segregated-picard-r3-evidence/native/current-direction-interference.json).

Commit `ff1c4fc5ac2ce75279f821dd944dcd8a97648bbc` replaces that same-base refresh
with a protected multi-secant update in the actual k coefficient-field metric.
Two-pass orthogonalization removes previously measured directions from the
new update covector, preserving all their residual, non-k-state and projected
trace actions. A relative independent component below 1e-4 is rejected without
changing the model. Every interpolation and retained-action defect is checked;
scratch measurements are cleared after an actual k update or a fresh model.
The response basis, exact independent k action, qualification cuts, actual
radius bound, tight two-size/current-branch checks, and physical acceptance
authority remain intact. Endpoint secants from different states are not treated
as exact protected derivatives.

The four-worker Linux build passes in **88.55s**, all six focused CTests pass
in **40.52s** at concurrency one, and a fresh standalone ASan/UBSan run passes
in **9.07s** after **25.88s** compilation. Its new known-answer test uses a
nonorthogonal field basis and correlated queries: the old update loses an
earlier independently known action, the protected update retains all actions,
and a dependent inconsistent query leaves the model unchanged. Sanitizer scope
is this helper and its tests, not the full solver or installed dependencies.
[Build, sanitizer commands and immutable identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-ff1c4fc.json).

Matched original-h0 native reruns use unique `protected-actions-control` and
`protected-actions-enabled` directories and the same producing image. Only
reference-weight reuse differs; both use the protected calibration. They keep
the original 2,165-unknown alpha-400/time-.0025 target and tolerance 2e-12.
`check_cls_epawf_protected_actions.py` additionally requires at least two actual
current-base measurements, retained-action checks, and an accepted qualified
original-start correction selected with those responses. Native outcomes are
pending; passing the helper is not a convergence result.
[Exact native commands and allocation](cls-epawf-segregated-picard-r3-evidence/native/protected-actions-native-launch.json).

The source-stencil recovery also now has independently read original-start
evidence: at numerical update five, five columns require eleven stencil
halvings, reducing the source derivative discrepancy from **0.05170 to 0.02507**
below the unchanged 0.03 limit. The resulting checkpoint residual is
**18.2593339487**, with complete qualification and bitwise fresh original norm
agreement. This remains an in-flight checkpoint, not a completed-stage history
receipt or physical acceptance.
[Original-start stencil recovery and separate-process audit](cls-epawf-segregated-picard-r3-evidence/native/source-precision-original-point5-review.json).

## Guarded reference-weight history bridge (e9bbda5)

The early original-start reuse control makes 8, 3 and 7 E/M profile-weight
switches over its first three numerical corrections. Its baseline B and
destruction P/C classifications do not change. At the first qualified-to-
qualified correction, the original linear prediction error is only **1.47374**
against the reuse policy's **12.7487** error limit, but the literal reference
branch comparison discards the response history. This explains why the
endpoint reuse improvement does not yet accelerate that original-start prefix.
It is a search-policy limitation, not evidence that H or the original equations
are solved or that a qualifying root cannot exist.
[Measured original-start switches and prediction errors](cls-epawf-segregated-picard-r3-evidence/native/original-start-reference-switches.json).

Commit `e9bbda5382ee45bc69739b8649493330b38b82a7` adds opt-in
`RANSWallReducedReferenceBranchReuse` /
`--reduced-reference-branch-reuse`, default Off. It permits only E/M switches
between **actual converged endpoints inside the qualification observation bounds** when forming a historical
secant model. It requires the measured-action refresh policy. Baseline and cap
switches still rebuild fresh responses, and the usual age, rank, signal and
actual-progress safeguards remain. It changes reuse policy; it does not change
the model's physical reference or wall-profile weights.

Before any larger trial, both fresh direction probes must remain on the current
selected branch and pass signed residual/individual-constraint checks. An
additional independent H solve reaches the identical prescribed k from the
other converged response. The residual action, all observation actions and
effective coefficients must agree within the declared limits; fixed k, attained
H errors and branch identity are recorded. Failure discards the historical
model while preserving the accepted state/history. The recorded old branch is
never presented as a current Jacobian certificate. Switch records include their
normalization/trace anchors.

The four-worker Linux build passes in **129.90s** and all six focused CTests
pass in **39.32s**, with test concurrency 1. Matched original-h0 runs request
three numerical updates, with three physical steps requested if the original
stage actually converges. They keep 2,165 unknowns, alpha 400, time .0025,
tolerance 2e-12, all original physics/history, and use the same executable and
numerical settings except this switch. The native branch-use test requires an
actually accepted qualified original-residual reduction using the bridged
direction, independent same-k/current-branch checks, unchanged numerical prefix
and fresh endpoint readback. Its name is
`check_cls_epawf_reference_branch_reuse.py`; native outcomes remain pending.
The longer original-start and complete-space attempts continue independently.
[Producing build, flags and identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-e9bbda5.json),
[matched original-start commands and allocation](cls-epawf-segregated-picard-r3-evidence/native/reference-branch-native-launch.json),
[Mac input preparation checks only](cls-epawf-segregated-picard-r3-evidence/native/reference-branch-input-checks.json).

## Current follow-up: measured action refresh and full-space progress

The original first stage is still unresolved. The earlier completed original-h0
eight-update run reaches residual **17.7529247432**, improving the preceding
18-update control's **19.1747990683**. Original qualification passes: gradient
departure **0.0009999992023**, epsilon lower ratio **0.9800000241**, and continuous
ground minimum **0.0138760891**. A fresh process reproduces the original norms;
accepted physical state/history remain unchanged. Later original-h0 searches
continue; an in-flight checkpoint is not a completed-stage receipt. All
physical acceptance counts remain **0 steps, time 0**. Three-step evolution,
same-history guess perturbations, and rejection after an accepted stage remain
outstanding.
[Completed original-start and matched follow-ups](cls-epawf-segregated-picard-r3-evidence/native/followup-completed-review.json).

The saved-guess full504 experiment now makes a separately audited correction
on the unchanged original operator. Its tightly equilibrated base residual
**4.1182787654** becomes **3.8960677838**. The final destruction-model prediction,
including one qualification correction, is **3.8948918678**; actual/predicted
squared-merit gain is **0.99488194**. Actual displacement **0.25842899** is within
radius **0.4**. The independent process reproduces all original residual norm
bits, checks the negative-ground guard, and rebuilds all dependent quantities.
The largest non-k norm is **1.80162e-12**, continuous ground minimum is
**6.67452e-7**, gradient departure is **0.0009999987966**, and epsilon lower ratio
is **0.9819776988**. Complete original qualification passes. This is a diagnostic
numerical correction with the original h0/history, not an accepted physical
stage or a replacement for original-start validation.
[Correction, prediction and independent readback](cls-epawf-segregated-picard-r3-evidence/native/qualified-full504-latest-9f-iteration0-review.json).

At that experiment's initial base, independent rank-504 first-order analysis
finds gradient norm **551.795** and near-active cone defect **519.182**, against
empirical two-size gradient uncertainty **2.530**. These observations and the
accepted correction do not establish a stationary obstruction. At the separate
4.85782 full-space base the corresponding defect is **1703.77**, with uncertainty
**4.107**. No nonexistence claim follows from either measured local model.
[Lower-base first-order audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-qualified-full504-latest-9f-iteration0.json),
[second base audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-qualified-full504-at5p182-iteration1.json).

The four-update 7ba50d3 response-reuse comparison has completed its focused
checks, full original readbacks, and failed-stage history checks. Reuse reduces
the controlled H evaluations from **807 to 470** and elapsed native target time
from **6244.89s to 4472.57s**, while giving the weaker terminal residual
**3.94084132 versus 3.93120535**. This is a work/progress tradeoff, not an
equal-accuracy speedup or successful M8 convergence. Initial state and all
numerical states before reuse match bit for bit; both endpoints pass original
qualification, accept zero physical steps, and preserve the physical history.
[Matched direction assertions](cls-epawf-segregated-picard-r3-evidence/native/response-chord-enabled-response-chord-checks.json),
[execution receipt](cls-epawf-segregated-picard-r3-evidence/native/response-chord-direction-comparison-receipt.json).

Tightening the common H reference for its first reused direction reduces the
two-size residual discrepancy from **0.0219014 to 0.0000238823**. The underlying
reference change in original k residual has mass-dual norm **3.44989e-7**,
whereas the half-probe signal is only **1.10908e-6**. This identifies an inner
accuracy sensitivity in small cancelling directions. Its historical residual
model discrepancy remains **0.0951773**, normalized by
`max(1, measured action norm)`; the measured action norm is **0.141164**. Thus
0.0951773 is not a 9.5% relative-error claim. Neither this model discrepancy nor
interpolation KKT accuracy certifies a current reduced derivative.

Inspection found a separate normalization error in the new observation
diagnostics: `OuterBound::Action` already divides by `scale*width`, but the
41ef144 diagnostic divided the action discrepancy by `width` again. Its reported
**1.51201/0.0257185** values are therefore superseded and must not be interpreted
as physical bound-width errors. The original qualification cuts, residuals and
acceptance audits were not affected. The initial optional 4e024ae refresh gate
used those overly strict diagnostic units; its numerical artifacts are retained
as a superseded variant, and its idle follow-up controller was stopped before
launching a new primary test. Its numerical processes were preserved.

Commit `7300df33eb5130415a6d615f0e61c49ef1aa3871` applies this scaling once and
records the raw field contractions, normalization scale and bound width at both
maximum-error witnesses. The automated check independently reconstructs each
reported error from those raw values. Protocol version 2 identifies corrected
units; the focused refresh test refuses the superseded protocol.
[Measured and historical failed checks](cls-epawf-segregated-picard-r3-evidence/native/followup-checks-20261107-0010.json).

The corrected update, introduced in `4e024ae` and repaired in `7300df3`, adds opt-in
`RANSWallReducedDirectionalRefreshUpdates` /
`--reduced-directional-refresh-updates`, default 0, maximum 4. A reliable
same-base two-size direction measurement now updates residual, effective-field
and projected-trace actions **before re-solving** the constrained direction.
Every individual signed qualification row is rebuilt from these actions;
independent k actions remain exact. The update does not change the accepted
nonlinear state or physical history. Historical directions must agree within
0.03 residual units normalized by `max(1, measured action norm)` and 0.1
qualification-bound widths. Measurement uncertainty must also pass the
unchanged two-size checks and be at most 0.05 bound widths before an action is
used. Exhausted or unreliable updates rebuild fresh responses. A refreshed
model remains explicitly historical; it cannot be promoted to a fresh
first-order or nonlinear stationarity certificate. Affine replay packets retain
distinct names across re-solves.

The corrected four-worker Linux build passes in **83.81s**, and all six focused
CTests pass in **43.37s**, with test concurrency 1. The precursor 4e024ae build
also passed its build and regressions; that does not validate its erroneous
diagnostic normalization. A corrected same-image two-update native
comparison uses refresh budgets 0 and 4 with identical physical inputs. The
focused checker requires an actually accepted correction using the re-solved
historical model, unchanged state/history during refresh, subsequent independent
signed action checks, and fresh full residual/qualification readback. It is
named `check_cls_epawf_direction_refresh.py` and does not assert M8 convergence.
The corrected same-image native comparison now **passes**. Both runs first
reach exactly the same original state/residual, **3.9711215458**. One measured
action refresh then reduces its direction-model discrepancy from **0.0951773
to 0.00746080** and its correctly normalized observation discrepancy from
**0.00241572 to 0.000803762**. A separately measured source-refined direction
has model discrepancy **0.00658842**, two-size discrepancy **4.65834e-6**, and
actual merit slopes **-0.09530548/-0.09530475**. The test independently rebuilds
the error normalization from raw field contractions; it does not accept the
superseded units or interpret interpolation defects as stationarity.

Its actual original residual reaches **3.9587753030**, versus **3.9637517592**
with refresh disabled. The source prediction after two qualification corrections
is **3.9540112155**. Original qualification passes: gradient departure
**0.0009999990455**, epsilon lower ratio **0.9829384300**, and continuous ground
minimum **6.75284e-6**. Actual displacement **0.59187820** stays within radius
**0.72**. All 214 H evaluations converge, with 936 separate flow/epsilon solves;
the control uses 212 H evaluations and 918 solves. Independent readback is exact
and both rejected physical stages preserve accepted state/history. Each accepts
zero physical steps. This is useful original-equation progress, not M8
convergence. The 96-direction historical model remains a measured search model,
not a freshly assembled 504-direction Jacobian certificate.

A fresh original-h0 three-step attempt has now started with corrected code,
and an additional eight-update full504 experiment starts from the separately
audited **3.8960677838** diagnostic point. The latter uses H tolerance 2e-12 and
keeps the original physical history. The matched corrected-code verification-Off
repeat now **passes**: intermediate numerical state/norm bits, native accepted
and rejected states, stage bookkeeping and all physical-history slots match.
The longer original-start/full-space tests have no completed-stage claim.
[Corrected build and immutable executable/library identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-7300df3.json),
[corrected matched native commands and allocation](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-native-launch.json),
[corrected native direction assertions](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-checks.json),
[assertion execution receipt](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-comparison-receipt.json),
[fresh original-h0 command](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-original160-3step-launch.json),
[additional full-space command](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-full504-at3p896-launch.json),
[verification-Off repeat command](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-unverified-launch.json),
[corrected bitwise verification/history assertions](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-verification-comparison.json),
[verification execution receipt](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v2-verification-comparison-receipt.json),
[superseded diagnostic scope](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-v1-superseded.json),
[Mac input preparation checks only](cls-epawf-segregated-picard-r3-evidence/native/direction-refresh-input-checks.json).

Commit `987bea8f4a9ea5444318661ced91d3100e2928eb` shares the actual pointwise
inner destruction law between the full diagnostic and a restricted value-only
evaluation. The opt-in `RANSWallReducedFastDestruction` selector refreshes all
wall/normalization dependencies and retains the exact native accumulation; it
avoids unrelated equation assembly when the source model asks only for this
component. The first three native samples agree bit for bit at all 504 free k
rows, independently of verification. Their restricted times are
0.174/0.179/0.177s versus full diagnostic 0.324/0.322/0.327s. These three timings
are component observations, not an end-to-end solver speedup. The Linux build
and six focused regressions pass. The matched enabled/disabled two-update
comparison now passes: all original state/norm bits agree and both original
readbacks and failed-stage physical-history checks pass. Native target elapsed
times are **1771.39s enabled versus 1969.29s control**, with the same 212 H
evaluations and 918 separate flow/epsilon solves. This is one matched fixture
under concurrent server load. A fresh original-h0
three-step attempt has started after
the guarded kernel produced an actual qualified correction. All actual trials
still use the complete original equations, complete nonlinear destruction,
qualification and physical acceptance gates. No new full-solver sanitizer claim
is made for either orchestration change; the small unchanged model helper
retains its earlier sanitizer evidence.
[Value-kernel build](cls-epawf-segregated-picard-r3-evidence/native/build-validation-987bea8.json),
[matched value-kernel test](cls-epawf-segregated-picard-r3-evidence/native/destruction-value-native-launch.json),
[bitwise component and full-trajectory assertions](cls-epawf-segregated-picard-r3-evidence/native/destruction-value-checks.json),
[check execution receipt](cls-epawf-segregated-picard-r3-evidence/native/destruction-value-comparison-receipt.json),
[original-h0 command](cls-epawf-segregated-picard-r3-evidence/native/destruction-value-original160-3step-launch.json).

The tighter 41ef144 verification-On/Off four-update native comparison also
passes: all numerical state and original norm bits agree. Its final original
residual is **3.9435906107**, original qualification passes, independent readback
is exact, and accepted physical steps/time remain 0/0. This confirms the
verification selector does not alter that trajectory; the superseded diagnostic
units are a separate issue.
[Tight-reference functional checks](cls-epawf-segregated-picard-r3-evidence/native/response-chord-tight-checks.json),
[verification execution receipt](cls-epawf-segregated-picard-r3-evidence/native/response-chord-tight-verification-comparison-receipt.json).

The older fa5ac96 eight-update warm run ends at residual **3.93839142**, but its
combined-direction check at iteration 6 has half-probe residual discrepancy
**0.0337936**, above 0.03. Its operator-check postprocessor therefore fails.
That failed check is preserved; the run is not counted as passing the focused
operator suite. Its converged H errors are 5.98e-10/1.81e-9 at a direction
fraction 1.42424e-5. Tight common-reference measurements and the new
bound-scaled direction gate address the demonstrated measurement/model issues;
they do not retroactively validate this historical check.

Heavy logs, all dense replay rows, checkpoints and failed-check logs remain under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard-r3/`.
The dated follow-up statements above supersede older in-flight snapshots below.

The fd8c034 full504 source-stencil experiment has also completed its focused
native check. One column requires one smaller stencil: its error changes from
**0.0382067 to 0.0178874**, within the unchanged 0.03 bound. All 504 surrogate
derivatives become available with 1,009 model evaluations, below the 2,080
budget. The resulting direction has predicted merit slope **-5.84426** versus
actual **-5.84520** at both probe sizes; H errors stay below 2e-12. It produces
an actual qualified correction from **5.1819703505 to 4.3905627178**, versus the
older same-input control's first **4.8578168219** endpoint. The final source
prediction is **4.3527817785**. Two qualification corrections reduce the
normalized violation **0.00303985 -> 7.66268e-6 -> 0**. Displacement
**0.47432403** is within radius **0.72**; ground minimum is **1.91595e-5**, gradient
departure **0.0009999990138**, and epsilon lower ratio **0.9818966470**.

All 1,019 controlled H evaluations converge, with 6,448 separate flow and
epsilon solves. Independent full residual readback is exact, original
qualification passes, and the one-update failed-stage receipt preserves h0 and
physical history, with zero accepted physical steps/time. This completes a
useful native recovery from the source accuracy failure; it does not resolve
the original stage or substitute its saved starting guess for original h0.
[Focused source-recovery assertions](cls-epawf-segregated-picard-r3-evidence/native/source-precision-checks.json),
[execution receipt](cls-epawf-segregated-picard-r3-evidence/native/source-precision-direction-check-receipt.json),
[completed histories and independent audits](cls-epawf-segregated-picard-r3-evidence/native/followup-completed-review.json).

## Common-reference accuracy for cancelling directions (41ef144)

The first historical-response correction is now actually accepted in the
7ba50d3 diagnostic: original residual 3.9711215458 becomes 3.9637524540,
qualification violation stays zero, and displacement 0.63675849 remains within
radius 0.72 after two qualification corrections. A second historical
correction reaches 3.9603257365. This is numerical progress on the original
equations, not a physical stage or a current per-column derivative certificate.
The matched four-update control and final readbacks are still running.
[Actual coarse-reference direction checks and trials](cls-epawf-segregated-picard-r3-evidence/native/response-chord-coarse-progress.json).

Its initial chosen-direction check has two-size residual discrepancy 0.02190
and measured merit slopes -0.08554 and -0.08182. A cancelling combination can
have a much smaller signal than its individual columns. These observations
motivate a tighter shared reference before claiming reliable small-step
responses, even though the larger actual accepted residual reduction is clear.

Commit `41ef14412554f336250f2569e03b5a22beb410ac` refines the fixed-k
reference and both chosen-direction probes to the smaller of the outer H
tolerance and stage tolerance. Every actual line-search trial retains the
same declared inexact H tolerance. Failure of this optional accuracy
refinement discards the historical model and rebuilds fresh responses without
changing the accepted state. The event records the original k-residual change
caused by tightening H. Signed observation records separately retain both
sides, both signs and actual/projected epsilon observations, with response
errors also expressed in bound-width units. All original acceptance gates
remain unchanged.

The four-worker Linux build passes in 80.96s and all six focused CTests pass
in 44.24s. The unchanged helper retains its prior sanitizer evidence; the new
orchestration has no full-solver sanitizer claim. Matched verification-On/Off
native tests have started from the same audited diagnostic input. A fresh
original-h0 follow-up is gated on an actual checked correction with the
tighter reference; the coarse-reference 7ba50d3 original-h0 attempt continues
as a separate control. Native accuracy and trajectory comparisons are pending.
[Build and producing identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-41ef144.json),
[Tighter matched native commands](cls-epawf-segregated-picard-r3-evidence/native/response-chord-tight-native-launch.json),
[Coarse-reference fresh original-h0 command](cls-epawf-segregated-picard-r3-evidence/native/response-chord-original160-3step-launch.json).

## Safeguarded short-lived response reuse (7ba50d3)

Reconstructing 96--504 converged, independently checked field responses at
every numerical correction is a substantial cost in the continuing tests.
Commit `7ba50d302cad7e7c0f676879a51029a4fcdc5416` adds the opt-in
`RANSWallReducedResponseReuseUpdates` / `--reduced-response-reuse-updates`
budget (default 0, maximum 4). The default retains fresh responses at every
update. The new native comparison uses budget 2 against budget 0.

After an actually accepted, qualified, converged-H endpoint, a field-metric
Broyden update interpolates the measured changes in original k residual,
effective independent fields and projected matching trace. Prescribed k
directions remain exact. Reuse is confined to the same profile branch,
qualified bases and a fixed response space. Age, basis enlargement, branch
changes, poor agreement or failed actual progress require fresh responses.
Below residual 1e-6 it always rebuilds fresh responses. No physical state or
history is advanced by retaining this search information.

These are explicitly labelled historical secant models, not newly verified
per-column Jacobians. Before any reused direction receives larger trials,
two fresh sequential H solves check its sign and size response. The two-size
residual/observation tolerances remain 0.03/0.05. The historical-model errors
must be at most 0.15/0.1, respectively, with residual normalization by
`max(1, measured mass-dual response norm)` and observation normalization by
`max(1, abs(measured response))`. Both measured merit slopes must be negative
and the probe branch must agree. A failed check rebuilds fresh responses at
the same accepted state. Source-refined directions receive the same mandatory
check. These checks operate independently of the observation/verification flag.

Actual trials still converge H, refresh all dependent data and scan the
original complete residual, ground polynomial and qualification inventory.
They retain the actual displacement radius. A model interpolation defect or
small quadratic KKT defect is not reported as original reduced stationarity:
the independent first-order auditor explicitly refuses historical secant
models. The focused native checker requires a freshly checked, actually
accepted original-residual reduction and tests that refusal.

The four-worker Linux build passes in 128.53s and all six focused CTests pass
in 45.79s with concurrency 1. The unchanged small model helper retains the
fd8c034 ASan/UBSan evidence; there is no new full-solver sanitizer claim.
An incorrect build-launch commit label was caught before compilation; the
corrected producing revision, failed guard log and successful build identities
are retained separately.
[Build and commands](cls-epawf-segregated-picard-r3-evidence/native/build-validation-7ba50d3.json).

Matched four-update enabled/disabled native tests are running from a freshly
audited **3.9946700988643817** diagnostic checkpoint. Original norms reproduce
bit for bit; complete qualification and the negative-ground guard pass.
They use the same executable, original h0/history, an initial 96-response
budget, 16 Krylov seeds, radius 0.4 and H tolerance 2e-9. Their four-update
comparison and cost remain pending. A fresh original-h0 follow-up has now
started after an actually accepted checked reuse correction. It uses 160
responses, 64 Krylov seeds, cap 4, three source refinements and requests three
physical steps; the fresh fd8c034 attempt continues as a separate control.
No speedup, root or physical acceptance is inferred.
[Matched native commands and initial independent audit](cls-epawf-segregated-picard-r3-evidence/native/response-chord-native-launch.json).

Mac input preparation emits the new reuse budget correctly and rejects a
negative budget, a budget above four and reuse without the reduced solver.
These are input checks, with no native numerical conclusion.
[Input preparation and guards](cls-epawf-segregated-picard-r3-evidence/native/response-chord-input-checks.json).

## Full-space curvature and source-stencil accuracy (fd8c034)

The full 504-direction experiment at the qualified 5.18197 checkpoint has
rank 504 and maximum column correlation 1.68e-16. Its original reduced merit
gradient has norm 1432.08. A separately computed cone defect using 58
near-active rows is 1348.71, compared with empirical two-size gradient
uncertainty 2.66. This base is not established to be a constrained stationary
point. These are measurements at this base, not a certificate about later
states or global existence.
[Independent first-order audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-qualified-full504-at5p182-iteration0.json).

The chosen linear-model direction predicts residual 4.00522. Its signed
original reduced-merit slope is -9.57888, versus measured half-step slope
-9.58109; both response solves meet H tolerance 2e-12. At the full proposed
length, however, the nonlinear destruction surrogate predicts ascent to
5.76306. The guarded ray chooses fraction 0.353927: its prediction is
4.83016 and the actual original residual is **4.857816821924814**. Qualification
passes and displacement 0.254314 is within radius 0.72. Thus locally useful
reduced feedback survives an actual original-equation trial; large finite-step
curvature, rather than a vanishing true reduced response, limits this update.

All 504 surrogate coordinate derivatives are available, but their maximum
two-size error 0.0382067 exceeds the unchanged 0.03 limit. The source-direction
refinement therefore aborts. This is a separate internal surrogate accuracy
failure; it is not a failure of the independently checked original reduced
actions or a proof of infeasibility.
[Direction, signed checks, failed surrogate refinement and actual trial](cls-epawf-segregated-picard-r3-evidence/native/qualified-full504-at5p182-iteration0-review.json).

Commit `fd8c0344bbbe3ea4a180df7358d0235ee8b675a4` retries an inaccurate
surrogate column with up to three smaller stencils. Each retry reuses the
previous half-step and needs one new guarded quarter-step evaluation. It
retains the 0.03 accuracy bound, the incumbent anchor and all ground/domain
guards. An unresolved derivative remains unavailable. The existing
`4*n+64` model-evaluation budget now applies with either anchor policy.
Per-column records retain sizes and initial/final errors. These refinements
affect the internal destruction model only; actual candidates still converge
H and evaluate the complete original nonlinear equations and qualification.

The four-worker Linux build passes in 86.30s, all six focused CTests pass in
43.10s with concurrency 1, and the changed standalone helper/tests pass
ASan/UBSan in 10.04s (30.39s compiler invocation). Independent quadratic
values test recovery of a resolvable truncation error; a discontinuous example
must remain unavailable after all three retries. The sanitizer scope is the
small helper and tests, not the full solver/dependency build.
[Producing image, libraries and commands](cls-epawf-segregated-picard-r3-evidence/native/build-validation-fd8c034.json).

Two new native tests are running. A one-update full504 diagnostic repeats the
same 5.18197 input, H tolerance 2e-12, probe sizes and original target, with
interior anchors disabled. Its focused assertion requires recovered accuracy
to produce an independently checked, actually accepted refined direction
that lowers the original full residual. A fresh original-h0 primary attempt
requests three physical steps, using 160 responses, 64 Krylov seeds, cap 4
and three source refinements. Neither pending run is counted as convergence.
[Exact native commands and allocations](cls-epawf-segregated-picard-r3-evidence/native/source-precision-native-launch.json).

Separately, a fresh fd8c034 probe reproduces every original norm of the
continuing fa5ac96 checkpoint **4.042756543242138** bit for bit. Complete
qualification passes: ground minimum 1.33977e-6, gradient departure
0.00099999888562 at x=-1/3, z=2, and epsilon lower ratio 0.9826997882.
Covariance/lambda and the deliberate negative-ground guard pass. This is an
in-flight numerical checkpoint, not a completed physical stage/history receipt.
A separate 4.11828 checkpoint has also started a tighter full504 diagnostic
using 9f65eed; its reduced first-order audit is pending.
[Fresh 4.04276 checkpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-previous-warm-current.json),
[Tighter full-space diagnostic launch](cls-epawf-segregated-picard-r3-evidence/native/qualified-full504-latest-9f-launch.json).

## Recovery after a failed radius-free model query (9f65eed)

The 6b1f28a original-h0 control reaches residual 24.6865, then repeatedly
exhausts 500 active-set iterations in its zero-multiplier quadratic query.
The 96-, 99- and 102-direction responses are accurate and full rank within
those spaces. Their model stationarity defects remain large; these are
unfinished model solves, not constrained stationarity of the original PDE.
The control eventually ends without accepting a physical stage.

Its independent terminal original residual is 24.68650069522063. All 2,002
non-k solves converge, with 8,821 flow and 8,821 epsilon solves. Accepted h0,
physical history and time remain unchanged. Tightening H after the last
numerical update changes the reported block norms by at most 8.41e-10; the
fresh terminal vector is retained separately as the authority rather than
claiming a bitwise match to that earlier logged point.

A standalone replay uses the production LAPACK factorization pattern and
identical libraries. At 500 iterations it reproduces the native failure
exactly, including stationarity defect 24.1012519835. A 2,000-iteration control
resolves the radius-free problem after the ball-certified cut filter; its
minimizer has displacement 10.8216, far outside radius 1. The filtered model
outside that ball is not a physical admissibility or root-existence test.
A positive-multiplier query resolves a better conditioned small problem.

Commit `9f65eedf2b59429ce3c04de0672064ad459baceb` therefore permits the bounded
solver to try positive multipliers after an unsolved zero query. An unsolved
query supplies neither a norm nor a lower bracket. Feasible-start failures
remain failures; exhausted positive queries remain explicit failures. The
existing 14-query bracket budget, original cuts and actual displacement bound
remain in force. New records distinguish an unavailable radius-free solution
from a solved bounded model. Governing equations, field responses and original
nonlinear acceptance are unchanged.

The repaired 500-iteration replay returns **exactly the same coefficient
vector** as the 2,000-iteration control. Independent checks of all 193,692
original cuts give primal defect 3.82e-18 and displacement 0.9998250752 within
radius 1. Model stationarity is 9.65e-13; its inexact radius complementarity is
0.00943. Total active-set steps are 697 versus 2843, and elapsed times 8.87s
versus 34.05s for this one helper comparison. This establishes bounded-model
recovery, not nonlinear progress, stationarity or M8 convergence.
[Matched model replays and independent checks](cls-epawf-segregated-picard-r3-evidence/native/radius-query-affine-review.json).

The four-worker Linux build passes in 81.80s, all six focused CTests pass in
42.73s with concurrency 1, and the changed standalone model passes ASan/UBSan
in 9.83s (28.85s compiler invocation). New known-answer tests exercise a failed
zero query, two failed positive queries, recovery of the independent bounded
answer, and explicit failure when all 14 positive queries fail. The sanitizer
scope remains the standalone model and tests.
[Build and immutable executable/library identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-9f65eed.json).

The eight-update cap-1 original-h0 fixture matching the failed control is now
complete. The cap-4, 48-update-per-stage three-step request continues with its
16-hour process limit.
[Native commands](cls-epawf-segregated-picard-r3-evidence/native/radius-query-native-launch.json).

The focused native recovery assertion passes. Initial state and all first six
numerical updates match the failed control bit for bit. At model 6, the zero
query again exhausts 500 iterations, but the positive-multiplier solve returns
a bounded direction with radius 1 and displacement 0.9998250752. Its measured
original reduced-merit half-step slope is -152.82904 versus prediction
-152.83393; relative half-action error is 0.0001155. After three qualification
corrections, actual displacement is 0.9998932922 and the original full residual
falls from 24.68650 to **23.1739120920**. The eighth update reaches
**21.102638703378357**. This is actual original-equation recovery beyond the
old stalled numerical solve, not a result inferred from its small quadratic.

Every original terminal norm reproduces exactly in a separate process.
Qualification passes: continuous ground minimum 0.0007594564, gradient
departure 0.00099999989355 and epsilon lower ratio 0.9800007655. All 1,588 H
evaluations converge, with 6,121 velocity/pressure and 6,121 epsilon solves;
target elapsed time is 9220.23s. The original 2e-12 gate still fails. Native
accepted steps/time remain **0/0**, and accepted state and physical history
are preserved. The passing test is explicitly named as bounded-model/operator
recovery, not M8 convergence.
[Focused native assertions](cls-epawf-segregated-picard-r3-evidence/native/radius-query-direction-checks.json),
[Check execution receipt](cls-epawf-segregated-picard-r3-evidence/native/radius-query-direction-check-receipt.json),
[Complete operator/history/readback evidence](cls-epawf-segregated-picard-r3-evidence/native/completed-original-and-expanded-review.json).

A further independently audited diagnostic checkpoint reaches original full
residual **4.253091674121031** in the continuing 160-direction previous-step
search. Every original norm matches the saved point bitwise. Complete
qualification passes: ground minimum 6.04969530e-6, gradient departure
0.00099999934768 at x=-2.625, z=4, and epsilon lower ratio 0.9823238740.
Covariance/lambda and the negative-ground guard pass. This is an in-flight
checkpoint, with no final native stage/history receipt or physical acceptance.
[Fresh checkpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-previous-warm-third.json).

## Configurable measured radius growth (be04686)

Several original-h0 models repeatedly hit the scaled k-displacement radius 1,
with positive radius multipliers and unconstrained proposals between 5.55 and
16.57. Eleven recorded updates have actual displacement above 0.8 and measured
actual/predicted gain above 0.75, yet retain radius 1 because of the hard cap.
These observations justify testing more growth; they do not establish that a
larger step will be useful or that a qualifying root exists.
[Measured cap-limited controls](cls-epawf-segregated-picard-r3-evidence/native/radius-cap-motivation.json).

Commit `be04686bb3dac8987b517f380c8a60d075f7b940` adds
`RANSWallReducedMaximumRadius` / `--reduced-maximum-radius`. Its default is 1,
retaining the historical control. A larger configured cap allows the existing
agreement test to grow the radius; poor agreement still shrinks it. The
configured radius continues to bound the actual proposed correction, including
constraint handling. Original residual, H accuracy, ground checks, complete
qualification, model and physical history are unchanged.

A Linux four-worker build passes in 125.52s, all six focused CTests pass in
43.40s with concurrency 1, and the changed standalone model test passes
ASan/UBSan in 9.42s (23.83s compiler invocation). A new identity-model test
checks the known constrained correction (2.5,0), ensuring that a radius above
1 remains a field-displacement bound. The Python native checker also asserts
the configured cap on every model radius and subsequent radius update. Input
preparation and three invalid-radius guards pass separately on the Mac; those
are input checks, with no solver execution.
[Build, commands and executable/library identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-be04686.json).

Native tests use the same be04686 executable and original h0: a four-update
cap-1 control and a cap-4 primary run with 48 numerical updates per stage and
three requested physical steps. Both start at radius 0.4, with 160 responses,
64 Krylov seeds, one source refinement and H tolerance 2e-9. The primary
process has a 16-hour limit, chosen from measured model costs; the control has
a four-hour limit. Run budgets differ, so matched claims concern the common
prefix only. The longer cap-4 test remains in progress; its four-update
checkpoint and the completed cap-1 control can now be compared as follows.
[Exact native commands](cls-epawf-segregated-picard-r3-evidence/native/radius-cap-native-launch.json).

| Original-h0 numerical update | Cap 1 control | Cap 4 experiment |
| --- | ---: | ---: |
| 1, qualification restoration | 63.74366608 | 63.74366608 |
| 2 | 49.09361202 | 49.09361202 |
| 3 | 43.31318662 | 40.81722751 |
| 4 | 37.52204608 | 29.99279194 |

Both use the same immutable executable/libraries. Their initial state and
first two updates reproduce bit for bit. The numerical XML differs only in
the declared maximum radius and update budget. The cap-1 control completes
with zero accepted physical steps/time, accepted state/history preserved and
an exact fresh original norm readback. All 1,310 H solves converge, with
4,804 velocity/pressure and 4,804 epsilon solves; target elapsed time is
7942.69s. The four-update cap is a control budget, not a convergence result.

The cap-4 fourth correction has actual displacement **1.2947276287**, within
radius **1.296** after three qualification corrections. A separate process
reproduces every original norm at that checkpoint bit for bit and passes
complete qualification and the negative-ground guard. Ground minimum is
0.0203734663, gradient departure 0.00099999900209, and epsilon lower ratio
0.9806180819. The local outer k ratios 0.99000001--1.00999999 also pass.
This establishes useful native use of a radius above one. It is an in-flight
numerical checkpoint, not a qualified algebraic root or completed native
physical/history receipt for the cap-4 run.
[Matched comparison and independent audit](cls-epawf-segregated-picard-r3-evidence/native/radius-cap-four-update-comparison.json),
[Completed cap-1 operator/history/readback review](cls-epawf-segregated-picard-r3-evidence/native/completed-original-and-expanded-review.json).

## Interior source-model anchors (900a9aa)

The continuing 160-direction search exposes a separate limitation in the
nonlinear destruction surrogate: at iteration 2 its best ray proposal has
zero available coordinate derivatives, and the refinement reports
`unresolved surrogate derivative`. The original reduced actions and actual
trial remain valid; this abort concerns the internal curvature model. A
constrained corner can forbid both signs of a coordinate perturbation.
The old rejection log does not distinguish the ground and scalar-inventory
families, so the exact rejection family is not inferred from that message.

Commit `900a9aabc58f368f7f61a0aa95e56f834551bc76` adds the opt-in
`RANSWallReducedInteriorSourceAnchor=On` selector, exposed by
`--reduced-interior-source-anchor`. The disabled selector retains the control.
After a blocked or inaccurate model derivative, it retries the linearization
at one half and one quarter of the incumbent coefficient vector. Every model
sample still passes the actual ground/domain guard. All derivatives must pass
the existing two-size error bound. There are at most two anchor reductions
and `4*n+64` additional model evaluations per refinement.

The affine constant uses the new anchor, while the incumbent proposal stays
unchanged. A replacement must improve its surrogate residual and regularized
objective, satisfy the original frozen cuts and displacement bound, and have
negative original reduced-merit slope. Every actual trial still equilibrates
H and scans the original complete residual, ground extrema and qualification.
The anchor is an internal model point; no physical state, history, source,
acceptance envelope or final 2e-12 tolerance is changed.

The focused C++ test uses two active affine positivity cuts for which neither
sign of the first coordinate probe is defined. It verifies rejection at the
corner, recovery at the interior anchor against independent analytic
derivatives, preservation of an already valid anchor and explicit exhausted
retry status. The Linux build with four workers passes in 117.91s, all six
focused CTests pass in 38.48s, and the changed standalone model passes
ASan/UBSan in 8.91s (24.33s compiler invocation). This sanitizer scope covers
the small model helper and tests, not the full dependency/solver build.
[Build, commands and immutable executable/library identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-900a9aa.json).

Matched enabled/disabled native tests and a verification-Off repeat use the
same 900a9aa executable, physical target and 7.13029 saved guess. Each permits
two numerical updates. The focused native assertion requires recovered
derivatives to produce an independently checked, actually accepted refined
direction, with the original operator and history assertions. A separate
original-h0 run enables the method with 160 responses/64 Krylov seeds, eight
numerical updates per stage, an eight-hour process limit, and three requested
physical steps. The original-h0 test remains in progress.
[Exact native commands](cls-epawf-segregated-picard-r3-evidence/native/source-anchor-native-launch.json).

The matched two-update native comparison is complete:

| Source-anchor policy | Update 1 | Update 2 | Qualification | Physical steps |
| --- | ---: | ---: | --- | ---: |
| Off | 5.7749372899 | 5.2418504922 | Pass | 0 |
| On | 5.7427345303 | 5.2018304001 | Pass | 0 |

Both models recover all 96 surrogate derivatives at the half-ray anchor, with
201 model evaluations against a 448-evaluation budget. Two-size derivative
errors are 0.00718 and 0.01195. The resulting directions have original reduced
merit slopes -20.7916 and -6.66947, versus measured half-step slopes -20.7795
and -6.67619. Both refined directions are actually accepted. The final source
prediction is 5.18600872 versus actual 5.20183040; vector error is 0.22738,
compared with 0.94178 for the original linear prediction. Qualification
curvature still requires two actual corrections per update. The final
correction has displacement 0.71464605 within radius 0.72. Passing the source
model checks alone is not its acceptance criterion.

All 403 H evaluations converge to their requested accuracies, no larger than
2e-9; 1690 velocity/pressure and 1690 epsilon solves are recorded. All 192
retained true response checks pass, with maximum residual/observation errors
0.002129 / 0.022650. The first 96 response records match the disabled control
bitwise apart from timing. Target-process elapsed times are 3146.63s enabled
and 2497.31s disabled under the recorded concurrent allocation; the modest
residual improvement is not a convergence or performance claim.

The enabled endpoint has ground minimum 1.34386009e-5, maximum gradient
departure 0.00099999911131 and epsilon lower ratio 0.9801652658. Fresh readback
reproduces every original block norm bitwise and passes complete
qualification, covariance/lambda and negative-ground checks. Native rejection
preserves accepted h0, every physical history slot and time. Verification On
and Off also produce bitwise identical initial and both numerical states,
original residual vectors, rejected endpoint and physical history.

`check_cls_epawf_source_anchor.py` passes its explicit requirement for an
independently checked, actually accepted direction from recovered derivatives;
its status is `passed_source_anchor_direction_checks_not_m8_convergence`.
The verification comparison passes separately. The original full residual
remains 5.201830400141313, so no physical stage is accepted.
[Focused direction assertions](cls-epawf-segregated-picard-r3-evidence/native/source-anchor-direction-comparison.json),
[verification equality](cls-epawf-segregated-picard-r3-evidence/native/source-anchor-verification-comparison.json),
[command receipts](cls-epawf-segregated-picard-r3-evidence/native/source-anchor-direction-comparison-receipt.json),
[fresh endpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-source-anchor-enabled.json),
[compact completed histories](cls-epawf-segregated-picard-r3-evidence/native/completed-original-and-expanded-review.json).

Separately, a full 504-response search starts at the independently audited
5.18197 point, with H tolerance 2e-12, probe step 5e-6 and the learned radius
0.72 (preceding actual/predicted gain ratio 0.9023). Its two-update budget has
an eight-hour process limit. It retains the original h0/history and physical
target. This tests the remaining response-space limitation; it has no original
root or physical-acceptance certificate yet.
[Full-space command and radius evidence](cls-epawf-segregated-picard-r3-evidence/native/qualified-full504-at5p182-launch.json).

The continuing canonical 160-direction run has reached numerical update 5,
with original full residual **4.396589811578443**. A new 900a9aa probe process
reproduces every original block norm bitwise. The complete qualification audit
passes: ground minimum 1.17275968e-5, gradient departure 0.00099999900310,
epsilon lower ratio 0.9816889339, and covariance/lambda checks. The gradient
witness has moved to x=-2, z=4. The negative physical-ground guard also passes.
This remains an in-flight checkpoint with no completed native stage/history
receipt and no physical acceptance.
[Fresh checkpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-canonical-warm-fifth.json).

A further 900a9aa saved-guess search starts from this checkpoint with three
source refinements, H tolerance 2e-10 and its learned radius 0.72. This is a
one-stage development diagnostic; the original-h0 tests remain separate. An
initial invocation requested multiple stages and was rejected by the saved-
guess argument guard before any native process or state was created. That
invocation is preserved and the corrected one-stage run is in progress.
[Corrected command](cls-epawf-segregated-picard-r3-evidence/native/source-anchor-warm160-three-refinements-launch.json),
[guarded pre-native invocation](cls-epawf-segregated-picard-r3-evidence/native/source-anchor-warm160-three-refinements-failed-invocation.json).

## Latest implementation and native checks

Code `6b1f28a9b316713ba1d06db5c794bd1e61eabb15` fixes another response-budget
problem: the previous successful k correction was queued after all504 free
modes. The96/160-direction searches exhausted their budgets before reaching
it. It is now queued after the conditional/Krylov and rate seeds. Its current
reduced response must still pass independent feasible two-size H solves and
rank checks before retention. This changes direction construction without
reusing an old action as a current derivative or changing any physical row.

The Linux build passes with8 compiler workers; all six focused CTests pass in
37.96s with test concurrency1. Bounds assertions remain enabled. The standalone
reduced-model header/test are unchanged from389d0a5, so their recorded ASan/UBSan
evidence is reused only for those files; it does not sanitize the changed
Picard driver or the full dependencies.
[Build, commands and binary/library identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-6b1f28a.json).

Diagnostic-only commit `fa5ac96b5b396a759002a104d2c0d7cf88bb6b3b` replaces the
stale protocol literal "through160" with the compiled `DirectFaceThreshold`
and records its numeric value64. The actual face threshold was already64 in
cea2a06 and all subsequent numerical runs; their historical text label was
incorrect. Linux compilation with2 workers passes in56.93s and the unchanged
reduced-model CTest passes in.97s. The five shared native regression results
are reused from6b1f28a; the recorded numerical-source diff changes only these
diagnostic fields. The independently tested numerical binaries below retain
their producing revisions and immutable identities.
[Diagnostic-only build and exact source diff](cls-epawf-segregated-picard-r3-evidence/native/build-validation-fa5ac96.json).

The matched two-update saved-guess tests are complete. They have identical
physical input and numerical settings; the initial state, first numerical
update and all 96 first-model response records agree bitwise (elapsed times
are excluded). Both second models have rank 96.

| Code / previous-direction policy | First update | Second update | Final qualification | Physical steps |
| --- | ---: | ---: | --- | ---: |
| 389d0a5 / queued after the free-mode pool | 5.7749372899 | 5.2426008180 | Pass | 0 |
| 6b1f28a / queued before the pool | 5.7749372899 | 5.2418504922 | Pass | 0 |

The modest difference supports retaining useful history in the finite response
budget; it does not establish uniform improvement or robust convergence. The
new previous direction has independent component 0.44359 and passes two-size
residual/observation checks with errors 0.0001122 / 0.01505. Its signed merit
response is +60.0394 versus measured half-step +60.0418, so it is information
for the constrained search rather than a blindly repeated forward step. The
constructed combined direction instead predicts slope -7.16691 versus measured
-7.17103. Every one of the 400 H evaluations converges to its requested
accuracy, at most 2e-9; 1656 flow and 1656 epsilon solves are recorded. All 192
retained response checks pass, with maximum residual/observation errors
0.001890 / 0.022650.

The final source prediction is 5.226529677 versus actual 5.241850492; vector
error 0.262808 is smaller than 1.062896 for the linear model. Three actual
qualification corrections are used. Final displacement 0.713299592 stays
within the adapted radius 0.72. Ground minimum is 2.13234591e-5, maximum outer
gradient departure 0.00099999904294 (x=-2.125, z=2.66667), and minimum epsilon
ratio 0.9801856994. Independent readback reproduces every original block norm
exactly and passes complete qualification, covariance/lambda and the negative
ground guard. Both rejected native stages preserve accepted h0, every physical
history slot and physical time exactly. They fail the original residual gate;
that delivery's best completed qualified saved-guess result is 5.241850492167935.
The later source-anchor comparison above improves it to 5.201830400141313.

`check_cls_epawf_previous_direction.py RUN --control CONTROL` passes its focused
assertions for retained current responses, identical physical input, identical
first model/update, original operator checks and rejected-stage history. This
is explicitly named an operator/comparison test, not M8 convergence.
[Matched focused test](cls-epawf-segregated-picard-r3-evidence/native/previous-direction-matched-review.json),
[complete compact histories](cls-epawf-segregated-picard-r3-evidence/native/completed-original-and-expanded-review.json),
[fresh original endpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-previous-pair-new.json).

A separate 6b1f28a native original-h0 run requests three physical steps with 48
numerical updates per stage and an 8h process limit. It remains in progress;
its prefixes are not completed native receipts.
[Exact native launch commands](cls-epawf-segregated-picard-r3-evidence/native/previous-direction-native-launch.json).
A matched verification-Off repeat completes with the same 6b1f28a binary,
physical data and two-update budget. The initial state and both accepted
numerical-update state hashes and original norm vectors agree bitwise. Accepted
physical state, every native history slot/time, rejected endpoint, accepted
state list and native bookkeeping are also identical. The existing
`test_cls_epawf_reduced_response.py VERIFIED --compare-unverified UNVERIFIED`
check passes. This verifies observation independence for this fixture and
retains the original failed-stage outcome of zero physical steps.
[Verification-Off launch](cls-epawf-segregated-picard-r3-evidence/native/previous-direction-observability-launch.json).
[Full On/Off comparison](cls-epawf-segregated-picard-r3-evidence/native/previous-direction-verification-comparison.json),
[successful comparison command](cls-epawf-segregated-picard-r3-evidence/native/previous-direction-observability-comparison-receipt.json).

The qualified 5.24185 endpoint is being followed with 160 responses and 64
Krylov seeds using fa5ac96, allowing eight numerical updates within a six-hour
process limit. The exact audited readback packet supplies only the nonlinear
guess. All physical history values, alpha 400, time 0.0025 and original
acceptance criteria remain fixed. A separate latest-code original-h0 run uses
the same measured response configuration, requests three physical steps, and
allows twelve numerical updates per stage within eight hours. Both are
pending; the saved-guess test does not replace the original-h0 test.
[Saved-guess extension command](cls-epawf-segregated-picard-r3-evidence/native/previous-warm160-extension-launch.json),
[latest-code original-h0 command](cls-epawf-segregated-picard-r3-evidence/native/current-original160-native-launch.json).

The saved native numerical packet and the accepted-h0 readback packet have
identical state bits and all 2165 history values match exactly. Their history
hashes differ at twelve signed zeros in the native callback's computed h0
combination. The runner checks the callback's normalized-value hash and checks
accepted physical state/history slots with raw IEEE hashes. The follow-up uses
the exact independently audited readback packet; its byte hash is recorded.
[Packet/history equivalence check](cls-epawf-segregated-picard-r3-evidence/native/previous-packet-history-equivalence.json).

The389d0a5 canonical160-direction warm search makes a useful qualified numerical
update from7.5004290706 to7.1302877208629. The actual continuous-ground minimum
increases from6.2733e-7 to5.32264357e-5; minimum epsilon ratio increases to
.9812197154. Maximum gradient departure is.0009999998290 with witness x=1,z=2.
The full scalar evaluation inventory, complete qualification, covariance and
lambda checks pass. A fresh process reproduces every original block norm
exactly and rejects the deliberately negative continuous-ground test before
scalar evaluation. This remains an ongoing numerical checkpoint with no final
native physical-history receipt, and is not an algebraic root or accepted step.

The160 field probes have maximum correlation9.98e-17; the mass-weighted
response-action rank is160. The
model has24 active constraints and11 added ground cuts. Its combined local
predicted merit slope-29.63178 agrees with measured half-step-29.60851. The
nonlinear source refinement uses all160 directional derivatives; its local
prediction-18.92105 agrees with measured-18.90670. After two actual qualification
corrections, the source prediction is7.328295228 versus actual7.130287721, with
residual-vector error.491810 versus3.212774 for the linear model. Actual final
displacement.212164569 is within radius.4; attained H error1.0344e-9 is below
the declared2e-9. The logged next-witness inventory is checked for canonical
ordering, one observation per family/spatial region, and membership in the
active physical observations.
[Compact checkpoint/model/trial review](cls-epawf-segregated-picard-r3-evidence/native/canonical-warm-first-checkpoint-review.json),
[independent original-operator audit](cls-epawf-segregated-picard-r3-evidence/native/audit-canonical-warm-first.json).

Its next numerical update reaches5.5688705520357455 with complete original
qualification still passing. Ground minimum is2.34885008e-5, maximum gradient
departure.00099999902235 (witness x=-2.60465,z=4), and minimum epsilon ratio
.9807155533. The source prediction5.562813859 differs from the actual norm by
.00606; vector error.035155 is smaller than the linear model's.488066. The
refined direction predicts local merit slope-17.85791 versus measured
half-step-17.85391. Actual displacement.291132591 stays within radius.4 and
attained H error1.5031e-9 is below2e-9. Independent6b1f28a probe readback is
bitwise identical to every original block norm and passes complete
qualification/ground/covariance/lambda and the negative-ground guard. This
lower-residual checkpoint also remains a qualified nonroot; its run continues
and has no final native physical-history receipt.
[Second checkpoint/model/trial review](cls-epawf-segregated-picard-r3-evidence/native/canonical-warm-second-checkpoint-review.json),
[fresh second-checkpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-canonical-warm-second.json).

The third canonical checkpoint reaches original residual 5.1819703502288155,
with ground minimum 5.85989889e-5, gradient departure 0.00099999921905 and
minimum epsilon ratio 0.9812599216. A fresh fa5ac96 probe reproduces every
original block norm exactly and passes complete qualification and the negative
ground guard. This is the best independently audited ongoing checkpoint; its
run still has no completed native physical-history receipt.
[Fresh third-checkpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-canonical-warm-third.json).

The dated 2026-11-06 18:24 UTC snapshot records ongoing original-h0 progress:
the older 90dd801 run has 15 numerical updates and residual 16.2513885503;
6b1f28a has four updates and residual 32.6605130307. Both checkpoints pass
qualification. The latest fa5ac96 original160 run is still constructing its
first restoration model, with its initial equilibrated-H residual 63.0021188124
outside qualification. None has a completed native first-stage receipt.
The latest completed original-h0 result remains 19.1747990683 below. The
snapshot also distinguishes the timed-out c02f965 run, which has no final
native history receipt, from the completed matched tests.
[Dated in-flight histories and executable identities](cls-epawf-segregated-picard-r3-evidence/native/progress-published-20261106-evening.json).

At the preceding7.5004290706 full504-response base, the true reduced gradient
norm is3172.69 and the near-active cone defect is2847.89, compared with empirical
two-size gradient uncertainty9.54. All504 response coordinates are independent.
This tested base is not shown stationary. It supports continued search and
does not establish global root existence, impossibility, or stationarity at
the later7.13029 checkpoint.
[Full-space original-merit first-order audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-qualified-full504-at-7p50043.json).

## Latest completed full-space search

The strict-accuracy `d710404` full504-response run completes two numerical
updates from the audited16.4201949282 diagnostic guess, reaching11.5210784187
then7.498542874129245. The physical alpha400/time.0025/h0/history and original
discrete problem remain fixed. The original complete residual still fails
2e-12. Its other block norms are[9.30e-13,1.36e-12,4.90e-14,6.43e-16,
8.06e-16,1.78e-17,0], leaving k dominant.

The final continuous-ground minimum is5.78454706e-7, maximum gradient
departure.0009999990012 and minimum epsilon ratio.9800289529. The independent
process reproduces the original norm vector exactly and passes complete
qualification, covariance/lambda and the negative-ground guard. Native
physical steps/time are0/0; accepted h0 and every history slot/time are
preserved. This completed saved-guess result does not replace original-h0
validation, whose latest completed run remains19.1747990683 as reported below.

All2033 controlled H evaluations converge to their declared accuracies, at most
2e-12;12272 velocity/pressure pairs and12272 epsilon solves are recorded.
All1008 retained two-size actions pass, with maximum residual/observation
relative differences.0004717/.009930. Both models have rank504. The last source
prediction is7.572041303 versus actual7.498542874, with residual-vector error
.272731 versus3.838011 for the linear model. The actual final normalized
displacement.399969604 lies within radius.4. Reduced-solver log time is20032.17s.
Further full-space and original-start tests continue; this is a qualified
nonroot, not a stationarity or root-existence certificate.

[Complete compact native history](cls-epawf-segregated-picard-r3-evidence/native/completed-original-and-expanded-review.json),
[fresh original endpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-filtered-search504.json).

## Completed matched Krylov follow-up

The `90dd801` matched searches start from the independently qualified numerical
checkpoint at16.4201949282, with the original alpha400/time.0025/h0/history.
They change only `RANSWallReducedKrylov` (16 versus64) in otherwise identical
session files. Executable/library identities, initial coefficients/norms and
the first17 response checks are identical; elapsed times are excluded from
the response comparison. Both models have96 independent directions.

| Krylov directions | First numerical update | Second numerical update | Final ground minimum | Final gradient departure | Final minimum epsilon ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 | 13.4942429683 | 9.3285649274 | 2.39380747e-4 | .000999999039 | .9817519306 |
| 64 | 11.7713161888 | 9.9029254238 | 1.55559381e-4 | .000999999049 | .9803458666 |

Every accepted numerical update passes the original qualification and ground
checks. Both final norm vectors reproduce exactly in a fresh process; the
complete covariance/lambda audit and negative-ground guard pass. Both native
stages fail the original2e-12 residual gate and preserve accepted h0 and all
physical-history slots/time. Physical steps/time remain0/0. More Krylov
directions help the first update but do not dominate after two updates.
These saved-guess experiments do not replace original-h0 validation.

All402 H evaluations in each search attain their declared tolerances, at most
2e-9. The16-seed run records1649 velocity/pressure pairs and1649 epsilon solves;
the64-seed run records1673 of each. All192 retained actions in each run pass
two-size checks. Maximum retained residual/observation relative differences
are.0008106/.017623 and.009291/.026122, respectively. The64-seed run rejects one
additional response. The independent current-operator branch/combined checks
remain enabled. These are inexact search responses; final acceptance retains
the original complete2e-12 gate.

The signed checks also show why feedback matters when choosing directions.
At this16.4202 base, the x=-2.6 sampled-K seed has logged conditional merit
slope+285.059, reduced prediction-1685.636 and measured half-step slope-1685.215.
The independent x=1.6 seed instead has reduced prediction+94.978 and measured
half-step+95.142. Each direction therefore supplies different residual and
constraint information. The epsilon-lower seed predicts a signed normalized
lower-bound change-3.14410560e-4 at x=-2,z=2 versus measured-3.14410224e-4.
These reported probes preserve the observed profile branch selections; field
and moving-panel quantities are rebuilt at both perturbation sizes.

Finite-step curvature remains substantial. On the second accepted update the
nonlinear destruction surrogate predicts9.240817519 versus actual9.328564927
for16 seeds, and9.975535343 versus9.902925424 for64. Corresponding residual-vector
errors are.641892 versus3.184512 for the linear model, and.186829 versus6.007674.
The latter correction has actual normalized displacement.719296630 within
radius.72. The former has displacement.399668275 within radius.4.

The full qualification scan is necessary: the16-seed first full proposal has
predicted violation near zero but actual normalized gradient violation.165945
at x=-2.6666667,z=2, together with an epsilon lower-bound violation. The
backtracked, re-equilibrated half proposal passes the complete inventory.
Thus the signed response model is useful locally, while finite-step
qualification predictions still require actual scans and safeguards.

[Matched-input, signed-action and witness checks](cls-epawf-segregated-picard-r3-evidence/native/response-krylov96-matched-review.json)
are reproduced by
`solvers/IncNavierStokesSolver/Tests/RANS/check_cls_epawf_krylov96_comparison.py RUN_ROOT OUTPUT`.
This specifically named check verifies two completed qualified nonroots; it
does not certify M8 convergence. The compact completed-history collection below
also includes both native receipts, commands and independent full audits.
Both endpoints are being followed with additional numerical updates, preserving
the original physical history and stage target;
[16-seed launch](cls-epawf-segregated-picard-r3-evidence/native/krylov16-extension-launch.json),
[64-seed launch](cls-epawf-segregated-picard-r3-evidence/native/krylov64-extension-launch.json).

Separately, the running `ce7e53a` full504-response follow-up reaches a second
numerical checkpoint at7.5004290706. An independent process reproduces its norm
vector exactly and passes complete original qualification: ground6.27326187e-7,
gradient.0009999990032, minimum epsilon ratio.9800480092. This ongoing run has
no final native stage/history receipt yet; the audited checkpoint remains a
qualified nonroot. It continues with further full-space updates.
[Independent checkpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-warm-full504-update2.json).
Independent full-space first-order audits at the earlier16.4202 and11.5216
bases include the original qualification cuts. Their near-active cone defects
are12411.90 and1947.16 versus empirical two-size gradient uncertainties12.82
and16.21, respectively. Both models have rank504. These earlier bases are not
shown stationary; this does not certify the later7.50043 checkpoint's
stationarity or global root existence.
[Qualified full-space model audits](cls-epawf-segregated-picard-r3-evidence/native/first-order-qualified-full504-20261106.json).

## Latest completed original-target validation

The0633d23 native original-h0 run completes18 numerical updates with original
block norms [9.34492e-13,1.02304e-12,19.1747990683,4.58916e-14,7.23219e-16,
7.42248e-16,4.74338e-18,0]. It fails the unchanged2e-12 complete residual test.
The final continuous-ground minimum is1.90615456e-5; maximum gradient departure
is.000999999004 and minimum outer-epsilon ratio.9846519511. A fresh process
reproduces the norm vector exactly and passes the original complete qualification,
covariance, lambda and negative-ground guard checks. The rejected stage preserves
accepted h0 and every physical-history slot/time. Physical steps/time are0/0.
Its1089 controlled H evaluations all attain their requested accuracies;7480
velocity/pressure pairs and7480 epsilon solves are recorded. All retained
responses pass the independent two-scale checks. This is a qualified nonroot.

The matched six-update extension from18.7645225166 reaches18.5778313750 with
source refinements versus18.4988245721 with ray-only selection. Both remain
qualified nonroots with exact fresh readback and preserved accepted history.
The ray-only policy expands its response space earlier: this counterexample
prevents a claim that the source refinements uniformly dominate adaptive
response enrichment.

A second matched pair starts at18.6652852921. Two updates using64 then80
responses reach18.6338419815;96 then96 reach18.6063142340. Both have full
numerical rank, pass all operator checks and original qualification, reproduce
their residuals independently, and accept zero physical steps. The final96-space
source prediction is18.60678920 versus actual18.60631423, vector error.01228865
versus.31720522 for its linear model. Its refined direction predicts merit slope
-1.21494328 versus measured half-step-1.21388370. These tests support useful
reduced directions and continued missing-space investigation, not stationarity
or successful M8 convergence.

[Compact completed histories and independent receipts](cls-epawf-segregated-picard-r3-evidence/native/completed-original-and-expanded-review.json)
retain commands, identities, H errors, signed active observations, displacement
radii, true residual histories and complete fresh audits. Larger raw payloads
remain in the linked Linux directories. The160-direction control also completes two updates, reaching17.7530797901
with original qualification passing, continuous ground1.27874280e-5 and minimum
outer-epsilon ratio.9805312098. Its second source prediction is17.75068755 versus
actual17.75307979; vector error.01986048 versus.88429921 linearly. The actual
normalized coefficient displacement.27286294 is within radius.4. Its local
predicted merit slope-20.38793030 agrees with the measured half-step-20.38184430.
All320 response columns pass two-scale checks, all658 H evaluations attain the
requested2e-12,3857 flow/epsilon pairs are recorded, and fresh readback has zero
norm discrepancy. This confirms that additional directions materially improve
this particular search; it remains a qualified nonroot. The504-space and new
large-face native controls continue and do not replace the original-h0 test.

## Latest full-space checks and execution improvements

`389d0a5` makes selection of subsequent response witnesses independent of the
active-row order returned by the quadratic solver. The retained and optimized
native controls have the same first physical active-observation set, with row
differences at most6.83e-16, but different order. The old first-seen-per-region
rule therefore picks different points and changes the second response space;
their second residuals8.61860502 and8.45610707 are not trajectory-equivalent.
Model-local cut indices also differ and are not used as cross-model physical
identities.

The new rule selects the lowest original observation index for each family
and spatial region, then orders the selected indices canonically. Successful
models record those indices and physical points. All original constraints,
response checks, physical equations and acceptance criteria remain. This
removes dependence on permutations of the same active observations; it does
not assert that different active sets or inventory topologies give identical
trajectories. A focused unit checks120 index permutations including duplicate
observations, distinct families/regions and conflicting-metadata rejection.
Six Linux regressions pass (38.46s), and the reduced-model unit passes
ASan/UBSan (8.87s).

New native tests use original h0 with96 responses/64 Krylov seeds, and a
separate160-response/64-seed diagnostic from the independently audited7.50043
checkpoint. Their outcomes remain pending. Every numerical rank has one OpenMP
and BLAS thread; the added two singleton ranks remain within the72-worker cap.
[Native ordering diagnosis](cls-epawf-segregated-picard-r3-evidence/native/active-order-native-diagnosis.json),
[build and checks](cls-epawf-segregated-picard-r3-evidence/native/build-validation-389d0a5.json),
[exact native commands](cls-epawf-segregated-picard-r3-evidence/native/canonical-native-launch.json).

`cea2a06` applies the existing checked primal-face solve, conservative
ball-inactive-cut test and feasible-face reuse above64 coordinates, including
the96-direction searches. The established path through64 coordinates remains.
This addresses measured12–28minute optimization intervals after constructing
each96-response model. New64/96 closed-answer tests exercise both paths,
obsolete-bound removal, original cut indices and the final displacement ball.
All six Linux regressions pass (37.85s); ASan/UBSan passes (8.87s).

Replaying the same native96-response model with all193747 original cuts at
radius.72 gives maximum coefficient difference1.79e-15 and model-objective
difference-8.53e-14. An independent scan of every original cut gives primal
defect1.56e-17; actual displacement.719328752 is inside the ball. The new affine
stationarity defect is2.11e-12. The replay takes134.34s with `90dd801` and12.30s
with `cea2a06`, recording17 warm queries and one successful cold retry. This is
one optimization-model comparison under concurrent jobs, not a nonlinear
stationarity or general performance result. New unchanged-h0 and matched
saved-guess native searches have launched with explicit singleton ranks and
thread limits; their outcomes remain pending.
[Build, commands and identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-cea2a06.json),
[identical-model comparison](cls-epawf-segregated-picard-r3-evidence/native/medium-face96-radius-replay-review.json),
[new native commands](cls-epawf-segregated-picard-r3-evidence/native/medium-face-native-launch.json).

The first optimized saved-guess numerical update independently reproduces its
original norm vector exactly and passes the complete qualification/ground
audit. Its residual8.83235687962 differs from the retained control8.83235687951
by1.05e-10; maximum coefficient difference is4.54e-13. All96 initial response
records match bitwise before optimization. The actual final displacement is
.344982331 within radius.72, with H error7.685e-10 below2e-9. The source model
predicts8.619766019 versus actual8.832356880; two qualification corrections
produce actual violation zero. The source/linear residual-vector errors are
.976475/2.231108. First-update log times are1448.12s and1108.16s under
concurrent jobs. This is a tolerance-consistent numerical checkpoint comparison;
both searches continue and neither has a final native physical-stage receipt.
[Actual native correction and fresh audit](cls-epawf-segregated-picard-r3-evidence/native/medium-native-first-update-comparison.json).

A separate original-h0 run now constructs all504 free-k responses from its
first model, keeping the original physical target/history. Its bounded budget
is six numerical updates per stage and eight hours, with the native three-step
horizon requested. This tests fresh-start subspace limitations directly;
its outcome remains pending.
[Full-space original-start command](cls-epawf-segregated-picard-r3-evidence/native/full504-original-launch.json).

`90dd801` extends the optional response-Krylov ceiling from16 to128; its
default remains0. The existing enrichment applies the measured non-k feedback
through the restricted conditional-k preconditioner, generating directions in
the space of `(I - Jkk^-1 Jred) q`. Each direction still needs a feasible
two-size, converged-H response; rank, action-error and feedback-above-noise
checks remain in force. This is motivated by the measured dominant pressure
feedback and the high cost of repeatedly sampling504 arbitrary field
directions. It neither assembles an all-field matrix nor changes a governing
equation. Six Linux regressions pass (39.67s) and ASan/UBSan passes (7.56s).
The matched96-dimensional two-update searches with64 versus16 Krylov seeds
complete as described above. Their extensions and an unchanged-h0 native run
with64 continue. Their outcomes are not inferred from build or operator success.
[Build and identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-90dd801.json),
[exact native commands](cls-epawf-segregated-picard-r3-evidence/native/response-krylov-launch.json).

`d710404` omits a frozen-model cut only when a conservative field-Gram bound
proves it cannot become active anywhere in the prescribed displacement ball.
The bound uses a Gershgorin lower eigenvalue estimate, reduced for roundoff.
It is disabled for restoration slack variables or a poorly conditioned field
Gram matrix. Every original frozen cut is still scanned at the final
inside-ball solution; actual candidates still receive the complete original
inventory and continuous-ground audit. The 192-coordinate model replay
has exactly zero coefficient difference, the same radius multiplier, primal
defect 1.96e-17 and affine stationarity 2.52e-12. Helper replay times are 89.83s
and 26.39s. These are optimization-model checks, not nonlinear stationarity or
native convergence. Six Linux regressions pass (37.48s), including a known
radius/active-cut answer and an ill-conditioned-field counterexample that
disables omission; ASan/UBSan passes (7.79s).
[Build](cls-epawf-segregated-picard-r3-evidence/native/build-validation-d710404.json),
[original-cut radius replay](cls-epawf-segregated-picard-r3-evidence/native/ball-filter-radius-replay-review.json).
The completed matched 192-direction native update also reproduces the entire
numerical endpoint and original norm vector bitwise, with accepted state and
physical history unchanged. It omits 151,975 of 193,721 frozen cuts within the
radius 0.4 ball, then passes the original full inventory and fresh residual
readback. This single update is not a trajectory or convergence certificate.
[Native cut-filter control](cls-epawf-segregated-picard-r3-evidence/native/ball-filter-native-review.json).

The first complete 504-response model at original residual 18.6652852921 has
rank 504. Its original-merit gradient norm is 1908.20. Near-active cone defects
at normalized windows 1e-10,1e-8,1e-6 are 1908.20,1804.91,1759.12, versus
empirical two-size gradient uncertainty 9.15. Thus this base is not shown
stationary even when every free k coefficient is represented. This is a
first-order audit of the original reduced operator at the model's base, not of
the later checkpoint and not a global root-existence test.
[Full-space first-order audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-full504-audit.json).

That bounded504-response run publishes one accepted **numerical** checkpoint
at original residual 16.4201949282. A separate process rebuilds its original
norm vector and complete qualification: continuous ground 2.63392640e-6,
gradient departure 0.0009999370393, minimum outer-epsilon ratio 0.9800205038.
The process then times out during the next response construction without a
final native stage report. Therefore its final physical-history receipt is
unavailable; the independently audited checkpoint is not physical acceptance.
The ordinary constrained follow-up keeps the original h0/history and stage
target, starting only its nonlinear guess from this qualified nonroot.
[Fresh checkpoint audit](cls-epawf-segregated-picard-r3-evidence/native/audit-wide504-checkpoint-d710404.json).

The completed matched 192-direction native control reaches 18.1656906795,
with original qualification passing and zero physical steps. Verification
On/Off reproduces every logged initial/accepted numerical state hash and norm
vector bitwise, together with accepted state, rejection and physical history.
The constant-row-norm cache produces the same complete native trajectory
bitwise on the same input. These controls do not establish that192 directions
outperform160; the160 comparison above follows two updates.
[Verification receipt](cls-epawf-segregated-picard-r3-evidence/native/fast-face192-verification-comparison.json),
[native cache identity](cls-epawf-segregated-picard-r3-evidence/native/cached-native-identity-review.json).

`c02f965` adds opt-in `RANSWallReducedResponseCache=On` (driver
`--reduced-response-cache`). Small response probes can reuse the separate base
velocity factor, pressure/gauge Schur factor and epsilon/trace factor. At every
iteration the original nonlinear H is rebuilt, including epsilon's actual
matching-trace residual after the flow correction. The original fixed-k,
accuracy, admissibility and rollback tests remain mandatory. A failed cached
subsolve retries from its incoming guess with current-operator field solves.
The canonical half-response branch check, combined-direction checks, actual
trials and final acceptance use the current operators. Every model rebuilds
these preconditioners at its equilibrated base; no all-field matrix is stored.
The default is Off, preserving the retained controls. Six focused Linux
regressions pass (38.48s) and the reduced-model ASan/UBSan check passes (7.78s).
The matched 32-direction native controls complete at 18.6379733298 with reuse
and 18.6380670320 without reuse, both qualified nonroots with bitwise fresh
residual readback and preserved accepted h0/history. All 32 responses pass
two-size checks; the cached path's maximum residual/observation relative errors
are 6.05e-5/0.001215. The independent current-operator canonical half response
agrees to 8.02e-6 in relative action and 2.32e-5 in normalized observation action.
All 74 logged H evaluations satisfy their requested 2e-12,386 cached flow/epsilon
applications are recorded, and no fresh retry is needed. Verification On/Off
reproduces every numerical update and all accepted physical history bitwise.
Reuse versus current-operator iterations are tolerance-consistent, not bitwise
identical algorithms. Native target invocations take677.1s versus755.3s
(reduced-solver log times656.1s versus734.1s), under concurrent experiments;
this single comparison is not a general performance claim. A new unchanged-h0
run continues. Build/operator progress must not be read as M8 convergence.
[Producing build and library identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-c02f965.json).
[Controlled native reuse comparison and independent audits](cls-epawf-segregated-picard-r3-evidence/native/response-cache32-review.json).

`ce7e53a` additionally reuses a feasible point and checked active face between
radius queries of the same large quadratic model. Obsolete constraints leave
through their multiplier signs; numerical failure retries the original feasible
point with an empty face. This is only an optimization-coordinate warm start.
The 192-coordinate model replay differs in coefficients by 4.30e-16,
retains the same radius multiplier and passes every original cut. It exercises
14 warm solves and one successful cold retry; the final primal defect is
2.18e-17 and affine stationarity2.16e-12. Replay time is 20.87s. Six Linux
regressions pass (37.77s), including a changed-gradient problem that removes an
obsolete active constraint and a deliberately inconsistent warm face;
ASan/UBSan passes. The ordinary full-response follow-up uses this implementation
with the original target and an explicitly declared inexact H tolerance that
tightens toward the unchanged final 2e-12 gate.
[Build and tests](cls-epawf-segregated-picard-r3-evidence/native/build-validation-ce7e53a.json),
[real-model warm-face replay](cls-epawf-segregated-picard-r3-evidence/native/warm-face-radius-replay-review.json).

The eta-zero R2-endpoint diagnostic also constructs all 504 orthogonal field
responses (rank 504, maximum column correlation3.11e-16). At its initial
original residual 2.9088342324 the reduced gradient norm is102.388, its
near-ground cone defect 3.397 and empirical two-size gradient uncertainty 0.716.
This is an original-alpha400 model, although its wrapper forbids physical
advance and observes qualification only. It is not shown stationary. The
direction predicts residual 2.90880695 and merit slope -1.48724843e-4; an
independent feasible half probe measures -1.48720552e-4 with relative vector
error 4.81e-5. The regularized model has eigenvalue range [0.001353,5.660e6],
condition 4.18e9. The completed two-update diagnostic reaches original residual
2.908809697352376, continuous ground minimum8.98767433e-15, gradient departure
.003821016357 and minimum epsilon ratio.9745118805. The original qualification
still fails. Fresh-process readback reproduces the original norm vector exactly;
the native receipt verifies zero accepted steps and unchanged accepted state
and every physical-history slot/time. The full run takes18502.93s and records
2029 controlled H evaluations and11524 velocity/pressure and epsilon solves.

The second full-rank model, based after the first update, has original gradient
norm102.604, near-active cone defect5.83457 and empirical two-size gradient
uncertainty.689805. Neither tested base is shown stationary; these models also
do not certify stationarity of the final numerical endpoint. A small affine KKT
defect is not substituted for original reduced stationarity. The useful actual
descent is small despite all504 field directions and the original ground cuts;
this diagnoses the local difficulty without proving global nonexistence.
[Operator-labelled full-space diagnostic](cls-epawf-segregated-picard-r3-evidence/native/first-order-r2-full504-audit.json).
[Completed two-update diagnostic and independent final audit](cls-epawf-segregated-picard-r3-evidence/native/r2-full504-completed-review.json),
[both full-space first-order models](cls-epawf-segregated-picard-r3-evidence/native/first-order-r2-full504-completed.json).

`b89f2be` caches immutable cut-row norms within each quadratic solve and computes
the direction norm once per QP iteration. It changes no matrix entry, cut or
acceptance bound. Six Linux regressions pass (37.80s), together with ASan/UBSan.
Replaying the actual192-coordinate native model before and after caching gives
identical output bytes, coefficients and affine KKT values (69 steps,68 augmented
face solves,0 fallbacks). The standalone replay takes13.15s versus6.83s; this is
a helper replay timing, not a native convergence or general performance claim.
[Build and checks](cls-epawf-segregated-picard-r3-evidence/native/build-validation-b89f2be.json),
[native-model identity replay](cls-epawf-segregated-picard-r3-evidence/native/cut-cache-native-replay-review.json).
The original-h0 four-update prefix on217e47d also matches0633d23 bitwise in every
state hash and original norm vector, with only run budget/ceiling changed.
[Exact prefix receipt](cls-epawf-segregated-picard-r3-evidence/native/original-start-prefix-217-review.json).
The matched cached-kernel native control has completed as recorded above;
full-response follow-ups continue.

`217e47d` adds a checked augmented solve for active faces of models above160
variables. It keeps the orthonormal constraint QR, solves a scaled saddle system
in response coordinates, performs two residual refinements and reprojects the
correction onto the active face. A projected linear-defect check gates this
path; numerical failure returns to the existing Euclidean null-space solve.
This addresses the measured roughly17-minute-per-cut cost in the504-direction
control; the speed and resulting native steps are still being measured. The
saddle matrix contains optimization coordinates and multipliers, never the
independent governing fields. Six Linux regressions pass (36.94s), as does
ASan/UBSan including a192-variable closed-form case and deliberate saddle-factor
failure/fallback. [Producing build and checks](cls-epawf-segregated-picard-r3-evidence/native/build-validation-217e47d.json)
record commands and executable/library hashes.

`cfec988` extends only the optional response ceiling to512, permitting a
benchmark-sized complete free-k response diagnostic (504 free coefficients).
The default remains32; no all-field matrix or governing-equation solve is
introduced. It passes all six focused Linux regressions (37.40s). Matched160-
and504-direction constrained native searches from the same qualified diagnostic
state are in progress; requesting504 columns is not yet a completed rank or
convergence result. [Build, identities and regressions](cls-epawf-segregated-picard-r3-evidence/native/build-validation-cfec988.json)
record the exact producing revision and unchanged installed libraries.

The current numerical direction implementation is `0633d23` (constrained nonlinear-source
direction refinements, opt-in and verified in matched native controls).
`04f275f` supplies the completed source-aware step-length comparison described
below. Both pass the six focused Linux regressions and standalone ASan/UBSan
checks. The earlier producing solver is `7715312` (optional larger response spaces,
retaining the point/trace seed prefix beyond64). `bdcf657` permits an explicit
ceiling up to160; the default remains32. `97ec1d7` supplies the expanded
free-mode seed pool.
`0926974` repairs the response-probe construction at continuous ground contact;
`8f420d1` supplies adaptive response-space growth and seed-history refresh;
`849c70f` supplies measured constraint-response secants and a configurable
64-column ceiling. Verification comparisons also cover diagnostic continuation
inputs (`0a9b095`, `e651396`). The earlier producing solver is
`cbbd486` (explicit reduced continuation diagnostic).
The ordinary constrained curvature-correction runs use `ba89d11`; orthogonal
field-response probes and combined-direction checks use `62201de`.
The radius-constrained model control is `9a5b58d`.
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

## Completed larger-space and nonlinear-source controls

The completed `7715312` enlargement used96,120,120 independent responses in
three numerical updates and reached original residual18.8837693663 with a
positive continuous-ground minimum3.32426e-5 and original qualification passing.
It remains a nonroot. The matched free-mode64 control reached18.9235364098
versus18.9453816015 for the point/trace-only59-dimensional control after two
updates. Both preserved accepted h0/history and reproduced their original
residuals independently. The first64 response records and original gradient
are byte-identical in the64/96 comparison; their original action-Gram prefix
differs relatively by1.81e-20 after removing each radius multiplier. These are
actual added directions, not a changed physical problem.

The120-dimensional original-merit gradient norm is2611.72. Near-active cut
cone defects for normalized windows1e-10,1e-8,1e-6 are2611.72,2514.18,2321.34,
against empirical two-size gradient uncertainty13.01. This is evidence of
remaining descent in the tested response space, not a constrained stationary
point or a global infeasibility result. The old affine-projection KKT defect
does not replace this original-gradient calculation.

At the last rejected full enlarged-space trial, vector prediction error is
9.07004; inner destruction accounts for9.06188, leaving.010613. The small
combined-direction check is accurate (half-step relative residual error.00215),
but the accepted nonlinear fraction is only1/128. Increasing the response space
alone therefore did not resolve finite-step curvature.

`04f275f` adds optional `RANSWallReducedSourceRay=On`. Its internal surrogate is

\[
M(z)=r_0+Bz+D(s_0+S'z)-D(s_0)-D'_0[S'z].
\]

Here B contains measured converged non-k responses, S' predicts their independent
field coefficients, and the independent k displacement is inserted exactly.
D is the complete native nonlinear inner destruction component, with profiles,
moving quadrature and normalization rebuilt. The final subtraction is calibrated
with a small directional component difference. This affine-field state is an
internal model, never labelled Y(k). Every actual trial still converges H,
rebuilds the original residual and scans the full qualification and ground sets.

The selector samples11 dyadic model fractions and at most one parabolic
interpolation, promotes only a model decrease, and retains the original dyadic
fallback. The declared policy also permits nonlinear qualification corrections
at smaller fractions and uses source-model gain for radius adjustment. Thus the
following comparison is of these complete numerical policies, not an isolated
one-line interpolation change. The target and acceptance conditions are identical.

| Completed matched warm-start run | Source-aware selection | Linear-ray control |
| --- | ---: | ---: |
| Initial original residual |18.8837693663|18.8837693663|
| After two numerical updates |18.8135344559|18.8761902618|
| Maximum outer-gradient departure |.000999999468|.000999999954|
| Minimum outer-epsilon ratio |.9841946074|.9841502701|
| Continuous-ground minimum |2.27570e-5|3.27262e-5|
| Controlled H evaluations |301|313|
| Velocity/pressure pairs; epsilon solves |1887;1887|2132;2132|
| Solver elapsed seconds |2898.00|3207.38|
| Physical steps/time |0/0|0/0|

Both use64 then80 columns with matching numerical rank, initial radius.4 and
H tolerance2e-12. The source-aware final actual displacement is.11826566 within
.4. Its final original block norms are
[1.0210e-12,1.5109e-12,18.8135344559,4.2990e-14,8.0064e-16,7.3629e-16,
6.2511e-17,0]. Fresh native complete audits report r2_qualified=1, positive
covariance and lambda1; the complete original residual still fails overwhelmingly.
Both failed native stages preserve accepted h0/history bitwise, with zero fresh
readback norm discrepancy. The negative-ground pre-evaluation guard also passes.

The first actual source-aware trial uses fraction.09664672: predicted norm
18.85479443 versus actual18.85497599 after qualification correction. Vector
error drops from the original linear model's.0933476 to.00380373. The second
uses fraction.29583513: predicted18.81466100 versus actual18.81353446, with
vector error.0223100 versus.224006 linearly. These are useful finite-step
predictions and actual reductions, not convergence certificates. The final
gradient witness moves to x=-1,z=4; the control's witness is x=-1,z=2. Complete
signed inventories continue to catch observations away from the initial active set.
The verified/unverified source-aware runs have bitwise identical terminal
coefficients, rejected state, accepted physical history and time. Verification
only adds observations and work.

`0633d23` adds an optional budget of up to three constrained nonlinear-source
direction refinements (`--reduced-source-refinements`; default0). Starting from
the source-ray candidate, each refinement builds a two-size-checked derivative
of this surrogate in the measured response coordinates and solves a small
Gauss-Newton model with the same signed cuts and actual displacement radius.
It requires surrogate improvement and descent against the original measured
reduced gradient. Undefined or inaccurate model derivatives, failed model solves
and unhelpful candidates retain the original ray. Actual trials and nonlinear
qualification corrections remain controlled H operations. This is a model
direction refinement in the existing free-k response space, not an all-field
Jacobian or an unequilibrated trial acceptance path.

The matched0633d23 native controls now complete two updates. Refinement reaches
18.7645225166 versus18.8136254786 for ray-only selection, from the same
18.8837693663 guess and physical history. The first update alone uses the same
64 response columns and reaches18.8264587558 versus18.8549759911. On update2
the refined method retains64 columns/rank64, while the adaptive ray-only control
uses80/rank80. This comparison includes the declared adaptive policy; it is
not a claim that adding source curvature changes only one coefficient.

Both refined full proposals are accepted after controlled H and nonlinear
qualification correction. Final actual displacement.25833959 is inside radius.4.
The second refined direction's original predicted merit slope is-2.5010258;
actual small-step slopes are-2.4988820 and-2.4989143. Half-step residual-action
error is.00128031. All64 surrogate derivatives pass their two-size checks
(maximum.00483935 on this refinement). At the final actual trial, the nonlinear
source prediction is18.76287582 versus18.76452252 observed; vector error.0187001
compares with.534987 for the original linear prediction.

The refined endpoint's full original block norms are
[9.9746e-13,1.1673e-12,18.7645225166,4.0851e-14,5.3893e-16,7.3850e-16,
2.0668e-17,0]. Independent complete qualification passes: gradient
.000999999106, epsilon lower ratio.9839749693, ground2.38331e-5, positive
covariance and lambda1. Verification On/Off matches every initial/accepted
numerical state hash and original norm vector, as well as terminal coefficients,
accepted h0/history and time. Neither run accepts a physical step.

The refined verified run uses273 converged H evaluations and1773 separate
flow/epsilon pairs, versus301 and1890 for its ray-only control. It additionally
spends549 internal source-model evaluations constructing its two directions.
Elapsed time is3525.21s versus2824.72s. The larger actual residual reduction
therefore comes with model-construction cost; this bounded comparison supplies
useful progress, not a convergence or general performance guarantee. The fresh
original-h0 run and six-update extension now complete as recorded above; both
fail the original residual tolerance.
At their common initial k, a disturbed/tighter non-k warm start changes the
reduced residual by3.37e-11, while H itself meets2e-12. H tolerance is therefore
not claimed to bound the propagated k-residual error. The signed response,
two-size and actual combined-direction checks supply the empirical action
accuracy evidence; the original complete residual remains the final authority.

The dated03:25UTC original-start snapshot records14 numerical updates and
residual20.32393218 for04f275f, and9 updates/residual40.95530307 for0633d23.
Both have entered original qualification, but neither has solved or accepted
its first physical stage. These are different progress points, not a matched
performance comparison or completed trajectory validation. The
[history-bound snapshot](cls-epawf-segregated-picard-r3-evidence/native/source-model-progress-20261106-0330.json)
records the original initial-guess selection, alpha400, t.0025, 2,165 unknowns,
physical-history hashes, complete block histories, ground minima and XML/log
identities. Running cases are not counted as completed failed-stage history tests.

The six focused CTests pass on04f275f (37.41s) and0633d23 (38.50s), and both
standalone model binaries pass ASan/UBSan. New known-answer checks exercise
nonlinear-ray curvature, an invalid full surrogate sample, no-descent fallback,
and two-scale one-sided differentiation where the opposite branch is undefined.
The120-response native affine packet also passes ASan/UBSan replay in31.87s:
41 active-set iterations, primal defect2.70e-18, affine stationarity3.00e-12,
complementarity1.12e-15. The initial replay was rejected by a legacy65-variable
input guard, not a sanitizer memory finding. Test-only commit`e65b0be` extends
that guard to161 (160 responses plus restoration slack); production solver code
is unchanged. These affine defects do not supersede the nonzero original
gradient audit above.

See [completed source controls and independent audits](cls-epawf-segregated-picard-r3-evidence/native/source-ray-completed-review.json),
[completed direction-refinement controls](cls-epawf-segregated-picard-r3-evidence/native/source-refine-completed-review.json),
[completed enlarged-space controls](cls-epawf-segregated-picard-r3-evidence/native/enlarged-completed-review.json),
[120-dimensional original-gradient audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-large-audit.json),
[nested response-prefix check](cls-epawf-segregated-picard-r3-evidence/native/nested-prefix-check.json),
[source-ray build identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-04f275f.json),
and [refinement build identities](cls-epawf-segregated-picard-r3-evidence/native/build-validation-0633d23.json).
The [large sanitized replay receipt](cls-epawf-segregated-picard-r3-evidence/native/replay-160-validation.json)
preserves the initial rejection and successful corrected replay.

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
start. From `0926974`, the first used response of every model is also solved
from the canonical warm start regardless of the optional verification setting:
inconsistent local responses reject a column, so this branch-selection check
must be algorithmic. The extra verification observations do not select columns.
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

## Orthogonal responses and nonlinear qualification curvature

The raw 18-column model could pass individual two-scale checks while amplifying
column errors: its saved field Gram matrix had condition number 18,746, maximum
column correlation .97494, and coefficient-l1/actual-displacement amplification
116.62. `62201de` probes the normalized two-pass QR directions themselves when
the actual ground and scalar inventory permit them; otherwise it retains the
feasible inward seed. Every retained column still receives its own two-size H,
original-residual and signed-observation checks. It also checks the **combined**
proposed direction at two small admissible fractions before larger trials.

At the qualified saved checkpoint, all 18 or 32 probes are orthogonal (maximum
correlations below 6.1e-17), with response rank 18 or 32. The combined 18-column
merit slope is predicted -2.93994496 versus measured -2.93983456/-2.93984242.
Its residual-action errors are 4.38e-5/5.14e-5 and normalized observation-action
errors 3.89e-6/8.08e-6; the full inventory remains qualified at both probe sizes.
The 32-column errors are 5.61e-5/4.37e-5 and 1.51e-5/1.41e-5, with no recorded
profile branch change. These are local empirical checks with attained H<=2e-12.

The 32-column full trial predicts qualification excess about 1.1e-14 but measures
.00190403 at x=-3.75,z2. Halving exposes another signed gradient witness at
x=-1,z4, excess .000118772; the quarter trial still exceeds by 1.116e-6. Only
one eighth is qualified, reaching residual 19.83066648. Thus accurate local
responses do not remove nonlinear qualification curvature or witness migration.

`ba89d11` adds optional `RANSWallReducedSecondOrder=On` / `--reduced-second-order`.
For a qualified base and residual-improving trial outside qualification, it
minimizes a correction in the checked response space using the **actual trial's**
complete signed inventory. A penalized elastic slack provides a feasible model
start. Domain and ground cuts remain hard; each correction is bounded by half
the initial trial length, and the **total actual displacement** must stay inside
the outer radius. Each assessed correction re-solves H to the same declared
accuracy, rebuilds dependent data and scans the actual full inventory/extrema.
There are at most three such H evaluations per model. Failures retain the caller
state and return to backtracking. Ordinary final acceptance is unchanged.
Residual actual/predicted gain controls the radius inside qualification; a
useful half step no longer automatically causes successive .6 radius shrinkage.
The previous control remains selectable with the option Off.

The first 18-column native correction restores the rejected full trial from
qualification excess .000107983 to 0 using correction norm 3.85259e-5. Its actual
residual is 19.76496439 and total displacement .25905542<=.2592. The control had
to halve the step and reached 19.81849948. This is useful nonlinear constrained
progress, **not an algebraic root or an accepted physical step**. A fresh process reproduces the residual to 3.6e-15, continuous ground minimum
.00131160, gradient .000999999522 and sampled-K maximum 1.0099999900; the
original complete qualification audit passes. The second update reaches
19.63663504 while remaining qualified.

The 32-column full trial is restored in three bounded corrections: normalized
qualification excess .00190403 -> 6.27958e-5 -> 6.59128e-7 -> 0. The final
residual is 19.59825833 and displacement .25915022<=.2592, versus the control's
19.83066648 at one eighth of the step. Corrections are 2.23961e-4, 2.51467e-6,
and 6.53987e-8. These are candidate-based correction predictions; the original
base's flat constraint model need not predict qualification after a nonlinear
correction. Actual rebuilt observations remain the acceptance authority.
Extended original-start and larger-space tests continue; completed results are
recorded separately as they finish.

The one-update original-h0 verification control passes the native operator
checker. Verification On/Off gives bitwise-identical rejected coefficients,
accepted h0 and native history; their XML differs only in that observation flag.
The rejected-state SHA256 is
`e3d5aee52a87fdd2fbbc88ba82bffca8d891d257af157b086b4e96e411bca1fb`.
The standalone curved-qualification known-answer model passes Mac and Linux
ASan/UBSan. The six focused native regressions pass for the correction source.

## Expanded responses, constraint secants, and contact repair

The following completed comparisons retain the original physical alpha400,
t=.0025, h0 and complete qualification inventory. They use saved guesses only.
All endpoints in this table are qualified **nonroots**, with zero accepted
physical steps and unchanged accepted h0/history. Independent final residual
readbacks agree with the native block norms exactly in these controls.

| Matched numerical search | Updates | Final original full residual | Ground minimum |
| --- | ---: | ---: | ---: |
| 32 responses, existing nonlinear qualification correction | 4 | 19.4017629319 | 4.52781210e-5 |
| 32 responses, measured constraint secants | 4 | 19.1478499169 | 6.57257635e-5 |
| Expanded budget64, measured constraint secants | 3 | 18.9970573458 | 6.13252399e-5 |

The first two start from exactly the same state, radius .2592, response step
5e-6, H tolerance2e-12, basis32 and Krylov budget8, using the same849c70f image.
The secant policy updates the constraint-response model with measured changes
in converged y and projected trace, restoring the exact independent k response
after each update. It leaves the original verified residual-action matrix B
unchanged. The finite-step secant is an approximate constraint model, not a new
analytic Jacobian. Its minimum physical field signal is1e-6. It permits six
qualification-correction evaluations rather than the control's three; the
comparison measures this complete policy, not an isolated secant coefficient.
At the third update the control fails to restore a full step in three
corrections and accepts1/32. The secant policy restores and accepts the full
step in four corrections. Its attained H errors remain below2e-12.

A separate one-update restart from the same later state gives19.2463018061
with secants and19.2463013182 without them: no material secant benefit there.
Restarting the witness/previous-direction seed history itself changes the
response space and can recover a much better direction. These comparisons
must not be conflated. The adaptive option grows a slow qualified search's
budget, and refreshes only numerical seed history after severe shortening at
the ceiling; original equations, observations and physical history remain.

The expanded run actually used59, then60 independent responses, rather than
64 in every model. Its final gradient departure is.0009999994783995714,
outer-epsilon minimum ratio.9832639559914378, local-k minimum ratio
.9905905192146565 and sampled-K range[.9990614789558008,1.0058610123511045].
All complete qualification checks pass. Its finite full-step model can still
be inaccurate: one prediction18.76032761 gives actual19.30544167 and a vector
error3.80071212; the qualified half step gives19.20896975. A later full trial
has vector error12.00285512 and is rejected. Local action accuracy and an
accurate affine QP do not guarantee a useful finite trial.

`97ec1d7` appends free-k coefficient modes ranked by the mass-inverse original
residual when an adaptive budget reaches40. They pass the same rank,
continuous-ground, controlled-H and two-scale action checks as every other
seed. Only the declared number of usable responses is constructed. This fills
directions missing from the point/trace seed span without assembling an
all-field matrix. Coefficient signs choose probe orientations, not positivity
certificates. The first matched native update gives18.9466404975 with64
directions versus18.9548244686 with the59-direction point space, both qualified.
The two-update comparisons continue. A96-direction comparison and a fresh
original-h0 run with adaptive ceiling160 are also in progress; this is a
numerical search budget, not an extension of serial/flat support.

The original reduced gradients were also audited separately from affine QP
KKT checks. In the59-dimensional model, the gradient norm is2337.58. Even
treating all normalized cuts with slack<=1e-6 as active leaves an NNLS
normal-cone defect1777.79, versus empirical two-size gradient uncertainty10.72.
The32-dimensional live and restarted models likewise leave defects607.71 and
832.67, with uncertainties4.77. These are coordinate-window diagnostics on the
measured finite spaces, not nonlinear stationarity or global feasibility
certificates. They do not support a stationary-obstruction claim.

An unchanged-cut radius replay of the59-column model predicts residual
18.73071675 at radius1, versus18.59878782 at radii4 and8. The latter minimizer
has displacement3.39334472 and inactive radius. Without cuts, its affine
prediction is.01604 at displacement5.99592. Constraint effects dominate this
local affine model; increasing the radius alone gives little benefit. These
are model-only results and do not establish that no nonlinear feasible root
exists.

At a separate eta-zero algebraic diagnostic's ground contact, the previous
one-orientation probing policy retained only12 orthogonal probes among32
columns; maximum correlation reached.9999999984191911, and the QP cycled.
`0926974` tries both signs and smaller admissible orthogonal probes before
using an inward fallback. The matched contact model retains22 orthogonal
probes, maximum correlation.9864776663, rank32, and solves its QP. Its full
step decreases the original residual18.4547287117 ->8.9157262358, with actual
displacement.9993335851 <= radius1 and continuous ground minimum6.38992e-15.
The remaining nonlinear error is3.83436665, of which3.82783217 is inner
destruction and.03444716 remains after that component. This decomposition
does not separate intrinsic destruction curvature from nonlinear y response.
The small-step original-merit slope is predicted-269.14804 and measured
-269.11826/-269.11878 at two sizes; corresponding residual-action errors are
2.82e-4, and signed-observation errors about1.08e-2.

That diagnostic observes qualification without enforcing it and is permanently
guarded against physical advancement. Its gradient departure is about.00467448,
so it is not a successful original stage. A fresh process rebuilt its residual,
continuous ground and original qualification; a deliberately negative ground
state was rejected before scalar constitutive evaluation. Verification On/Off
has bitwise identical rejected coefficients, accepted h0, history slots, time
and bookkeeping, with guess files compared by exact content. The extended
diagnostic continues; a decreasing constrained residual is not a root certificate.

An independent component study at this new contact state keeps its continuous
ground polynomial unchanged. At h=1e-5/H tolerance2e-11 and h=5e-6/H
tolerance2e-12, the conditional original-merit slopes are+261.80639/+261.77797,
while the reduced slopes are-889.82835/-889.84907. The feedback is essential
to that direction's descent. The result is additional derivative evidence,
not convergence or a statement about other branches.

The independent readback observer reconstructs actual Gauss-panel endpoints,
checks all ordered Gauss nodes, and records limiter selections and evaluation
elements. Both small reduced probes retain2481 panels on112 rays, with no
observed limiter, cap, profile-category, panel-count or evaluation-element
switch. Panel endpoints nevertheless move by6.08222e-6 and3.04112e-6, so a
fixed-quadrature interpretation would be wrong. The finite full contact step
changes2487 ->2481 panels,73 ray inventories, six profile categories and at
least220 comparable limiter selections, with no cap switches. It observes
counts, geometry and ordered element membership; it does not infer origin
labels of coincident internal breakpoints or prove differentiability there.

See [completed space controls](cls-epawf-segregated-picard-r3-evidence/native/completed-space-controls-review.json),
[contact response and independent audit](cls-epawf-segregated-picard-r3-evidence/native/ground-contact-response-review.json),
[original-gradient audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-completed-audit.json),
and [unchanged-cut radius replay](cls-epawf-segregated-picard-r3-evidence/native/radius-model-audit.json).
Large per-response receipts, saved independent coefficients and matrix packets
remain in the listed Linux artifact directories.

[Contact component responses](cls-epawf-segregated-picard-r3-evidence/native/contact-feedback-components.json),
[small-step panel/limiter observations](cls-epawf-segregated-picard-r3-evidence/native/profile-branch-contact-response-audit.json),
[finite-step branch changes](cls-epawf-segregated-picard-r3-evidence/native/profile-branch-contact-full-audit.json),
and [dated original-start/large-space progress](cls-epawf-segregated-picard-r3-evidence/native/large-space-progress.json)
retain the scope and producing runtime identities. The earlier response's
[panel observer](cls-epawf-segregated-picard-r3-evidence/native/profile-branch-feedback-audit.json)
also finds unchanged observed selections at both probe sizes, with moving
endpoints5.96200e-6/2.98101e-6.

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
| Completed original-h0 run (`0633d23`,18 numerical updates) |19.1747990683|.000999999004| qualification passes; original residual fails;0 accepted physical steps |

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

The original-start qualification failure occurs before a k update. A fresh
probe at the actual target alpha400/time.0025 confirms that original h0 is
qualified (gradient.000500194924, ground.5105080131). Keeping every k
coefficient at h0 and solving H raises the gradient to.002857296287, with the
negative witness at x=1.625,z=2. The k residual changes62.60959507 ->63.00211881
while the non-k blocks meet their declared2e-9 inner tolerance. Thus a k update
must account for flow/epsilon feedback and qualification even from the original
initial state.

The full inventory tracks witness migration during restoration: x=.59170018,
-.75,-3.66666667,.33333333, then-3.66666667 at z=2. Some intermediate epsilon
and local/sample-k bounds are also violated; evaluation-domain and continuous
ground admissibility remain enforced. The seventh numerical update first enters
full original qualification, with gradient.000999999001 and residual63.34723722.
The next two source-refined updates retain qualification and reduce residual to
49.11413627 and40.95530307. This is restoration/search progress, not acceptance
of a physical stage. The [fresh qualification trace](cls-epawf-segregated-picard-r3-evidence/native/qualification-first-violation.json)
binds the baseline, unchanged k during H, h0 history and target parameters.
Its first strict raw-bit h0 guard rejected serialized zero signs; the repeat
asserts exact equality of every coefficient value and the original serialized
history hash. No incoming state or history coefficient was replaced.

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

The reduced-response continuation option (`cbbd486`, driver `--method reduced`)
observes original qualification while solving the labelled auxiliary algebraic
problem. It retains actual domain and continuous-ground cuts, separate field
solves, full original residual readback and the wrapper's unconditional rejection
of physical advance. Ordinary reduced stages still enforce qualification during
direction construction and final acceptance. There is no unqualified ordinary
stage path. The eta-zero guard test deliberately requests auxiliary tolerance
1e-8; the native target retains 2e-12 and preserves accepted state/history.
A reproduced eta500 root also returns through explicit native no-advance
rejection, with original-target residual 25.54458179, gradient .00160673758 and
ground .07065555. A stricter reduced repeat reaches auxiliary residual
1.92635e-11 after one k correction. These are positive-eta diagnostic roots.
A separate demand for H<=1e-13 fails explicitly at 1.62051e-12 with
`non_k_line_search`; it is not silently used as Y(k) and does not establish a
fundamental precision floor. A second branch experiment declares positive-eta
auxiliary tolerance 2e-10; eta zero remains 2e-12. The followed lower-eta reduced
branches are still in progress and do not yet certify ground contact or an
eta-zero root.

The controlled reduced branch has since certified eta450,400 and300. At
eta300 the continued norm is6.86352505669e-12 (declared positive-eta
tolerance2e-11), original-target norm18.9335635180, gradient departure
.00195261177912, epsilon lower ratio.982196382536, and continuous ground
minimum.0281343569353. Both continued and original residuals were rebuilt
independently. Original qualification fails the gradient bound, and accepted
physical state/history remain unchanged. Eta200 subsequently certified continued
norm7.54681124186e-12, ground.0110053572728 and gradient.00218914012184; it
also fails original qualification and advances no physical history.

A secant guess extrapolated from the eta400/450 roots to eta300 had ground
minimum-.235017 at x=-2.748466; the native guard rejected it before scalar
evaluation. Yet the safeguarded eta300 solve above obtained an admissible
auxiliary root. This is a concrete example of why a failed initial predictor
must not be interpreted as branch nonexistence. See the
[dated reduced-continuation receipt](cls-epawf-segregated-picard-r3-evidence/native/reduced-continuation-branch-review.json).

## Continued roots near ground contact

An explicitly separate positive-eta sequence declares auxiliary tolerance1e-10;
the original eta-zero gate remains2e-12. It certifies roots at150,140,135,132.5,
131.25,131,130.98,130.96,130.955 and130.952, with both original and continued
residuals independently rebuilt and physical history unchanged. These are not
solutions of the original alpha400 equations.

| Eta | Continued norm | Continuous-ground minimum | Gradient departure |
| ---: | ---: | ---: | ---: |
|131.25|6.83683e-11|4.63368342e-5|.00238854189|
|131|1.16219e-11|7.56300042e-6|.00238933420|
|130.98|2.27125e-12|4.46192363e-6|.00238939761|
|130.96|3.11374e-11|1.36082225e-6|.00238946102|
|130.955|8.18585e-11|5.85543060e-7|.00238947687|
|130.952|9.38468e-11|1.20374808e-7|.00238948638|

All violate the original.001 gradient bound. The three closest certified-root
secants predict ground contact at eta130.95122362..130.95122367, with ground
slope about1.55056e-4 per eta. This is a local branch estimate, not a certified
root/nonexistence bracket. At130.951 the tested search reaches an arithmetic
ground margin and residual about4.93e-6. That bounded process has no final native
history receipt; it must not be counted as a verified rejected stage. A separate
64-response restart has completed with a native no-advance receipt, a fresh
continued residual4.92248169e-6 and original-target residual9.55312366.
Its actual continuous ground minimum is5.20104e-17 and gradient
departure.00238948989. The original and continued residuals are distinct:
this state equilibrates non-k equations at operator alpha530.951, not400.
Accepted state and every physical-history slot/time are preserved. At the
restart's base the tested64-dimensional continued-merit cone defect9.46e-6
still exceeds empirical action uncertainty1.50e-6, so this experiment does not
establish a local stationary obstruction even in that space.
[Completed contact receipt and both readbacks](cls-epawf-segregated-picard-r3-evidence/native/contact64-completed-review.json),
[continued first-order audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-contact64-audit.json).
The branch qualification exit was already observed much earlier, at
eta1168.75..1209.375 with positive ground. No statement excludes other roots or
other H branches.

[Compact root/readback receipts and secants](cls-epawf-segregated-picard-r3-evidence/native/continued-ground-collar-review.json)
retain the declared auxiliary versus original tolerances and incomplete-case
statuses. The first-order auditor now labels the actual operator explicitly;
[its positive-eta check](cls-epawf-segregated-picard-r3-evidence/native/first-order-continued-root-audit.json)
says `continued_operator`, with operator alpha530.952. It is not original-target
stationarity information. The
[original160/96-space audit](cls-epawf-segregated-picard-r3-evidence/native/first-order-160-audit.json)
retains large original-merit cone defects (2345.02 and1794.46 at window1e-6)
versus empirical action uncertainties12.01 and11.32. These states are not shown
stationary in their tested spaces.

The ordinary504-response search enforces qualification. An additional explicitly
labelled eta-zero diagnostic starts from the preserved R2 endpoint with the same
physical alpha400/time.0025/history and2e-12 original tolerance. It observes
qualification while enforcing ground/evaluation validity and unconditionally
rejects physical advance. Its purpose is to distinguish a ground obstruction
from qualification, and it cannot replace the original-h0 native test.

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
All six pass on0633d23 (38.50s),04f275f (37.41s),7715312 (41.99s),97ec1d7 (38.24s),0926974 (41.65s) and849c70f
(39.80s). The larger120-response native packet subsequently passes ASan/UBSan
with the e65b0be test-only input-guard repair. Pressure-verification source is unchanged from
the0926974 bounds-checked two-sweep comparison.

Heavy artifacts are under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard-r3`.
The shared Linux checkout and installation were preserved; the shared revision
is `f4b70baec10016dd098005bc5a5507b87e4a57bf` and its pre-existing symlink content
status is recorded in the build manifest. Source/build trees and immutable copied
executables are isolated. Eight explicit compiler workers, one test at a time,
one rank and one thread per BLAS/OpenMP runtime were used, within the72-worker
limit including the simultaneous bounded experiments. No support claim extends
beyond the existing serial/flat tensor-cubic ground implementation.

The completed18-update original-start command can be reproduced in the isolated Linux tree
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
  --solver "$R/binaries-0633d23/IncNavierStokesSolver" \
  --resolution-probe "$R/binaries-0633d23/RANSWallResolutionProbe" \
  --reduced --reduced-verify --reduced-basis 24 --reduced-krylov 8 --reduced-radius .4 \
  --reduced-restoration-weight 0 --reduced-second-order --reduced-constraint-secant \
  --reduced-adaptive-space --reduced-maximum-basis 160 --steps 3 \
  --reduced-source-ray --reduced-source-refinements 2 \
  --sweeps 18 --timeout 21600 --output "$R/source-refine-original-reproduction"
```

The harness deliberately returns nonzero unless the original stage completes;
run the operator checker separately on bounded failed-stage controls. Each
native numerical command is one rank. Logs and manifests distinguish solver,
readback, expected rejected-stage and startup return codes.

[Orthogonal verification controls](cls-epawf-segregated-picard-r3-evidence/native/orthogonal-controls.json),
[larger-space source build and sanitizers](cls-epawf-segregated-picard-r3-evidence/native/build-validation-7715312.json),
[free-mode source regressions](cls-epawf-segregated-picard-r3-evidence/native/build-validation-97ec1d7.json),
[contact repair and pressure verification](cls-epawf-segregated-picard-r3-evidence/native/build-validation-0926974.json),
[constraint-secant source validation](cls-epawf-segregated-picard-r3-evidence/native/build-validation-849c70f.json),
[verification observability](cls-epawf-segregated-picard-r3-evidence/native/orthogonal-verification.json),
[raw response-basis conditioning](cls-epawf-segregated-picard-r3-evidence/native/basis-conditioning.json),
[curvature-source regressions and sanitizers](cls-epawf-segregated-picard-r3-evidence/native/build-validation-ba89d11.json),
[Current-code regressions and continuation guards](cls-epawf-segregated-picard-r3-evidence/native/build-validation-cbbd486.json),
[fresh corrected-state audit](cls-epawf-segregated-picard-r3-evidence/native/audit-curvature-checkpoint18-1.json),
[Completed orthogonal and radius controls](cls-epawf-segregated-picard-r3-evidence/native/orthogonal-checkpoint-completed.json),
[fresh 32-direction corrected-state audit](cls-epawf-segregated-picard-r3-evidence/native/audit-curvature-expanded32-1.json),
[explicit interim curvature histories](cls-epawf-segregated-picard-r3-evidence/native/curvature-progress.json),
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
Rejection after an already accepted physical step also remains outstanding;
the completed rejection checks exercise startup history.
