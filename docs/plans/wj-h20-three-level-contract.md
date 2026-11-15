# H20 frozen three-level target-domain study

Starting/fetched HEAD: fd009626de46e164c0f645b9879a68d2102ae11f. The WJ worktree
is clean, no H20 root/job exists, and healthy ext4 has about84.7 GB free.
Source/Git is Mac; authoritative execution is Linux. Prior ledgers stay closed.

## Fixed physics and identities

Retain the full-height neutral6 mm flat control, physical interface9.5 mm,
d0=3.5 mm, h_s=4.341801417014184 mm, z0,e=.227050812848672 mm and
K=.008809015145515192. All700 Q3/Q1 elements/24688 coefficients, IDs/signs/maps,
q12/normal16, Functions, original full3D WJ/root/raw work/Cd, finite-nu Budget AWF,
current-K cutoff/epsilon, nonconstant inner k, pressure convention and complete
boundaries remain fixed. IMEX1, eta0 and one flow-k-bound-epsilon sweep stay.
Use CheckpointAuditV1/BoundedV1/ExactStepV1, ImmutableRowReuse Off, StopOnBound On.
No new native feature/build is planned. Reuse H19's b2365e52 capsule/ELFs/private
closure only after actual input/dependency/loaded-library verification.

Original XML SHA256137566f8704fed63d881dfbb5b804a02776a09c2081b9f29c42cadd91e794511;
original coupled qh0 SHA256631caa32a59ff43253f1fcac4effd955eb68b8529d9c28f29e615cf14ca0f077.
The literature hashes1959e96e…dbd9124 and0ab52e0d…fbc47c and their accepted
quad/\quad and `(d+a*h)/a^2` errata remain unchanged. No height/mesh/model/BC/
regularization/forcing/refit/pressure-regauge/EP or physical-experiment work.

## Immutable matrix and order

| Order | Level | dt | Restored step/time | Additional attempts | End | Scientific stride |
|---:|---|---:|---|---:|---|---:|
|1|M|.0025|H18 step4/.01|4|8/.02|2|
|2|C|.005|H19 step2/.01|2|4/.02|1|
|3|F|.00125|one compatible CP0|16|16/.02|4|

Exactly3 maximum governing blocks/22 new attempts,1 new target CP0 and at most3
final strict reads, serially after each process exits. Recovery stride1; audit
restored startup and final states. XML Time remains0, NumSteps means additional
work; TimeStep remains in the signature. Explicit ceilings are M8/C4/F16.
No rerun of old C/M evolution and no other timestep or horizon beyond.02.

Common times and absolute steps: .005=(C1,M2,F4), .01=(C2,M4,F8),
.015=(C3,M6,F12), .02=(C4,M8,F16). Each retains its own producer/history/time
bits/signature/callback, alpha200/400/800; initial F is alpha0. Fine-old M2 keeps
its H18 A producer and separate H18 B restored audit. CP suffixes are not steps.
Match recorded times within128*eps*max(1,abs(t)); no interpolation/regauging.

F's existing native step-zero helper strictly imports original fine CP0 under
its own signature/dt, preserves coupled fields/pressures/reference/defect/eta
bitwise and calls native IMEX1 Initialise once at.00125. Zero spatial Project,
Lift, initial velocity transform, reference reconstruction or governing solve.
No history copy/rescale. Run-local anchors are each level's own dt-compatible
CP0; common spatial lineage is separately pinned. Reuse H19 refusal fixtures.

## Phase roles and support

M receives physical qualification; comparisons are NOT_YET_AVAILABLE as planned.
C receives physical qualification and C-M absolute screens against existing M
partners. F receives complete triplet/all-pair screens and contraction diagnostics.
A partner missing when required is an integrity failure, not planned absence.
Restored-startup screen stops are labelled separately from physical failure and
preserve the old qualified anchor; completed REVIEW_STOP commits current physical
verification before stopping. No prepared acknowledgement promotes a cursor.

Use fixed H16/H18 support11456 points/17521 owners plus44800 physical Gauss8 rows.
Add all current participating native witnesses with established +/-1e-6 physical
neighborhoods, each state's current cutoffs +/-1e-6 at fixed matching-outer x,
and all legitimate owners. C-M uses its two-state union; F replays all three on
the final current triplet union. Reassessed old-state failure retains history and
halts new work. Immutable cache bindings prove dt independence, not changed
signature validity. Recompute every state/profile/root/history-derived decision.
Retain fixed probes (15,0,.25607517888672315), geometry251 upper/upper, and
(-7.5,0,.1), geometry57 lower-x/upper-normal, with actual adjacent owners.
No averaging or lowest-ID selection. Required current rows cannot be evicted
from coverage. Variable-row catalog512; current view/state/action arrays bounded
to this finite support, with estimated extra comparison working memory300 MiB.

## Frozen screens and contraction definitions

Unchanged native hard/scalar/root/covariance/transverse/viscosity/cutoff/outlet
checks; quality>=.002068763923454 and all255+525 Q1 moments<=2e-15. Strong
divergence is separate; evolved states do not inherit CP0 pointwise-zero demand.
Every successful sweep has one U/P/k/epsilon solve, one bound call, zero activation/
material correction and equal raw/corrected k hashes.

All C-M, M-F AND direct C-F pairs use H19 limits per regional physical RMS/max:
U1e-4(scale1), k1e-3(scaleK), R/pi1e-3(scaleu_tau^2), G1e-2(scale
max(||G_eq(r)||_F,.051/19.9125)), outer epsilon1e-2(scalepositive epsilon_eq(r)).
Physical region measures497.8125/3.75/494.0625; vector/tensor Euclidean/Frobenius
norms. Floor4e-11 scaled, with sqrt(measure) for unnormalized L2; action replay
accuracy2e-11. These are provisional numerical-review limits, not application
uncertainties/production budgets. Derived inner destruction is separate.

For pair coarser-finer and SAME outer owner, B=.1*(a_finer-a_gate). Compare
minimized actions and both states at both member minimizing directions against
that same B. C-M uses M; M-F/C-F use F. Retain dimensional/normalized actions,
nu_eff, directions and signed differences. Unresolved allowance/action sign or
resolved absolute exceedance halts target work without constitutive relabelling.

On one triplet support, E1=qC-qM, E2=qM-qF, D1=norm(E1), D2=norm(E2), rho=D2/D1.
Report regional RMS/max and signed component/fixed-owner cancellation. D1<=floor
means UNRESOLVED; D2<=floor alone gives an upper bound, not zero/exact order.
Resolved p_obs=log2(D1/D2) is descriptive only. rho>=1 or unresolved ratio ALONE
does not halt if physical/absolute gates pass. Do not impose .5 or any ratio gate.
At each thetaC/thetaM/thetaF evaluate ALL THREE actual states at that same owner
and direction before forming principal differences; retain minimized differences
separately. No manufactured constitutive hybrids, Richardson bound/state or
global stability claim. H11's different-domain findings remain historical.

## New finite allocation and exit

New root /home/jay/build/wj-h20-three-level-20261115 and7200s ledger: setup60,
closure300 fully charged; protect3*120s final reads. Combined component/bootstrap/
target/strict<=3600; optional affected build<=900 (none planned). Reserve every
operation/failure/unknown, require25% prospective measured remaining-cost headroom
before each branch/scientific barrier, no native retry/reset/allowance transfer.
One native plus bounded observer, all6 numerical threads1, host<=72/build<=8.
Science RSS<=8 GiB/build<=16; source/runtime<=8 GiB/output<=2 GiB. Estimated
base native+observer6.874 GB plus300 MiB triplet work fits the cap, measured guard
remains decisive. Keep every short-study recovery and first-failure evidence.

Frozen prior: ordinary compute about40s/attempt, audit44s, pair observer21s,
triplet observer conservatively45s, startup110s, helper150s, final read120s.
Update from M's retained cost; protect all remaining strict/closure reserves.
Separate inclusive waits from nested audits/observer work. No future fine oracle
is assumed for production forecasts. Hypothetical recovery/science/terminal
physical intervals.025/.25/.125 imply C5/50/25, M10/100/50, F20/200/100, not this
qualification schedule. Preserve unexecuted horizon505.02466521278404/window
[447.8483832973273,505.02466521278404] and ceiling-step distinction.

Publish integrity, physical gates, absolute screens, resolved/unresolved
contraction, reference/drift, terminal strict residuals, actual costs/storage and
closed accounting separately. Missing required work is NOT_RUN/unavailable.
Recommend one finite next decision, execute none. No production dt, long flat
preservation, curved terrain/hill/experiment or physical-height admission follows;
7 mm observations remain inside the9.5 mm physical interface.
