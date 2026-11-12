# WJ-H3: corrected pressure seeds and bounded principal diagnosis

Baseline: reviewed/fetched `b83ec8979cb029ab39136ef89ff14f88a5031a75`,
tree `83f9bfd6bd88cad496f62a45a9b38256f596b8eb`. Existing Mac WJ worktree,
normal commits/pushes; no reset, rebase or unrelated source change.
H2's original results, ledger and T=8 prelaunch denial remain immutable.
The refused continuation driver/flag will not be invoked or modified.

This NEW allocation permits only fresh corrected controls through T=.08.
Original WJ/Budget equations, coefficients, raw-work convention, full3D tensor,
positive wall nu, current endpoint-K cutoff/trace, uncapped source, Q3/Q1 spaces,
one-sweep accounting and hard predicates remain unchanged. No cap, clipping,
epsilon alteration, pressure stabilization, state reset or hill run is included.

## Storage and runtime

The server rebooted and auto-mounted the same NTFS volume
UUID `5AB8A669B8A642FD` at `/media/jay/新加卷1`. The requested old path is an
empty root-owned directory. No remount/unmount/installation change is made.
Use explicit read-only path translation for original H2 identities; do not
edit their report, XML, checkpoint or signature bytes. Actual new root:
`/media/jay/新加卷1/nektar-v2026-runs/wj-h3-20261112`.

Reuse governing ELF `9ae37a1aba62755b10ad99c9484c498cef9d927a75fb950f00344f4ffc5d06fe`
from scientific producer `c344dc2209f030e426e2ee5501d5e5b72a6d2c46`, and observer
`6f35b837f0dd13eb155e0fe34f0dc670c609ee6cbe739609a5de288dd49ed06d`
from `e62b5f4d265f068dd5af3f090bfb0e385c3ae53e`, after native-input/library
equivalence checks. RANS_BUILD_SOLVER=ON, Release C++17,
-ffp-contract=off -fno-fast-math retained. Python generator/observer changes
have their own source identities and do not imply a new governing build.
An isolated source export is used; no native rebuild is planned unless a
concrete changed native behavior requires one. Separate build ceiling1200s.

One rank; OMP, OpenBLAS, MKL, BLIS, NumExpr and VECLIB threads1. Inspect current
scientific workers before each launch and keep the combined total<=72.
Post-reboot initial load~.3, available memory~264GB; inspect actual mounted
storage before work. Hard limits:44 scientific attempts, plus4 reserved ONLY
for a genuine changed-code regression;3600s numerical/observer/closure charge;
8GiB process-tree RSS and invocation output,16GiB total root. Native timeout
300s initially (bounded at900s if measured cost warrants), passive300s, closure
reserve300s. Persist reservations before launch; unknown attempts/time remain
charged. No reset/increase. Passive setup/analysis/plots are charged too.

## Pressure and completion correction

Before edits verify the continuous old expression and retained base/refined
CP0. For parallel initialization, T_w,nn=-pi_w and
T_o,nn=2K/3-pi_o-R_o,nn, hence pi_w-pi_o=R_o,nn-2K/3.
The H2 physical expression has the opposite sign. The analytic neutral seed
already has the correct sign and must remain unchanged.

Version the minimal physical seed `SignCorrectedV1`: pi_o=0,
pi_w=Rzz(H)-2K/3. Independently form both complete physical stress tensors;
the old formula must fail, corrected continuous normal balance must pass near
roundoff. Represented Q3/Q1 CP0 has separate projection/lift errors and no
extra strong interface constraint.

Version a distinct `VerticalBalanceV1` seed: C=Rzz(H)-2K/3,
pi_w=C, pi_o(z)=2k(z)/3-Rzz(z)+C. Then continuous Tzz=-C everywhere.
R uses the SAME declared U/k/epsilon and actual physical gradients. This is
initial data only, not an equilibrium source or an ongoing pressure condition.
Record expression approximation, Q1 projection, normal traction and vertical
spatial residual separately. Keep initial shear/scalar/streamwise defects.

Correct labels only: input L=f*delta differs from actual asymptote
L+kappa*H. Record both; do not change the historical completion formula/bytes.

## Passive retained-state and integration design

Read/hash original base/refined CP0, refined accepted .04/rejected .06, and
refined half-dt accepted .04/.06/.08. Use genuine report/audit callback and
checkpoint validation; rejected state is a candidate, never an accepted CP.
Reconstruct element-owned gradients on both normal-face owners, x=8 and an
inward band (x=7.999999,7.99,7.9,7.5,6,4), including H and the failed face.
Local maps retain R spectrum, root, principal direction/margin, Gxx/Gxz/Gzx/Gzz,
k/epsilon, divergence, actual/applied outlet traction and H traction defects.
Compare weak pressure moments with strong divergence; use explicit offline
pressure-space projection/complement, not a governing correction.

Independently replay failed, prescribed-profile and swapped-gradient/ratio
states. Check isotropic-gradient removal leaves WJ stress/principal unchanged;
tracefree/pure-shear cases are diagnostics only. Do not infer cause from bulk
L2, positive covariance or a universal q threshold.

Outer load integration: native8/12/16 versus independent24/32, at most48
if reference convergence is not resolved. Hold coefficients, physical measures,
element/physical breakpoints and three test modes fixed. Report element/mode/
flux-component and regional sums. Preserve original component limits
2e-9+2e-8*abs(value), regional4e-9+2e-8*abs(value). The native outer quadrature
uses TangentialOrder; NormalOrder stays16. Check retained fields AND corrected
CP0 before selecting smallest adequate q_star. Use an explicit passive
integration view, never changed-signature checkpoint readback.
No qualified q_star => all B trajectories NOT_RUN. Check decisive B terminal
accepted/rejected fields against the converged higher rule within the observer
budget. A later integration failure stays UNRESOLVED; no adaptive rerun.

## Fixed new case matrix

A: original refined43-cell normal mesh, Nx1 (1966 states), zero complete
traction, original quadrature8/16, original U/k/epsilon/h_s/H. Only pressure
offset changes. Fresh CP0; dt=.02,<=4 attempts to .08, first hard/correction
exit stops. One fresh dt=.01,<=8-attempt companion is allowed only to resolve
the retained timestep sensitivity. No further halving. A max12 attempts.

B: vertically balanced seed, original base23-cell normal mesh, dt=.01,
fixed qualified q_star/16,<=8 attempts per case to .08:

| Order | Outlet | x knots | Expected state size |
|---|---|---|---:|
| B1 | zero full traction |0,8|1062|
| B2 | frozen continuous seed full traction |0,8|1062|
| B3 | zero full traction |0,8/3,16/3,8|2630|
| B4 | frozen continuous seed full traction |0,8/3,16/3,8|2630|

All four are needed to separate the paired outlet/x effects if resources and
q_star permit. Confirm Nx2 base1846 with zero steps; refined Nx2 theoretical
3418 is outside the3000 evolution guard and will not be evolved.
Frozen outlet T_seed*n comes from the same continuous completed fields and
balanced pressure, separately in each region. No pressure Dirichlet/gauge
addition. Fits retain errors; this is a startup/boundary diagnostic, not a
validated production outlet. A/B differ in several inputs and are not a
single-factor comparison.

All cases stop immediately on their first hard/nonfinite/correction exit;
continue only a DIFFERENT predeclared control. Each complete attempt retains
one U, pressure, k, epsilon solve and one bound call; preserve partial counts.
Do not restart a failed state, retry unchanged data, tune coefficients, change
q_star adaptively, or advance any case past .08. Compare available common
times/physical points and local gradients/principal margins, not only bulk
field percentages. A short positive sampled result is not developed flow or
global stability. The later precursor-development design is a proposal only.

## Publication

Compact H3 code/tests/contract/evidence, scientific tables/figures and receipts
index immutable heavy data. Distinguish initializer sign, expression/projection,
integration, outlet/x sensitivity and unresolved constitutive effects.
Keep H2 historical discrepancies/approval records visible. Commit and push
normally after focused checks; no broad historical replay or new SVD campaign.
