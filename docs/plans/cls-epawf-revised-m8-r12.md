# M8-R12: verified frozen Lf direction and resumed comparison

R12 completes its primary goal: native Lf and an independent reference pass
the full original separation oracle and the unchanged numerical certificate
checks. The corrected native solve needs 123 admitted cuts and 124 scans.
Optimizing at the trial scale improves actual residual progress. Three Nf
tangents produce a distinct direction with a further, small residual gain
and substantially better trace behavior than Lf. Neither policy passes every
practical promotion gate. No continuation or physical timestep is run. M8
remains open.

The reviewed head was 6ab80ffe1895811e958e703c2a4a51e3e0496999; no newer
remote commits existed before editing. The numerical producing revision is
cf6a4bddcb8a733abc8f6be0374548060d96d438. Its source, binary, runtime, inputs
and large payloads remain at
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r12/epoch-004`.
The [compact ledger](cls-epawf-revised-evidence/m8-r12/decision-ledger.json)
pins exact commands and payloads. Historical producing builds, all M0–R11
evidence, shared Linux source/installations and the separate Picard worktree
and jobs were preserved.

## Frozen contract and scope

| Quantity | SHA256 |
|---|---|
| R9 initial / physical R8 iterate-10 state | eb063e15f16b2610f1c5a832280d87894f567e8dab321560202e6dd9af95194c |
| Physical J at that state | 2194965dee9c5960ba4e1c1604f76e8e85372cabef2d5ced7c8a394f593c8dc6 |
| Original h0 | b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71 |

E(q)=2.911641510023811, continuous ground minimum is
6.733016805109465e-6 and scaled actual matching discrepancy is
4.2429689461870394e-4. The bound residual, constraints, mass metric, source
value/action and XML remain R11's. The experiment retains 2,165 unknowns,
outer6 p3/p1 spaces and quadrature, alpha400, time .0025, eta0, zero forcing,
original h0, SaturatedR2 rho .8, beta .4, GammaScale1/GammaUniform0, original
outer C(q), viscosity/diffusion conventions, samplers, sensors, profiles and
projected trace. Destruction remains Capped with Linear smoothing, theta .2
and delta .1. Radius .125, 5% reserves, true-contact handling, the disabled
positive-value inward rule, original box and 3,000-unknown guard remain.
Root, action and original-row gates are 2e-12, 1e-7 and 5e-12. Existing KKT
gates remain 1e-9.

`RANSWallR12Experiment=Replay|Complete|Verify|Compare` is opt-in.
`RANSWallR12CutLimit` is a numerical SOLVERINFO setting. Prior modes retain
their original exchange schedule, including historical failed limits.

For s=.125 the frozen affine solve uses r_D=R_D/s on free differential rows,
r_C=R_C on all overwritten/algebraic rows, and A=J. Its origin is always
-J^-1 r_f. L1 is supplied separately as a witness; it never replaces that
origin. Full total displacements satisfy C+J_C d=0. Fractional native
predictions still contain (1-s)C. No pressure-continuity condition is added.

## Exact limit-stage diagnosis and repair

R11's archived `exchange_limit_candidate` is recorded before admitting the
last incoming cut. The old loop processes that cut, then exits without
reconstructing or checking the updated working set.

R12 adds one final reconstruction and complete oracle/KKT check after the
last allowed admission. That check cannot admit another cut. It retains
the pre-update snapshot separately. All original rows are scanned in native
and separately accumulated extended arithmetic; wall-reserve extrema are
recomputed for every candidate. RefinedR6 forward/transpose refinement and
direct reconstruction remain active.

| Frozen replay point | Admissions | Scans | Worst original-row violation | Dual lower bound in Psi units |
|---|---:|---:|---:|---:|
| Historical/pre-update limit snapshot | 95 processed, incoming 96th | 96 | .127584602606 | 271.748641736 |
| Corrected post-update limit check | 96 | 97 | .119452598684 | 271.769863064 |
| Completed native Lf | 123 | 124 | 4.03849e-12 | 271.957191723 |

Thus the missing completion check was real, but its correction alone does
not qualify the 96-admission result. The original failure was an unfinished
numerical construction, not nonlinear stagnation or physical infeasibility.

The instrumented replay discovers 73 distinct original affine rows, with
23 exact alias readmissions and no repeated physical working sets. Its
nonnegative dual bound increases throughout: .0302643 at scan24,
176.6477 at scan48, 181.0272 at scan72 and 271.7699 at scan97.
The worst violation need not decrease during dual cut discovery.
These observations support continued discovery rather than a cycling claim.

Before the completion attempt, a single extension to 256 admissions was
declared. This uses the same dual algorithm and original constraints.
It finishes at 123, preserving distinct nearby wall rows. No constraint
coalescing, relaxed accuracy or second production strategy is used.
The final active system has 31 rows; its minimum QR diagonal is
1.39325e-5. Ninety-two releases occur. Full original/row/transform payloads,
all discovered and released row identities, per-pass candidates, bounds and
exchange histories are retained outside Git. Selected active coefficient
and transformed gaps agree within 1.66e-14.

## Feasible witness and independent certificate

L1 passes all Lf scalar/coefficient/outer/continuous-wall restrictions, the
native and extended scans, full affine algebraic restoration, native slope
and full source-domain checks. Its worst original violation is 4.63410e-12;
its scaled algebraic defect is 5.83e-20 and slope is -.0149347535214,
below the -.000582328302005 bound.

The native witness gives Phi_Lf=4.264985002760582 and
Psi_Lf=272.95904017667726. The independent metric reconstruction gives
4.264985002760593 and 272.95904017667795. This supplies a feasible upper
bound for the frozen QP, not a physical-root existence result.

The independent reference uses an isolated Python environment with SciPy
SLSQP in a pivoted-QR subspace, followed by 80-digit active-system KKT
polishing. Its feasible starting coordinate is
B_D^-1(r_f+J d_L1)_D. Equivalently, recentering around L1 requires the
objective shift .5||y+c||² as well as shifted inequalities. The coordinate
regression checks both, and checks unchanged algebraic rows.

The 31-row finite pool is a relaxation. Both its independently reconstructed
candidate and the native candidate are subsequently checked with the
COMPLETE original native/extended oracle, including moving continuous-wall
reserve extrema. Both pass, with full anchor source domains valid.

| Check | Native | Independent reference |
|---|---:|---:|
| Psi_Lf | 271.9571917228855 | 271.95719172288557 |
| Phi_Lf | 4.249331120670079 | 4.249331120670087 |
| Stationarity norm | 5.43e-15 | 1.48e-15 |
| Complementarity | 8.07e-16 | 1.33e-13 |
| Worst full original-row violation | 4.03849e-12 | 4.03870e-12 |
| Scaled algebraic defect | 1.87e-19 | 3.03e-19 |

SLSQP succeeds in 17 iterations. Its objective differs from the polished
reference by 4.72e-11. The polished nonnegative dual lower bound is
271.95719172288557; the displayed double-precision primal/dual gap is zero.
The minimum normalized multiplier is .000563708. Original/transformed
reference gaps differ by at most 9.99e-14. The predeclared objective-agreement
and primal-dual gates are 2e-8, with a 2e-9 feasible-upper-bound allowance.
These are numerical certificates, not interval proofs of exact arithmetic
optimality.

The independent wall-reserve polynomial has a sign-resolved minimum
-2.06177e-12 at x=1.56383864737, between two distinct active wall locations.
Normalized by K0 it is -4.03867e-12, within the unchanged 5e-12 cut gate.
This small direction-model defect is reported; it is not a strictly
nonnegative exact polynomial. Full-anchor and selected actual-state k
remain positive. Nothing is floored, clipped or reprojected.

The final-cut schedule regression accepts a solution completed by its last
allowed cut, rejects a still-infeasible post-limit candidate, and verifies
that the final scan admits no extra cut. Read-only replay of the actual
completed trace also gives success for limit123/scan124 and rejection for
limit122/scan123 without another native solve.

## Matched native comparison

All three rays use radius .125, the unchanged R8 fitting/held-out/proposal
and fallback selector, and native decrease 1e-4*f. The selector does not test
full f=1. Each policy evaluates two native states; all chosen fractions are
.125. Selection is bound by the explicit chosen name, fraction, state and
residual hashes, including multiple accepted trials.

| Policy | Actual E | E reduction | Ground minimum | Actual trace block | Charged policy cost |
|---|---:|---:|---:|---:|---:|
| L1 | 2.910253517943 | .001387992081 | 5.9334752793e-6 | .000377244015026 | 27.72 s |
| Lf | 2.905663388314 | .005978121710 | 5.9334717079e-6 | .000564403040731 | 146.88 s |
| Nf | 2.905421256019 | .006220254005 | 5.9334724882e-6 | .000372403166812 | 531.10 s |

Lf's reduction is 4.307 times L1's. Nf adds 4.05% to Lf's reduction and
reduces its trace block by about 34.0%. All selected rays deplete the base
ground minimum by about 11.875%. Nf is a distinct constructed
direction, not an unchanged fallback seed.

Nf completes three safeguarded tangents, using 146, 119 and 135 scans.
Accepted optimizer fractions are .25, .25 and .0625. Phi_.125 decreases
monotonically through 4.2527163169, 4.2525597506, 4.2522018447 and
4.2520162205. Moved source-action checks pass at 5.92e-9 and 2.56e-9;
the second check follows ordered profile-branch changes. No redundant
third check is required after unchanged branches. There are 20 new source
values, four cache hits and three source-action assemblies. All tangent
original-row and KKT gates pass; all accepted anchors retain full affine
algebraic restoration, the original reserves and native-slope screen.
The model stops at its declared three-tangent budget with continuing
improvement, not a stationarity certificate.

At the selected Lf state, affine model E is 2.904497607194, whereas native E
is 2.905663388314. At the selected Nf state, destruction-model E is
2.905420778717, within 4.77e-7 of native E. Its blockwise vector errors are
6.41e-5, 5.40e-7, 1.98e-5, 2.01e-7 and 1.14e-6 for the two velocities,
k, epsilon and trace respectively. The governing norm error is distinct
from these vector errors. Model trace remains .000371259782791, retaining
the fractional algebraic residual. Actual trace is always reevaluated.

Lf fails the ground, trace and efficiency promotion gates against L1.
Nf passes the gain/trace gates against L1 but fails ground and efficiency;
against Lf it fails the required 25% extra gain and efficiency tests.
The ground losses relative to L1 are only 3.57e-12 and 2.79e-12, but the
predeclared no-lower-ground criterion is retained exactly. No threshold is
moved. No policy is promoted and no conditional continuation is run.
The lowest-E admissible research trial, Nf, is saved separately from
production-policy selection and physical history.

## Costs, verification and provenance

The frozen-QP execution budget is 900 s. All earlier failed entry-command
costs, 64.69 s, are charged. Successful verification finishes with 127.32 s
remaining. Builds, source integration and desktop analysis are recorded
separately. The earlier setup/development clock is also preserved.
The matched comparison takes 402.23 s; the native process, including
initialization/final reporting, exits normally after 438.64 s. Its
cooperative/native-supervisor limits are 600/650 s. Independent audits are
separately charged.

The 96-admission replay takes 102.11 s inside the frozen entry, and native
Lf completion takes 119.09 s. The independent QP reference takes 4.61 s.
Nf's tangent construction/factorization times are 125.27, 108.78 and
132.44 s; source values take 8.95 s and source actions 2.06 s. Its own
model/trial work totals 389.78 s. Direction construction dominates cost.

Policy charges include required seed construction. L1 uses its bound,
unchanged producing construction charge of 22.392846542 s plus current
validation/trials. Lf adds current completion work, and Nf adds Lf seed and
its own model/trial work. Shared work is reused in the study; these
counterfactual per-policy costs do not sum to comparison wall time.
Independent reference, failed development attempts and audits are separate.

The isolated build uses 16 compiler workers. Numerical runs use one rank,
allocated CPUs and explicitly one OpenMP/OpenBLAS/MKL/BLIS/vecLib/NumExpr
thread. Concurrent Picard solvers and MPI/runtime threads are accounted for
within the total 72-worker ceiling. Shared jobs and installations are not
modified.

Three native CTests, nine predictor-binding tests, five selection tests,
two coordinate/signature tests and the compiled cut-budget regression pass.
The mathematical producing revision remains cf6a4bd. A final metadata-only
correction names the row family `family`, avoiding collision with the
common event `kind`. The archived schema is read by its stable identity and
row pins. The corrected C++ file passes the Linux warning-enabled compiler
syntax check; the native event-writer fixture passes strict duplicate-key
parsing. No numerical equation or returned direction changes in that fix.

Initial integration failures are retained: a failed compile, the missing
R12 predictor-role registration and the use of a numerical PARAMETER that
changed the physical signature. Both entry omissions were fixed before
any QP iteration. An initial PMIx warning and a rejected runtime-fingerprint
attempt are retained; the private MPI environment is explicitly bound.
An audit helper's initial sibling-module import error is corrected to use
the tested source selector, with failed read-only logs preserved. These
failures are distinct from QP or physical failures. Exact wall time for the
three failed read-only audit import attempts was not recorded; successful
audit timings are reported separately.

The new Lf and Nf states are independently audited in 249.36 and 245.09 s
on two separate single-threaded CPUs. Independent E values are
2.9056633883142866 and 2.9054212560186072. Scalar and fixed-box checks pass;
80/110-digit wall reconstruction confirms the positive minima in the table.
Source identity errors are 1.93e-15 and -1.65e-15; same-state energy identity
errors are 9.68e-14 and -7.11e-15. L1 reuses the exact previously audited
state/residual/history tuple. Exact runtime/command bindings are in the ledger. At the retained Nf state, weak inner/outer continuity block norms are
1.58e-15 and 2.10e-15, and the gauge norm is 1.91e-17. Native strong
divergence L2 norms are .0339426 and 9.30810e-5; they are reported separately.
Independent normal/tangential traction RMS defects are .000383262 and
.000206310, and k-flux RMS defect is .000653314. These nonzero interface
observations are not relabelled as balance at a root. Mechanical-pressure
continuity is not imposed. Native history checks verify original h0/startup
slots and zero accepted physical steps. Accounting identities at a non-root do
not establish balance at a converged RANS stage.

## Reproduction and decision

Use the committed `m8-r12/experiment.py` to prepare fresh external inputs
for Replay (96), Complete (256), Verify and Compare. The committed
`reference.py` independently solves/polishes the frozen finite-pool QP;
`collect.py` resolves exact chosen candidates and history; `audit-selected.py`
audits or reuses identical tuples; `review.py` builds the compact ledger.
Recorded phase commands have this form (archived paths must be replaced
by fresh external outputs for a rerun):

```bash
r12_epoch='/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r12/epoch-004'
r12_ev="$r12_epoch/source/docs/plans/cls-epawf-revised-evidence"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python3 -B "$r12_ev/m8-r12/experiment.py" --epoch "$r12_epoch" \
  --out "$r12_epoch/replay-bound" --phase Replay --cuts 96
python3 -B "$r12_ev/m8-r3/run-stage.py" --epoch "$r12_epoch" \
  --session "$r12_epoch/replay-bound/case.xml" \
  --manifest "$r12_epoch/replay-bound/predictor.json" \
  --output "$r12_epoch/replay-run" --phase r5-local --seconds 240 \
  --expected-source cf6a4bddcb8a733abc8f6be0374548060d96d438 \
  --cpu 0 --metadata-cpu 50
taskset -c 48 "$r12_epoch/reference-env-clean/bin/python" -B \
  "$r12_ev/m8-r12/reference.py" --epoch "$r12_epoch" \
  --case "$r12_epoch/complete-bound" --output "$r12_epoch/reference-result.json"
```

The private MPI environment is the pinned `private-mpi-launch.json`; the
runtime fingerprint must be prepared with that same environment. Complete,
Verify and Compare commands and exact direction arguments are retained in
their manifests/receipts and comparison declaration. The ledger pins each
R3 `run-stage.py` invocation, source/runtime/XML,
reference command, budget declaration, regression commands and result hashes.
A rerun requires a fresh output directory and the same historical payloads;
never overwrite the archived epochs.

| Question | Result |
|---|---|
| A. Why original construction stopped | Missing post-update qualification plus insufficient 96-admission work budget; productive discovery, not demonstrated cycling |
| B. Independent frozen Lf | Solved and qualified against the complete original oracle and numerical certificate gates |
| C. Native Lf construction | Qualified at 123 cuts / 124 scans |
| D. Actual direction progress | Lf improves E reduction 4.31x; Nf adds 4.05% and improves trace, at much higher cost |
| E. Promotion / root / physical advance | None |

Next decision: retain the existing production policy and benchmark a
proper feasible-start construction for successive Nf tangents, transporting
physical incumbents/row identities and rebuilding transformations for each
current A. Construction consumes almost all additional work. Keep the
same objective, original-row, native ground and trace gates; any remaining
no-lower-ground failure still prevents promotion. Do not increase model
iterations or launch a new physical campaign on these results.

Equilibrium preservation, numerical construction, constitutive coverage,
nonlinear convergence and physical-time advancement remain separate claims.
All historical broad-box constitutive counterexamples remain unchanged.
Curved-wall MMS, active-cap/variable-height coupling, full-element and fully
discrete positivity, refinement, long-time robustness, scalability and
physical floor/hill validation remain outstanding. M8 remains open.
