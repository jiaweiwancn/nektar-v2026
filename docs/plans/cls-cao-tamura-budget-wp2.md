# WP2: neutral open approach-flow qualification

2026-11-16–17, Asia/Shanghai. Checkpoints 1 and 2 pass for both deployed regularized
CLS and calibrated standard k–epsilon. Checkpoint 3 qualifies the short timestep
pair and measured cost. Both models have completed the unperturbed physical-transit,
half-step and finer-mesh runs. WP2 remains **OPEN** pending complete perturbation
recovery, including the final-window and trend criteria.

The physical calibration, full-height analytical inlet, local transported-K
cutoff budget and one-sweep equations are unchanged. The new
`NeutralEquilibriumOpenV1` policy supplies full roof shear and epsilon conormal
flux, stress-consistent regional outlet traction and absolute pressure, with
immutable analytical backflow data. Top U/k/epsilon remain free; W is essential.
The old boundary policies retain their meanings.

The stored pressure is pi=p/rho+2k/3. At Pout=.003 the independent selected-kernel
check gives CLS pi_inner=-.004039556572351707 and
pi_outer=-.0014583858291560801; standard gives both .003. CLS's retained
mechanical pressure jump is -.0025811707431956273. No separate roof/outlet pressure
Dirichlet datum or closed-box mean gauge is added.

All native face modes, including nonconstant modes and both outflow/backflow
signs, agree with independent GL32 integrations within 1.2e-15 relative error.
With configured GL24 face/source quadrature, analytical free-block weak loads are
at most 1.58e-13 in native normalization. Wrong-top and omitted-outlet-shear
controls respectively give U=.04529 and W=.04033 defects. Both models pass new
production-policy direct/split restart, post-assignment rejection/rollback and
changed-boundary restart rejection. These tiny fixtures qualify implementation,
not a transit interval. Full receipts, actual source hashes and immutable solver,
probe and library identities are in [checkpoint1.json](cls-cao-tamura-budget-wp2-evidence/checkpoint1.json).

Radial energy projection is an optional representation operation, including
outer epsilon and model-correct regional pressure. It supplies no balancing
volume load and is not a stationary PDE solve. The common full six-column
progress is:

| Cubic slots | CLS U / W / k / epsilon residual | Standard U / W / k / epsilon residual |
|---:|---|---|
| 2804, L2 | .13704 / .004816 / .04859 / .001826 | .13686 / .0001048 / .04651 / .002950 |
| 2804, energy | .02261 / .002042 / .005188 / .0006376 | .02199 / .0001243 / .004970 / .0007684 |
| 5339, energy | .003470 / .0004199 / .001045 / .0001060 | .003343 / .00001984 / .0009998 / .0001234 |
| 8002, energy | .001272 / .0004486 / .0001386 / .00006463 | .001056 / .00001863 / .0001302 / .0001034 |
| 9882, upper refinement | .001236 / .0004280 / .0001343 / .00006288 | .001026 / .000003697 / .0001265 / .0001025 |

The upper-only refinement was ineffective. Signed panel contractions against
the full inverse free-row mass located large cancellation around the matching
plane, with shared CG rows and cross terms retained. Splitting the last inner
and first two outer panels gives the selected **8650-slot** common case:

| Case | CLS U / W / k / epsilon | Standard U / W / k / epsilon |
|---|---|---|
| Selected | .000469076 / .000166622 / .000126918 / .000036735 | .000412545 / .000018553 / .000120549 / .000029991 |
| 11014-slot confirmation | .000329744 / .000137559 / .000072366 / .000033939 | .000277433 / .000018548 / .000067656 / .000026954 |

Both cases use the complete six-column topology, cubic fields and regional P1
pressures. Increasing source/face quadrature from GL16 to GL24 changes the
selected block norms by at most 2.2e-11. Independent GL28 plus the frozen
physical stations and exact native side limits pass both models' full physical
budgets. Independent covariance minima (inner, outer) are
(.00631381050,.00609524100) for CLS and (.00683838645,.00671295638) for standard,
far above each model's unchanged 1e-4 K0 target. Independent NumPy eigenvalue
replays agree within 1e-12. No stationary PDE correction or artificial supply
is used.

Selected versus confirmation physical differences, integrated on the union of
their interfaces, use at most 8.41% of the allowed quarter-budget limits.
Absolute regional pressure offsets are retained. One-ulp ReferenceK integration
roundoff between native panel partitions is recorded explicitly, without
changing normalization or a physical tolerance.

The historical default evolution limit remains 3000. The explicit opt-in
`WP2MeasuredDenseV1` permits at most 12000 only for this serial affine one-sweep
problem, following measured 5339–9882-slot audit peaks below 3.73 GiB. Physical
consistent masses, pressure inverse and all acceptance inventories are
unchanged. Actual larger first/continued-step RSS and cost must still pass
checkpoint 3 before longer evolution is adopted. The policy and exact free rows,
case/runtime identities and physical gates are in
[checkpoint2.json](cls-cao-tamura-budget-wp2-evidence/checkpoint2.json) and
[numerical-adaptation.json](cls-cao-tamura-budget-wp2-evidence/numerical-adaptation.json).

The exact composite inlet flux is 1.050195996901104 m²/s per unit span, U_bulk is
1.31274499612638 m/s and T_req is 1.1942387693600194 s. Both selected dt=.02 runs
have now completed 60 accepted steps to actual T=1.2000000000000006 s, with zero
rejection or correction supply and passing per-state native physical checks.
Both half-step runs completed 120 accepted steps to actual T=1.2000000000000008 s;
both 11014-slot confirmation runs completed 60 accepted steps to
T=1.2000000000000006 s. All 61 executed common native times pass timestep and mesh
comparisons. Complete dynamic closure has not yet been issued. H2k's hill startup remains preserved;
no hill run is authorized here.

The selected-space dt=.02/.01 s trials reach the same actual T=.04 s with 2/4
accepted steps per model. Every equation and logical bound operation occurs once
per accepted step; solved epsilon is retained, and correction activity/supply are
zero. Differences use at most 3.29% (CLS) and 1.03% (standard) of the allowed
quarter budgets. Continued native steps cost 26–28 s; cold two-step driver runs
cost 120–123 s with peak RSS below 4 GiB. The confirmation mesh's two steps cost
179–180 s total, continued steps about 39–40 s, peak 5.54 GiB. Early selected/finer
differences remain within 8.41% of the quarter budgets.

The final diagnostic runtime explicitly checks prescribed coefficients before
free-row masking and reports the post-k dependent epsilon refresh lag separately.
It changes no evolution arithmetic: both models reproduce old-runtime governing
stage hashes, states, native histories, reference/defect and audit values bitwise
from the same accepted checkpoint. All earlier trial prescribed rows independently
replay with zero error. Frozen runtimes are retained, and successors name their
exact parent checkpoints. The epsilon refresh lag is not removed by a second
solve or Lift.

Both predeclared recovery seeds pass the native admissibility checks at zero
physical steps. A constrained consistent-mass projection adds only outer k on
.012<z<.05, with an actual maximum 4.0492% departure, zero endpoint-K change and
unchanged U, W, pressures, epsilon, reference, defect and boundary data. The TKE
increments are 5.35438e-6 (CLS) and 5.59378e-6 (standard), per unit span/density.
These are setup increments, not positivity-correction supply. At T=1.2 both
recovery legs have 120 accepted steps and zero bound supply. CLS passes its
native final profile budgets and the independent full GL28 profile check, but
its full homogeneity/window/trend assessment still fails. Standard retains a
6.14% native pointwise k error and pressure/W errors. The same accepted
trajectories were extended through T=2.4 and now have a final bounded T=4.2 target; no seed, physical parameter, dt or
acceptance threshold is changed. Recovery is not yet demonstrated.

Use dt=.02 and .01 to T=1.2 s, with observations at .3/.6/.9/1.2; the last .3 s
exceeds T_req/4. Continue the same accepted trials instead of repeating them.
The finer mesh has completed the same physical interval; recovery uses .01. The checkpoint-3
remaining-program forecast was about 4–4.5 h with two serial native cases, including
setup/reporting, restart/rollback and headroom, inside the original 12 h envelope.
The [checkpoint-3 receipt](cls-cao-tamura-budget-wp2-evidence/checkpoint3.json)
separates inclusive step costs, exclusive CPU categories, actual splitting
defects from their algebraic identity, and declared versus executed duration.

The first implementation attempt exposed an immutable-reference dispatch using
candidate coverage assertions. That setup failure is retained on Linux; the
reference now calls the selected kernel directly, and complete hard audits remain
active. No historical witness is reinterpreted by this repair.

Reproduction: `run_cao_tamura_budget_wp2.py` builds common fresh model cases and
charges the aggregate work ledger; `wp2_native.py` runs fail-fast production
restart/rollback fixtures. Heavy artifacts and exact command logs are under
`/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-wp2-20261116`.
See the [predeclared plan](cls-cao-tamura-budget-wp2-plan.md),
[policy](cls-cao-tamura-budget-wp2-evidence/policy.json) and
[literature/source crosswalk](cls-cao-tamura-budget-wp2-evidence/source-crosswalk.json).

The full coefficient-derived, all-state readback checks full 3x3 Reynolds and
mean stresses, complete gradients, homogeneity, absolute pressure, source work
and covariance at 254016 GL28 points per state on the union physical partition.
Independent AWF ReferenceK primitives reproduce native values within 1.18e-16;
source-grid primitive fields/gradients replay native data bitwise. Inner dependent
epsilon and viscosity are recomputed from each actual state's endpoint K, not
a report's terminal transient sampler. Native exact sided/interface inventories
remain complementary to this finite volume sampling. All original raw reports
and earlier analysis receipts are retained.

The unperturbed dt/half-step comparison uses at most 20.15% (CLS) and 5.12%
(standard) of the allowed quarter-budget limits. The maximum selected/finer
fraction over the whole interval is 51.40%, at the initial inner sampled shear.
Full-stress regional RMS errors reach at most .002211 (CLS) and .001525
(standard), normalized by u_tau^2, below .01. The measured final .3 s window
and transit-scaled least-squares trend pass every unperturbed leg.

The actual preserved upstream profiles, with frozen observations and weights,
give full 15-point mean-velocity RMSE .0912971 Uinf for both models, stress
RMSE about .001822 Uinf^2, and streamwise-RMS RMSE .90110 (CLS) versus .77668
(standard), in sigma_u/u_tau units. This is the selected full-height surface-layer
idealization, without a finite-layer tail. Numerical preservation does not
remove its experimental outer-profile discrepancy or calibrate the standard
model's isotropic normal stresses. Fig.3d and hill observations remain held out.

The complete unperturbed final-window check now bounds **every pair of executed
times**, using the componentwise GL28 ranges; it does not rely only on equal
endpoints or a linear trend. All six cases pass. Denser GL28 timestep differences
use at most 21.85% (CLS) and 5.87% (standard) of the allowed quarter budgets;
mesh differences use 56.65%. Native wall/interface/continuity comparison bounds
also pass, using the physical triangle inequality rather than adding component
norms as error percentages. Independently aggregating the two regional GL16
strong-divergence L2 norms reproduces the reported physical RMS within 9.3e-16.

At T=2.4, CLS's native physical profiles pass except the unforced block norm
(.002268); its complete recovery window does not qualify. Standard has k point
error 16.88%, strong-divergence ratio .0015904, wall-shear error 1.264%, and
mechanical-pressure RMS error .13102 in u_tau^2 units. These are failed accuracy
criteria, despite positive scalar/covariance states and zero correction supply.
The read-only assessor now records failed transient accuracy checks and writes
a complete negative receipt; it still requires the same final limits and never
bypasses scalar or covariance requirements. The original aborted assessment log
is preserved.

One conditional standard dt=.005 control completes 60 accepted steps to actual
T=.30000000000000016 from the exact same pulse as dt=.01, with independently
formed histories. First callback history hashes equal the perturbed seed for
both models and both tested timesteps. At T=.3, the dt difference is about 2.82
times the quarter-budget point-pressure limit; field differences are 2.8–10.8%
of the fine trajectory's response. This warns that the transient is time
sensitive. It neither establishes temporal order nor proves a physical model
instability. No extra long timestep trajectory is launched.

The final recovery milestone remains within the original 12-hour envelope, with
no automatic extension beyond T=4.2. Floating-point time comparisons after T=1.2
use a predeclared global-step addition bound, at most (steps+2) ulps of the
expected time; saved times and histories are never edited or interpolated.
Ordinary incompatible-dt restarts remain rejected. Exact completed trajectory
snapshots bind subsequent readbacks, preventing a later continuation from
changing the assessment's input manifest.

The production native source chain matches the final Mac production C++ files.
The one C++ difference is the separately compiled, geometry-only row-export
probe, frozen under its own identity. Python workflow changes are recorded
separately from the unchanged production executable. The closure collector
checks executable and every frozen linked-library hash, input receipt hashes
and exact final state/history/reference bindings; missing or failed recovery
evidence cannot close WP2. Source and runtime records are indexed in
[native-source-chain.json](cls-cao-tamura-budget-wp2-evidence/native-source-chain.json).

Native turbulence budgets retain actual production, cutoff sink and its separate
counterfactual reference production, outer epsilon production/destruction/Yap,
physical storage from saved histories, all physical faces, constrained reactions
and signed free/interface row projections. Regional weak identities close within
1.5e-16 at the retained unperturbed endpoints. The raw face column labelled
`outward_diffusion` is the conormal D*grad(q) dot n; the outward physical
diffusive transport is its negative. Sampled face-gradient balances are distinct
from exact weak row identities. CG interface row projections are not unique
control-volume fluxes, and constrained rows are not unconstrained evolution
equations. The dependent inner epsilon has no invented transport equation/flux.

Useful reproduction entry points, with exact argument arrays retained under the
Linux root, are `run_cao_tamura_budget_wp2.py`, `run_wp2_preservation.py`,
`verify_wp2_long_transaction.py`, `analyze_wp2_approach.py`,
`audit_wp2_dense_profiles.py`, `audit_wp2_dense_comparisons.py`,
`audit_wp2_window_envelope.py`, `audit_wp2_native_boundary_bounds.py` and
`collect_wp2_closure.py`. Every analysis is read-only. Existing output owners
reject duplicate launches; use their recorded unused-root setup for reproduction.
