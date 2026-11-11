# Directional polynomial order and independent A93 resolution qualification

The directional capability passes its declared native tests, and the previously blocked A95 streamwise experiment now has eight fresh trajectories through T=.016. Increasing streamwise velocity/scalar degree from five to six and compatible pressure degree from three to four reduces all four **perturbed total** scalar-budget defects by 67–80%. Complete budget accuracy still fails: the neutral baseline retains its boundary-normal defects. Several response observables also remain spatially unresolved. No longer response is admitted.

The independent A93 mesh candidate passes its native spatial-residual reductions at .016, including common vertical momentum, but fails the unchanged outer-epsilon field-refinement gate. It stops there. Its signed common-w load comparison is also temporally unresolved. The historical A93 T=.5 failure therefore remains open.

The actual starting revision was `d0182aa8df9b3ebb03577811ce658ee93a594f45`, clean on `codex/a93-a95-e94e95-awf`, with no intervening commits. All historical evidence is unchanged. The adjacent [evidence directory](a93-a95-e94e95-awf-directional-resolution-qualification-evidence) contains the versioned contract, complete result bodies, immutable requests/completion seals, native accepted/stage receipts, primitive classifications, runtime identities, failures, commands and hashes. Heavy data remain on Linux at the destinations below.

The contract is `DirectionalResolutionExperimentV1`, SHA256 `bab9ace03e97a3d4f7b39e06f89f76dfed14a9aa0ed4bbb9a771bd0d6fbc26a2`. `DirectionalInitializationPhaseV1` was recorded before either A93 fine trajectory advanced. Initialization requires construction, physical/field admissibility, analytical weak balance and common-space inclusion; initial .85 ratios are diagnostic. The original .85 gates still apply at evolved milestones. The first zero-step attempt that incorrectly applied an evolved reduction gate at initialization is retained separately.

The closures, coefficients, natural/prescribed kappa resolution, wall laws, high-Re outer molecular policy, finite inner viscosity, total-epsilon semantics, signed production, augmented pressure, one-sweep order and admissibility floors are unchanged. [source-scope.json](a93-a95-e94e95-awf-directional-resolution-qualification-evidence/source-scope.json) records identical hashes for the model/wall kernels, main PIMPLE ordering, physics signature, resolution-term implementation, initialization solver and both specifications. `GeneratedVolumeOwnerV1` remains the operator. Eta and corrections/bounds remain zero/disabled; no reference cancellation, reset, trace repair, clipping or regularization was introduced.

The supported new layout is uniformly directional Modified_A on canonical affine serial flat 2D, velocity/k/epsilon (px,pz)=(6,5), pressure (4,3), under the existing new-model one-sweep route. Arbitrary elementwise p, anisotropic 3D, curved/MPI, legacy action routes and systems above 3000 unknowns remain guarded.

`RANSWallGlobalField` carries per-axis degrees, counts and strides for native packing, signed assembly, vertices, values, first and mixed/radial second derivatives, faces, restriction and export. Fixed-x faces use the z basis; fixed-z faces use the x basis. Ground traces use px; radial projection/injection uses pz. Diagnostic inventories use the relevant directional degree. Row caches retain field, generating geometry, native element, side and coordinate identity. Corrected common tests require containment and unchanged physical coordinates; the old tolerance path is explicitly legacy read-only analysis.

Directional preservation identity adds each axis's basis/modes independently of quadrature. Native histories remain strict about each axis's basis **and points**. Isotropic identity strings and arithmetic are retained. The operator version was not changed solely to admit this separately guarded layout.

The continuous ground-k test converts Modified_A traces through degree six to polynomials and checks endpoints and every relevant real derivative root. Recursive derivative-root isolation includes repeated/nearly repeated stationary points. The arithmetic guard is 512 long-double eps times the modal L1 bound, covering conversion/evaluation error; uncertain tiny interior margins reject without altering k. Native conversion checks, positive/negative sextics, a narrow negative pocket missed by 81 ordinary points, endpoints, repeated roots and tiny margins pass. Existing p3/p4/p5 behavior is retained. This is a continuous ground trace check, not a global positivity or realizability result.

Native tests compare custom rows and signed restriction with independent Nektar BwdTrans/PhysDeriv and mass integration: every basis mode, constants, affine/mixed polynomials, dx/dz/dxz/dzz, four faces and exported outer coefficients. Maximum primitive defects are 2.85e-13 or smaller against per-quantity limits (3e-12 values/first derivatives/mass, 2e-10 second derivatives, 5e-12 polynomials/faces/restriction). P5 and X6Z5 radial projections are bitwise equal. The small two-region test uses its own admissible manufactured reference and boundaries; it does not claim a coarse logarithm is an equilibrium.

| Small test layout | Actual unknowns | Gauged mixed rank | Divergence / pressure-action rank | Ungauged PIMPLE pressure rank |
|---|---:|---:|---:|---:|
| P3/P1 | 109 | 29 | 7 / 7 | 7 |
| P4/P2 | 179 | 61 | 17 / 17 | 17 |
| P5/P3 | 267 | 105 | 31 / 31 | 31 |
| X6Z5 / X4Z3 | 314 | 131 | 39 / 39 | 39 |

Mass-Gram whitened SVDs use relative threshold 1e-11 and the actual constrained PIMPLE pressure action D M_free^-1 G-C. There is only the expected common pressure constant before gauging, with no extra null mode in these finite cases. This is not a general inf-sup theorem. The production ownership regression passes all 11 controls, including permuted IDs, tiny/moving-cutoff panels, exact-H sides, nearby interiors, shared owners, partition/gradient/volume/Stokes and free-row constant-flux/gauge identities.

Old/new isotropic p3 and p5 frozen rows, derivatives, masses, unanchored residuals, signed coefficients, exports, represented references and identities compare exactly. Directional 2769-unknown fresh initialization, uninterrupted/restarted evolution, cache on/off and deliberately rejected-step rollback give bitwise equal states and complete native histories where equality is expected. Focused ordinary CLS cache/restart/rollback and pinned isotropic checkpoint restoration pass, as do affected admissibility, interface and ownership regressions. Required-scope and budget negative controls pass, including each shared prerequisite failing while all A95-specific labels pass. Reconstruction cannot substitute for accuracy; missing/foreign/stale receipts and forged primitive pass labels reject, and required numerical failures return exit 1.

Unequal-point restart negatives initially crashed the SDK's automatic MatrixFree interpolation timing before checkpoint reading. Those failed jobs remain failed. Separate negative fixtures pin the existing NoCollection path with both session .opt sidecars; changes to either axis then reject for the intended checkpoint-signature reason. Scientific layouts retain their qualified 8x8 native points. Unequal-point automatic MatrixFree startup is unqualified.

The frozen scalar reconciliation uses sealed old Accepted baseline/response states and the actual one-sweep stages: epsilon uses corrected flow with old k/epsilon coefficients; k uses corrected flow, old k and newly solved epsilon. Both transport newly solved scalar values/gradients, including dependent epsilon-boundary storage. Backward storage, frozen RHS/reaction, conservative advection and diffusion reproduce the actual unmasked native load to 5.56e-15; free rows are below 2.69e-15. Advective integration by parts closes to 5.11e-15. Independent physical-flux derivative h/h2 controls close diffusion integration by parts to 7.36e-15 (limit 2e-10). No implementation inconsistency was found. Scalar times div(U) appears only in the conservative-IBP diagnostic, never as a budget correction. Prescribed rows remain algebraic reactions, not independent conservation evidence.

All eight A95 histories are new and bound to the rebuilt runtime. Both layouts share B5 vertical nodes and the physical boundary specification, with their own projected references. Same-space dt companions have identical actual initial coefficients; perturbed hashes are distinct from neutral reference hashes. Raw/effective and source/target boundary differences are zero; explicit/implicit live boundaries agree. Actual StreamfunctionPolyV1 initial divergence is at roundoff. The common physical grid resolves shared partitions/profile breaks, with explicit interface sides and one global pressure gauge.

| A95 layout | Actual unknowns | Fields / pressure | dt values | Endpoint |
|---|---:|---|---|---|
| B5-IsotropicP5 | 2355 | (5,5) / (3,3) | .002, .001 | .016 |
| B5-UniformX6Z5 | **2769** | (6,5) / (4,3) | .002, .001 | .016 |

The legacy member key “B11” here means B5-UniformX6Z5, **not** Balanced11's vertical mesh. Neutral and amplitude 2e-5 trajectories are run for each layout/timestep: eight fresh histories, 8 or 16 steps each. Wall-normal quadrature is 48, tangential/mass/outer/interface orders 12; native basis points are 8 per axis. Frozen controls compare volume orders 48/64 and physical-face orders 16/32, with all other actual orders recorded. T=.016 means theta=.0032 for H=2 and u_star=.4.

All eight trajectories pass absolute field/traction accuracy, domain, source/storage and no-correction gates through .016. Across accepted states and audited stages, k stays above .6115196, outer epsilon above .03959798, covariance minima above .1772641, actual A95 radicand above .6917658 and required cutoff ratio above 5.89870. Bound calls/activations are zero. Maximum refreshed dependent-trace lag is 4.79e-8 at .002 and 2.57e-8 at .001; it is recorded as splitting and not repaired. Source/storage defects stay below 2.25e-16 and corrected weak continuity/gauge below 9.10e-17. Locations, sides, times/stages and separate final extrema are retained.

Neutral absolute fields and timestep uncertainty pass. Changing only streamwise orders is not a balanced radial-refinement qualification: neutral value errors remain limited by the unchanged vertical representation. Full/native and fixed-common alpha-zero residuals remain visible alongside stage residuals. Both neutral and perturbed native audit sets pass: analytical weak balance, residual reconstruction, shared/interface rows, scalar accepted-state coverage, frozen quadrature and model-specific constitutive/traction comparisons.

The 17 response observables use perturbed minus matching unperturbed numerical fields, separately from analytical-reference error and actual-initial drift. The table gives scaled weighted physical RMS differences. Spatial limits are max(5e-7,.1*maximum response RMS); temporal limits are max(1e-7,.05*maximum response RMS), unchanged. “Floor” is not a measured order.

| Observable | Spatial difference | Spatial status | Largest dt difference | Temporal status |
|---|---:|---|---:|---|
| U | 9.3454e-7 | pass | 6.3192e-8 | pass |
| w | 3.3058e-6 | pass | 9.7292e-8 | pass |
| k | 1.2558e-7 | pass; signal below floor | 5.9346e-9 | pass |
| outer epsilon | 1.1297e-7 | pass; signal below floor | 1.0170e-8 | pass |
| Pi | 9.8854e-5 | unresolved | 3.0036e-6 | pass |
| mechanical pressure | 9.8876e-5 | unresolved | 3.0068e-6 | pass |
| Tzz | 1.0061e-4 | unresolved | 3.1636e-6 | pass |
| Txz | 2.4256e-5 | unresolved | 7.7459e-7 | pass |
| Rxx | 1.5058e-5 | unresolved | 6.6905e-7 | pass |
| Ryy | 8.0635e-6 | unresolved | 3.9050e-7 | unresolved |
| Rzz | 1.7113e-5 | unresolved | 7.9227e-7 | pass |
| Rxz | 3.1810e-5 | unresolved | 1.0674e-6 | pass |
| Qxx | 1.5210e-5 | unresolved | 5.2179e-7 | pass |
| Qyy | 0 | pass; signal below floor | 0 | pass |
| Qzz | 1.5210e-5 | unresolved | 5.2179e-7 | pass |
| Qxz | 1.0285e-5 | pass | 7.5767e-7 | pass |
| full signed P | 5.4916e-5 | unresolved | 2.0844e-6 | pass |

Read-only decomposition reproduces native Ryy and production to 5.56e-17. For this A95 flat-2D law, Qyy=0 and Ryy=2k/3+2 nu_t tr(A)/3 under DeviatoricStrainV1. Its response's divergence part has RMS 7.300e-6 / 4.270e-6 on the primary layouts, versus scalar parts 7.017e-7 / 7.453e-7; divergence-part spatial difference is 8.068e-6. Weak continuity at roundoff does not imply pointwise tr(A)=0. This explains the measured Ryy sensitivity without modifying the extension. Full production keeps isotropic work -2k tr(A)/3. A95's commutator Q has Q:A=0 analytically; its observed contraction is below 1.3e-20 in response units despite resolved perturbed Qxz. Neutral Qxz=0 and neutral Rxz=-u_star^2 remain distinct identities.

All 16 physical-budget reconstruction leaves complete. Volume48/64 and face16/32 uncertainty checks pass. The signed frozen-stage defect is backward storage minus frozen source minus independently integrated outward (D grad(f)-U f) flux, divided by region volume times K0*u_star/H or epsilon(H)*u_star/H. These are integrated rate diagnostics, not field RMS. Whole-patch external totals exclude both internal-interface sides.

| Kind / equation | B5 defect | Directional defect | Absolute ratio | Spatial / temporal |
|---|---:|---:|---:|---|
| Perturbed inner k | 1.22755e-5 | 4.00221e-6 | .3260 | pass / pass |
| Perturbed outer k | 1.42630e-5 | 2.79557e-6 | .1960 | pass / pass |
| Perturbed outer epsilon | 1.38970e-5 | 3.22639e-6 | .2322 | pass / pass |
| Perturbed whole-patch k | 1.32304e-5 | 3.36008e-6 | .2540 | pass / pass |
| Neutral inner k | -1.70299e-6 | -1.74120e-6 | 1.0224 | failed / unresolved |
| Neutral outer k | -3.10247e-9 | -1.88424e-9 | .6073 | unresolved floor / pass |
| Neutral outer epsilon | 1.46403e-7 | 1.52834e-7 | 1.0439 | failed / unresolved |
| Neutral whole-patch k | -8.22448e-7 | -8.40351e-7 | 1.0218 | failed / pass |

The .85 spatial ratio, 5e-9 residual floor, quadrature uncertainty and dt<=max(1e-9,.1*spatial difference) criteria are unchanged. Accuracy classification and exact configuration/amplitude/horizon/four-member binding both exit 1. Baseline subtraction is additional evidence and does not replace these total-budget gates. A fresh .001/.0005 pair is not triggered: required neutral spatial reductions already fail, in addition to response spatial uncertainty.

A separate probe-only runtime exports advective/diffusive face components without changing any original total. All eight highest-order replays reproduce the scientific runtime's original budget object **exactly**, including coordinates/weights. Component integration agrees to 1.97e-14 scaled; pointwise net-flux reconstruction is exact. In neutral inner k, physical ground diffusive flux is 7.85207e-7 / 8.11249e-7 dimensionally; it is not forced to zero. Streamwise advective faces cancel but are each about 6.0182. In the perturbed state the two inner streamwise diffusive faces sum to 1.01563e-6 / 2.36771e-6 while the ground term remains near 7.85e-7 / 8.10e-7. Storage/source are separate. Together with the old endpoint control, this supports distinct streamwise-response and boundary-normal baseline limitations. Velocity, scalar and compatible pressure discretization change together, so this does not establish a scalar-only cause.

A93 uses only RoughnessTransferV1: remove B11-only node 2.75 from (2.5,3), add .1375 in (.1,.175), preserve every B5 node/active physical breakpoint, and retain p5/p3. Prior assembled signed common-w supports and varying-flux terms motivated the target, rather than large cancelling isolated element norms. Actual count is **2967**. Initial nesting, gradients, mass congruence, fixed rows and field/domain/weak gates pass. Exact sealed old B5 requests/histories are reused; both fine histories are new on the pinned GeneratedVolumeOwnerV1 runtime.

At .016 common-w falls from 8.03578e-6 to 2.24990e-7 (ratio .0280). Full-w ratio is .0506; common-u .00428; common-k .00306; inner-only common-k .00302. Required residual reductions, native audits, quadrature and field timestep bounds pass. But outer-epsilon reference RMS changes from 6.4376344e-9 to 6.2771666e-9, ratio **.9750735**, above .85 and above its unchanged 1e-10 field-refinement floor. This is an accuracy/refinement failure, not a physical step rejection. Removing the outer 2.75 split leaves that outer approximation equal to B5 and gives almost no epsilon refinement benefit. No .1/.5 candidate advancement or second mesh search was performed.

The assembled signed-load audit retains actual row supports and pressure/stress/advection/constant-versus-varying flux terms. Common-w spatial load difference is 8.04001e-6; B5/fine timestep differences are 1.14766e-6 / 1.18833e-7, against limit 8.04001e-7: **temporally unresolved**. Current accepted stage-history storage is 8.16577e-6 / 8.83974e-8 and current stage residual 8.08928e-7 / 1.94348e-7, with load identities below 2e-23. Free constant-flux norms are about 1.3–1.7e-13. This distinguishes finite-time storage/splitting from alpha-zero spatial load; no isolated stress norm is presented as an uncancelled governing error. Promising .016 w reduction does not close the original .5 question.

Linux is authoritative. The Mac task worktree is `/Users/jay/Documents/nektar-v2026/.worktrees/a93-a95-e94e95-awf`. New heavy work uses `/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-directional-20261111` (R), with isolated source/build/private libraries/immutable runtimes/tmp/runs. The old campaign is `/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3` (OLD). Existing CLS checkouts, installations, controllers and checkpoints were preserved. GCC 11.4 Release uses -fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS. The SDK is MPI-enabled but every qualifying run is serial with numerical-library thread counts one.

| Runtime purpose | Compiled source | Executable SHA256 |
|---|---|---|
| New scientific solver | 11ace15e28684cc9314c3192dec0d1c2d85a8b2d | db28f6f91264d918f2d552d7a2fd032a1914481009613bbf0c0df8ed7660e798 |
| New scientific preservation probe | same | c3aa1c23e5f211d32ba18cdd60abd6f0cb14668d77ce92d2466a9eca90b9e913 |
| Pinned A93 solver | 7cb6757e06e6b1a619b01a98a2475dddb40a293d | 99828f96dcc95e0d529f315a136e72cec63b62f8e500e6c9017aa71712d4f914 |
| Extra physical-face probe | 79cfa8dc, probe-only changes | bda4a990715de5fbdd0d23200469f6006f7864e47703a8cb0257e1b13f9018d4 |

The signed-export auxiliary probe was built separately at 40e653d; its identity and sole test-header difference are recorded. Later commits change diagnostic exports/Python/reporting, not the frozen scientific runtime. Complete compiled source/header hashes, raw linker identities and actual /proc loaded maps are included. The unchanged model probe SHA is fd7d3346e68fd67138f4755b3f5f6c0513927a2545651d63213d98d1030d85b0. Historical T=.5 and response evidence remain attached to OLD. The eight unresolved historical epsilon-order observations are unchanged; timestep pairs here bound uncertainty, not order.

The owner-aware guard reserves 16 startup workers, counts helpers/diagnostic threads and pauses only owned families. Maximum observed aggregate was 65, below 72; no debugger was used. One owned family's child processes were suspended for 4.27 seconds and its controller for 122.77 seconds; raw events are retained. Solver CPU is separate from elapsed time and wrapper/sampling/audit work.

| Work | Process CPU seconds | Elapsed seconds |
|---|---:|---:|
| A95 neutral isotropic, 8 steps | 58.86 | 61.36 |
| A95 neutral directional, 8 steps | 67.95 | 70.08 |
| A95 perturbed isotropic, 8 steps | 60.04 | 63.57 |
| A95 perturbed directional, 8 steps | 69.33 | 73.29 |
| A95 perturbed directional, 16 steps | 122.97 | 128.48 |
| A93 new fine, 8 steps | 67.46 | 72.73 |
| A93 new fine, 16 steps | 130.20 | 154.07 |
| A95 finishing audits/evolution/budgets controller | 1563.39 | 1724.50 |
| Eight read-only physical-face replays | 302.30 | 312.20 |

Controller/child rows overlap and must not be summed. Face replays additionally waited 148.46 seconds for sealed inputs, recorded separately. Final scientific build CPU was 81.40 seconds, wall 106.91. One successful probe compile had a failed packaging controller (wrong build/bin path); its lost captured CPU cost is marked unavailable. Build/link attempts, coarse fixture inadmissibility, MatrixFree startup failures, runtime-alias preflight failure and derived-diagnostic subtraction-order failure remain archived. The alias fix canonicalizes the identical executable/library path. A separately sealed supplement reuses neutral endpoints and advances only already initialized perturbed histories; no failed receipt was overwritten.

Exact command arrays are in execution receipts and copied controllers. For the retained Linux environment set R and OLD to the paths above and S=$R/source/solvers/IncNavierStokesSolver/Tests/RANS. Use fresh output labels and start the owner-aware resource guard before native work.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$R/tmp" LD_LIBRARY_PATH="$R/runtime/lib"
ctest --test-dir "$R/build" -j 1 --output-on-failure \
  -R '^(RANSNativeDiagnosticChecks|RANSDirectionalRequiredScopeContract|RANSOwnershipRequiredScopeContract|RANSBudgetAccuracyContract)$'
python3 -B "$S/check_rans_directional_layout.py" --runtime "$R/runtime" \
  --output "$R/reproduction-layout-new-label"
python3 -B "$S/audit_rans_physical_scalar_budget.py" --root "$R" \
  --label A95-N-Accepted-directional-dt-0.002-response-B11 \
  --time .016 --quadrature-index 2 --face-order 32 --runtime-name runtime \
  --output-label reproduction-directional-Q64-F32-new-label
python3 -B "$S/classify_rans_physical_budget_accuracy.py" \
  --source "$R/directional-resolution-supplement-v1/budget-analysis/results.json" \
  --output "$R/reproduction-budget-classification-new-label"
```

The last command must exit 1 for these primitives. Fresh A95 evolution uses `run_rans_budget_axis_control.py` with the exact sealed mesh/hypothesis, binding, layout, amplitude and dt command arrays. A93 uses `run_rans_a93_roughness_transfer.py` with both contracts and OLD; its first evolved failure is reproduced by classifying the .016 endpoint, not resuming a failed gate. Read-only probes may evaluate old coefficients, but no history's dt, signature, boundary or state is edited. `collect_rans_directional_resolution.py --root "$R" --old-root "$OLD" --output <new-directory>` also exits 1 while preserving the complete failed required-scope export.

Demonstrated: guarded directional capability, continuous sextic ground checking, isotropic compatibility, ownership/history/cache/rollback behavior, frozen scalar reconciliation and the declared bounded native experiments. Failed or unresolved: response field uncertainty, complete physical-budget accuracy, A93 epsilon refinement, A93 signed-load temporal uncertainty and historical A93 .5 reduction. Unattempted: second mesh search, additional dt cohort, A93 .1/.5 candidate continuation, A93 perturbation, longer response, prescribed-kappa directional campaign and all excluded capabilities.

The streamwise hypothesis is supported for perturbed budgets, but is insufficient for neutral baseline defects. A93 target supports matter at .016, while the donor tradeoff prevents duration admission and does not establish a .5 cure. The next bounded package should first screen combined directional/endpoint-normal redistribution while preserving important outer representation: adding both endpoint splits directly predicts 3044 unknowns and is forbidden, so a measured nonphysical donor is needed below 3000. Keep horizon .016 and total-budget/field gates. Separately, an A93 donor must preserve outer-epsilon improvement before another .5 attempt. Measured Ryy sensitivity makes pressure/pointwise-divergence resolution explicit before response extension. No longer response or open-boundary work is admitted; a future multi-column neutral inflow/outflow capability still precedes terrain validation.
