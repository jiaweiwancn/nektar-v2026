# CT-Accepted-12: actual-target pressure response and campaign feasibility

The actual target has a different leading pressure driver from the CT10 strip: represented inner linear and molecular shear dominate its initial momentum-load pressure image. The first inner radial slab reinforces that image with signed alignment **1.04718**; nearby slabs partly cancel it. All six captured finite-step corrections and absolute-pressure telescopes reconstruct. This diagnoses a numerical departure; it does not pass the unchanged CT11 uncertainty criterion or admit longer preservation.

## Contract, source and retained scope

The starting/fetched remote revision was `b62e59365a82f60d987c51ba1db478ff7314fe0b`, with a clean dedicated Mac worktree and no intervening commits. Source/Git work uses that implementation worktree. Linux numerical work uses the verified `/media/jay/新加卷1` mount and the isolated archive `/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-12-20261115`. Shared installations, unrelated worktrees and all CT4–CT11 seals remain unchanged.

The [prospective contract](cao-tamura-accepted-12-evidence/contract.json) fixes the exact CT8 physical manifest, 60,907 stored slots, 930 affine elements, P4/P2, x=[−25,30], roof=20 and H=.1875. Original/Natural A95, the displaced Accepted/composite policies, finite inner viscosity, zero outer molecular diffusion, current sampled K, transported inner k, dependent epsilon, full raw-gradient signed production, covariance checks, RN top and absolute nonlinear outlet traction are retained. No C++ governing or probe source was changed or rebuilt. The CT9 governing executable and CT11 probe stay pinned by executable and actually loaded-library hashes.

The reconstruction uses helper revision `7bb625021e0605de2c584054fee0d538d214726f`; the completed replay/retention proof uses `60a60ef94316bd98d97a471c90de22d7c9b15197`. The static preflight uses `dd6efdf473ffe03b336cc1342811f22376f074d3`; the bounded axis solve uses `fd6efc7b08eb4ae49763aec3d6520a01156b9527`; its external-traction correction, face control and final collector use `e74b4c0a2f267a161ce71511edd2d8fad65a63be`. Documentation/publication revisions are separate. The three native replay operations are exactly two accepted steps and one deliberate rejection; wrapper recovery repeated none of them.

## Target pressure reconstruction

The [method](cao-tamura-accepted-12-evidence/method.md) uses captured free restrictions, A/G/D/C and the consistent free-velocity mass. C is exactly zero and the unused gauge slot is excluded; all physical pressure tests retain the absolute traction-determined level. Native signs are reproduced by `du*=−A⁻¹rm`, `S dpi=D du*+rc`, `S=D M⁻¹G/alpha−C`. No dense Schur matrix or mass-inverse pressure columns are formed. Signed band loads are streamed.

All-row frozen action differences are at most **1.49e−16**, below 2e−10. Native pressure-increment reconstruction is below **5.64e−16** scaled by u_tau², and signed component sums below **4.05e−15**, against 1e−8. Tight repeats and native telescopes preserve whole/inner/outer means, variations and interface quantities. Work-limit flags on tiny diagnostic right-hand sides remain visible; they are not labelled relative convergence. True residuals and physical repeat/sum bounds determine the reconstruction status.

| Common time | Whole RMS ΔPi/u_tau² | Signed mean | Variation RMS | Inner RMS | Outer RMS |
|---:|---:|---:|---:|---:|---:|
| 0.001 | 9.70030709e-06 | 8.37092166e-06 | 4.90139045e-06 | 1.08088431e-05 | 9.68921043e-06 |
| 0.002 | 3.82461413e-06 | 3.28426733e-06 | 1.95991361e-06 | 4.32810684e-06 | 3.8195322e-06 |

No pressure component is removed for admission. The T=.002 spread remains **3.8246141343e−6 > 5.7111742990e−7**. It is about 0.3825% of the absolute .001 preservation tolerance; that observation does not replace the frozen relative uncertainty test.

| Signed contribution to the T=.002 pair difference | Alignment with total ΔPi |
|---|---:|
| representation/inner/volume/linear_shear | -18.6877943 |
| representation/inner/volume/molecular_shear | -9.09763233 |
| evolution/inner/volume/linear_shear | 13.0473182 |
| evolution/inner/volume/molecular_shear | 6.56908404 |
| evolution/outer/volume/augmented_pressure | 9.11157346 |

Alignment is a signed projection onto the measured difference, not a positive error fraction. Large reinforcement and cancellation are essential: the aggregate initial-representation and evolving-load alignments are −27.7772 and +28.7771. Continuous-reference imbalance and backward-storage/lifting images are negligible at this scale. Initial advection volume/outlet terms cancel with their correct signs. Prescribed reactions are retained separately and excluded from the free-row pressure map. One-sided interface flux observations are not added again to shared assembly.

The independent axis reconstruction gives streamwise/normal/continuity alignments .950085/.049794/.00012034 at T=.002. Whole/outer pressure means are predominantly streamwise-driven, while the normal image remains important for the inner/interface variation. An additional initial-shear check maps its streamwise load to alignment .998590 with the first native increment; the corresponding normal face-load alignment is only .00009781. Thus neither a large normal RHS nor a large interface response identifies the dominant whole-target pressure driver.

This target-specific result does **not** transfer the strip’s outer-normal-stress ranking. The full target has both a different radial grid and a much greater streamwise extent; no isolated domain-length experiment was performed. The evidence supports a first-inner-slab representation hypothesis, not a unique explanation based solely on model stress norms.

The small-step consistent-mass image and the finite-alpha predictor/mass-limit remainder have signed alignments .05395/.945995 at T=.001 and .22668/.773199 at T=.002. The latter is not identical to PIMPLE-versus-coupled splitting error: it includes finite-alpha predictor diffusion. The first exact frozen coupled diagnostic still gives a nonzero pressure correction (whole RMS about 2.4681e−5 scaled); its difference from the actual split correction is about 1.6042e−7 in whole RMS. This supports a represented initial momentum imbalance, with finite-step response and feedback, rather than a free pressure gauge or continuity-only startup defect. These additional solves are read-only.

Continuity’s endpoint alignment is .00012034, and its scaled maximum pressure image is below 8e−10. Actual stage data show that stored Pi changes at the flow update and remains unchanged through epsilon and k solves. Final fine-step fixed-operator images of flow/epsilon/k load changes have physical RMS 1.4811e−8, 1.7322e−11 and 5.2785e−9 respectively; these images describe changed loads, not additional native pressure updates. Imposed dependent-trace defects and the refreshed post-k lag remain separate. No trace repair occurs.

The held-state q48 governing-quadrature control changes the scaled pressure image by at most **4.56e−12**. The q64 aggregate probe fails with `std::bad_alloc` under the unchanged 16-GiB address-space bound. Its request, state and failure remain archived. Reconstruction passes; complete evolved governing-quadrature coverage remains blocked. This is separate from the four passing streamed physical-scalar controls.

## Replay, physical diagnostics and output retention

The two detailed replay steps reproduce CT11’s accepted states, all numeric flow/k stage primitives, CSR matrix bytes and genuine final native history exactly. The sole deliberate rejection retains the actual previous accepted history and candidate/stage witnesses. Recovery fixed directory/empty-checkpoint bookkeeping only. All native physical checks remain enabled, extra governing sweeps/corrections remain zero, and the actual physical horizon stays **T=.002**.

All native/q48f16/q64f16/q64f32 scalar controls pass, each containing inner_k, outer_k, outer_epsilon and whole_k plus both k/outer-epsilon reconstructions. They consume the actual consecutive, corrected-flow, epsilon and pre-k states. Conservative advection, dependent-boundary storage, signed production and actual physical ground/interface/external gradients are retained. No scalar*div(U) correction or reaction substitution is used.

Independent physical-face/volume Stokes reconstruction of the initial normal inner-shear difference agrees on all rows to **6.15e−19** and in pressure image to **6.95e−16** scaled. Reversing outlet orientation gives a resolved 1.1788e−5 maximum pressure-image defect. Actual face gradients are used; reactions never replace them.

The raw CT11 fine endpoint ground error is about 9.12076e−5 and full interface-traction RMS .000807998, leaving approximately .000192002 of the .001 interface allowance. The last fine step changes mean ground shear by −1.63068e−8 through the flow-gradient update and −3.25878e−9 through the post-k endpoint-viscosity update, with zero cross term. Epsilon transport changes neither of these ground quantities. The stage product identity closes below1.56e−16 scaled. Corrected total outlet reconstruction RMS decreases from8.44223e−5 initially to7.31824e−5 at the endpoint; top tangential mismatch is1.98295e−7. The inner-outlet-only RMS is.000751380 at the endpoint. These are physical constitutive reconstructions, separate from prescribed natural data. The first new external-traction postprocessor wrongly used D(0) at finite height; its inner-outlet values are explicitly invalidated and preserved in the correction notice. Correct radial D(r;K_H) agrees with the pre-existing independent traction implementation within9.49e−19 normalized. Axis/ground/interface results were reused unchanged. Top normal traction is diagnostic only; w remains essential and no top normal-pressure datum is added.

`LosslessSegmentCASV1` is an external retention controller, not a solver cache. Its 47 newly produced payloads map to 33 unique lossless blobs. It measures 555,786,201 raw bytes, 266,659,550 retained bytes, and 647,791,478 bytes at raw-plus-archive peak. Compression/validation/pruning costs 10.42 process-CPU seconds and verified restoration 2.95 seconds. Missing blobs, altered plaintext hashes and foreign native bindings reject before any unlink. Every restored byte matches. Required native reports, histories, configurations and audit receipts remain directly readable. The demonstration restores the new raw payloads afterward; historical archives are never pruned. Online rolling-window eviction and compressed native-checkpoint reading are not qualified.

## Measured performance and admission

| Exclusive accepted-step category | Step 1 CPU s | Step 2 CPU s | Calls per step, approximately |
|---|---:|---:|---:|
| Element-row construction | 88.6604 | 89.3899 | 34.606 million |
| Domain computation, excluding nested timers | 72.7030 | 73.6320 | 5 |
| Basis lookup/evaluation | 50.9432 | 52.1682 | 69.211 million |
| Evaluation-row work | 35.4161 | 35.4937 | 25.933 million |
| Flow assembly, excluding nested timers | 24.1008 | 25.0664 | 1 |
| Geometry location | 22.3539 | 22.3503 | 11.381 million |
| Owner catalog lookup | 17.4803 | 17.6138 | 21.907 million |

These are exclusive categories; inclusive totals are never added. Detailed accepted-step totals are 359.1161/364.1897 seconds versus CT11’s uninstrumented 143.8562/137.9906. The 2.566-fold instrumentation cost prevents treating detailed geometry fractions as achievable clock-off savings. The second native replay uses 528.164 process-CPU seconds including 163.697 exclusive startup seconds, 531.050 elapsed seconds and 3,762,260 KiB peak RSS. The first native run’s temporary RSS/library-map object was lost in a wrapper failure; its detailed times and exact state/stage proof survive. Second-step and audit actual library maps are retained.

No governing reuse optimization was selected. The observed call counts identify immutable row/owner work for a future controlled clock-off comparison, but no safe measured change has demonstrated the roughly 10.83-fold complete-campaign CPU improvement required by the existing plan. ExactStepV1 remains guarded and all state-dependent gradients, cutoff/endpoint values and mandatory validity checks remain live.

## Single static contrast and precise stopping gate

The [prospective selection](cao-tamura-accepted-12-evidence/static-first-inner-v1.json) bisects [0,H_A/64] at H_A/128=.00146484375 H_h (0.05859375mm). It preserves every original x/r node, P4/P2, the physical wall/interface, roof, projection, quadrature, physical functions and fit. The same logical amendment is recorded for the matched terrain definition; no terrain field count is claimed.

Native construction confirms **62,533 slots**, **961 elements**, and block sizes
[15625,15625,15625,11625,1071,2961] plus the unused gauge slot. Fresh zero-step
initialization/history succeeds. The first frozen predictor has30,628 rows,
2,152,435 matrix nonzeros and **40,497,972 factor entries**, exceeding40M.
The native guard stops this one contrast. Its exception has the inherited
“EP-AWF patch” prefix; the actual configuration is displaced A95 Accepted/RN,
and no EP route was enabled.

Later detailed reference/operator, alpha2000 and quadrature controls remain
unattempted after that first guard. Initial raw native residual norms are
retained, but they are not a common-test convergence qualification. No evolved
contrast, second mesh, altered pressure initialization or increased cap is used.
Whether this refinement resolves the pressure uncertainty remains untested.

## Closed resource and feasibility decisions

The [evaluated decision](cao-tamura-accepted-12-evidence/evaluated-v1/decision.json)
keeps reconstruction, quadrature coverage, physical qualification, retention,
static capacity and campaign admission separate. The closed-census maximum is
**14 project /18 aggregate workers**, including helpers, below16/72. Numerical
libraries use one thread, native execution is one MPI-linked serial rank, and
no compilation was required. The finite CT12 budget is12,000 process-CPU
seconds,14,400 active-wall seconds and12GiB output. The additive
[closed publication accounting](cao-tamura-accepted-12-evidence/closed-publication-accounting.json)
includes all47 closed owned families: **9,598.85 process-CPU seconds** and
**7,650.42 seconds active-wall upper bound**. Measured output before packaging
is2,682,848,976 bytes; a64MiB publication reserve bounds the final package below
2.750GB, well within12GiB. This reserve is labelled separately from measured
native file sizes. Mac source/Git work and transfer are not native CPU benchmarks.
Suspension events and conservative active-wall accounting are retained.

The unchanged evaluated receipt precedes the collector/publication closure and
records9,575.43 process-CPU seconds,7,623.35 seconds summed active-wall upper
bound and2.679GB new run output. It is not overwritten by the final supplement.
The bounded axis lifecycle repair has3,793,256KiB RSS /5,508,716KiB VmPeak;
the earlier diagnostic allocation failure remains archived. The independent
q64 aggregate failure remains a capacity block, not a model-domain failure.
No historical resource failure is relabelled.

| Prospective .407 campaign quantity | Current conservative forecast | Ceiling | Decision |
|---|---:|---:|---|
| Process CPU, spent plus remaining | 269,570s, including closed publication | 24,000s | no-go |
| Active-wall upper bound | 266,511s, including closed publication | 28,800s | no-go |
| Output with measured CAS policy | 176.31GB planning upper bound | 24GiB | not admitted |

The forecasts retain1215 additional pair steps,39 restart segments and25%
headroom, including startup, controllers/helpers and required audits. The output
estimate conservatively applies measured one-step native retention; amortizing
native-report headers over32-step segments or online eviction is not qualified
by this two-step replay. Retention works byte-for-byte but does not supply a
positive campaign storage plan. Actual restored replay output contains158.69MB
native reports,193.37MB sparse matrices,60.45MB stage snapshots and133.84MB audit
packets. None of these checks may be removed to claim savings.

Even the optimistic calculation that removes all measured domain work from the
old CPU forecast leaves about122,520s, far above24,000s. No achieved governing
CPU saving is claimed. Bounded reconstruction, output compression or a different
pressure interpretation cannot authorize an unaffordable campaign.

| Scope | CT12 outcome |
|---|---|
| Six target corrections, signed loads, telescopes, consistent-mass/G actions | pass |
| Independent initial-shear face/volume identity and orientation negative | pass |
| Two accepted replay states, stage/CSR equality, genuine restart and rollback | pass |
| Four-control/four-equation scalar audit and both reconstructions | pass |
| Limited lossless payload retention and restoration | pass |
| Evolved q48 governing-quadrature sensitivity | pass |
| Complete evolved q64 aggregate diagnostic | blocked at16GiB |
| CT11 absolute-pressure timestep criterion | historical fail, unchanged |
| First-inner-slab static contrast | blocked by40M factor cap |
| Fresh resource and finite package budget | pass, with final closed accounting |
| Longer full-domain preservation | unattempted / not admitted |
| Terrain, A93/E94E95 displaced transfer, EP | unattempted / outside scope |

The actual physical paired horizon remainsT=.002, approximately4.92348e−5 of a
transit. CT12 replay does not add physical horizon. The .407 first stage and
later .1/1-transit stages remain unstarted. No new uncertainty threshold is
adopted: a future criterion must allocate an accuracy budget prospectively and
use fresh qualification evidence; CT11 is not relabelled.

## Next bounded implementation

The smallest concrete capability now needed for the selected causal test is a
native **no-solve frozen CSR export before factorization**, followed by a bounded
fill-ordering assessment on that exact62,533-slot state. Keep70k/40M/16GiB and the
old CSR component unchanged. A separately versioned predictor ordering may be
qualified only after complete action/update equivalence; it must preserve
consistent mass, the full pressure action, absolute traction and one sweep.
This is a narrow factor-capability task, not a general backend rewrite.

Only then can the already selected first-inner-slab contrast finish static
admission and receive one newly contracted short paired causal trial. No further
mesh or timestep search is justified. The independent evolved q64 observer
allocation also needs a bounded phase/streaming correction. Long-campaign CPU
and output readiness remain separate gates; reduced factor fill alone cannot
supply a10.8-fold throughput improvement.

Before a rough-hill pilot, target flat pressure/preservation and affordability
must pass, together with mapped region-aware numbering, joined-map Cartesian
gradients/cofactors, conservative transmission, and the outstanding curved
spatial/temporal qualification. Roughness-sublayer/interface adequacy,
thin-layer accuracy, smooth preparation and approach-layer fidelity remain
physical decisions. Numerical flat preservation is not experimental hill
validation. See the [ready-to-use next prompt](cao-tamura-accepted-12-evidence/next-package-prompt.md)
and [archive/reproduction index](cao-tamura-accepted-12-evidence/README.md).
