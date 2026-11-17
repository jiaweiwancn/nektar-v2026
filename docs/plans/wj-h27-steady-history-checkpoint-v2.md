# Approved-v2 steady, MMS and history checkpoint

This checkpoint preserves substantive numerical repairs and their Linux evidence.
It does not close H27 or authorize Task3. The user-approved tensor remains
`FiniteNuIsotropicWallFullWJCollarV2`; inner k diffusion remains
`nu+Cd*k*k/(sigma_k*Ew)`. No further physical amendment is active.
`TASK3_READY=NO`.

## Tested source and runtime

The source-v11 numerical identity is
`dfdd5d8fb2384d08a5a54e5fd8d76f0ec7e37334fbebeb3b16e439553d8c080b`.
The immutable Linux source capsule records parent revision plus file hashes;
all 105 source-v11 overlay files match this worktree except two test-only
serialization fixes, identified below. The production C++/headers match exactly.
This commit binds the tested overlay; it is not a claim that the build was made
from this as-yet-uncreated commit. The previous source-v9 checkpoint is
`0802c605b54e1d89485b62e688b542a56c2f5114`.

The evidence directory is
[`wj-h27-linux/native-steady-history-v2c2`](wj-h27-linux/native-steady-history-v2c2/inventory.json).
It contains source-v10/v11 manifests, runtime/binary/shared-library identities,
compile/link commands, prospective inputs, independent reports, compact native
results, all operation receipts, failed-operation logs and an OPEN account.
Every compact result points to its complete Linux raw artifact and SHA256.
Large fields remain under
`/home/jay/build/wj-h27-coupling-r1-20261117`; no raw scientific state was deleted.
The archive is 2,598,016 bytes, SHA256
`1f087acebd8493b33c2ea1fbdcc8d24b0fbcd8ea176ba61e9e761ba5f1f92899`.

Build parallelism was four; source-v11 build elapsed 217.823 s. All numerical
runs used one MPI rank with OMP/OpenBLAS/MKL/BLIS/VECLIB/NumExpr threads one.
The shared-host watchdog enforced 72 actual workers, science RSS 8 GiB, storage
caps and free-space reserve. The installed Nektar libraries and unrelated Linux
checkout were preserved. Mac checks were hashes, source/AST and Git only.

## Numerical repairs

* The positive directional shear split retains the complete original stress
  remainder and current equations, including weak divergence. The former
  isotropic-mu split remains a control. Full-row direct scalar/stress actions
  are checked against assembled actions, including essential reactions and
  the nonzero augmented pressure multiplier. BDF residuals form state
  differences before mass action to reduce cancellation.
* Scalar Krylov and final acceptance now use the same explicitly evaluated
  current-matrix residual. This fixes a reproducible rounding disagreement;
  it does not loosen a physical residual gate. A separate micro-dt case still
  hits a true linear roundoff floor and is preserved as a failure.
* Full current principal validation precedes steady nonlinear-state evaluation
  and retention. Derivative probes retain root/covariance/domain guards; fatal
  base or unrelated exceptions remain fatal. Exact-state data reuse is confined
  to one fixed steady context.
* The steady inverse approximation uses the complete frozen-scalar planar
  velocity-gradient tangent, all signed rows and cross blocks. Only the
  preconditioner changes. Actual equations and finite-difference actions retain
  the complete current nonlinear stress.

The tangent probe initially serialized nested finite-difference arrays with six
digits, producing a false 7.73e-7 check failure. The test-only max_digits10 repair
was compiled against the immutable source-v11 libraries. Native tangent FD
errors are 1.878e-10/3.752e-10; the independent assembled matrix control is
2.117e-9/5.288e-10 against the unchanged 2e-8 limit. Full signed matrix action
error is 8.66e-16. Both failed and corrected results are retained. Repaired probe
SHA is `3e80f02efbb99187f790e3d0aee95615407eafaba6dfde854be740b89ca01568`;
checker SHA is `0747e5a160d677adc3eb7f3a60163fc6c20108af3442e6b832f2aed69d813baa`.

## Independent verification

The refined evolving temporal MMS study passes its frozen mesh/quadrature/
tolerance separation and owner-gradient budget controls. The two observed orders
are U 2.11355/2.05316, k 2.01494/2.01092, epsilon 2.01863/2.01267 and Pi
1.98973/1.98865. The earlier study with insufficient U mesh separation remains
a failed prerequisite, not a successful result.

The Historical advancing spatial matrix shows convergence, but its final dt
halving changes U and k by 33.3% and 57.6% of their finest spatial errors,
exceeding the fixed 10% subordination gate. It is NOT qualified. The separately
admitted frozen continuous Historical steady MMS matrix removes analytic time
loads and advances no physical time. Its finest pair passes spatial orders:
U 3.83031, k 5.16615, epsilon 3.75542, gradU 3.68280, gradk 4.34908 and
gradepsilon 2.77005. Higher quadrature and tighter tolerances are subordinate.
This verifies spatial convergence separately from the evolving BDF2/EX2 study;
no steady iteration is described as physical spinup.

The seven-step polynomial MMS history foundation passes independent final
equations. Cold continuation and the supported same-dt retry reproduce current,
previous and explicit histories exactly. Cold continuations from retained
recovery writes 4 and 6 differ by at most 2.317e-12 and 7.402e-13 in the fixed
coefficient reference norms. Three actual writes occur at steps 2/4/6; step 2
is actually removed and 4/6 remain independently usable. Scientific outputs
2/4/6 remain retained. Changed-dt and checksum-corrupted restarts reject.
The first inventory checker rejected affine geometry ulps; the frozen v2
checker uses 64-epsilon prescribed-axis classification without altering owners,
fields, original bounds or scientific limits. This is an MMS foundation only:
physical cold restart, terminal handling, monitoring/storage and H27-D stay open.

## Physical reference and timing

The same-inlet 700- and 2800-cell source-v11 steady states independently meet
their current discrete equations. Native covariance and planar principal
checks pass. These facts do not establish an accurate developed reference.

| Calculation | Whole Linux operation s | Native Newton/Krylov | Independent momentum / k / epsilon |
| --- | ---: | --- | --- |
| 700-cell steady from labelled earlier candidate | 184.915 | 8 / 432 | 3.214e-12 / 3.916e-10 / 2.155e-11 |
| 2800-cell steady, every x/r interval bisected | 417.867 | 6 / 186 | 1.590e-12 / 1.584e-10 / 5.138e-14 |
| 700-cell first physical BE step, positive split, dt=.005 | 41.32 | 36.437 s native advancement | Current coupled equations pass |
| Same step, former isotropic split control | 85.39 | 79.846 s native advancement | Same-time independent fields pass <=1e-8 |

The first-step timings are diagnostic. They are not warmed multi-block costs,
an accuracy-qualified dt, or a one-day production forecast.

| Fixed physical or numerical gate | 700 cells | 2800 cells | Limit |
| --- | ---: | ---: | ---: |
| Inner k departure from retained inlet, RMS/Kin | .132331 | .132795 | .05 |
| Inner k departure, sampled MAX/Kin | .459519 | .440730 | .10 |
| Coarse/fine inner k difference, RMS/Kin | .0138243 | same comparison | .001 |
| Coarse/fine inner k difference, MAX/Kin | .129202 | same comparison | .001 |

U/W, outer-epsilon and ground-shear physical envelopes pass, but several
numerical field and owner-gradient flux budgets also fail. Thus these results
show a persistent material inner-k concern, not a mesh-resolved continuum
obstruction or a qualified physical reference. Preserve all regions, stations
and failed budgets. An inactive inner-diffusion proposal and independent column
diagnostic are supplied for review; neither changes active physics nor proves
that the proposed law would remedy the actual 2D departure.

## Remaining work and resources

The completed OPEN account snapshot is 6947.866607507807 aggregate seconds in
248 entries, including previous closed feasibility 703.6786409680499 s and
120 s bootstrap/control. Both failed publication attempts and every failed
numerical operation are charged. The 43200 s allowance and protected 3600 s
reserve have not reset. No operation is unresolved at this snapshot.

H27 model/reference remains NOT_MET because the physical reference and its
development/accuracy gates are unresolved. Coupled numerical foundations now
have independent temporal and steady-spatial evidence, but this checkpoint
does not declare all H27-B conditions closed. Physical baseline/startup/finite
perturbation pilots, selected h/p/q/dt, actual physical restart and rotations,
terminal scientific retention, denser monitoring, full measured cost attribution
and the Task3 launch manifest remain open. The original T=505.02466521278404
and terminal window [447.8483832973273,505.02466521278404] remain binding.
No full Task3 trajectory has been launched.
