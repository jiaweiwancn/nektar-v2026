# M8-R3: nonpressure actions and a durable exact target

The exact SaturatedR2 stage did not converge. The new actions pass the
specified fixed-state checks, and the target finished with a complete durable
terminal record. Scalar step limits led to stagnation at residual `16.0695`,
well above `2e-12`. The terminal outer gradient also exits the fixed box.
No physical step was accepted. M8 remains **open**.

The branch started exactly at reviewed baseline
`e4c9d44758f08109ba5ea8c2b909790a02b3acc8`, without intervening commits or
tracked changes. Unrelated user literature and HIT material is preserved.
Source/Git work uses the Mac checkout. Numerical work uses the designated
Linux server through Tailscale. The canonical Linux checkout remains on
`codex/hit-inflow-generator` at
`f4b70baec10016dd098005bc5a5507b87e4a57bf`. Historical source trees, binaries
and M0–M8-R2 evidence are preserved.

## Frozen problem

The physical law is SaturatedR2, rho=.8, with the Linear inner-law control;
beta=.4, GammaScale=1, GammaUniform=0 and original outer C(q). Inner molecular
viscosity remains present and outer molecular viscosity remains zero. Scalar
diffusivities, capped source definitions (theta=.2, delta=.1), H, roughness,
sampler/sensor positions, transition and the dependent epsilon constraint
retain their R2 definitions. The target uses the identical outer6 geometry,
p3/p1 expansions, GLL8 and source quadrature8/16, with 2165 unknowns and the
3000-unknown guard. The existing finite inventory, constitutive box, absolute
action threshold1e-7 and coupled residual threshold2e-12 remain fixed.

The sole expensive coupled target is alpha400, physical time .0025, original
h0, eta0 and zero prescribed body/scalar forcing. The failed SaturatedR2
alpha3200 terminal candidate is an explicitly authorized initial guess. It
retains its historical nonconvergence and supplies no physical time sample,
history advancement or changed target equation. No new auxiliary sequence is
required. A Linear target retry requires a specific unresolved causal question.

`ScalarFeasibleR3` and its paired `VerifiedR3Action` selector identify numerical
research changes. Legacy and ScalarFeasibleR2 remain reproducible. The R3
action name is an implementation identifier; independent error measurements
determine qualification. Curved/3D extensions are outside this focused action
implementation. Actual stress, production, traction, conservation identities
and scalar evaluation-domain definitions retain their physical meanings.

## Fixed-state inputs and reproduced baseline

Heavy inputs and results live under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r3`.
`inputs-bitwise-003/input-ledger.json` pins the producing R2 source/binary,
original XML, complete native records and lossless state/history/residual
payloads. The original native initialization h0 has float64-LE SHA256
`b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71`.
The current production zero-step preflight reproduces both the initialized
state and the frozen reference bit for bit. The latter has SHA256
`4ad26ed0167bb1765bd0be2c150928280b69dbd82ace8a27ea5f9c8c7f4d0409`.
Evaluating the retained zero-amplitude initial functions changes 46 zero sign
bits between these vectors, without changing any coefficient value. The R2
auxiliary histories and the R3 target use the actual initialized h0. Native
JSON `-0` tokens are parsed losslessly. The earlier extraction and unsuccessful
setup receipts are retained.
The exact saved Linear and SaturatedR2 alpha3200 terminal states remain separate
from the frozen-scalar path. Saved path iterations2 and3 have native minimum k
`1.2762700326467293e-5` and `6.38135016469755e-8`. Iteration 2 supplies the
verified positive boundary-path action base. The terminal value near8e-15
remains numerically unresolved and is not a positivity proof.

The retained initial outer-k component discrepancy was reproduced on Linux
without a native call or trajectory, using the immutable R2 analytic diagnosis
source `a268962a5dcc9adfe59fb412981b57242d671654` and the actual retained matrix.
All63 diagnostic controls complete; the original action gate remains failed:

| Output block | Stored Linear J minus independent analytic action |
| --- | ---: |
| streamwise momentum | 7.15996593112029e-7 |
| normal momentum | 5.607262291281462e-9 |
| k | 6.410811162823151e-8 |
| outer epsilon | 4.859374905274791e-8 |
| pressure/gauge/constraints | 0 |

This reproduction is in `reproduce-initial-k-001`. It establishes the retained
Linear witness again; it supplies no SaturatedR2 action result or nonlinear
causality claim. Lambda remains one in the relevant R2 completed records.

The first source epoch, `b19299eafe572401032b9922b13fbcbcc2847283`, stopped at
configuration because prefix-only package discovery did not locate the installed
lib64 CMake package. Its receipt/logs remain in `epoch-001`; no C++ or numerical
run completed there. The next helper supplies the documented Nektar++_DIR.
Static review also aligned invariant, limiter and cap selection with the native
pre-cast long-double branches before compiling the next source epoch.
Epoch002 configured successfully, then exposed a C++ name ambiguity between
the prototype dense matrix and Nektar's matrix type in the new resolution
diagnostic. Qualifying that type changes no arithmetic or equation. The failed
build receipt and complete compiler log are preserved.

## Producing build and action verification

The producing epoch is `epoch-005` at
`ce408c874ff12aaf6133cd8660e8c5b171f17a11`. Its build receipt and runtime
fingerprint pin the images, toolchain and libraries. Epochs 003 and 004
retained probe compilation errors involving private mass access and a
temporary filename vector; the probe now uses the flat field's consistent
mass and a named filename vector. The first fixed-state attempts also exposed
missing explicit numerical selectors and a rank-local session discrepancy.
The subsequent probes use the rank-local session settings and both explicit
R3 selectors. A focused probe-only guard now rejects a missing selector before
assembly; the valid initialized-h0 matrix is bitwise unchanged. No physical
parameter or production solver code was changed by this last probe correction.
Compiler parallelism is16. Native numerical runs use one rank, and
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, BLIS_NUM_THREADS,
VECLIB_MAXIMUM_THREADS and NUMEXPR_NUM_THREADS are all explicitly one.
Independent analyses use separately allocated single CPUs. The three focused
legacy time/evolution/stage tests passed, the current affine and active
pressure-action regressions passed, and the four compiled durability checks
and nine predictor metadata checks passed.

The component assembly differentiates moving quadrature panels, weights,
evaluation locations, local polynomial rows and gradients, the discrete
cutoff budget, normalization, sampled K, transition and actual projected
trace. The existing exact time-mass and pressure/gauge actions are retained.
The independent value oracle uses indexed constitutive tensors and a separate
complex evaluation of the discrete operator, with both `1e-20` and `1e-30`
steps and checks of the cutoff slope. It checks the selected strain/rotation
and Yap branches, rather than differencing across a switch. The R3
nonpressure cap selection is uncapped at the exact knot; the existing
pressure column selection at that knot is directional. A common mixed
generalized action at that exact knot is outside the verified states. The
moving-point action assumes an interior basis evaluation in the selected
element; it is not qualified for a point displaced outside its selected
element and clamped by the native element-search tolerance.

| Fixed state | Maximum absolute physical action error | Result vs `1e-7` |
| --- | ---: | --- |
| Linear h0, including retained direction | `1.891e-12` | Pass |
| SaturatedR2 h0, including retained direction | `1.898e-12` | Pass |
| SaturatedR2 saved alpha3200 | `1.294e-11` | Pass |
| SaturatedR2 positive boundary-path iteration 2 | `1.782e-12` | Pass |
| Limiter and active transition fixture | `1.906e-11` | Pass |
| Active cap and projected trace fixture | `1.169e-11` | Pass |

The repaired action of the retained outer-k direction is `1.503e-12` in both
laws at h0. This result applies to the tested states and directions. It does
not revise the historic failed Linear or SaturatedR2 gates or establish the
cause of any earlier nonlinear stagnation. The first constructed positive
fixture activates the limiter and transition but has minimum b=`.11869` and
lambda=1. A second, targeted cap fixture reduces only the declared fixture's
low polynomial background from `.01 K0` to `.001 K0`; its measured minimum
b=`.07630`, lambda=`.76297`, and its maximum action error is `1.169e-11`.
The sampled-K derivative of the actual projected trace has a nonzero
coefficient action of `.58481`. The admissibility interpolation gate is on its
unit plateau in the tested states. These are operator inputs outside the fixed
physical box, not field floors, clipping, or physical retuning.

## Fixed-state residual resolution

Each original native residual is reproduced bit for bit. These comparisons
keep coefficients, histories, time, alpha and quadrature fixed. The maxima
below are all in the existing scaled physical norms and are dominated by
streamwise momentum. The root criterion remains `2e-12`.

| Evaluation | Linear alpha3200 | SaturatedR2 alpha3200 |
| --- | ---: | ---: |
| Actual native | `5.33424e-12` | `4.59464e-12` |
| Independent term assembly | `6.08434e-12` | `5.33409e-12` |
| Native term split, compensated sums | `5.56173e-12` | `4.79480e-12` |
| Native rounded primitives, extended products/sums | `5.55832e-12` | `4.79064e-12` |
| Native minus extended, norm of the difference | `1.88998e-12` | `1.68776e-12` |
| Native minus independent, norm of the difference | `3.13178e-12` | `2.93534e-12` |
| Extra outer isotropic cancellation, norm of the difference | `1.84e-15` | `1.99e-15` |

The exact-polynomial mass certificate encloses the norm of each **fixed**
load and agrees with the corresponding native mass-dual norm at the reported
precision. The existing code subtracts state minus history before applying
mass and accumulates its dot products in long double. Using the same delta,
the largest observed stored-mass action discrepancy is `9.99e-16` in
coefficient units. The projected constraint differences are below `7.1e-15`
in the scaled norms. For the frozen combined spatial increments, the
summation-only bound is `4.66e-12` in maximum coefficient units, with up to
3,794 additions. It does not enclose errors in the basis, closure, other
pointwise inputs, mass entries, or projection. Spatial momentum assembly is the
material resolution concern. These fixed-state comparisons do not establish
a `5e-12` floor or demonstrate an improved attainable root from an arithmetic
change. The exact target therefore uses the original production residual.

## Scalar contact

The contact diagnostic uses conservative scalar coefficient loads and the
consistent free/boundary mass partition. A physical rate is a basis-row
evaluation of the mass solve, with dependent-trace coupling stated explicitly.
The single constructed polynomial is nonnegative by its exact tensor
coefficients, keeps U and stored reduced Pi fixed, and reconstructs the actual
dependent trace (relative discrepancy `5.31e-15`). It passes the declared
finite scalar, covariance and outer-box checks. The verified trace branch
has lambda=1 and zero epsilon feedback in its derivative. The independently
projected sampled-K trace rate includes its induced epsilon boundary mass.
The free consistent-mass equations balance to `8.89e-16` for k and
`4.78e-18` for epsilon in coefficient units.

| Row in the constructed state | k | Physical k rate | Largest contributions |
| --- | ---: | ---: | --- |
| Constructed zero at `(1,0,0)` | `0` | `+34.9138` | Diffusion `+28.4932`, production `+6.67569`, destruction `-.254739` |
| Original diagnostic row at `(1.6046496089512394,0,0)` | `.186642` | `-633.545` | Destruction `-677.698`, production `+26.0784`, diffusion `+18.0769` |

The assembled component loads recover the unchanged scalar load, with net
destruction counted once. The constructed zero row has an inward rate. The
negative rate at the other row is at positive k; it proves no discrete
inward-invariance defect. This diagnostic neither explains the full O(10)
frozen-path residual nor proves absence of an admissible implicit root.

## Complete assembly and exact target

The first benchmark invocation stopped in predictor binding, before any
Jacobian or nonlinear update. The debugger localized the failure to the
signature/history check. The corrected manifest selects the actual initialized
state rather than the equal-valued frozen reference. The one complete benchmark
then assembled and factored the actual alpha400 matrix from the unchanged
failed auxiliary predictor, without a nonlinear update. Its actions agree
with the independently qualified saved-predictor actions after the exact
mass-coefficient change to `8.97e-12` in the existing norms. The current
predictor has minimum covariance `.08200` inner and `.17692` outer and
outer gradient departure `.00053225`; it passes the fixed finite box. It
remains a numerical initial guess, not a converged auxiliary or time sample.

| Operation | Complete benchmark calls / seconds | Exact target calls / seconds |
| --- | ---: | ---: |
| Scalar domain | 1 / `.250` | 13 / `3.134` |
| Constraints | 1 / `.062` | 7 / `.432` |
| Native residual | 1 / `.150` | 7 / `1.602` |
| Complete action assembly | 1 / `9.193` | 6 / `56.758` |
| Factorization | 1 / `.161` | 6 / `1.042` |
| Line search | 0 / `0` | 6 / `5.860` |

These are measured solver phase times; line search includes its domain and
residual calls. The 600-second target cap uses the measured assembly,
factorization and trial costs, up to 12 updates and 32 trials per update,
with extra allowance for setup, durable output and terminal audit. The one
target invocation finished in `116.79 s` of native wall time. It did not
repeat any auxiliary stage and did not run a Linear target.

| Completed iterate | Scaled residual max | Minimum native inner k | Accepted fraction |
| ---: | ---: | ---: | ---: |
| 0, saved predictor | `47.38339197` | `.3446858` | — |
| 1 | `17.32380644` | `1.72343e-3` | `.582500` |
| 2 | `16.07587660` | `8.61715e-6` | `.0692303` |
| 3 | `16.06953786` | `4.30857e-8` | `3.94294e-4` |
| 4 | `16.06950616` | `2.15429e-10` | `1.97240e-6` |
| 5 | `16.06950600` | `1.07714e-12` | `9.86201e-9` |
| 6, terminal | `16.06950600` | `5.35357e-15` | `4.93098e-11` |

All six line-search trials were feasible accepted descents. All six
factorizations and solves completed. At each update the limiting fixed
location was the native reconstruction row `(1.6046496089512394,0,0)`.
The fraction-to-boundary steps shrank while the residual stayed nonzero,
ending at `feasible_residual_stagnation`. This identifies the barrier on the
observed Newton path; it does not establish nonexistence of an admissible
implicit root. The last small positive native k value does not certify its
sign independently. No resource limit or residual exception stopped the run.

The native terminal maximum comes from k (`16.06950600`); streamwise
momentum is `4.64580`, normal momentum `.0141459`, epsilon `.0681894`, and
the scaled constraint error `8.52433e-5`. Weak continuity and gauge are
small separately. The native covariance minima are about `1.04e-15` and
`.17711`; the near-zero inner value requires the same numerical caution.
The terminal outer gradient entry departure is `.00196589 > .001`, while
the k and epsilon ratios remain inside their intervals. Lambda stayed one
and no cap was active; cap-induced trace attenuation is not the barrier.
The outer-box exit is a separate physical failure, not the globalization
stopping criterion. No native startup step was authorized by this outcome.

The independent reconstruction gives `16.06950600086231` against native
`16.06950600086237`. The largest block norm of their difference is `3.11e-12`.
The common finite inventory has one sign-unresolved value at the limiting
row: independent k is `5.42e-15`, within its `9.16e-15` arithmetic scale.
The audit confirms the outer gradient departure `.00196588655`, with the
scalar ratios inside their specified intervals. It also verifies the native
report and durable residual bytes and the unadvanced integrator history slots.

| Terminal observation | Independent value |
| --- | ---: |
| Weak continuity, inner / outer, scaled | `2.02e-15` / `2.19e-15` |
| Gauge, scaled | `4.34e-17` |
| Strong divergence L2, inner / outer | `.00411030` / `2.33256e-5` |
| Normal / tangential total-traction jump RMS | `9.73167e-5` / `2.06061e-4` |
| k-flux jump RMS | `1.74611e-5` |
| E minus net D minus redistribution | `1.74e-16` |
| Integrated redistribution | `-1.74e-16` |
| Actual-state energy weak work minus direct sum | `-8.53e-14` |
| Free energy work at this non-root | `.802665` |

Pressure uses the existing gauge and reduced-pressure representation. The
traction checks retain the full normal and tangential stress; no mechanical
pressure continuity is imposed. Source and energy identities on the actual
state are separate from balance at a root. The native energy **report** fails
its same-state crosscheck by `3.53e-8` at this unconverged iterate. It computes
direct terms with the dependent epsilon trace substituted, whereas the weak
residual uses the retained matching coefficients. A read-only substitution of
the 15 differing trace coefficients reproduces the native direct stress/source
term to `1.0e-13` and accounts for the discrepancy. The terminal trace constraint
is not satisfied at the root tolerance, and this diagnostic discrepancy has
not been promoted to a physical energy defect.

The [decision ledger](cls-epawf-revised-evidence/m8-r3/decision-ledger.json)
binds the exact source, executable, XML, predictor, h0, action errors by block
and direction, and all small receipts. The complete target is under
`epoch-005/target-run-001` and `epoch-005/target-case-002`. The atomic
`case.r3-progress.jsonl.last.json` and its immutable terminal snapshot bind
the full state, history, residual and constraints. Each completed iterate
was published before starting its next matrix. The initial benchmark stop,
failed setup attempts and historical evidence are preserved. The matrices,
field vectors, and detailed logs remain on Linux.

## Next decision

Before another expensive stage solve, use near-contact iterate 4 and its
retained current matrix, after checking the limiting sign, to test for a coupled residual
descent direction tangent to the limiting native k row. Keep the same
residual and trace constraints. This isolates whether a different feasible
direction can continue past the observed Newton step restriction. A nonzero
constrained residual minimum would remain a failed stage. The present
contact test supplies no verified positivity defect that would justify a
discretization change.

Curved MMS, active-cap curved/variable-H coupling, fully discrete positivity,
long-time/refinement, scalable feasibility and developed physical floor/hill
validation remain open. The broad constitutive counterexamples are unchanged.
