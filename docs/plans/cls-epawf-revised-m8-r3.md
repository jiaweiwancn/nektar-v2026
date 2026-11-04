# M8-R3: nonpressure actions and a durable exact target

Status: implementation and fixed-state verification in progress. M8 remains
**open**. No new target solve is qualified by this working report.

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
`inputs/input-ledger.json` pins the producing R2 source/binary, original XML,
complete native records, state/history/residual binaries and extracted packets.
Original history h0 has float64-LE SHA256
`f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca`.
The exact saved Linear and SaturatedR2 alpha3200 terminal states remain separate
from the frozen-scalar path. Saved path iterations2 and3 have native minimum k
`1.2762700326467293e-5` and `6.38135016469755e-8`; their signs and the full
current inventory must be checked before directional verification. The terminal
value near8e-15 remains numerically unresolved and is not a positivity proof.

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

## Execution gates

The ordered work is component action implementation and independent checks;
saved-state residual-resolution audit and scalar-contact accounting; durable
iterate/timeout verification; one complete repaired Jacobian/factorization
benchmark; then a budgeted exact SaturatedR2 target from the revalidated saved
predictor. The target is not launched before its independent action gates.
Compiler parallelism is16. Native numerical runs use one rank, and
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, BLIS_NUM_THREADS,
VECLIB_MAXIMUM_THREADS and NUMEXPR_NUM_THREADS are all explicitly one.
Parent-owned allocations include concurrent analyses and remain below72.

The revised action must include moving quadrature panels, weights, evaluation
locations, local polynomial rows, cutoff budget, normalization, sampled K,
transition and actual projected trace derivatives. A closure value evaluator
with a selected invariant branch is not a derivative API. Nonsmooth choices and
appropriate one-sided comparisons are recorded. Small relative errors do not
replace the absolute gate.

The resolution audit keeps coefficients fixed while comparing arithmetic.
The existing code already forms state minus history before applying mass, and
its matrix-vector products already accumulate in long double. These are not
invented defects. Any stable equivalent evaluation requires measured evidence;
the2e-12 threshold is unchanged.

The contact diagnostic uses conservative scalar coefficient loads and the
consistent free/boundary mass partition. A physical rate is a basis-row
evaluation of the mass solve, with dependent-trace coupling stated explicitly.
It cannot be inferred from a load coefficient's sign. A constructed contact
state is a separate diagnostic. Any discrete inward-invariance conclusion
must retain its finite-state/operator scope and does not prove implicit-root
absence or continuum instability.

## Open decisions

The final report will distinguish arithmetic resolution, action accuracy,
valid non-descent, scalar contact, physical range and resource limits using
complete retained states. Accounting identities and balance at a root remain
separate. Only a target passing the original residual, independent reconstruction,
scalar/covariance/applicability/trace and accounting gates permits one actual
native startup step with proper history bookkeeping.

Curved MMS, active-cap curved/variable-H coupling, fully discrete positivity,
long-time/refinement, scalable feasibility and developed physical floor/hill
validation remain open. The broad constitutive counterexamples are unchanged.
