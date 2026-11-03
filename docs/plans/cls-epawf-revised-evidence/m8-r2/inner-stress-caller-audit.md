# Selected inner-stress caller audit

This records source inspection and the intended R2 verification scope. It is
not a Linux build or numerical pass receipt. The compiled selector is
`RANSWallInnerStress=Linear|SaturatedR2`; Linear remains the default. The R2
continuation has fixed rho=.8 and no field floor or eigenvalue clipping.

| Actual caller | Selected tensor/action | Deliberately retained |
|---|---|---|
| `RANSWallEvolutionPatch::Residual`, inner volume | `InnerStress(k,nuA,G)` supplies full R and actual `-R:G`; `RANSWallInnerMeanStress` supplies the momentum tensor | Scalar `nu+nuA/sigmaK`, analytical auxiliary profiles, roughness-force definition, and capped source normalization |
| `EnergyBudget`, inner volume | Uses the identical selected production and mean stress | Physical molecular strain dissipation, actual pressure work, storage, convection, prescribed force work and essential reactions |
| `Audit`, finite covariance | Uses selected R, including exact local k=0; outer singular-scalar checks remain separate | The existing declared finite sets and their ownership/reductions; no polynomial-everywhere inference |
| `Audit`, interface | Uses the selected full inner total traction in every Cartesian component | Outer full R and zero resolved molecular viscosity, actual scalar k flux, geometry and normals |
| `Audit`, wall/source budgets | Passes selected actual P to the unchanged source map and exposes the R2 active volume/a/c/tangent statistics | Pre-E/pre-D budgets, lambda, actual E+=lambda Epre, and the actual outer matching trace |
| `RoughnessBudget` | Uses actual represented velocity and the same physical force definition as Residual | It contains no substitute `nuA*S` stress or analytical shear-power approximation to replace |
| `RANSCurvedWallProbe::Interface` | Calls the selected mean-stress helper | Scalar diffusion and the existing outer tensor; this caller change is not a new curved MMS |
| Native directional actions and velocity condensation | Use the current full residual, which re-evaluates the selected inner stress | Local k, sampled K/profile viscosity, chi, cutoffs, normalization and dependent trace columns remain coupled; affine pressure treatment is verified separately |

The auxiliary `Original().EddyViscosity/TotalViscosity` calls in profiles,
reference reconstruction, cutoff construction and scalar diffusivities are
intentional. They must not be replaced indiscriminately by the saturated
strain factor. The original analytical off-reference velocity columns are
not used as an exact active-R2 reconstruction: the actual evolution patch
already solves for its full velocity/pressure/k fields with the new nonlinear
volume tensor.

The older `RANSWallPatch`/`RANSWallInterface` prototypes have no R2 selector and
remain Linear. The selected research implementation and active fixture use
the actual coupled `RANSWallEvolutionPatch`, rather than those older local
prototypes.

Both actual pressure unknowns retain their existing reduced convention.
Writing c=2k/3, the implementation uses
`Ti=nu*S_i+F-pi_i*I` and `To=-R_o+(c-pi_o)*I`. Thus
`p_i/rho_f=pi_i-c`. At the shared continuous-k interface the two c shifts
cancel in the mechanical-pressure jump. Actual production retains
`P_i=F:G-c*div(U)`. There is no outer molecular term and no change of pressure
unknown, gauge or physical pressure-to-sensor map. The outer dimensional
increment remains `GammaScale*nu_t0(x)`, with Gamma prescribed independently
of the evolving state.

The physical/checkpoint signature receives
`:M8-R2:InnerStress:SaturatedR2:v1:rho=0.8` only when SaturatedR2 is selected.
The Linear and inactive branches keep the original `(nu+nuA)*S` arithmetic;
the inactive selected kernel returns identical R and P to Linear. The active
branch uses `nu*S+F`. Exact k=0 is R=F=0 for R2, while Linear keeps its original
possibly nonrealizable tensor. Outward scalar directions at exact zero are
rejected. The k/strain and k/nuA joint-origin limits are directional actions,
not a claimed full Frechet Jacobian.

The focused verification sources cover the compiled kernel versus independent
indexed and frozen-R1 formulas, all nine gradient directions, mixed local-k
and viscosity actions, full production/traction and zero inputs. Eight fixed
flat native fixtures evaluate actual unanchored alpha=0 weak loads: active
Linear/R2, chi=0/1, inactive identity, one legal depleted-but-uncapped chi=0
case, plus Linear/R2 exact-local-zero covariance countertests. The zero
fixture analytically prescribes the represented zero ground trace before any
operator evaluation; it is not a post-step scalar repair and takes no outward
central k perturbation. `Audit[inner_zero_covariance]` reports the number of
exact-zero inner evaluations and their minimum/maximum covariance, after
actually evaluating the selected kernel. The originally proposed cap fixture
used rK=H=2 and was rejected by the existing strict sample-height guard before
operator evaluation. Its failed XML and log are retained. The corrected
cubic-depletion fixture keeps the physical rK=1 and H=2: an independent Linux
calculation gives b=0.1486821878953737 and lambda=1, so it is explicitly
uncapped. Compatible external EPS corner
data are fixed once, while interior lambda*Epre(H) modes remain dependent and
are rebuilt at every perturbed state. Two mixed direction levels capture
actual Lift/constraints/residual actions. Their independent all-row audits
do not replace the stricter pressure/chi action test or the actual implicit
stage result. All active fixtures are explicitly outside physical outer-box
qualification. Genuine curved MMS and the remaining M8 gates stay open.

The already executed pressure-active fixture at the same producing binary
epoch separately exercises the cap at legal rK=1. Analysis-only replay checks
its complete constrained residual against independent unforced alpha=0 weak
loads plus the exact time-mass term, the complete normalization rays and
budgets, and the actual baseline-plus-correction projected lambda*Epre(H)
trace. This is a chi-active capped test, not a cap-active chi=0 test. The
chi=0 depleted control independently checks the unchanged budget mechanism
with selected actual production. Algebraically the pre-E/pre-D
normalization and dependent trace do not contain P; replacing P by the same
selected stress contraction in the k equation retains their cancellation.
No extra native call, CPP change, hard-guard relaxation or model-coefficient
retuning is introduced to activate a cap.

The first B verifier attempt compared diagnostics from distinct states. The
native probe preserves `Audit(state)` on its supplied projected state, while
its physical weak loads, EnergyBudget and Interface evaluate the lifted state.
The corrected postprocessor retains both inputs: supplied-state diagnostics
are checked against that same projected state, and actual lifted-state
interface traction is checked independently at its own order-12 rule and
reported samples. No production stress or Audit result is changed to make
these comparisons agree. Exact retained-report reanalysis and execution of
missing fixed fixtures are separate modes, with original failures and
producing probe/source hashes preserved.

The first legal depleted fixture produced a native report, then exposed an
overflow in the historical independent M1 nested complex exponential at its
large outer q. A nonfinite error subsequently prevented that attempt's JSON
serialization; its wrapper, log and native bytes remain preserved. R2 now
uses a separately indexed outer oracle with the same dimensional polynomial
tensors and actual q. For q>8 it rounds the factor
`1-exp(-.36*exp(.75*q))` to 1. On the positive real axis, the omitted factor
is below 1e-63 and its q derivative below 1.1e-61; hence the C value error is
below 3e-64 and its q-derivative error below 4e-62. These are analytic omission
bounds, separate from ordinary floating-point contraction error. Per-state
tensor omission bounds are recorded for the actual finite oracle calls.
The historical M1/R1/M8 oracle source is unchanged: R2 reuses its quadrature
assembly through copied function globals with explicit R2 stress bindings.
This changes no field, invariant, model parameter or production equation.
The verifier now records nonfinite comparison diagnostics before computing
an error and preserves any remaining invalid numeric serialization value as
an explicit failed marker, so a second JSON error cannot hide the first
scientific failure. A separately pinned additional-retained-run option
reuses the successfully produced sixth native report rather than executing
it again.

The exact-zero R2 fixture subsequently exposed a separate diagnostic
conditioning limit. At its nearest positive source point r is approximately
4.85e-7 and k is approximately 1.20e-13, while the reported maximum a/c is
approximately 9.17e10. The coefficient-derived maximum differs by about
2.28e-5 relative, so that comparison remains explicitly **unresolved at the
original 3e-10 tolerance**. The replay retains the strict full weak, energy,
traction, covariance and exact-zero checks and the five well-conditioned
saturation statistics. It separately checks the exported worst-point
location/primitive inputs against the independent coefficients, recomputes
the selected law using those checked actual inputs at the original tolerance,
and records a directed-rounding, a-posteriori same-point quotient sensitivity
bound. That conditional arithmetic bound does not certify all native basis
evaluation errors or the global maximum. Nonzero fixtures retain the full
original statistic comparison. The report flags that all saturation
statistics have **not** been independently certified. Neither the physical
model nor its fields or evaluation domain changes for this diagnostic issue.
The separate Linux 80-digit polynomial diagnosis finds a 3.0237e-9 relative
native k discrepancy at that point and a modal cancellation condition of
4.12558e6; the older direct binary64 independent basis evaluation is less
accurate there. The native exported-input ratio recomputation agrees within
2e-15 relative. These are finite reference calculations, not rigorous native
basis-error or continuum certificates; their exact inputs and results are
retained in `inner-zero-ratio-diagnosis.json`.
