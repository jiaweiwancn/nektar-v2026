# M8-R2 operational policy, version 1

This policy supplements the immutable v2 candidate contract committed at
`13a95d66abdbc7b29e95fe6a3b617cfe4e4d77e6`. Its Markdown and JSON hashes
remain `77e20f90cf931f3572d9a9eed3bd9658eb3d2f44ca177b7e05e735a6a81c3bdd`
and `7cbc782674178dfd7366df5e31a9152f0d2fa6c139fd25c4002f66a8f8524e80`.
This new record specifies implementation restrictions and bounded solve
constants before numerical experiments. It changes no physical coefficient,
residual, applicability bound, tolerance or historical result.

The opt-in `ScalarFeasibleR2` prototype supports one MPI rank, a flat 2D
upstream-transition patch and capped destruction. It uses the full dense
coupled solve, without Schur condensation, and retains the 3000-unknown guard.
These restrictions are narrower than the mathematical inner constitutive
domain. Unsupported combinations reject explicitly. Linear inner stress
and Legacy nonlinear policy remain the defaults.

Each accepted Newton iterate receives a fresh current Jacobian. Direct
pressure/gauge actions use the exact affine contribution plus the complete
nonlinear sensor/source/trace action. The stage-mass contribution `alpha*M`
is exact. Other columns use dimensionless, scalar-domain-aware differences,
including admissible one-sided differences at a boundary. Independent
directional-action error must satisfy the unchanged `1e-7` threshold.
The original native convergence tolerance remains `2e-12` in the maximum
reported dimensionless block norm, including algebraic constraints.

The starting step is bounded on the represented fixed scalar inventory.
For decreasing scalar values v with dv<0, define
`maximum=min(1,min(v/(-dv)))`. Apply fraction-to-boundary tau=.995 only if
maximum<1; otherwise start at one. A strict boundary reached at a bound of
exactly one is screened before its residual and may require backtracking.
Exact boundary directions require the stated admissible one-sided treatment;
if no coordinate one-sided action exists, terminate with a classified failure
instead of substituting positive k. This estimate
does not replace trial screening. Rebuild and check all fixed and moving
volume/trace normalization rays, cutoffs, actual lift and dependent epsilon
before each trial residual and before internal acceptance. Native-GLL or
ray negativity is not hidden by positive volume-Gauss values.

Use strict Armijo decrease with coefficient `1e-4` in the maximum
dimensionless true-residual block norm, at most32 backtracking trials,
and classify a step below `2^-40` as collapse. Three accepted improvements
smaller than `1e-4` relative improvement constitute stagnation. These are
bounded-search failures, not accepted physical stages or root-absence
proofs. Save the limiting set/location, margin, step limit, branches and
true block residuals for every attempted/accepted update.

The fixed stage family uses the actual native equation-system entry point
`RANSWallStageTest=MatchedR2` (`Off` is the default). It has t*=.0025,
actual startup dt=.0025, eta=0, the
same represented initial history h0 at every level, and alpha levels
`[3200,1600,800,400]`. Allow at most12 Newton updates per level and no
extra retry or continuation level. If an auxiliary level fails, skip the
remaining auxiliary levels and attempt alpha400 exactly once, using the
last converged scalar-feasible iterate or the initial state if none converged.
A failed alpha400 is never retried. Each law starts from the same admissible
initial field; the preceding auxiliary solution may provide the next
initial guess while h0, time, forcing and constraints remain unchanged.
Neither an auxiliary solve nor a comparison probe calls Advance() or advances
native time/history; `stage_test_levels` is a separate report record. If both
target roots are unresolved, at most one `RANSWallStageTest=FrozenScalarR2`
alpha400 solve may freeze velocity/pressure for scalar localization. This is
a diagnostic rather than a solution of the coupled stage.
The producing execution protocol records the actual arguments, step family,
source/runtime/input hashes and all failures before any measurement.

`ScalarDomain` checks scalar evaluation admissibility throughout the
internal solve. `R2LocalV1` adds the selected covariance, unchanged fixed
outer box, positive sampled K, actual dependent epsilon and lambda>=.5 at
converged-stage acceptance. A comparison probe may preserve a converged
but unqualified state and its reasons; no native accepted history is
advanced by that observation. Complete finite sets remain distinct from a
polynomial-everywhere or invariant-region certificate.

Representation order stays outer6 then, only if coverage requires it,
outer9, with the unchanged inner mesh and base/high independent quadrature.
All Linux numerical work awaits root's verified-operator/action gate and
explicit allocation. Numerical threads are one in all six thread variables;
the reserved analysis worker is CPU49. Parent owns all transfer, execution,
source freezing and conditional short-stage/trajectory allocation.
