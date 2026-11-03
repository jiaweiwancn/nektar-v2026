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
Linear/R2, chi=0/1, inactive identity, and one separately labelled capped
chi=0 case, plus Linear/R2 exact-local-zero covariance countertests. The zero
fixture analytically prescribes the represented zero ground trace before any
operator evaluation; it is not a post-step scalar repair and takes no outward
central k perturbation. `Audit[inner_zero_covariance]` reports the number of
exact-zero inner evaluations and their minimum/maximum covariance, after
actually evaluating the selected kernel. The cap fixture uses rK=H=2 solely for operator verification;
it is not a changed physical E/F control. Compatible external EPS corner
data are fixed once, while interior lambda*Epre(H) modes remain dependent and
are rebuilt at every perturbed state. Two mixed direction levels capture
actual Lift/constraints/residual actions. Their independent all-row audits
do not replace the stricter pressure/chi action test or the actual implicit
stage result. All active fixtures are explicitly outside physical outer-box
qualification. Genuine curved MMS and the remaining M8 gates stay open.
