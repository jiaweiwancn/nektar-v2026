# R2 scalar boundary and derivative evidence, addendum 3

This clarifies the immutable contract v2 and policy v1 without changing
coefficients, equations, tolerances, the fixed box or the running stage
allocation. M8 remains open.

The physical restricted contract requires lambda>=.5 and actual outer
epsilon/epsilon0 in [.98,1.02], including the projected dependent trace.
The unchanged capped-source relation is epsilon(H)=lambda*Epre(H), also
at chi=0. As lambda approaches zero this trace approaches zero, while the
original outer high-Re equations require positive epsilon. This prototype
adopts no high-Re continuation through epsilon=0 and adds no epsilon floor.
Internal scalar screening rebuilds normalization and projected targets and
rejects a nonpositive target or actual outer epsilon before a residual or
internal acceptance. An accepted physical stage must additionally satisfy
the stricter fixed applicability guards. A continuation outside them is an
unresolved model choice, not an automatic consequence of the inner law.

The mathematical inner law permits exact local k=0. Compiled admissible
one-sided and joint-origin actions and the zero-covariance definition are
verified separately. The solver uses admissible coordinate differences;
if neither side of any required coordinate action is evaluable after its
fixed shrink ladder, it terminates as a classified derivative/domain
failure. This is an implementation restriction, not a proof that no
admissible combined direction or stage root exists. Fields are never
replaced to obtain an evaluable derivative.

The native weak continuity residual is the linear form
`B*U + gauge*pressure_mean`, evaluated with adequate polynomial quadrature.
The exact pressure-gradient/roughness/gauge directional assembly and the
complete nonlinear pressure/sensor/source/projected-trace action are
separately verified. Other Jacobian columns in ScalarFeasibleR2 use scaled
residual differences; this includes their continuity rows. They are not
described as analytically assembled columns. The independent stage analyzer
assembles the polynomial B and pressure means separately, checks actual
continuity/gauge residual values, and compares the stored current matrix's
columns/actions after external lifting against that affine operator using
the unchanged absolute 1e-7 diagnostic gate. Agreement of a cached linear
solve or a small pressure residual is not derivative evidence.

Internal Newton iterates satisfy only the declared finite scalar evaluation
domain. Their covariance and outer-box membership are separate observations.
Auxiliary mass-parameter roots are not physical time samples. Only a true
original-target residual below the unchanged 2e-12 norm tolerance, followed
by all required scalar, constitutive and accounting checks, can qualify
the target. No residual minimizer with nonzero residual is accepted.
