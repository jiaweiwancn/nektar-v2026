# Approved-v2 radial reference refinement admission

The latest x-only5600-cell reference still misses inner-k physical and
numerical gates. Its refinement retained100 radial intervals. It is not a
resolved continuum barrier proof. Add the complementary radial refinement,
without changing any physics, frozen norm/scale, station, boundary or inlet.

Start with the actual2800-cell XML axes28x100. Bisect every physical-elevation
interval and retain every x coordinate:28x200=5600 cells. Keep Q3/Q2,
native split quadrature8x12, the approved-v2 tensor/transport/reference,
physical molecular viscosity and all current scalar/root/covariance/principal
guards. No v4 diffusion/cutoff or equilibrium forcing is admitted.

Use source15/ELF `be3dd8fa87236e3c25a5e6b890077090587dd6f85fa83cc309a31a46415b1a05`;
no source change or build is needed. The existing source11 2800-cell field is
only an explicitly source-bound steady initial guess, allowed by the existing
SteadyReference import contract. It is projected to the new signed spaces;
dependent traces, explicit data and fresh zero-time histories are rebuilt.
It is not an imported accepted trajectory or an accuracy-qualified reference.
Source15 must solve its own current equations from that guess. Keep
FullGradientV1 steady preconditioning, nonlinear1e-10, linear1e-12,
Krylov dimension80 and the1800s internal steady deadline. Physical time is0.

Independently assess the full steady momentum/continuity/k/epsilon equations,
weak accounts, boundaries, all original-owner interfaces and every frozen
physical region/station. Compare2800 to radial5600 using the unchanged
numerical RMS/MAX budgets, and compare radial5600 to the prior x5600 at common
physical owner points. The latter meshes are complementary anisotropic
refinements; neither is described as jointly finer. Do not average C0
gradients or omit a failed region. Source11 versus source15 identities remain
explicit in these comparisons; only their unchanged-v2 equations are shared.

Measured prior x5600 native time was708.726278s/RSS6,621,138,944bytes;
its independent equation audit101.018735s and2800/5600 comparison169.910764s.
Forecast the radial native700--1200s, equation audit100--180s and the two
comparisons180--400s each. The differently shaped sparse factors may consume
more memory; the hard8GiB science process-tree watchdog remains binding.

Admit fixtures30s, native2000s (1800s internal), independent equations400s,
2800/radial5600 comparison450s, complementary5600 comparison600s,
compact localization60s and evidence90s: <=3630 aggregate additional
operation seconds. This shares the same43200s account, last verified
OPEN12204.860811439692s; it leaves the3600s protected evidence allowance.
No allowance reset. Keep actual shared-host workers<=72, one science rank,
all six library thread limits1, science RSS8GiB, source/build8GiB,
science storage8GiB and free-space reserve15GiB. Use a new isolated namespace
on the authorized healthy ext4 path. Preserve failed/unknown outcomes before
any retry. No production pilot or full Task3 horizon is admitted.

The first source15 radial run fails construction after14.572948s at the
unchanged2e-13 native H velocity-trace/gradient check; fixtures cost1.604624s.
No steady or physical-time progress occurred. Preserve that namespace and all
receipts. Before changing production code, add a standalone diagnostic linked
to the exact existing source15 static libraries: compare physical-H point
rules with native lower-face rules and the independently written Modified_A
edge polynomial, retaining each owner and the same2e-13 threshold. Admit
configure30s/build100s/two diagnostics80s/assessment10s, adding<=220s to this phase
(<=3850s total). It does not disable a production check, alter a coefficient,
or qualify a reference. A repair, if supported by the witness, needs a new
immutable source/build and prospective qualification before reuse.

The exact-source15 diagnostic confirms radial-H point gradient error
3.258504577274837e-13 versus2e-13; the native lower-face value/gradient errors
are6.439293542825908e-15/4.996003610813208e-14. The failing requested material
point maps to eta=-.9999999999999641. The2800-cell control still passes.
This supports explicit one-sided native reference coordinates for an H trace,
not loosening the check or changing a physical state. Introduce
RuleFactory::HorizontalTrace with exact eta=-1/+1 and endpoint-relative x.
Use it for coupled H rules, the unchanged2e-13 constructor comparison and
the steady dependent-Gamma derivative. Preserve arbitrary Points, field
ownership, model equations, tolerances and all full current acceptance guards.

Materialize source18 from the immutable source15 capsule plus exact Mac
changes, preserving older files. Build the normal solver and trace diagnostic
with four compiler workers. Before the radial retry, require the new native
edge-polynomial comparisons on700/2800/radial5600, source18 continuous
Polynomial/Historical advancing MMS with independent current equations, and
a short700-cell physical source18/source15 comparison at the unchanged1e-8
scales plus a separate current-equation audit. Reuse source15 inputs only as
fresh session definitions; no old accepted native history is imported.

Add materialization30s/build360s, three trace diagnostics120s, two MMS natives
60s/two independent MMS audits70s,12-step physical native120s, physical
current audit150s, cross-source assessment90s and preparation30s. This adds
<=1130s; the complete radial/diagnostic/repair/refinement phase is bounded by
4980s under the same campaign ceiling, including earlier failure charges.
Any failed prerequisite preserves evidence and stops dependent execution.
No gate or physical model changes. The post-storage OPEN12204.860811439692s
charge carries forward; evidence protection and all resource caps remain.
