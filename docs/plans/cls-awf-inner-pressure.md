# H=.5 cutoff-v1 AWF: controlled inner-pressure P2 test

The preceding package qualified a P1-inner/P1-outer selected discrete steady
candidate through three complete half-unit windows, ending at T=4.02. Its
remaining traction deficit is .2558%, and its inner strong-divergence L2 is
.00940571. This experiment changes **only inner pressure from P1 to P2**.
It tests pressure-space sensitivity; it does not assume pressure degree is
the sole source of reference error.

The [predeclared policy](cls-awf-inner-pressure-evidence/policy.json) fixes
the matched problem, pilot endpoints, norms, numerical comparison floors,
resource limits and bounded continuation decision before execution. Heavy
artifacts are under
`/media/jay/新加卷/nektar-v2026-runs/cls-awf-inner-pressure-p2-20261111`.

## Implemented control

`--pressure-space inner-p2-v1` is restricted to `Accepted-AWF-05G` with
`cutoff-v1`. The existing default and the separate H2 outer-P2 control retain
their previous generated XML. The new control has native field sizes
`[784,784,784,400,187,54,1]`: 2,994 coupled unknowns, within the unchanged
3,000-unknown pilot limit. Inner pressure has quadratic tensor expansions;
outer pressure remains linear. Velocity, k and outer epsilon remain cubic.
Both regional pressures and the single common gauge remain separate.

The physical mesh, H=.5, endpoint K=k(H), z_K=z_m=H, h=.1, nu=1.5e-5,
u_star=.4, calibrated budget cutoff, actual P_h-E_C source and zero forcing
are unchanged. The selected companion is Regularized/Linear. Native IMEX1
uses dt=.000625, unit relaxation, one U/p/epsilon/k solve and one logical
final bound call. ExactStepV1, NativeElementSideV2, AcceptedBudgetV1,
RetainedPool, floor 1e-8 and controls 2000/20000/2000/900 remain selected.

The native signature already binds every field's dimension, regional
geometry, basis modes and quadrature. Thus the changed pressure space has
a strict restart identity without adding an unused physical parameter.
No checkpoint migration, pressure coefficient injection, trace restoration
or extra transport solve is used.

## Focused Linux qualification

The generated default H=.5 and existing H2 outer-P2 XML match the archived
3714cf3 generator exactly after normalizing transient output paths. The
new case has identical geometry, physical parameters, source quadrature,
controls and non-pressure expansions. Unsupported combinations are rejected.

The exact-control-closure native checks in `verification-matched` pass:

- Actual 2,994-unknown native initialization and inactive startup bound.
- Bitwise identical U/k/epsilon seed coefficients, differential history
  slots and time labels against P1. All initial physical differences,
  including commonly aligned pressure, are exactly zero.
- Two uninterrupted steps equal one step plus restart bitwise in coupled
  state, native histories/time labels and complete checkpoint payload.
- A deliberate post-solve hard rejection retains the actual candidate and
  restores the accepted coupled state, physical time and native histories.
- P1-to-P2 and P2-to-P1 ordinary restarts both reject signature mismatch.

Initialization uses the established native zero-perturbation projection
and two initialization Lift iterations. These are separate from physical
steps. Verification performs five solves of each governing equation in
independent test paths; production is a separate fresh h0 trajectory.

## Runtime and comparison provenance

The solver is the unchanged Linux Release image used by the P1 pilot,
SHA256 `177ce5e54e7de514fa019122e46d3a5d709455cd091bd3023fbfa745c48c9bfc`.
No C++ change or rebuild is required for this native regional-expansion
control. The strict read-only probe is
`67e9d5093d8bdf164ca7fdf7dc97f5a5cb70a7012dce1096be399e45363f2a43`.
Source/configuration and focused tests are `b6c65b4`; matched readback
orchestration is `711eff8`. Both executable and linked-library closures
are immutable in the new run root. One rank and all six computational
thread controls one are explicit.

An initial healthy exploratory segment in `pilot` was stopped at step64,
T=.04, after checking dependency identities. It inherited the preceding
long-run library closure, which has a different `libSpatialDomains` image
from the original short P1 pilot. It is preserved as a separate dataset.
The definitive `pilot-matched` experiment starts freshly from h0 with
the **exact P1 pilot dependency bytes**; no preliminary checkpoint is used.
The native qualification was repeated under that exact closure before
launch. This avoids treating dependency provenance as an uncontrolled
pressure-space effect.

The analysis reuses hash-verified P1 readbacks and the existing strict
checkpoint observer. New P2 readbacks do zero governing solves, Lift calls,
physical commits or history advances. Physical integration panels are common;
epsilon covers the complete .5<=z<=4 region with the fixed E0*sqrt(12)
normalization. One pressure-mean shift across both regions preserves their
relative offset. Reference error, drift from actual h0, pressure-space
difference, selected residual and correction supply are separate metrics.

Matched native observations are requested at T=.02, .05, .1 and .5.
An observation is used only after its exact native checkpoint and completed
block report are published. The long P1 control remains immutable; it is
not rerun. The completed H2 dataset and other experiments remain untouched.
