# Calibrated budget AWF: independent reference and curved-column checks

The authoritative specification is
[`accepted_awf_appendix.tex`](../literature/accepted_awf_appendix.tex), copied
verbatim from the supplied Mac file. Its SHA256 is
`9998548d33bd2ca856bb6cc85551c82b86d3b499fd54ae8dc9ba26453450fdbf`.
The accompanying revised note has SHA256
`d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f`.

The independent Linux oracle is
[`validate_cls_accepted_awf_kernel.py`](../../solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_accepted_awf_kernel.py).
It solves the original CLS coefficient/Yap roots rather than prescribing the
rounded reference constants. It constructs the calibrated wall coefficients,
retains the original intercept correlation, evaluates closed metric-weighted
production and dissipation primitives, and solves the monotone cutoff equation
independently of the native source quadrature. Its GL64/GL96 factor-three
panels differ from the native panels. A real-dilogarithm primitive is a second
evaluation of the flat rough k profile.

## Independent Linux results

The standalone `--self-check` passed 770 checks on 2026-11-08. The complete
compiled-native comparison then passed **10,827 checks**: seven cases at native
quadrature orders 16 and 32, and three expected applicability rejections.
Neither check advances a native governing equation. The frozen probe SHA256 is
`138d7f53178785662113782edc1062c487044953389c0d2a3e4629c9866c5c0e`.
Its 17 invocations took 1.41 seconds in total, excluding Python orchestration.
The maximum native/reference scaled difference is 3.52e-13 (smooth-reference
k gradient). All native reference weak residuals are below 2.57e-16 absolute.
Native order-16/order-32 cutoff values are identical; weak loads differ by at
most 9.21e-17 scaled by the production budget, and ground k by at most
4.33e-16 scaled by endpoint K.

For h=.1, nu=1.5e-5 and reference u_star=.4:

| Quantity | H=2 | H=.5 |
|---|---:|---:|
| Endpoint K | .510508013059886 | .510508013059886 |
| Reference production/A | 5.909280150953166 | 4.531094243588622 |
| Cutoff | .014755594641345 | .014636433344370 |
| Matching velocity | 5.379706230389028 | 4.128270816644428 |
| Aerodynamic z0 | .005342939598684 | .005300459545886 |
| Reference k at ground | 1.022910439492602 | .994298029088087 |
| Reference k at z=.1 | .544798281032111 | .530743075243479 |

The flat rough integral and closed k evaluations agree within 3e-15 scaled by
endpoint K. The independent weak residuals are below 1.13e-16 absolute for
shifted Legendre modes 0--5 and for every assembled ordinary cubic modified-CG
row on the declared graded wall-normal partitions. Diffusion and source loads
are retained separately; their nonzero higher modes cancel. This is continuum
reference balance, not a claim that a projected native polynomial has zero
residual.

The native forced-momentum check uses the same calibrated viscosity and
cutoff. Its accepted and calibrated-EP source comparison retains identical
matching velocity, forcing and K. The whole-layer sources agree while their
higher modes differ. For H=2, selected shifted-Legendre loads are:

| Degree | Accepted source | Calibrated-EP source |
|---:|---:|---:|
| 0 | .054691928357523 | .054691928357523 |
| 1 | .019999933562998 | .023600353680282 |
| 2 | -.005333142558859 | -.013808577243728 |
| 3 | -.004129702550364 | .009775541144865 |

These are source-only kernel counterfactuals on the same forced analytical
momentum field. They are not separately evolved coupled solutions or a claim
that the old Upstream/Uniform campaigns are matched controls.

The smooth h=0, H=2 control has an admissible cutoff
.000628974425719 but analytical ground k=-.304067694181446. It is reported
without clipping as a stationary-reference positivity qualification failure.
Cutoff applicability and positivity of this auxiliary reference are distinct.
The oracle also rejects a budget below its minimum (h=0, nu=.05, H=2) without
clipping the cutoff. No general evolving-state qualification is inferred from
the sign of an associated auxiliary stationary reference.

## Curved scalar scope

The plus-curvature convention is

    J(r)=(1+kappa1*r)(1+kappa2*r),   x=X+r*n,   B=grad_Gamma(n).

For Q_J=integral_0^r J(E_C-p_A0), the frozen-column relation is
Q_J=J*D_k*k'_column. The physical upward diffusive flux is -Q_J/J.
The independent metric budget is

    B(c)/A = 1+log(H/c)
             +(kappa1+kappa2)*(H-c/2)
             +kappa1*kappa2*(H^2/2-c^2/6).

Its derivative is -A*integral_0^c J/c^2. The oracle checks the derivative,
strict monotonicity, minimum budget, cutoff uniqueness, J=1 limit, and
normal-column endpoint fluxes. It checks principal stretches individually;
positive J alone cannot validate two negative stretches.

At H=2 with the same physical inputs:

| kappa1 | kappa2 | Cutoff | Frozen-column k(0) |
|---:|---:|---:|---:|
| 0 | 0 | .014755594641345 | 1.022910439492602 |
| .1 | 0 | .014792183273770 | 1.026967008714417 |
| -.1 | 0 | .014719042683243 | 1.018769827433485 |
| .12 | -.06 | .014776408356599 | 1.025122739031065 |

These are reference scalar columns. They are not terrain-flow equilibria or
curved PIMPLE evolutions. Positive local stretches do not prove global normal
projection uniqueness. Curved native evolution still needs full mapped
momentum/continuity; curvature and normal-motion traction; tangential scalar
transport; actual matching-surface normals; conservative flux matching;
geometry admissibility; and pressure/positivity coverage consistent with that
geometry. With varying H, endpoint K must sample the actual matching face
using local physical H.

## Reproduction

Heavy results are under
`/media/jay/新加卷/nektar-v2026-runs/cls-accepted-budget-awf-20261108`.
The standalone receipt is `kernel-independent-self/validation.json`, including
source/specification hashes, physical values, weak loads and numerical errors.
The full native receipt and per-case JSON/logs are in
`kernel-native-v2/validation.json`. `kernel-native-v1` is retained: its harness
initially scaled a zero budget residual by 1e-14 and therefore imposed an
inappropriate absolute tolerance near 1e-24. The native curved residual was
1.11e-16. The corrected oracle scales this residual by the physical production
budget and requires 2e-12 relative accuracy. No kernel or acceptance tolerance
was changed for that correction.
Use a fresh output directory when reproducing:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
python3 -B "$R/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_accepted_awf_kernel.py" \
  --self-check --output "$R/reproduce-independent"
```

With a frozen `RANSAcceptedAWFProbe`, replace `--self-check` by `--probe PATH`.
Use its frozen linked-library closure in `LD_LIBRARY_PATH`, as recorded by
`kernel-runtime/identity.json`.
The complete native comparison checks orders 16 and 32, calibrated forced
momentum, all source moments for identical-field accepted/EP counterfactuals,
and explicit invalid-budget/stretch rejection. `--quick` uses native order 16
while retaining independent quadrature refinement. Neither command advances a
production trajectory.
