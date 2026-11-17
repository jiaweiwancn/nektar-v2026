# Refined MMS study v2, prospective recipe

The v1 temporal orders passed but the coarse 8-cell base changed U by 13.27%
of the finest temporal error on a 32-cell comparison, beyond the fixed 10%
separation requirement. Keep all thresholds and the continuous equations.

Use source-v10's registered solver and PolynomialScalarsV1 at terminal time
0.4 with dt=0.08,0.04,0.02 on nx=nr=4 (32 cells). Compare the dt=0.02 result
with nx=nr=8 (128 cells), 12x18 instead of 8x12 integration, and tighter
linear/nonlinear tolerances 3e-13/3e-12 instead of 1e-12/3e-11. Native fields
remain Q3/Q2. The existing independent study's fixed RMS order >=1.8,
10% separation, 20x24/24x32 rule stability and improving physical owner-gradient
budgets apply unchanged. Current residuals use fixed physical reference norms.

For genuine scalar spatial error, use the Historical MMS (quartic k and
rational outer epsilon) at T=0.02, dt=0.0025, with nx=nr=2,4,8. Add the fine
mesh at dt=0.00125. Report U/k/epsilon/Pi/gradient rates and require temporal
contamination <=10% of total fine-grid field error before interpreting spatial
rates. Analytic smooth scalar targets support asymptotic L2 order four and
H1 order three; field coupling and small-error saturation are reported
explicitly, without asserting a rate from a representable target.

The independent helper samples original owners and split material/cutoff
panels; it does not average C0 gradients. The spatial-only report can establish
endpoint/current/budget evidence but cannot set the temporal study gate.

Admit at most 1200 aggregate Linux operation seconds for this study tranche,
including generation, all native variants and analysis. No dependent launch
until the source-v10 build and focused current-action checks pass. Use one
science rank and all numerical-library threads one under the shared-host
watchdog. Retain failed results and do not change limits after seeing them.
