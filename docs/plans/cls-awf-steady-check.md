# Steady-convergence assessment of the six campaigns

2026-11-09, Asia/Shanghai; native evidence captured 20:48:50–20:49:14.
**None yet meets the complete declared persistence screen.** The two H=.5
budget-AWF runs have successfully completed T=5 and are the closest: their
latest half-unit windows pass the common-panel screen, but only one of the
three inspected windows passes. Four other campaigns remain running toward
their existing T=5 targets. No simulation was restarted or reconfigured.

The [compact evidence](cls-accepted-budget-awf-evidence/steady-check-20261109/native.json)
and [artifact index](cls-accepted-budget-awf-evidence/steady-check-20261109/artifact-index.json)
retain the exact selected steps, time labels, reports, checkpoint and native
history identities. Heavy immutable copies are under
`/media/jay/新加卷/nektar-v2026-runs/cls-awf-steady-check-20261109-final`.

## Method and unchanged criteria

The existing [AcceptedAWFExploratoryWindowsV1 policy](cls-accepted-budget-awf-evidence/covariance/convergence-policy.json)
requires three consecutive windows of duration .5. Normalized U/W/k changes
must be <=1e-4, epsilon/pressure <=1e-3, selected spatial residual blocks
<=1e-3, relative traction change <=5e-4, relative matching mismatch <=1e-5,
and negligible correction supply. The same fixed common normalization is
used for all six cases: U=.4, K=.5105080131, epsilon=.0726358333,
pressure=.16 and volume=24. H-dependent native scales are retained separately.

The analysis selects the latest completed endpoint with three preceding
endpoints spaced by exactly .5, using availability alone. It additionally
reads the latest completed report. The Linear H=.5 windows end at 4.99
because its handoff blocks have that offset; its separate T=5 final report is
also checked. No time interpolation is used.

All spatial residuals are the already computed **alpha=0 selected uncancelled
operator** on the retained fields. They are separate from small linear-solve
residuals or finite-alpha stage defects. The native `HeightComparison` path
calls `Residual(...,alpha=0,anchored=false)` and verifies its input state hash.
Report bytes, coupled state, reference/defect, physical signature and every
native GLM slot/time label match the copied checkpoint and accepted ledger.
Signed JSON zeros are preserved. This check performs zero governing solves,
zero Lift/projection calls and zero physical-time/history advancement.

## Actual results

| Campaign | Latest analyzed T | Latest .5-window normalized k change | Latest normalized spatial k residual | Passing windows / required 3 | Assessment |
|---|---:|---:|---:|---:|---|
| Budget H2 refined Regularized/Linear | 3.78 | 2.51e-4 | 2.47e-3 | 0/3 | Continuing drift |
| Budget H=.5 graded Regularized/Linear | 5 | 8.33e-5 | 6.24e-4 | 1/3 | Approaching; persistence unproven |
| Budget H=.5 graded Regularized/Saturated | 5 | 8.24e-5 | 6.24e-4 | 1/3 | Approaching; persistence unproven |
| Existing Uniform G | 3.655 | 1.35e-4 | 8.90e-4 | 0/3 | Improving; k drift remains above screen |
| Existing Uniform H | 4.665 | 2.32e-5 | 2.556 | 0/3 | Small recent field drift with sustained repair |
| Existing Upstream one element | 3.575 | 1.62e-2 | 8.612 | 0/3 | Substantial bounded transient |

For the H=.5 budget cases, the three consecutive k-window changes are:

| Case | First window | Second window | Third window |
|---|---:|---:|---:|
| Linear: 3.49–3.99 / 3.99–4.49 / 4.49–4.99 | 2.276e-4 | 1.411e-4 | 8.325e-5 |
| Saturated: 3.5–4 / 4–4.5 / 4.5–5 | 2.257e-4 | 1.398e-4 | 8.243e-5 |

The first two windows fail the 1e-4 k-drift criterion. Their endpoint k
residuals also remain above 1e-3 until the final window. At T=5 all four
common-normalized U/W/k/epsilon residuals are approximately
1.14e-4 / 6.23e-5 / 6.24e-4 / 2.98e-4. Matching discrepancies are about
1.18e-7, and no correction is active in the inspected windows. This is
evidence of approach to a selected-equation steady state, not three-window
confirmation. The reported mean wall traction is .152657 versus the .16
neutral reference, still 4.59% low; settling and reference preservation are
distinct assessments.

H2 still fails both k drift and k spatial residual. Uniform G's most recent
spatial residuals pass the 1e-3 screen, but its k drift is about 35% above
the stricter common screen and earlier windows fail residual/matching checks.

Uniform H's small field change cannot establish original-equation convergence:
its k residual is 2.556 and the bound remains active. In its last .5 window,
signed/positive/negative/absolute repair rates are
.5044 / .9613 / -.4569 / 1.4182 integrated-k units per unit span per physical
time. Upstream-one rates are 1.3775 / 2.5818 / -1.2043 / 3.7861.
Both exhibit substantial redistribution as well as net supply. Their absolute
rates normalized by K0*V/(2/.4) are .579 and 1.545. These are rates relative
to a fixed field/time scale, not percentages of total turbulent energy.
The observation supports investigating a bounded fixed point for Uniform H;
three-window original-equation or bounded fixed-point convergence is not
established. No statistical-stationarity claim is made for Upstream-one.

## Scope and next check

Window differences use the archived common physical quadrature and one mean
pressure alignment across both regional fields, preserving their relative
offset. Epsilon changes cover the common z>=2 region; the spatial epsilon
residual covers the entire actual outer transport space. A positive steady
classification would also require full-outer-epsilon field confirmation.
The refined H2 native reporting panels do not split every fine inner interface;
its negative verdict independently follows from the fully assembled k
residual, without relying on that field-norm quadrature detail.

These exploratory thresholds do not certify continuum/spatial accuracy or
accumulated long-time temporal accuracy. The older, broader 1e-3 long-time
screen is not silently substituted for the supplied common screen.

The smallest useful continuation assessment is two more half-unit windows for
the H=.5 graded budget pair, with full-outer-epsilon changes. If those windows
also pass, the existing persistence requirement could be met around T=6;
this is a conditional target, not an observed result or a newly launched run.
Both completed cases remain at their authorized T=5 endpoints. The other four
supervisors continue unchanged.

Re-run the read-only checker on Linux with the same policy and explicit
campaign paths, using an unused output directory:

```bash
S=/media/jay/新加卷/nektar-v2026-runs/cls-awf-audit-cache-temporal-20261109/source/solvers/IncNavierStokesSolver/Tests/RANS
python3 -B "$S/check_cls_awf_steady_windows.py" \
  --policy "$S/../../../../docs/plans/cls-accepted-budget-awf-evidence/covariance/convergence-policy.json" \
  --campaign LABEL=/absolute/campaign/trajectory \
  --output /media/jay/新加卷/nektar-v2026-runs/UNUSED-STEADY-CHECK
```

Repeat `--campaign` for each case. A later live recheck intentionally selects
later available windows. The published values can be independently reproduced
from the indexed immutable report/checkpoint copies and the script's existing
physical-norm functions. Authoritative execution used one Python process with
all six computational-library thread controls set to one; no build was needed.
