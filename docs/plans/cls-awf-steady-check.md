# Steady-convergence assessment of the six campaigns

**2026-11-10 update:** the [new exact-checkpoint assessment](cls-awf-steady-check-20261110.md)
confirms three passing windows through T=6 for both H=.5 budget cases,
including the full outer epsilon region. The dated findings below remain
historical; the new report separates candidate discrete steady states from
reference accuracy and correction-supported settling.

**Later continuation:** the H=.5 graded pair has now
[restarted from T=5 toward next target T=6](cls-awf-05-continuation-to-six.md).
This does not change the dated assessment below or establish that the added
windows pass the unchanged persistence criteria.

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

## Why the one-element k residuals are large

These are the earlier capped CLS–EP-AWF cases, with column-average K. Their
selected spatial residual excludes the post-transport positivity correction.
The corrector preserves admissibility of the declared inventory while changing
the numerical step map; it does not make the corrected field solve the original
k equation.

The residual measure is the actual consistent-mass dual norm on free k rows:

\[
\widehat R_k=\frac{\sqrt{r_f^T M_{ff}^{-1}r_f}}
 {K_0(u_*/H_{\rm ref})\sqrt V},\qquad
H_{\rm ref}=2,\quad V=24.
\]

The denominator is .5001936566, fixed independently of current k. Uniform H
has dimensional norm 1.2785437 and normalized norm 2.5560974; Upstream-one has
4.3074891 and 8.6116428. No division by the approximately 5.1e-9 contact value
causes these large numbers. The H-dependent native norms are .6390 and 2.1529;
the common scaling multiplies both by four, preserving the large imbalance.

A zero-PDE same-state decomposition uses the archived **PicardKMass** M_ff,
not the frozen transport matrix. Both mass exports have 378 free rows and the
same offsets as their corresponding current native reports. Component loads
sum exactly to the retained free residual, mass-inverse errors are below
1.2e-15, and reconstructed norms match the independent native values at roundoff.
The [component evidence](cls-accepted-budget-awf-evidence/steady-check-20261109/k-residual-components.json)
contains all source/mass hashes and component norms.

To expose cancellation across modes, let
\(c_j=r_j^TM_{ff}^{-1}r/\|r\|_{M^{-1}}\). These signed projections sum to the
dimensional residual norm; separate component norms cannot simply be added.

| Contribution along the total residual | Uniform H | Upstream one |
|---|---:|---:|
| Inner destruction | +4.5272 | +5.6163 |
| Inner production | -1.5319 | -.7362 |
| k diffusion | -1.6858 | -.5859 |
| Advection and outer source terms combined | -.03093 | +.01333 |
| Total | 1.27854 | 4.30749 |

Thus inner destruction is the main uncancelled contribution, partly opposed
by production and diffusion. Even approximate column-integral balance does
not balance every free polynomial test mode. Both cases use one ordinary cubic
over [0,.5]. The controlled initial state already has ground gradient 32.2468
against analytical 181.9194, and production .89897 against reference 28.61069.
Saturation is inactive at that initial point. With chi=1 the point source is
about -27.7117; with chi=0 it is about -1335.092. The previous 16-versus-32
source-quadrature test changed the initial free source load by only about
2.36e-6 relatively, separating this representation error from initial integration
error. These initial point-source values are distinct from integrated late-state
rates and from the residual norm.

At T=4.665, Uniform H still has chi=1 at x=-3.5,-1,1.6. At T=3.575,
Upstream-one has chi=0 at x=-3.5,-2.625,-1,1.6: EP is inactive even at the
upstream sample. The source there uses the original-AWF branch,
which has the much stronger near-wall destruction. Stress regularization and
SaturatedR2 do not enforce balance between actual production and the selected
destruction law.

At the saved steps the raw k linear residuals are 3.89e-16 and 4.44e-16, but
the raw ground minima are -.0030643 and -.0094988. The corrector raises them
to the prescribed floor. If b is its coefficient correction, the exact frozen
linear residual changes by A_k b. The measured post-correction maximum-row
residuals are .167934 and .460129, agreeing with the added A_k b values;
the distinct storage-only contributions are .158559 and .356565. These are
different metrics from the mass-dual spatial norm. A successful raw linear
solve therefore coexists with a large corrected-equation imbalance.

The evidence establishes current destruction-dominated weak imbalance and
sustained correction support. It does not apportion every long-time change
between temporal splitting and spatial representation, or turn an almost
stationary corrected field into a steady solution of the unmodified equations.

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
this is conditional, not an observed steady result. Following authorization,
both cases have now [restarted from their exact T=5 checkpoints toward T=6](cls-awf-05-continuation-to-six.md),
with unchanged histories, runtimes and criteria. The completed parent campaigns
are preserved; the other four supervisors continue unchanged.

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
