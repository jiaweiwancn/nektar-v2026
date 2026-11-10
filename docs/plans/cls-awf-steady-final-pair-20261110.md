# Final status of Uniform G and Upstream-one

Observed **2026-11-10 17:00:23 Asia/Shanghai** on the authoritative Linux
server. Both formerly live campaigns completed their authorized T=5 targets
normally. Their supervisors, drivers and native solvers have exited, with no
new numerical rejection or infrastructure interruption. Neither was restarted,
retargeted or reconfigured during this inspection.

| Campaign | Accepted step / T | Terminal progress recorded | Steady assessment |
|---|---|---|---|
| Existing Uniform G, H=.5 graded | 8000 / 5 | 16:14:10, Nov 10 | Passes 3/3 declared windows, including full outer epsilon |
| Existing Upstream-one, H=.5 | 8000 / 5 | 14:49:23, Nov 10 | Passes 0/3; substantial correction-supported transient |

The original six campaigns have now completed their respective authorized
targets. Three separate native processes from the newer
`a93-a95-e94e95-balanced-20261110` experiment were observed; those are outside
this pair's assessment and were left untouched.

The [compact evidence](cls-accepted-budget-awf-evidence/steady-final-pair-20261110/native.json)
binds actual terminal progress, exited process identities, checkpoint/state/all
history hashes, runtime identity, immutable native reports and diagnostic
outputs. The [artifact index](cls-accepted-budget-awf-evidence/steady-final-pair-20261110/artifact-index.json)
locates heavy files under
`/media/jay/新加卷/nektar-v2026-runs/cls-awf-steady-final-pair-20261110`.

## Unchanged convergence screen

The same `AcceptedAWFExploratoryWindowsV1` policy requires three consecutive
.5-unit windows: normalized U/W/k changes <=1e-4, epsilon/pressure <=1e-3,
selected alpha=0 spatial residual blocks <=1e-3, relative traction change
<=5e-4, relative matching mismatch <=1e-5 and negligible absolute bound supply.
All normalizations and thresholds remain unchanged. Numerical acceptance,
steady convergence and qualification inside the older narrow reference box
remain distinct.

The existing checker reads exact native report/checkpoint pairs at
T=3.5,4,4.5,5, without interpolation. It checks the saved reference, coupled
coefficients, all GLM slots/time labels, physical signature and accepted ledger.
Residuals are the retained selected uncancelled alpha=0 operator; no Lift or
epsilon-trace replacement is applied. No governing solve or history/time
advancement was added.

## Uniform G: a selected discrete steady candidate

| Window | Normalized k change | Full-outer epsilon change | Complete screen |
|---|---:|---:|---|
| 3.5–4 | 7.653e-5 | 8.383e-5 | pass |
| 4–4.5 | 2.864e-5 | 3.726e-5 | pass |
| 4.5–5 | 9.806e-6 | 1.435e-5 | pass |

At T=5 the common-normalized U/W/k/epsilon spatial residuals are
2.716e-5 / 1.519e-5 / 6.095e-5 / 7.678e-5. Refreshed matching mismatch is
1.533e-8. Correction activation and signed/absolute supply are zero in all
three windows and the cumulative recorded correction integral is zero.
This supports **a steady candidate of the selected discrete legacy equations**
at the declared tolerances, rather than a fixed point maintained by bounding.

The missing full-epsilon confirmation now covers the entire transported
.5<=z<=4 region, not just the previous z>=2 samples. A new read-only probe
operation exports the actual full physical consistent epsilon mass matrix,
including prescribed and matching-trace coefficients. It retains the existing
strict native-checkpoint reader: exact state, reference, physical signature,
native differential/history slots, time labels and checksum are validated.
It performs no constitutive substitution, transport solve, state correction or
initialization Lift.

The matrix is symmetric positive definite, its constant-field measure matches
physical area 21, and the full-window mass norms agree with independent native
GL4 element-row integration to 1.4e-20. The fixed common epsilon normalization
remains E0*sqrt(12), so the enlarged reporting region does not relax the test.
All four coefficient vectors match the checkpoint bits exactly. Corrupted
state/history/signature/checksum and nonzero-alpha requests are rejected.

This is not an exact-root or spatial/temporal-accuracy certificate. Mean wall
traction is .15280236 versus .16, **4.50% low**. Inner/outer strong-divergence
L2 norms remain .06509/.003235 despite near-roundoff weak continuity integrals.
The covariance margins remain positive, .08201 inner/.15971 outer. Settling
therefore does not establish neutral-reference preservation or pointwise
incompressibility.

## Upstream-one: completed, but not steady

| Window | Normalized k change | Common screen |
|---|---:|---|
| 3.5–4 | .01897 | fail |
| 4–4.5 | .02062 | fail |
| 4.5–5 | .02234 | fail |

The k drift is increasing across these windows. At T=5 the selected
U/W/k/epsilon residuals are .2305 / .1589 / 8.5176 / .4679, and the refreshed
matching mismatch is .0007447. These are far outside the scientific steady
screen, even though every scheduled step passed its retained hard numerical
acceptance checks.

The bound activates on all 800 steps of the last half-unit window. Its signed,
positive, negative and absolute correction **amounts over that window** are
.68915 / 1.29737 / -.60822 / 1.90559. Dividing by the .5 duration gives the
corresponding rates; the absolute rate normalized by K0*24/(2/.4) is 1.5553.
The cumulative signed correction integral is 7.21957. These are integral/rate
budgets, not percentages of total turbulent energy or correction/raw-step
norm ratios.

Mean wall traction is .0054965, about 96.56% below reference. This remains a
substantial bounded transient, not an uncorrected-equation steady state;
statistical stationarity is not demonstrated. Full-outer-epsilon confirmation
is unnecessary for this negative verdict because several independent required
conditions already fail by large margins.

## Diagnostic build and reproduction

Source base was `ba0f306`; read-only metric/confirmation implementation is
`6b4d554`, with exact executed source hashes retained in the evidence. Linux
Release flags retain `-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS`.
Only a separately isolated diagnostic probe was built, with parallelism four.
Its SHA256 is
`0df6a85daa4d5c88274c952bb1babb0a68c0e7dfe5e03bbdbbdf2694785b55d6`;
its linked-library closure exactly matches the previously validated frozen
probe. No production executable or library was rebuilt or overwritten.
Native probes use one rank and all computational-library thread controls one;
concurrent independent experiments were preserved within the 72-worker cap.

```bash
R=/media/jay/新加卷/nektar-v2026-runs/cls-awf-steady-final-pair-20261110
S="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
python3 -B "$S/confirm_cls_awf_outer_epsilon.py" \
  --common-panel "$R/common-panel" --case UniformG \
  --probe "$R/runtime/probe/RANSWallResolutionProbe" \
  --output "$R/UNUSED-full-epsilon-confirmation" --contract-tests
```

Use an unused output directory. The common-panel snapshots are reproduced
with the existing `check_cls_awf_steady_windows.py`, the unchanged policy and
the two original campaign paths. The explicitly corrupted checksum fixture
is a rejected test copy, never a restart or a modification of a production
checkpoint. This assessment creates no new evolving campaign.
