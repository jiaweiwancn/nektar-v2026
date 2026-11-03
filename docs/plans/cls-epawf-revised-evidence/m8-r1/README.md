# M8-R1 evidence

This directory retains compact evidence for the bounded first-failure and
joint constitutive-feasibility task. **M8 remains open.** A completed audit
or an expected failed-run reproduction is not physical qualification.
The [report](../../cls-epawf-revised-m8-r1.md),
[first-failure ledger](first-failure-ledger.json) and
[event table](event-table.md) are the main entry points.

| Evidence | Contents |
|---|---|
| [Archive ledger](archive-ledger.json) / [index](archive-curation-index.json) | Complete twelve-epoch reconstruction, source history, pressure/continuity/traction/flux/energy/source observations and historical producing binary. |
| [All offending points](archive-points-index.json) | Exact external 18,180-point export with full tensors and same-point initial comparisons; collected as an explicitly later archive-only supplement. |
| [B1](nonlinear-b1-ledger.json) / [B2](nonlinear-b2-ledger.json) / [index](nonlinear-curation-index.json) | Full failed425-coefficient state/history/residual, all original/fresh line-search outcomes, invalid-point reconstruction, full directional ladder and retained pressure/gauge accuracy failure. |
| [Preflight](preflight-summary.json) / [notes](preflight-notes.md) / [index](preflight-curation-index.json) | Actual original425 and improved1571 unknown initializations, two quadrature levels, raw unanchored loads, positive improved finite covariance and failed applicability coverage. |
| [Controlled trajectories](controls-ledger.json) / [index](controls-curation-index.json) | Final two independently initialized dt controls; both fail at their first nonlinear stage. |
| [Runtime](runtime-summary.json) / [focused checks](focused-checks-receipt.json) | Exact source/executable/library/tool/package/environment pins, before/after identity and the passing three native probes plus two zero-step guards. |
| [Kernel](kernel-summary.json) / [proof](kernel-outer-proof.json) / [retained witnesses](kernel-retained-witnesses.json) | One fixed bounded suite, rigorous restricted outer interval bound and all historical counterexamples, with production adoption false. |
| [Candidate contract](candidate-contract.md) / [preregistered bytes](kernel-contract-preregistered.md) | One restricted next research contract, full inner continuation/dependencies/budgets/trace/traction changes and unresolved lambda→0 behavior. |

`run-*.py` are the exact root-owned Linux supervisors. The committed
`cls_epawf_m8_r1_*.py` helpers live under
`solvers/IncNavierStokesSolver/Tests/RANS`. Historical executed helper
snapshots are retained where a later corrected helper has different bytes.
Failed postprocessing/build attempts are retained and identified separately
from native scientific failures. Pure `*-curate.py`/`curate-*.py` helpers and
`assemble-ledger.py` select/hash already computed data; they run no solver.

The large native results, dense Jacobians, full field histories and exact
unabridged fingerprints remain under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r1` on the designated
Linux server. Component indexes name their exact paths, file-byte hashes,
canonical vector hashes and any omitted arrays. No omitted point population
or successful stage is inferred from a minimum or a residual norm.

Original, CubicScaled and Regularized model identities and all historical
M0–M8 receipts remain intact. No new constitutive law has been compiled or
selected. The six-slot matrix is exhausted: four native trajectory starts and
two conservatively counted same-stage diagnostic attempts. Four preflight
and two focused guard calls have NumSteps=0; direct unit probes are separate.
