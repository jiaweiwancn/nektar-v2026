# H=.5 graded budget-AWF continuation from T=5 to T=6

2026-11-09, Asia/Shanghai. Both requested campaigns have restarted natively
from their completed step-8000 checkpoints. Their **next production target is
T=6, global step 9600**. Intermediate T=5.5 is a retained block endpoint.
The completed T=5 campaigns, frozen runtimes and checkpoint bytes are preserved.

The source/guard implementation is `211ab2d` on
`codex/cls-epawf-pimple-one-sweep`. Linux authoritative execution uses the same
previously validated solver SHA256
`177ce5e54e7de514fa019122e46d3a5d709455cd091bd3023fbfa745c48c9bfc`
and an exact copy of its linked-library closure. No C++ rebuild or numerical
step-map change was required. Native dt remains .000625, IMEX1, with one
U/p/outer-epsilon/k solve and one final logical bound call per accepted step.
The accepted budget source, H=z_K=z_m=.5, endpoint K, model companions,
floor/reserve, work controls, pressure correction and hard checks are unchanged.

## Verified parent checkpoints

| Model | Parent step | Native parent time | Parent State.bin SHA256 |
|---|---:|---:|---|
| Regularized/Linear | 8000 | 5.000000000000251 | 54f28fea7ac1f8fa532849cfae926bdbadad7c69ebc542e3489dc4fc5078ec6b |
| Regularized/Saturated | 8000 | 5.000000000000251 | 1e7912e26bfdcfd3219b7c96626c361b49b7910023207e28eda0d8487200a61d |

The controller validates successful parent completion, exact report/checkpoint
signature, coupled state, reference/defect, every native history slot and time
label, model/tail/inventory settings and executable/library hashes. It excludes
a live parent supervisor, driver or orphan native child and positively acquires
the parent lock. New output roots retain that lock through their supervisors;
the parent seal records the successor instead of changing the old progress
or checkpoint payload.

One isolated restart step per model is **verification-only**. Startup bounding
is inactive and the initial committed state is unchanged. Each first production
successor checkpoint equals its isolated restart checkpoint byte for byte,
including all native histories. The actual production passes retain exactly
one governing solve of each equation, one final bound call and no additional
governing solves. Focused Linux endpoint/history/signature, terminal-rejection,
model/tail/inventory and target guards passed.

## Observed continuation

At **21:38:24 Asia/Shanghai**, both have completed and validated the first
32-step production block after restart and remain running:

| Model | Verified completed block step | Verified completed T | Status | Next target |
|---|---:|---:|---|---:|
| H=.5 graded Regularized/Linear | 8032 | 5.020000 | running | T=6 / step 9600 |
| H=.5 graded Regularized/Saturated | 8032 | 5.020000 | running | T=6 / step 9600 |

Each completed block has exactly 32 U/p/epsilon/k solves, 32 implicit callbacks
and 32 logical bound calls, with zero additional governing solves. The common
normalized alpha=0 U/W/k/epsilon residuals are
.00011175/.00006097/.00061063/.00029293 for Linear and
.00011185/.00006102/.00061110/.00029326 for Saturated. These small endpoint
residuals do not replace the required window persistence test. The first block
receipts cost 263.29/246.42 wall seconds; at similar load, the remaining
extension takes approximately 3.6/3.4 hours, conditional on admissibility.

Both have zero post-restart bound activation and zero cumulative correction
supply. The other four campaigns remain under their existing supervisors and
physical settings. The observed native/compiler worker count is 18, below the
hard 72 limit; all six computational-thread controls are one and each native
pilot uses one rank. Run-volume free space exceeds 1.4 TiB.

The [compact evidence](cls-accepted-budget-awf-evidence/to-T6-20261109/native.json)
retains current supervisor/driver/native start tokens, immutable checkpoint
copies, state/history hashes, solve counts, actual launch commands and the
independent restart work. The [artifact index](cls-accepted-budget-awf-evidence/to-T6-20261109/artifact-index.json)
locates heavy inputs, native reports, verification packets and frozen libraries
under `/media/jay/新加卷/nektar-v2026-runs/cls-awf-05-to-T6-20261109`.

## Scientific interpretation and next assessment

The unchanged `AcceptedAWFExploratoryWindowsV1` policy requires three
consecutive .5 windows: normalized U/W/k changes <=1e-4, epsilon/pressure
<=1e-3, selected spatial residuals <=1e-3, relative traction change <=5e-4,
relative matching <=1e-5 and negligible correction supply. T=6 supplies two
additional half-unit windows after T=5; it does not guarantee their passage.

The previous [steady check](cls-awf-steady-check.md) found only the most recent
window qualifying. Its Linear windows used the native block offset, ending
at T=4.99; an exact 4.5/5/5.5/6 comparison must use the retained exact native
checkpoints, not interpolate or relabel that earlier window. Positive steady
classification also needs full-outer-epsilon field confirmation, as already
documented. No criterion is weakened or automatically declared satisfied by
reaching the target. Reference-traction error remains separate from settling.

The controllers use the existing durable supervisor's native checkpoints,
32-step blocks and explicit 7200-second timeout. Numerical rejection remains
terminal; infrastructure recovery preserves the committed history. The runs
are **running**, not completed T=6 assessments.

The already executed launch commands use this committed helper and separate
output roots:

```bash
S=/media/jay/新加卷/nektar-v2026-runs/cls-awf-05-to-T6-20261109/source/solvers/IncNavierStokesSolver/Tests/RANS
O=/media/jay/新加卷/nektar-v2026-runs/cls-awf-audit-cache-temporal-20261109
N=/media/jay/新加卷/nektar-v2026-runs/cls-awf-05-to-T6-20261109
python3 -B "$S/extend_cls_awf_completed.py" \
  --campaign "$O/handoff-budget-05-linear-v2/trajectory" \
  --output "$N/regularized-linear" --block-steps 32 --timeout 7200
python3 -B "$S/extend_cls_awf_completed.py" \
  --campaign "$O/handoff-budget-05-saturated-v2/trajectory" \
  --output "$N/regularized-saturated" --block-steps 32 --timeout 7200
```

These parents now have recorded live successors; the guard rejects a duplicate
launch. The commands explain the recorded operation, not permission to bypass
ownership. Existing successor recovery uses its own supervisor and manifest.
