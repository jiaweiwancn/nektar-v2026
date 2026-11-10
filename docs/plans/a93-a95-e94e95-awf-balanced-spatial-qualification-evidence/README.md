# Balanced spatial qualification evidence

New supplemented export from the isolated authoritative Linux experiment:
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-balanced-20261110/evidence-final`.
The collector at `4a0382dfc6b8cba10da1d9dc8d410e5daeceb471` returned **0**.
All six aggregate gates pass. Historical evidence is unchanged.

`summary.json` is the index. `complete-results.json.gz` contains the complete
new comparison bodies, candidate screens, physical-term/localization records,
full native and common-test norms, twelve primary evolved audits, trajectory
and history identities, all accepted-state/stage margins, final-state receipts,
cutoff motion and correction budgets. Its SHA256 is in the index.
`previous-twelve-complete-comparisons.json` separately contains all twelve old
required bodies with endpoint/request/runtime/consumed-input provenance.
The native reports, full basis/load packets, checkpoints and logs remain at
their hashed Linux paths; large numerical artifacts are not copied into Git.
`support-energy-EP.json` and `support-energy-Accepted.json` preserve the signed
dual-energy localization on the initial screening spaces. Element contributions
sum to the assembled squared norm; isolated stress loads are not additive norms.

Six natural model/wall combinations have independent four-member spatial/dt
comparisons at T=.0016 and .016. E94E95 Natural Accepted additionally has
unchanged-method restart continuation through .1 to .5, with both spaces and
their dt/2 companions at the final horizon. A95 Accepted and A93 EP have
selected prescribed-kappa=.41 short controls. Coverage is explicit; the new
duration does not qualify other models/walls/policies.

The failed .004 pressure-uncertainty comparison remains in `duration` in the
complete results. Exploratory/incomplete-endpoint failures remain separate in
the previous supplement. Negative controls succeed only for their intended
rejection reason. The old overall unresolved receipt and its eight epsilon
temporal-order observations retain their original interpretation. New dt pairs
bound uncertainty and do not establish temporal order.

Runtime manifests distinguish the 577b5f1 solver snapshot from the revised
0a5496f diagnostic probe. Actual loaded-library maps accompany both. The
governing-source reuse proof, targeted CTest receipts, native diagnostic
regressions, complete final quadrature body and final audit commands are
included. `SHA256SUMS` verifies every other file in this directory.

```bash
sha256sum -c SHA256SUMS
python3 - <<'PY'
import gzip, json
with gzip.open('complete-results.json.gz', 'rt') as stream:
    result = json.load(stream)
assert result['status'] == 'pass'
assert all(value == 'pass' for value in result['gates'].values())
assert len(result['required_short_comparisons']) == 12
print(result['gates'])
PY
```

The portable common-test sampling witness can be replayed without the solver:

```bash
python3 /path/to/source/solvers/IncNavierStokesSolver/Tests/RANS/reproduce_rans_common_sampling.py \
  --witness portable-sampling-witness.json --output /new/path/witness-result.json
```

The full report and Linux reproduction commands are in
[`../a93-a95-e94e95-awf-balanced-spatial-qualification.md`](../a93-a95-e94e95-awf-balanced-spatial-qualification.md).
