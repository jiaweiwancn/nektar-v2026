M5 evidence uses the authoritative Linux run root
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised` and Tailscale only.
The small reports and figures here accompany the
[M5 closure report](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-m5.md).
Full fields, meshes, checkpoints, compiler/test logs and failed attempts remain
on Linux, indexed by `raw-artifact-index.json` there.

`verification.json` records the complete native-CG qualification at
42b8ba83897cca1202b91032a97492c7d463c503. `note-verification.json`
records the separately registered exact-note checks. `regression.xml` records
composed coverage of the other 26 registered tests, including the full M4
qualification. The original run passed 25; its optimizer-cache probe failure
and the corrected repeated execution are preserved separately.
Together these cover all 28 registered tests without repeating the passing
full M5 solve. Source/runtime fingerprints justify reuse of that qualification
when only the supplementary test and CMake registration were added. A subsequent filename-only probe fix is
validated by `session-replay-verification.json`; the solver and libraries stay
identical.

`integration.json` binds these results, actual final-stage/physical-output
checks, all ten frozen M5 equation-owner IDs, original physical failure gates,
and preservation of the unrelated HIT checkout. The approved contract and
revised note retain their original hashes and planned baseline records;
current delivery is recorded separately here.

The collection and plotting scripts are reproducible tools. For example,
with all four BLAS/OpenMP thread counts set to one on Linux:

```bash
python3 collect-integration.py \
  --root /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised \
  --qualification m5-attempt-005 --note m5-attempt-006 \
  --regression m5-attempt-007 --replay m5-attempt-008 --output /path/to/new-audit-directory
python3 plot-verification.py verification.json /path/to/figures
```

Do not use a pre-existing output directory for collection. The collector
requires passing receipts, identical qualification inputs/runtime, actual
converged stage constraints, physical Cartesian VTU points/values, and unchanged
Linux source/HIT checkout states before it emits a closed milestone receipt.
