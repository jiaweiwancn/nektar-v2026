# CT9 evidence and reproduction

The immutable parent is CT8 at `0c28610e685aef95f7e6f21b9a4fdcd305dbf319`.
CT8 physical, extraction, observation and curved-history decisions are unchanged.
`contract.json` precedes the new trajectories. Later numerical-component or
bounded-memory amendments are separate files; they do not alter physical gates.

The authoritative archive is
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-9-20261114` on
`jay@100.105.247.110`. The actual mount was verified as `/dev/sda2`, fuseblk.
The older `/media/jay/新加卷` paths in sealed history are not silently rewritten.
Full matrices, CSR snapshots, vectors and histories stay in that archive.

Run every numerical command through the ownership-aware launcher, with the
same thread environment **before** Python imports or solver launch:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
root='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-9-20261114'
old='/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-8-20261114'
export TMPDIR="$root/tmp" OMPI_MCA_orte_tmpdir_base=/tmp/nektar-ct9-ipc
tests="$root/source/solvers/IncNavierStokesSolver/Tests/RANS"
# Every CONTROL and OUTPUT below must be new. Never overwrite a sealed run.
python3 -B "$tests/ct_accepted3_command.py" --control "$CONTROL" -- \
  timeout -k 30s 7200 python3 -B "$tests/ct_accepted9_build.py" \
  --root "$root" --ct8 "$old" --output "$OUTPUT"
```

The other owned commands use these scripts and arguments:

| Purpose | Script and arguments after the launcher/timeout |
|---|---|
| Sole new dense history | `ct_accepted9_pressure_run.py --ct8 "$old" --runtime "$old/build-v3/runtime" --contract "$root/source/docs/plans/cao-tamura-accepted-9-evidence/contract.json" --output "$OUTPUT"` |
| Correct read-only pair aggregation | `ct_accepted9_pressure_compare.py --ct8 "$old" --finer "$root/pressure-finer-v1" --output "$OUTPUT"` |
| Frozen pressure decomposition | `ct_accepted9_pressure_diagnose.py --ct8 "$old" --runtime "$old/build-v3/runtime" --output "$OUTPUT"` |
| Small native controls | `ct_accepted9_controls.py --ct8 "$old" --runtime "$RUNTIME" --old-backend-controls "$root/controls-v3" --output "$OUTPUT"` |
| Fresh sparse physical equivalence | `ct_accepted9_equivalence.py --ct8 "$old" --runtime "$RUNTIME" --controls "$CONTROLS" --output "$OUTPUT"` |
| Increasing flat capacity | `ct_accepted9_capacity.py --ct8 "$old" --runtime "$RUNTIME" --output "$OUTPUT"` |
| Exact-coefficient numerical recheck | `ct_accepted9_frozen_recheck.py --runtime "$RUNTIME" --parent "$PARENT_CASE" --output "$OUTPUT"` |
| Primitive and process-exit tests | `test_ct_accepted9_gate.py` |

The executed command arrays, source revisions, executable/library identities,
elapsed/process CPU, RSS and closed resource high-water records are the source
of truth for the final report. `request-manifest.json` binds each immutable
request; `completed-artifacts.json` binds its results. A changed binary, solver
component or history gets a new destination and, where applicable, fresh
initialization. Cross-backend equivalence is not cross-backend restart.

The compact final export must bind the **actual** build, controls, equivalence
and capacity labels through the collector's required options. A successful
collector does not admit target evolution. The physical pressure gate is
recomputed from complete primitives independently of implementation completion.

Preserve failed build attempts, the initial incomplete temporal export,
`capacity-v1` and the exact-state numerical rechecks. A corrected numerical
component cannot retroactively qualify an earlier runtime or resume its evolved
history under a modified signature.
