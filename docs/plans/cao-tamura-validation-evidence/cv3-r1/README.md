# CV3-R1 native evidence

The [execution record](../../cao-tamura-cv3-r1-reference-feasibility.md) explains
the results and limitations. The [protocol](../../cao-tamura-cv3-calibration-protocol.md)
defines the gates. The generated `evidence.json` records the exit decision;
short, failed or deliberately retired controls cannot close R1.

`record_cao_tamura_r1.py` exports reference receipts, source/binary bindings,
native diagnostics, checkpoint comparisons, exact producer versions and the
historical compiled-source snapshots needed to reproduce earlier controls.
Restore archived producers at their recorded source paths for reproduction.
The native CSV and figure include accepted steps only; forced startup rows
are retained in the CSV and excluded from the figure. Recomputed Newton-action
diagnostics do not replace the independent physical stage audits.

Large native fields, histories, meshes and full logs remain under
`/media/jay/新加卷/nektar-v2026-runs/cao-tamura-calibrated-validation/calibration/cv3/ep-transition-v2/r1/`
on Linux. Artifact hashes and producing paths are retained in the receipts.
Run the exporter there with the immutable build, reference, guard, history,
matched-field and Release CTest paths recorded in its command receipt.
Supplying `--developments` requires all four qualifying native cases.
The [terrain example](../../../../solvers/IncNavierStokesSolver/Examples/TerrainHill/README.md)
contains the isolated Release build recipe against the existing MPI installation.

These are unfitted feasibility results. They contain no calibration freeze,
reserved Fig. 3d or hill-validation scores. Experimental agreement and CV4
spatial, temporal, domain and sampling uncertainty remain separate gates.
