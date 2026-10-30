# Milestone 6 preservation evidence

Generated from the authoritative Linux `m6-preservation/final-study` outputs
with `Tests/RANS/plot_m6.py`. Large fields, logs and history payloads remain
under the Linux run directory recorded in the parent report.

- `preservation.png/pdf`: five-flow 2D/full-3D profile errors for standard,
  original, fixed-correction and switched-correction policies. Crosses mark
  stopped runs; the dashed line is the 0.1% engineering target.
- `refinement.png/pdf`: fixed-H vertical and polynomial refinement.
- `profiles.png/pdf`: common physical-height errors at four downstream stations,
  using the corrected element-local profile sampler.
- `wall-response.png/pdf`: production/dissipation mismatch and switch response.
- `metrics.json/csv`: final, peak and first/upper-layer error norms, flow times,
  and stage/wall summaries, including stopped and inconsistent controls.
- `evidence.json`: closure gates, source/binary/report hashes, environment,
  continuation and backend comparisons, perturbations, and numerical-column
  results. Its retained robustness trial records the earlier degree-five stop.

Interpretation and limits are in [the Milestone 6 report](../neutral-abl-rans-awf-milestone-6.md).
The numerical-column copy runs one flow; the independent analytical comparison
runs five flows. The engineering target does not assert exact discrete
preservation or calibrated terrain performance.
