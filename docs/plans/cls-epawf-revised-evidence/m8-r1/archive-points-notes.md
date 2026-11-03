# Complete offending inner-point export — later archive supplement

The original phase-A ledger remains unchanged. This supplement exports every negative covariance evaluation entry from its existing `source_volume_inner`, `physical_sample_K`, and `interface_inner` constitutive sets over all 12 accepted epochs. It introduces no new sampling set, solver run, trajectory or constitutive modification. These are finite reconstruction observations, not a continuum positivity/negativity certificate.

The [summary](archive-points-summary.json), [runner](archive-points-runner.json), [exact log](archive-points-run.log), and [integrity index](archive-points-index.json) retain the source, command and hashes. The unchanged helper is bound by SHA256 `87fbb11cf86163035499bd55da355ad568a63854439b086745144c169e453501`; it subclasses the completed original archive driver and imports only the clean ab5 numerical oracles. The Linux job used CPU 50, one worker and all six numerical thread settings equal to one. It passed 612 reconstruction/export checks in 230.98 seconds.

All **18,180** negative entries belong to `source_volume_inner`; the two other declared constitutive sets have zero negative entries at every accepted epoch. Per-set population counts, negative counts, state/history hashes and scalar minima agree with the original ledger. Repeated coordinates and evaluation entries are deliberately retained. Counts are not negative physical volumes or unique geometric-point counts.

| Accepted epoch | Negative source-volume entries | Physical K sampler | Inner interface |
|---:|---:|---:|---:|
| 0 | 1240 | 0 | 0 |
| 1 | 1240 | 0 | 0 |
| 2 | 1240 | 0 | 0 |
| 3 | 1241 | 0 | 0 |
| 4 | 1247 | 0 | 0 |
| 5 | 1272 | 0 | 0 |
| 6 | 1304 | 0 | 0 |
| 7 | 1461 | 0 | 0 |
| 8 | 1528 | 0 | 0 |
| 9 | 1629 | 0 | 0 |
| 10 | 2004 | 0 | 0 |
| 11 | 2774 | 0 | 0 |

The full `negative-inner-points.ndjson` stays on Linux at:

`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r1/archive-points-001/negative-inner-points.ndjson`

Its exact byte size is `42238882` and SHA256 is `4259f6dd08247a6a9ffc3136d81751738a2dc71ed927928eca15a85ea3f21e5a`. Each row contains time/epoch/set, location and element/side identity; full G, S and R; minimum eigenvalue/eigenvector; local k, physical sampled K, k/K, nuA and rR; the current source controls; and a full comparison with the actual initial coefficient state at that identical point. The initial sampler/nuA uses its own initial-state ray cache. The analytic neutral shear is separately labelled. `epsilon_profile` denotes the reconstructed inner profile, not an outer transported-field measurement.

The original first-event table uses the minimum over the original t=0 finite population: eigenvalue `-0.0976947699300747` at `[-0.9602898564975362,0,0.5092722755356407]`. That remains authoritative. NDJSON row zero instead follows export enumeration order and has eigenvalue `-0.005748304448373792` at a different t=0 point. It is not a replacement first-event minimum.

The full accepted-prefix minimum is `-0.09983101642016129` at epoch 11. Its same-point initial value `-0.09774042097023938` is reconstructed at the later adaptive-rule location; that location need not belong to the original t=0 rule. It must not be substituted for the original initial finite-set minimum. Every row exposes k/K and strain changes, but this export makes no source-cap or Newton-causality inference.

The original phase-A receipt finished at `2026-11-03T09:04:03.800222+00:00`. This supplemental extraction ran from `2026-11-03T11:26:02.869957+00:00` to `2026-11-03T11:29:53.854686+00:00`, after that gate and after subsequent diagnostic work had begun. It closes the all-offending-points export gap without claiming the full NDJSON existed before those later trajectories. Original gate timing and historical receipts are preserved. No extra trajectory or numerical rerun is required for this supplement. M8 remains open.
