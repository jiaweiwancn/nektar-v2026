# Flat WJ–AWF outlet B qualification

The baseline is `aabffc99f67a2dc941d02b356c32e1c496ad9cfb` on `codex/wj-earsm-awf`. B is the explicit flat candidate under `ZeroTransverseV1`. It imposes zero vertical velocity and the frozen normal traction at x=15, leaves axial velocity free away from the shared ground vertex and reports the transverse reaction. A is the `FullTractionV1` default. It missed 16 rows in the original two-mesh comparison and 180 in the isolated collar comparison. B passed the original 690 rows and missed 122 in the latest radial comparison. [Prior decision and exact B2 evidence](wj-h27-task2-outlet-decision.md)

Use the existing ext4 run root on `jay@100.88.70.60`. Keep the verified B1, B1c and B2 fields at zero simulated time. Their sessions use Q3 velocity and scalars, Q2 pressure, 6×8 native panel quadrature, 1e−12 linear tolerance and 1e−10 coupled tolerance. The pinned runtime uses `ExactInnerFiniteNuV1`, `ScalarCdV2`, finite molecular viscosity, the complete positive stress split, sampled AWF force, current-K cutoff and full production. H=.15 code is 6.0 mm above the modeled wall. Every solve uses a typed zero-time field guess with native lifts and fresh history. [Runtime, state and grid inventory](wj-h27-linux/task2-outlet-qual/cp0/outlet-qual-plan.json) · [Model note](../literature/wj_earsm_neutral_abl_ep_awf_note.tex) · [Cutoff appendix](../literature/wj_earsm_budget_cutoff_awf_appendix.tex)

## Spatial ladder

Replay the B2−B1c owner comparison with the original common subcells, endpoint sampling and signed fixed scales. Retain all 690 rows, the 29 stations, both owners at shared edges and the fixed experimental window x∈[−2.5,8]. The flat Fig. 3 reference is at the virtual hill station x=0. The rough-case 7 mm observation is 3.5 mm above the modeled wall. [Published measurement conventions and paper provenance](wj-h27-boundary-audit.md#5-prospective-physical-development-envelope-distinct-from-numerics)

B3 bisects the eight B2 descendants in the original [.08125,.10640877656902] radial patch. Its 111×77 grid preserves all 112 x edges and every existing radial edge. A 111×73 fallback splits only the four B2 lower children in [.08125,.10] and starts a separate refinement sequence. The 12 GiB admission uses a 10.41 GiB forecast for the full mesh from the B2 peak, its element-count ratio and a 10% margin. The fallback forecast is 9.87 GiB. Neither estimate certifies factor memory. [Exact proposed axes and measured anchors](wj-h27-linux/task2-outlet-qual/cp0/outlet-qual-plan.json)

The lowest-cost candidate is the coarsest state that passes each direct finer comparison. The independently changed radial witness bisects [.05,.06], [.06,.06875], [.06875,.08125], [.10640877656902,.1125] and [.1125,.11541061001735] on the selected mesh. These cover the previously shared blend intervals and both adjacent flanks. A failure in a newly split interval calls for a repeated split on its support. Keep separate comparisons for separate refinement supports. The prior x control split the terminal tip only. The proposed independent x control holds the radial grid fixed and bisects the last two x columns and the x owners touching retained experimental stations, subject to measured admission on the selected radial grid. A reduced x witness requires an explicit owner-coverage assessment.

Apply the following RMS and MAX limits to the original rows and the fixed region. The independent audit assembles the current residual, checks native positivity, root, covariance, principal and current-K guards, verifies weak budgets and tracks shared and one-sided interface fluxes. Raw gradients enter all acceptance comparisons.

| Gate | RMS / MAX or residual limit |
| --- | ---: |
| Spatial U and W | 1e−4 / 1e−4 |
| Spatial k/Kin, R/tau0, Pi/tau0 and shear/tau0 | 1e−3 / 1e−3 |
| Spatial velocity gradients and local outer epsilon | 1e−2 / 1e−2 |
| Current momentum and continuity | 1e−9 |
| Current k and outer epsilon | 1e−8 |
| Physical U and W | .01 / .02 |
| Physical k/Kin | .05 / .10 |
| Physical local outer epsilon | .10 / .20 |
| Physical ground shear relative to tau0 | .05 / .10 |
| Inlet projection U, k/Kin and velocity gradient | 1e−5, 1e−4 and 1e−3 for both norms |

For final controls, use an independent 8×12 native quadrature solve at the base tolerances and a separate 6×8 solve at 3e−13 linear and 3e−12 coupled tolerance, both on the selected mesh. Each also applies the 1e−4 inner-k RMS/MAX limit. Ground shear against Tref(Kin) is reported separately.

## Resources and conditional extent test

The 64,800 s ceiling includes all earlier charges, and the 3,600 s protected reserve carries forward. Stage the complete B3 experiment first. It includes replay, admission, setup, two possible continuations, current and owner audits, checkpoint publication, seven seconds of reservation grace per run and 300 s for final evidence. The subsequent stage estimates are provisional until measured admission. The watchdog uses one native rank on a quiet permitted CPU on the launcher's NUMA node and six one-thread settings. Use the 80-worker shared-host ceiling, with a preemptive pause at 78 and resume at a projected census of at most 76 held for four seconds. Retain the four-job build ceiling, two 8 GiB disk caps, 15 GiB free space and 16 GiB build cap. [Idempotent amendment and opening ledger](wj-h27-linux/task2-outlet-qual/cp0/outlet-qual-amendment.json)

If a spatially qualified recipe and a complete forecast admit the extent test, move the exit from x=15 to 17.5. Keep the upstream grid through x=14.375, move the terminal grading and retain x=15 as an interior station. The extent test compares both extents at matched p/q/tolerances after x=15 native parity. Report each state relative to the inlet reference, then its signed change between extents. Apply the physical limits independently to the signed U, W, k, outer epsilon and ground shear changes. Report stress, Pi and traction beside within-B numerical indicators. [Frozen B2 extent design and thresholds](wj-h27-linux/task2-outlet/cp4/outlet-location-eligibility.json)

**TASK2_COMPLETE=NO. TASK3_READY=NO.** The one-day criterion is `1.25*C_total ≤ 86,400 s` per advancing trajectory.
