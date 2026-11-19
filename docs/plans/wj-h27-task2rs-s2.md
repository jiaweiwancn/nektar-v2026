# Task 2R-S S2 resolution tests

**S2=INCONCLUSIVE. ORIGINAL_PROBLEM_SPATIAL_QUALIFICATION=NOT_MET.** The independent tests separate the collar's radial resolution sensitivity from the smaller outlet corner cells. A converged, strictly nested refinement increases the raw pressure, shear and velocity-gradient differences at the vertex. The coupled angular amplitude remains unresolved.

The original boundary problem uses the `ExactInnerFiniteNuV1` loads and reference with the registered coupled solver. These are unforced, zero-time Q3 velocity and scalar, Q2 pressure solves with native 6×8 panel quadrature. Every original observation and its fixed engineering tolerance remains in each comparison. The numerical results below use raw physical-owner fields.

## Radial and joint controls

The 114×62 mesh splits only the first modeled-wall interval twice at fixed x. Its last x width is .00244140625 and its first r width is .000125. The native solve from the registered 114×60 zero-time seed takes 489.957 wall seconds and peaks at 8,391,999,488 bytes RSS. The native norm is 4.53e−13. The independent current residuals are 3.60e−11 momentum, 1.22e−14 continuity, 5.78e−11 k and 5.44e−11 epsilon. All physical envelopes pass, including inner k RMS/MAX .0196460/.0402985 and ground shear .0121769/.0219842 against tau0. The parent 109×60 comparison still misses 28 rows. The same-x 114×60 to 114×62 comparison misses 12 rows, all in Pi, raw shear traction, and velocity gradient at the ground–outlet vertex. The largest k difference is 2.68e−7 in Kin units and the fourfold reduction in first radial width changes vertex D from .255468 to .254598. The H owner-gradient k diffusion jump is 5.09e−7 RMS, with shared mass and advective k flux cancellation.

The first joint design halves the outlet x and ground r widths to .001220703125 and .0000625 and splits both radial owner intervals adjacent to r=.10. Nine upstream merges allow 106×65 elements under the measured process-tree limit. Their screened source-field interpolation differences are at most 4.85% of each field budget before the solve. The native solution takes 770.183 seconds, peaks at 8,216,969,216 bytes and converges to 1.12e−12. Independent momentum, continuity, k and epsilon residuals are 1.02e−10, 1.79e−14, 2.43e−10 and 2.36e−10. It passes all physical envelopes.

The solved joint comparison misses 180 of 690 original rows and 54 of 207 supplementary rows on the merged strips. These failures span stations whose x intervals also remain identical, with peaks at the refined r=.10 interface between collar owners. For example, the Rrr difference at x=5 is .007223 in tau0 units with the same x interval on both meshes. This shows why an interpolation screen for upstream x resolution cannot certify the actual solution after changing the global radial collar. The joint physical strong-divergence RMS falls from 1.9492e−4 to 1.2362e−4. At the upper r=.10 outlet owner, D falls from .03104 to .003716. The Ryy difference there is −.0221801. Its trace-free-planar diagnostic is 1.42e−5, whereas the raw acceptance row uses the full gradient. At the ground–outlet vertex D grows to .337600. The corner and collar require separate evidence.

| Common-owner pair | U RMS / MAX | W RMS / MAX | k/Kin RMS / MAX | Stress RMS / MAX | Pi RMS / MAX | Traction RMS / MAX | Gradient RMS / MAX | Local outer epsilon RMS / MAX |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 114×60 to 114×62 radial | 3.03e−8 / 7.49e−7 | 1.43e−8 / 5.11e−7 | 1.09e−8 / 2.68e−7 | 1.40e−5 / .00067465 | 2.30e−5 / .00103615 | 2.86e−5 / .00144079 | .00059967 / .0310243 | 1.71e−10 / 1.02e−9 |
| 114×62 radial to 106×65 joint | 6.90e−6 / 2.65e−5 | 4.15e−6 / 2.06e−5 | 1.22e−5 / 4.37e−5 | .002166 / .0221801 | .001565 / .0125842 | .0007855 / .0120532 | .005623 / .2595401 | 2.77e−7 / 2.14e−6 |
| 114×62 radial to 110×62 isolated | 3.44e−7 / 1.02e−5 | 7.40e−7 / 1.37e−5 | 2.28e−6 / 1.01e−5 | .0001887 / .0052815 | .0002675 / .0112929 | .0003468 / .0112792 | .005517 / .2574163 | 1.63e−8 / 5.58e−7 |
| 110×62 isolated to 111×63 nested | 1.81e−7 / 6.75e−6 | 3.46e−7 / 9.01e−6 | 4.59e−7 / 4.27e−6 | .0001507 / .0079441 | .0002620 / .0180547 | .0003109 / .0169656 | .004939 / .3653162 | 2.52e−9 / 1.41e−8 |
| Fixed for both norms | .0001 | .0001 | .001 | .001 | .001 | .001 | .01 | .01 |

The independent wedge calculation uses exactly the S1 mode and alternatives. The three 114×60, 114×62 and 106×65 states supply 50 common resolved ray samples. Across their admissible fits the shear amplitudes span −.0297573 to −.0109906, with largest measured amplitude uncertainty .00135534 from cross-mesh differences. Their consistent sign is a diagnostic. The spread across annulus and smooth-background choices exceeds the prospective stability criterion, which leaves a nonzero asymptotic amplitude unresolved. No sampled loss of positive covariance or principal margin follows from the observed trend.

The admitted isolated corner test moves the last interior x edge and first interior radial edge of the 114×62 mesh. It halves those two terminal widths at the same aspect ratio, with the same counts and original upstream and collar edges. The adjacent intervals grow and are included with all retained observations. This tests the corner response without the global collar change from the joint case. Its admission binds the prior process-tree RSS, both comparison budgets, the matched angular analysis and the evidence allowance.

The isolated native process hit the 8 GiB watchdog after 222.978 seconds, at a process-tree peak of 8,592,818,176 bytes. Retrying the same corner widths with 110×62 cells releases the predicted factor-memory pressure by merging four screened upstream x pairs. Its collar widths stay those of the 114×62 radial case. Retain the four extra strip observations to measure the actual effect of that reallocation.

The additional near-corner fit probes radii .0025, .003, .0035, .004, .0045, .005, .0055 and .006 at 12°, 18°, 24°, 30°, 36°, 42° and 48°. Hold out the .0035 and .005 radii and 18° and 36° rays. Use windows [.0025,.005] and [.0035,.006], the same angular mode, absolute pressure anchor and smooth-background controls. This limits the change of the auxiliary shear coefficient across the observations. The broad-annulus results and all original 690 spatial rows remain in the assessment.

The 110×62 retry stays within RSS at 8,131,780,608 bytes and times out after 1,065.577 process seconds. It retains a typed, unqualified same-policy zero-time candidate. The solved Newton base contracts from 4.37e−6 to 1.85e−8, followed by a long restarted inner solve. The continuation uses the same exact geometry and current equations with Krylov dimension 64. Its forecast adds the two dense bases' maximum incremental storage to the last measured peak, and reserves a separate memory margin. The continuation rebuilds lifts and histories and qualifies only after the original residual tolerance is met.

The continued 110×62 state converges to 1.03e−11 in 425.651 process seconds with 8,078,295,040 bytes peak RSS. Its pressure-tested current equations, all physical envelopes and inlet representation pass. All 92 supplementary observations in the four x strips pass. The 16 failed rows against 114×62 are confined to the corner region. The upper r=.10 collar owner's D is .0312112 versus .0310404 on the radial comparison, with Ryy difference .0001425 in tau0 units. The joint-to-isolated comparison again finds the distributed collar differences. This isolates a distinct radial collar issue from the effect of the smaller corner cells.

The common broad fit on 50 points retains negative shear coefficients spanning −.0309817 to −.0109416, with maximum observed mesh uncertainty .00337365 and sensitivity to smooth background. The shorter annuli retain 43 common points and limit the largest shear-coefficient change to .716 of its frozen value. They leave independent owner holdouts unresolved and have uncertainty as large as .0341833. The next controlled level bisects only the two terminal intervals of this converged isolated geometry and keeps every existing edge. Its predicted footprint, both angular controls, full common-owner comparison, independent current audit and evidence allowance are reserved together.

The native process on the nested 111×63 mesh reaches its timed limit after 1,015.919 process seconds with 8,400,969,728 bytes peak RSS. The accepted norm contracts from 6.03e−6 to 4.84e−8 and an unqualified zero-time iterate is sealed. A bounded continuation starts from that iterate on the same geometry and retains the Krylov dimension of 64, using the measured memory headroom and the original acceptance norm.

The first nested continuation takes 770.353 seconds and peaks at 8,424,374,272 bytes. It reduces the accepted residual by a further factor of 25.4 to 1.91e−9. One final same-geometry continuation is admitted from that measured reduction, with the same basis, zero time, 670-second complete-process cap and the full independent comparison and evidence allowance.

The final continuation converges to 5.95e−11 in 267.504 seconds and peaks at 8,374,493,184 bytes. The independent mesh reader initially split 126 coincident face halves into different decimal keys at x≈14.9993896484375. Adjacent native copies differ by 3.55e−15. The repaired reader groups edge keys within 1e−12 while retaining every original owner coordinate for polynomial and field evaluation. It verifies 6,993 elements and 13,812 interior faces. The repaired comparator reproduces all 690 baseline rows with the original 28 failures. The charged audit retries include an invocation that resolved the old evaluator from its script directory. The final run explicitly imports the paired-topology evaluator.

## Final nested comparison and measured gates

The last geometry halves the previous terminal widths to .0006103515625 and .00003125 and keeps the corner aspect ratio. The same 16 rows miss on the latest two corner comparisons. They concern Pi, traction, stress and gradients. The last MAX entries for all four occur at x=15, r=0 on the original one-sided owners. Twelve of the previously flagged R3 rows are within budget on the last pair, and 16 still miss. The released evidence records all 28 original rows on every available comparison, including the earlier 180-row collar sensitivity. [Original-row dispositions, per-case and per-field tables](wj-h27-linux/task2rs-s2/summary.json)

The nested vertex D grows to .453255 from .350710. The independent pressure-tested continuity dual norm is 1.20e−14. The strong domain RMS is 1.94913e−4, with an interior Gauss peak .282316. The upper r=.10 collar D changes from .0312112 to .0311781 and its Ryy difference stays within 2.72e−5 in tau0 units. The corner-only response reproduces the earlier localized shape. The collar split in the independent joint case remains evidence of a distinct radial-resolution issue.

The final common broad wedge fit uses 50 points and gives anchored shear amplitudes −.0311361 to −.0109323, with cross-mesh uncertainty up to .00337384. The spread across regular backgrounds and windows remains large. The short-annulus fit uses 43 common points, with uncertainty up to .0341966. Some holdout partitions have zero unseen original owners. Neither fit supplies the predeclared independent amplitude evidence. The bounded wedge result gives a permitted local exponent .577445. A nonzero coupled coefficient and any eventual loss of covariance remain unresolved. The diagnostic W=0 outlet retains the frozen normal traction and replaces the vertical natural condition. Its frozen wedge has no real root in (0,1). The diagnostic is a local calculation. Every coupled solve in S2 uses the prescribed original total traction.

| Finest 111×63 physical or projection measure | RMS | MAX | Gate |
| --- | ---: | ---: | --- |
| Inner U departure | .00173127 | .00259542 | PASS |
| Inner W departure | 2.16824e−5 | .000254260 | PASS |
| Inner k/Kin departure | .0196461 | .04029853 | PASS |
| Local outer epsilon departure | .00135583 | .01277464 | PASS |
| Ground shear against tau0 | .01217729 | .04223895 | PASS |
| Ground shear against Tref in tau0 units | .01010803 | .06422353 | Diagnostic |
| Inlet U projection | 1.29486e−6 | 3.18171e−6 | PASS |
| Inlet k/Kin projection | 5.96875e−8 | 2.47974e−6 | PASS |
| Inlet gradient projection | 2.23986e−5 | .000621352 | PASS |

The physical evidence covers the whole, inner and outer regions and the 29 stations. The reference Tref offset is 2.19845847% relative to tau0. The final independent residuals for momentum, continuity, k and transported epsilon are 2.75e−10, 2.68e−13, 3.98e−10 and 3.67e−10. Weak identity error is 5.14e−14. The native sampled minimum covariance/k and principal/mu are .2809062 and .00406887. Shared mass and k advection at H cancel. Independent owner-gradient k diffusion jump RMS/MAX is 5.09e−7/1.17e−6 and total-traction jump RMS/MAX is 5.44e−7/1.17e−4. The steady configuration uses Q3/Q2, Gauss 6×8 per governing panel, 1e−12 linear tolerance and 1e−10 coupled tolerance. Mesh-independent qualification, higher native quadrature and tighter algebraic controls on a selected final recipe remain open.

## Admission and handoff

The measured 111×63 solver attempts plus the final valid independent analysis take 2,301.032 wall seconds before layout preparation, topology-repair retries and publication. The same-parent nonprotected ledger had 2,279.477 seconds before S2 publication. A further nested 112×64 grid scales the largest measured 111×63 peak to about 8.64 GB before contingency, above the 8 GiB cap. The next refinement with the same resources and tensor-grid layout has no measured complete admission. The two measured issues need independent resolution in the collar and at the vertex, including multiple resolved owner blocks in lower annuli. A prospective lower-memory or local-refinement path would require its affected native spaces, operators, projection and current-equation checks. [Full measured limits and receipts](wj-h27-linux/task2rs-s2/summary.json)

Reviewed native source `061c48fe6805b1b1bb66a8effffc3bfc28c5f86777c5e669f52afd9f8c0d861f`, application `7a94420aaa31f41d9d8b6754d20c16f11d1398e894f243eaa6beaa3af7924819` and table `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1` bind all the tested states. This package changes the independent evaluator's face matching, the S0 metadata label and the fixture and evidence tooling. The sealed raw cases, checks and hashes are in the original Task 2R ext4 root on `jay@100.88.70.60`. The bundle SHA-256 `bd6f39fe7249c76dc3f4875b66bfe63aaa356a6eb16554a0daaac58d1150506c` and all 276 member identities verified after transfer.

All solves use one rank/process and six single-thread library limits. The sampled shared-worker peak is 64. The settled total is 37,329.075597 wall seconds, including 6,589.628524 charged during 2R-S. Remaining capacity is 2,270.924403 nonprotected seconds and 1,219.215123 staged seconds, with 3,600 protected and zero active reservations. The final storage readings are 1,401,743,582 science bytes, 138,442,549 source/build/runtime bytes and 3,578,158,313,472 free bytes. The recorded retries and admitted allocations share the original 43,200-second ceiling.

| Checkpoint or launch gate | Status |
| --- | --- |
| S0 contract and identity | PASS |
| S1 replay and causal diagnosis | PASS |
| S2 original-problem investigation | INCONCLUSIVE |
| S3 final recipe controls | NOT_RUN |
| S4 physical three-dt, cold restart, stock comparison and cost | NOT_RUN |

**DIAGNOSIS=INCONCLUSIVE. ORIGINAL_PROBLEM_SPATIAL_QUALIFICATION=NOT_MET. RECALIBRATED_NATIVE_REFERENCE=NOT_MET. SHORT_SCREEN=NOT_RUN. COMPLETE_COST=NOT_RUN. TASK2_COMPLETE=NO. TASK3_READY=NO.**

[Settled account](wj-h27-linux/task2rs-s2/accounting.json) · [Sealed raw bundle](wj-h27-linux/task2rs-s2/bundle.json) · [Task 3 launch hold](wj-h27-linux/task2rs-s2/task3-launch-hold.json)
