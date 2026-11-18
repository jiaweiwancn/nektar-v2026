# Task 2R R3 physical and spatial evidence

**R3=NOT_MET.** The unforced flat references pass the physical envelopes, independently assembled current equations, and native guards. The common-owner mesh comparisons exceed several prescribed RMS/MAX budgets. The developed reference remains unqualified for the short timestep, restart, and cost screen.

The unforced flat solves use the qualified `ExactInnerFiniteNuV1` reference, the retained `ScalarCdV2` transport, and only the declared sampled AWF source. The same physical map applies to every remesh. The reference length is 40 mm and the modeled wall is 3.5 mm above the original floor. `H_code=0.15` gives a 6.0 mm matching height and a matching-plane elevation of 9.5 mm. The steady solves use the zero-time-mass equations and advance zero simulated time.

The first 28×112 solve reached 3.65e−9 before its wall deadline. A second solve with the typed same-policy candidate reached 1.20e−10. The next solve met the unchanged 1e−10 target at 3.14e−11. Each continuation reconstructs the boundary lifts, current data, and fresh histories from the saved candidate. The independent baseline residuals are 1.42e−12 for momentum, 1.47e−15 for continuity, 3.15e−11 for k, and 6.73e−12 for outer epsilon. The weak accounting identities and boundary conditions pass. Increasing quadrature from 8×12 to 12×16 returns the same coefficients. Tightening the linear and nonlinear tolerances to 3e−13 and 3e−12 changes k by at most 1.50e−8 in Kin units over the original common-owner observations.

| Physical mesh | Inner U RMS / MAX | Inner W RMS / MAX | Inner k/Kin RMS / MAX | Outer epsilon RMS / MAX | Ground shear in tau0 units RMS / MAX |
| --- | ---: | ---: | ---: | ---: | ---: |
| 28×112 | .0017443 / .0032225 | .0000402 / .0002725 | .0197938 / .0554613 | .0013722 / .0187376 | .0122324 / .0219813 |
| 39×168 | .0017432 / .0028587 | .0000261 / .0001924 | .0197788 / .0402940 | .0013582 / .0146280 | .0121296 / .0219842 |
| 56×112 | .0017418 / .0028428 | .0000248 / .0001913 | .0197581 / .0403578 | .0013581 / .0134499 | .0121291 / .0219813 |
| 106×60 | .0017393 / .0027898 | .0000230 / .0002059 | .0197272 / .0402818 | .0013582 / .0128886 | .0121151 / .0219597 |
| 109×60 | .0017338 / .0026569 | .0000219 / .0002190 | .0196720 / .0402993 | .0013566 / .0127933 | .0121410 / .0219597 |
| 114×60 | .0017313 / .0025954 | .0000217 / .0002522 | .0196460 / .0402985 | .0013558 / .0127746 | .0121542 / .0219597 |
| Limits | .01 / .02 | .01 / .02 | .05 / .10 | .10 / .20 | .05 / .10 |

Physical envelopes pass in the whole, inner and outer regions and at all 29 fixed stations on these meshes. The acceptance comparator uses the frozen recalibrated inner columns and original outer inlet. The historical inner comparator and ground shear relative to the independent equilibrium shear remain diagnostic. The latter uses tau0 in the denominator and preserves the reference offset of 2.19845847%.

At 114×60 the worst station RMS/MAX for inner k/Kin is .02419/.03985. The corresponding historical-inner diagnostic is .14235/.54390. The ground shear RMS/MAX relative to the new equilibrium shear is .01013/.03454 in tau0 units. Independent momentum, continuity, k, and epsilon residuals are 3.60e−11, 7.50e−14, 5.81e−11, and 5.45e−11. Native minimum covariance/k is .28091 and the sampled minimum principal coefficient divided by mu is .004069. The maximum weak accounting identity error is 4.23e−14.

The shared assembled interface mass and k advection cancel. The original-owner k diffusion jump RMS falls from 4.56e−6 on 28×112 to 5.05e−7 on 39×168 and 5.07e−7 on 106×60. Streamwise refinement alone leaves the radial owner flux RMS close to 4.54e−6 on 56×112. The 106×60 regrade retains the small matching-plane owner intervals while trimming smooth radial regions for streamwise resolution. Its measured inlet projection meets the U, k and velocity-gradient RMS/MAX gates. The full 106×60 solve peaks at 7.29 GiB RSS, and its 6×8 quadrature control independently preserves the equations and all compared field rows.

| Common-owner pair | U RMS / MAX | W RMS / MAX | k/Kin RMS / MAX | Stress RMS / MAX | Pi RMS / MAX | Traction RMS / MAX | Gradient RMS / MAX | Local outer epsilon RMS / MAX |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 28×112 → 39×168 | .0006774 / .0016571 | .0001449 / .0003822 | .018745 / .034603 | .042832 / .078129 | .008664 / .011157 | .003296 / .012781 | .015836 / .171585 | .0002924 / .0085831 |
| 28×112 → 56×112 | .0007289 / .0014849 | .0001395 / .0003270 | .017606 / .030659 | .040251 / .069223 | .008834 / .010990 | .003683 / .011411 | .016642 / .174389 | .0003309 / .0074459 |
| 39×168 ↔ 56×112 | .0002034 / .0007142 | .0000309 / .0001624 | .002375 / .012466 | .005577 / .028147 | .000927 / .004041 | .0016439 / .0066479 | .005931 / .103561 | .0000606 / .0022952 |
| 56×112 ↔ 106×60 | .0001659 / .0004058 | .0000412 / .0002240 | .002651 / .008623 | .006206 / .019469 | .003180 / .007720 | .001403 / .005041 | .005524 / .086389 | .0000586 / .0010357 |
| 106×60 ↔ 109×60 | .0001064 / .0002062 | .0000438 / .0001657 | .001199 / .004350 | .003145 / .009823 | .007541 / .014276 | .003295 / .006296 | .008992 / .048633 | .0000182 / .0001759 |
| 109×60 → 114×60 | .0000527 / .0000816 | .0000213 / .0001180 | .000565 / .002254 | .005567 / .020317 | .002795 / .034085 | .006195 / .021812 | .035880 / .469487 | .0000087 / .0000601 |
| Limits for both norms | .0001 | .0001 | .001 | .001 | .001 | .001 | .01 | .01 |

Common-owner differences localize first near x=7–8.75, then near x=13.125–14.76 after additional streamwise resolution. The largest remaining 56×112 to 106×60 k and radial velocity-gradient differences occur at x=14.375 and r=0. The 122×60 outlet-band bisection triggered the 8 GiB watchdog with both the full-gradient inverse and the positive frozen inverse with a shorter Krylov basis. These runs settle at measured peaks of 8,593,043,456 and 8,634,007,552 bytes. Reallocating cells from the measured quiet upstream bands retains the entire 56-cell parent x grid and puts more of the admissible resolution downstream.

All entries in the numerical table are maxima over the corresponding RMS or MAX rows, including the fixed station and owner rows. The 109×60 to 114×60 pair uses identical radial edges and 6×8 native quadrature. Five nested bisections shorten the last x cell from .078125 to .00244140625. Every U and outer epsilon row meets its budget across 690 checks. W, k, stress, Pi, traction, and gradient leave 28 rows open. The worst pressure, shear traction, and radial U-gradient differences occur on the ground–outlet corner owner. Reallocating the outlet x cells reduced the k MAX, with its remaining largest difference at x=15 and r=.0002625. The stress MAX occurs at x=15 and r=.10. The k diffusion jump RMS stays at 5.09e−7 at H. The traction jump MAX at H reaches 1.08e−4 after terminal refinement. The shared mass and k advection cancel on all original H faces.

## Corner evidence and operation limits

The separate one-sided trace check gives zero U, W, and streamwise U derivative at x=15 and r=0, with radial W derivatives .04001, .07774, and .25547 on the 106×60, 109×60, and 114×60 meshes. The local modified pressure in tau0 units moves from −.82447 to −.83875 to −.87283. These are pointwise results from the original corner owners. The independently assembled continuity residual remains within its fixed tolerance.

Below .4H, the frozen lower-wall v2 operator has normal viscosity nu and shear coefficient mu=nu+nuA. Eliminating pressure with a streamfunction gives the fourth-order principal operator `mu*(d_x^4+d_r^4)+(4*nu-2*mu)*d_x^2*d_r^2`. The no-slip and homogeneous traction wedge pencil at the local sampled K admits a real velocity exponent near .577445, with relative boundary null residual below 1.7e−16 across the three corner witnesses. The associated gradient and pressure exponent is −.422555. Its amplitude in the variable-coefficient, sampled-force coupled solution remains unverified. The one-sided trace and pencil witnesses support testing simultaneous local x and radial refinement under the same exterior data.

The 114×62 trial preserves the 114×60 terminal x grid and splits the first radial interval twice, from .0005 to .000125. It uses the accepted 114×60 state as a zero-time seed. The 8 GiB RSS forecast uses the measured 8,005,373,952-byte peak of the 114×60 solve. The initial trial preflight recorded 69 shared workers, leaving insufficient headroom for science admission. The runner deferred the trial before fixture generation. The ext4 checkpoint records the deferred operation, the nested radial recipe, and the digest-bound admission receipt.

The 114×60 native process takes 376.658 wall seconds. Construction and initialization take 23.668 and 17.647 seconds. The 322.656-second Krylov scope includes 96.534 seconds for residual callbacks, 157.705 seconds for preconditioner preparations, 19.728 seconds for preconditioner applications, 45.357 seconds for nonlinear validation, and 1.755 seconds for trial validity. The final field/history export takes .320 seconds. The remaining process time includes outer solver work, report serialization and launcher overhead. The independent current analysis takes 124.254 seconds and the mesh analysis takes 94.252 seconds. These are steady-process measurements. Advancing trajectory costs and a cost target remain unmeasured.

All runs use one rank/process with the six numerical thread limits at one on the ext4 server at `jay@100.88.70.60`. The runner keeps the 72-worker watchdog. An earlier worker census stopped a projection at 74, including parked MPI event listeners. Inspection showed sleeping launcher helpers. The subsequent census excludes parked MPI event listeners while counting every main rank and runnable helper. The retried projection passes. The serial pre-radial publication reached a sampled peak of 72. Both 8 GiB storage caps, the 8/16 GiB science/build memory caps, and the 15 GiB free reserve remain part of admissibility.

The sealed pre-radial checkpoint charges 30,739.447 seconds total, with 8,860.553 nonprotected seconds and 2,643.844 staged seconds remaining. The protected 3,600 seconds remain reserved for evidence. The staged additions of 1,500, 1,750, 1,000, and 850 seconds draw from the same parent. There are no active reservations at publication.

## Disposition

| Gate | Status | Evidence |
| --- | --- | --- |
| R0 prospective contract and retained geometry | PASS | [Contract](wj-h27-task2r-contract.md) |
| R1 calibration and independent cutoff family | PASS | [R1 evidence](wj-h27-task2r-r1.md) |
| R2 native coupling and inlet representation | PASS | [R2 evidence](wj-h27-task2r-r2.md) |
| R3 current equations, guards and physical envelopes | PASS | [Current and physical evidence](wj-h27-linux/task2r-r3/summary.json) |
| R3 independent quadrature and algebraic comparisons | PASS | Measured 28×112 controls and the 106×60 quadrature check |
| R3 common-owner spatial qualification | NOT_MET | [Common-owner comparison and exact failed rows](wj-h27-linux/task2r-r3/summary.json) |
| R3 three-dt accuracy, stock comparison, restart and cost screen | NOT_RUN | Depends on a spatially qualified developed reference |
| R4 full Task 3 entry matrix | NOT_RUN | [Launch hold](wj-h27-linux/task2r-r3/task3-launch-hold.json) |

**RECALIBRATED_NATIVE_REFERENCE=NOT_MET. SHORT_SCREEN=NOT_RUN. TASK2_COMPLETE=NO. TASK3_READY=NO.**

The fixed scientific identifiers bind source `061c48fe6805b1b1bb66a8effffc3bfc28c5f86777c5e669f52afd9f8c0d861f`, runtime SHA-256 `7a94420aaa31f41d9d8b6754d20c16f11d1398e894f243eaa6beaa3af7924819`, and reference table SHA-256 `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`. The pre-radial Linux archive hashes to `f4b005465b2ec339fa8af810c6c4b438564cfefefd92608db4686023c1ddcffd`. Its manifest binds the scripts, immutable fixtures, complete comparison records, receipts, and logs to their ext4 source paths.

[Settled account](wj-h27-linux/task2r-r3/accounting.json) · [Bundle identity](wj-h27-linux/task2r-r3/bundle.json) · [Corner trace](wj-h27-linux/task2r-r3/corner-witness.json) · [Frozen wedge pencil](wj-h27-linux/task2r-r3/frozen-pencil-witness.json)
