# Bh2 and one-sided H-interface convergence

Bh2 uses the exact original H support and the fixed Bf streamwise edges. It quarters both original H-adjacent widths on a 111×98 grid. Two accepted updates take 691.522 s, converge to 1.12273e−11 and peak at 11.93 GiB against the 13 GiB allowance. The independent audit and direct comparison take 217.583 and 115.666 s. Bh2 passes momentum at 2.64287e−10, continuity at 1.69777e−14, k at 3.90719e−10 and outer epsilon at 3.61763e−10. The remaining independent gates, all 690 signed rows and every fixed-window row pass. [Measured commands and costs](wj-h27-linux/task2-outlet-campaign/cp2-bh2/operations.json) · [Independent audit](wj-h27-linux/task2-outlet-campaign/cp2-bh2/outlet-Bh2-audit.json)

The Bh2−Bf maximum Rrr and Pi changes are +9.30915e−5 and −9.69822e−5 tau0. They occur at x≈14.99511719, r=H on the reported one-sided owners, with zero RMS weight. The largest shear change is −3.77900e−5 tau0 near x=15. These whole-domain changes are small compared with the available margins in the Bx3 stress and pressure results. [Signed 690 rows, fixed window and both owners](wj-h27-linux/task2-outlet-campaign/cp2-bh2/outlet-qual-Bh2-Bf.json)

The common upward normal and frozen H sampling rule give the following maximum jumps.

| Grid | k diffusion / A | x traction / tau0 | r traction / tau0 |
| --- | ---: | ---: | ---: |
| Bf | 1.171e−6 | 2.987e−5 | 6.926e−5 |
| Bh1 | 2.056e−7 | 1.657e−5 | 3.764e−5 |
| Bh2 | 5.006e−8 | 7.083e−6 | 1.753e−5 |

The Bf-to-Bh2 direct differences pass on both sides. The largest direct MAX is 3.85125e−5 for lower normal traction. Eleven of the twelve successive RMS and MAX ratios are below .8, with ratios from .119 to .644. The upper normal-traction MAX rises from 7.58857e−6 on the first difference to 9.19478e−6 on the second. Both lie within the frozen 1e−5 resolution floor. Shared mass and k advection agree exactly on all three levels. The proof costs 13.613 s. [Frozen test, signed peak owners and measured contraction](wj-h27-linux/task2-outlet-campaign/cp2-bh2/outlet-campaign-H-convergence.json)

Bq1 now uses the exact Bf mesh and a typed Bf initial guess with native 8×12 quadrature and the original tolerances. The full forecast is 4,357 s, with 13.01 GiB expected peak and an allowance of 14 GiB. The numerical controls require inner-k RMS and MAX of at most 1e−4 on the full inventory and the original limits on all other fields. [Bq1 forecast](wj-h27-linux/task2-outlet-campaign/cp2-bh2/outlet-qual-bq1-forecast.json)

**Campaign OPEN. TASK2_COMPLETE=NO. TASK3_READY=NO.**
