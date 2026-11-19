# Task 2 isolated collar control

M1c keeps all 112 M1 x edges, including the terminal width .0006103515625, and the .00003125 first radial interval. Its only new radial edges are .090625 and .10320438828451. They bisect [ .08125, .10 ] and [ .10, .10640877656902 ], respectively. This gives 111×65 cells. Both laws use the same M1c XML, the same zero-time Q3/Q2 discretization, 6×8 quadrature and the inherited equation and acceptance settings. The 106×65 historical joint mesh also changed x and the tip. [Frozen axes and geometry](wj-h27-linux/task2-outlet/baseline.json) · [Matched-mesh screen](wj-h27-task2-outlet-cp2.md)

The first A1c and B1c solves each accepted two Newton updates and reached 1.13e−9 and 1.14e−9 at their wall limits. Their typed continuations used the same M1c XML, reconstructed the current coupled data and converged to 3.65e−11 and 3.73e−11, respectively. All current, boundary, interface, physical and projection gates passed for both. The matching geometry SHA-256 is `72f5bb7b7ec02b62945a5bf6cc421754d46e47c3d6b45ce2989529ce68c65269`. [A1c audit](wj-h27-linux/task2-outlet/cp3/outlet-A1c-audit.json) · [B1c audit](wj-h27-linux/task2-outlet/cp3/outlet-B1c-audit.json)

| Policy | First solve and continuation charged wall seconds | Peak RSS across attempts | Peak shared workers |
| --- | ---: | ---: | ---: |
| A | 1,682.890 + 416.458 | 8.24 GiB | 64 |
| B | 1,680.404 + 411.286 | 8.10 GiB | 71 |

B1c was admitted after A1c using the measured costs and both policies' existing physical and numerical evidence. Independent assembly, common-owner comparison, raw station witnesses, publication, a continuation and reservation grace were included in its forecast. [Admission ledger and inputs](wj-h27-linux/task2-outlet/cp3/outlet-b1c-admission.json) · [Run receipts and archive](wj-h27-linux/task2-outlet/cp3/operations.json)

## Radial sensitivity under each outlet

Both original-to-collar differences miss 180 of the 690 rows. A's largest raw change is Ryy at the upper r=.10 outlet owner. At x=5, A's stress difference reaches .007223 and B's interior stress difference follows the same radial pattern. The two unchanged x grids and preserved corner widths isolate this influence to the collar refinement. [A1c−A1](wj-h27-linux/task2-outlet/cp3/outlet-a-collar.json) · [B1c−B1](wj-h27-linux/task2-outlet/cp3/outlet-b-collar.json)

| Full-domain scaled MAX | A1c−A1 | B1c−B1 | Limit |
| --- | ---: | ---: | ---: |
| U | 2.64e−5 | 1.48e−5 | .0001 |
| W | 2.06e−5 | 1.84e−5 | .0001 |
| k / Kin | 4.37e−5 | 2.85e−5 | .001 |
| R / tau0 | .022295 | .012257 | .001 |
| Pi / tau0 | .009157 | .009157 | .001 |
| Vxr / tau0 | .009890 | .008530 | .001 |
| Velocity gradient in code units | .020419 | .014687 | .01 |
| Outer epsilon / local reference | 2.19e−6 | 2.19e−6 | .01 |

At x=5 and r=.10, D on the upper owner falls from .003459 to .000403 under both laws. At x=15 the A upper owner falls from .031178 to .003716, and B falls from .010072 to .000959 on both owners. B's zero transverse outlet trace gives Wᵣ=0 there. The Ryy change at that upper outlet owner is −.022295 tau0 under A and −.007538 tau0 under B. The samples retain both raw owners, including k, dependent E, sampled K, production, full R and the original velocity gradients. [A1c raw profiles](wj-h27-linux/task2-outlet/cp3/outlet-A1c-stations.json.gz) · [B1c raw profiles](wj-h27-linux/task2-outlet/cp3/outlet-B1c-stations.json.gz)

## Boundary effect with the radial indicator

The full-height B−A estimates are stable when both laws move from M1 to M1c. Interior U, W and k/Kin maxima on M1c are 1.264e−5, 3.664e−6 and 1.926e−4. Pressure reaches .002903 tau0. The radial differences are much larger than the two-mesh tip indicators for peak stress, pressure and gradients. Using the sum of the A and B radial differences, the interior U, W, Pi and shear RMS estimates remain resolved. Their MAX estimates are unresolved with this indicator. Both norms of k and outer epsilon remain resolved. The gradient and stress flags are component-specific. [Matched M1c effect and radial indicators](wj-h27-linux/task2-outlet/cp3/outlet-collar-effect.json)

The same near-collar response at upstream locations identifies a separate spatial qualification obstacle even after B clears the two-mesh outlet screen. A uniform second bisection of the original collar patch under B isolates whether that error contracts. The next comparison can hold the outlet and the experimental-region x grid fixed.

**TASK2_COMPLETE=NO. TASK3_READY=NO.**
