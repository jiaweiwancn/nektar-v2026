# Steady flat WJ–AWF qualification and outlet decision

**Campaign CLOSED. Select Bf at x=15 with `ZeroTransverseV1`.** Bf is the least expensive fully qualified recipe tested. Its 111×92 mesh uses Q3/Q2, native 6×8 quadrature, a 1e−12 linear stop and a 1e−10 coupled stop. The matched x=17.5 solve confirms that x=15 is the shortest tested outlet that meets the frozen physical criteria over [−2.5,8]. Both cases pass the independent current, physical, boundary, weak balance, inlet projection and native admissibility gates. `FullTractionV1` remains the solver default. [CP2 gate matrix and exact recipe](wj-h27-task2-outlet-campaign-cp2.md) · [Matched extent decision](wj-h27-linux/task2-outlet-campaign/final/outlet-campaign-be1-assessment.json)

## Measured outlet effect

Be1 retains the 92 radial cells and all x edges through 14.375, then translates the terminal grading to 17.5 on a 118×92 grid. The rebuilt runtime reproduces the x=15 coefficients, native and independent residuals, and all 690 signed rows exactly before the extension. Be1 keeps x=15 as an interior station and measures its new outlet at 17.5. The signed differences below are Be1 minus Bf. Each physical RMS and MAX passes after adding the sum of the measured direct within-B indicators. [Exact x=15 parity](wj-h27-task2-outlet-campaign-cp3-parity.md) · [Both extents, owners and outlet profiles](wj-h27-linux/task2-outlet-campaign/final/outlet-campaign-be1-assessment.json)

| Field | RMS | Signed MAX | Within-B indicator RMS / MAX | Physical limit RMS / MAX |
| --- | ---: | ---: | ---: | ---: |
| U | 6.368e−7 | +8.376e−6 | 1.751e−8 / 1.752e−6 | .01 / .02 |
| W | 6.493e−7 | −1.984e−6 | 3.511e−8 / 3.378e−6 | .01 / .02 |
| k / Kin | 3.074e−6 | +1.292e−4 | 1.573e−7 / 1.732e−5 | .05 / .10 |
| Outer epsilon / (A/r) | 1.821e−6 | +2.121e−5 | 5.862e−8 / 1.777e−6 | .10 / .20 |
| Ground shear / tau0 | 1.605e−5 | +3.783e−5 | 1.771e−6 / 1.427e−5 | .05 / .10 |

The W maximum lies inside the measured numerical indicator, with its effect-plus-indicator below 5.37e−6. The largest effect-plus-indicator as a fraction of any physical limit is 0.147%, for the k maximum. All original station rows and the fixed-window rows are preserved. The signed inventory records 122 exceedances at the mesh comparison thresholds, including outlet-local peaks and three fixed-window Pi rows. The physical extent decision applies the five limits above. [Full signed rows and failing component owners](wj-h27-linux/task2-outlet-campaign/final/outlet-qual-Be1-Bf.json)

| Additional quantity in tau0 units | RMS | Signed MAX and location | Within-B indicator RMS / MAX |
| --- | ---: | --- | ---: |
| R | 1.723e−5 | Rxx +2.918e−4 at x=7.99936, r=.0166111 | 3.793e−6 / 5.119e−4 |
| Modified pressure Pi | .003225 | +.003621 at x=8, r=19.9125 | 7.473e−6 / .0006220 |
| Shear traction Vxr | 1.663e−5 | +3.783e−5 at x=8, r=0 | 2.411e−6 / 2.578e−4 |

The pressure shift is resolved in both norms, with a signed mean of +.003218 tau0 in the fixed window. The stress and shear-traction RMS effects exceed their indicators, while their maxima are within the respective indicators. The ground-only shear maximum exceeds its smaller ground indicator. [Signed means, both physical owners, and resolution record](wj-h27-linux/task2-outlet-campaign/final/outlet-campaign-be1-assessment.json)

## Departures from the inlet in the fixed window

| Field | Bf at x=15 RMS / MAX | Be1 at x=17.5 RMS / MAX |
| --- | ---: | ---: |
| U | .00026135 / .00227230 | .00026106 / .00227121 |
| W | .00003530 / .00006108 | .00003475 / .00006094 |
| k / Kin | .00176135 / .03222030 | .00175901 / .03220856 |
| Outer epsilon / (A/r) | .00131734 / .01219058 | .00131667 / .01216937 |
| Ground shear departure from tau0 | .01154628 / .01163992 | .01155956 / .01165618 |
| Ground shear departure from Tref, scaled by tau0 | .01043906 / .01064971 | .01042586 / .01064524 |

The separate Tref value is .002658181904845143 for both runs. The complete audit also records whole-domain departures and all retained stations. [Bf independent audit](wj-h27-linux/task2-outlet-qual/cp2/outlet-Bf-audit.json) · [Be1 independent audit](wj-h27-linux/task2-outlet-campaign/final/outlet-Be1-audit.json)

## Measured resource use

| Case | Native wall s | Native peak GiB | Independent audit s | Signed comparison s |
| --- | ---: | ---: | ---: | ---: |
| Bf | 1059.481 | 11.36 | 207.663 | 102.498 |
| Bfl | 1092.621 | 11.32 | 204.647 | 103.033 |
| Br3 | 1552.720 | 13.83 | 248.630 | 124.384 |
| Br4 | 1255.361 | 13.31 | 263.268 | 120.274 |
| Bx2 | 866.065 | 11.58 | 222.981 | 104.359 |
| Bx3 | 2763.131 | 14.49 | 255.183 | 122.941 |
| Bh1 | 685.491 | 11.61 | 223.602 | 102.241 |
| Bh2 | 691.522 | 11.93 | 217.583 | 115.666 |
| Bq1 | 347.804 | 10.74 | 207.908 | 100.302 |
| Bt2 | 1315.889 | 11.39 | 207.957 | 104.080 |
| Bp1 | 100.881 | 3.96 | 208.722 | 105.429 |
| Be1 | 2248.592 | 12.10 | 230.771 | 107.886 |

Native totals include the successful typed Newton chains from their recorded seeds. The isolated Release build took 317.559 s with four jobs and 4.26 GiB peak. The matched extent assessment took 75.855 s. All native science runs used one MPI rank with the six recorded thread settings set to one. The wider radial tests hit their watchdogs, with the second stop sampled at 21.06 GiB against its 20 GiB allowance. The combined Bx1 streamwise test missed independent momentum and its tighter continuation timed out. The passing staged witnesses supply the final spatial proof. The final Be1 recovery losslessly compacted three earlier failure reports and preserved their exact decoded hashes. [Radial resource stop](wj-h27-task2-outlet-campaign-br2.md) · [Bx1 and Bt1 diagnosis](wj-h27-task2-outlet-campaign-bt1-time-stop.md) · [Storage and recovery receipts](wj-h27-task2-outlet-campaign-be1-iterate.md) · [Final measured commands](wj-h27-linux/task2-outlet-campaign/final/operations.json)

The settled cumulative ledger charges 84,896.795338 s of the 108,000 s ceiling, including 25,555.465955 s since the opening settlement. It carries zero outstanding reservations, 19,503.204662 s of nonprotected headroom, and the full 3,600 s protected reserve. Science storage settles at 8,463,861,449 bytes and the build inventory at 394,931,219 bytes. The final ext4 evidence archive SHA-256 is `e5aba1a3fbf39d593350b6f9fbd277a59c04e04058a3f08a8acf44724b9b24eb`. [Settled ledger](wj-h27-linux/task2-outlet-campaign/final/accounting.json) · [Publication receipt, archive and exact operations](wj-h27-linux/task2-outlet-campaign/final/operations.json)

## Advancing-work handoff

Use the verified Bf geometry, exact-inner table and Bf typed zero-time packet on `jay@100.88.70.60`. The explicit outlet is `ZeroTransverseV1`. The pinned Bf mesh, native binary, coefficient packet, source identities, reference table, session, options and commands are in the machine-readable handoff. The exact x=15 replay also binds the diagnostic runtime built from `b5b6e9267bb389e6d334f495d5ed1314f9aea446`. The verified CP2 selection and CP3 parity checkpoints are `3e24a97cf29f04fecc5b44223b5a1caf25db4309` and `d132917162dbb00fc81e7cc4e300eaaa63d0bfde`. [Exact input hashes and environment](wj-h27-linux/task2-outlet-campaign/final/outlet-campaign-final-handoff.json)

Physical time advanced is zero. Physical timestep accuracy, cold restart, matched stock timing and full advancing trajectory cost remain downstream. The production target is `1.25*C_total <= 86400` wall seconds per trajectory.

**TASK2_COMPLETE=NO. TASK3_READY=NO.**
