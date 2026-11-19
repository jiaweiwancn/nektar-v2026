# Task 2 matched-mesh outlet comparison

The comparison uses archived A0 and A1 and the accepted B0 and B1 states on the same two native meshes. Each is a zero-time steady state with Q3 velocity and scalars, Q2 pressure, 6×8 native panel quadrature and the inherited linear and coupled tolerances. Both A state packets match their rebuilt native replays exactly. B0 and B1 use the opt-in zero transverse trace with the same frozen normal load. [Archive replay and action evidence](wj-h27-task2-outlet-cp1.md) · [Four-case matrix](wj-h27-linux/task2-outlet/cp2/outlet-matrix.json)

The direct archive-to-B0 guess failed admissibility when the outlet lift changed W across the terminal column. It made Wₓ=−.03831 near (x=14.999959, r=5.51834), with principal margin −.02038 in local viscosity units. The retry tapered just the starting W polynomial with a quintic across the last code unit. It passed the existing guard and converged in four accepted Newton steps with residual 5.14e−11. B1 started from the converged B0 fields on the nested mesh. Its first bounded solve reached 6.91e−10. The registered best iterate transferred as a typed steady continuation and completed at 1.21e−11. Both imports rebuilt their current lifts and histories. The charged wall times of the successful executions were 594.452 s and 354.044 s, respectively. The two unsuccessful attempts took 42.158 s and 1,328.959 s. [Inputs, receipts and Linux archive](wj-h27-linux/task2-outlet/cp2/operations.json)

## Current equations, physical bounds and tracings

The audits bind each native result to the corresponding executable, table, signed owner inventory, mesh, typed seed and settled run. They independently assemble the complete current balances and evaluate the physical bounds, inlet projection and shared owner fluxes. The pressure-test Riesz norm divides by the square root of area. The companion local trace output retains the earlier area-normalized statistic under its original key so that the S2 records remain reproducible. The corrected RMS key and the independent current residual use the inherited Riesz normalization. [A0 audit](wj-h27-linux/task2-outlet/cp2/outlet-A0-audit.json) · [A1 audit](wj-h27-linux/task2-outlet/cp2/outlet-A1-audit.json) · [B0 audit](wj-h27-linux/task2-outlet/cp2/outlet-B0-audit.json) · [B1 audit](wj-h27-linux/task2-outlet/cp2/outlet-B1-audit.json)

All four states meet the current equations, weak balance accounting, native and independent trace checks, shared interface mass and advective k flux, physical envelopes and inlet projection limits. The raw one-sided gradients and their diffusive flux differences remain in the evidence. B's weak pressure-tested continuity remains below 3e−13. Ground shear for B1 is .01209 RMS and .02198 MAX in tau0-normalized units against tau0. The audit also reports its deviation from Tref(Kin).

| Case | Ground–outlet vertex D | Whole-domain strong D RMS | Near annulus D RMS | Broad annulus D RMS | Outlet trace D sample RMS | Ground trace D sample RMS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A0 | .350710 | 1.94915e−4 | 1.456e−3 | 1.188e−3 | 2.264e−3 | 4.553e−3 |
| A1 | .453255 | 1.94913e−4 | 1.437e−3 | 1.188e−3 | 2.156e−3 | 1.752e−3 |
| B0 | 0 | 1.94864e−4 | 5.838e−5 | 1.105e−3 | 1.042e−3 | 8.479e−5 |
| B1 | 0 | 1.94864e−4 | 5.863e−5 | 1.105e−3 | 1.111e−3 | 8.460e−5 |

The ground no-slip trace gives Uₓ=0 and the B outlet trace gives Wᵣ=0 at their shared vertex. The near annulus spans 11 owners on each mesh and its D RMS is about 25 times smaller under B. The broad annulus spans 20 owners and its change is about 7%. Whole-domain strong D changes by about 0.025% on M1.

The station witnesses retain the signed owner values for U, W, k, dependent E or transported epsilon, sampled K, production, Pi, R, shear and the raw gradient. They include face endpoints and both sides at material and discretization interfaces. Outlet work records distinguish the prescribed normal action, the discarded A transverse load, the raw Cauchy traction and B's signed reaction. Both B runs carry outlet mass flux 19.0487942. The B1 integrated raw transverse traction is .05118909, while its outward reaction on a unit transverse displacement is .04631204 from outlet-only rows plus .00487706 from the two shared corner rows. Transverse prescribed work and vertical reaction work both vanish in B. [A0 raw profiles](wj-h27-linux/task2-outlet/cp2/outlet-A0-stations.json.gz) · [A1 raw profiles](wj-h27-linux/task2-outlet/cp2/outlet-A1-stations.json.gz) · [B0 raw profiles](wj-h27-linux/task2-outlet/cp2/outlet-B0-stations.json.gz) · [B1 raw profiles](wj-h27-linux/task2-outlet/cp2/outlet-B1-stations.json.gz)

## Spatial sensitivity and effect of changing the outlet

The comparison uses common physical owner cells, panel quadrature and the inherited 690 rows. The spatial gates apply within each policy. At all 29 stations and in the fixed interval [−2.5,8], the files give both within-policy differences beside each same-mesh boundary difference, with signed RMS or maximum effects flagged using the prospective sum of A and B spatial indicators.

**B meets the two-mesh screen.** The A pair reproduces the 16 archived failures. The B pair passes all 690 rows. A's four outlet maxima occur at the ground–outlet vertex, on owners 109 and 110. B's maxima move into the collar at x=15. These are two-mesh sensitivities.

| Outlet station maximum | A1−A0 | B1−B0 | B peak radius and component |
| --- | ---: | ---: | --- |
| Stress / tau0 | .00794412 | .000150600 | r=.06000, Rxr |
| Pi / tau0 | .01805471 | .0000418068 | r=.06875, Pi |
| Shear traction / tau0 | .01696561 | .000158529 | r=.06000, Vxr |
| Velocity gradient in code units | .36531620 | .00148177 | r=.066610, Uₓ |

At x=15 the M1 boundary change reaches .04333 for stress, .06697 for Pi, .05035 for shear and 1.08663 for gradients. Within the experimental window the changes are much smaller. The main U, W, k, Pi and Vxr effects exceed both spatial indicators combined at both levels. The signed Wr RMS and Rxy and Ryr rows carry the unresolved label.

| Full-height interior observable | B1−A1 RMS / MAX | A1−A0 MAX | B1−B0 MAX |
| --- | ---: | ---: | ---: |
| U | 1.30e−6 / 1.26e−5 | 6.08e−10 | 6.83e−9 |
| W | 1.39e−6 / 3.66e−6 | 1.39e−10 | 1.81e−9 |
| k / Kin | 5.19e−6 / 1.93e−4 | 8.70e−9 | 1.08e−7 |
| R / tau0 | 4.19e−5 / 4.35e−4 | 2.23e−8 | 2.45e−7 |
| Pi / tau0 | .001556 / .002903 | 1.29e−7 | 1.65e−6 |
| Vxr / tau0 | 4.19e−5 / 9.97e−5 | 4.68e−8 | 7.46e−8 |
| Velocity gradient in code units | 7.68e−6 / 9.98e−4 | 9.94e−7 | 1.46e−6 |
| Outer epsilon / local reference | 3.62e−6 / 3.67e−5 | 1.52e−9 | 1.96e−8 |

At the virtual crest x=0 and the retained 7 mm elevation from the original floor, M1 changes by −2.47e−6 in U, +3.25e−8 in W, −1.50e−5 in k/Kin and +.001944 in Pi/tau0. The two x owners agree on these continuous values.

The independent collar control remains necessary. The M0–M1 pair changes the terminal x and wall r widths, so the passing B screen leaves global radial collar sensitivity open.

**TASK2_COMPLETE=NO. TASK3_READY=NO.**
