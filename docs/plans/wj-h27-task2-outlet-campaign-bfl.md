# Bfl lower blend repeat

Bfl adds r=.0621875 and .0665625 to the Bf radial grid. It holds all 112 x edges and the upper companion interval above r=.06875. The 111×94 experiment tests whether the newly exposed lower band accounts for the outlet stress and shear change at r=.06 and the axial derivative at r=.06875. It retains the lower owner of the recorded near-inlet Wr peak. [Verified campaign admission plan](wj-h27-task2-outlet-campaign-cp0.md)

The first native segment uses one rank, a named 13 GiB cap and one Newton update. It writes a typed accepted iterate at residual 8.316e−8 and charges 299.928 s. Its continuation uses that zero-time packet with the exact mesh and reference and fresh lifts and histories. The physical model, Q3/Q2 orders, 6×8 quadrature and 1e−12 linear and 1e−10 coupled tolerances are held.

At the outlet and r=.06, χ=0 gives $V_{xr}=(\nu+\nu_A)(U_r+W_x)$. The signed difference uses the mean of the two effective viscosities for each gradient change and the mean of the two shear gradients for the viscosity change. The independent Q3 verification expands all four Jacobi-normalized Modified_A modes in unit coordinates and applies each element's signed coefficient map at 45 decimal digits. Both one-sided values and the unweighted face endpoints enter the sample record.

## Measured gates

The continuation converges at 3.585e−12 in 792.693 s after three updates. The peak native process-tree RSS is 11.32 GiB. The independent audit charges 204.647 s, the comparison 103.033 s and the derivative probe 8.041 s. The independent momentum, continuity, k and outer-epsilon residuals are 2.645e−10, 4.066e−14, 3.907e−10 and 3.613e−10. The weak accounting, boundary, shared interface, physical, inlet projection and native guard checks all pass. [Audit and state hashes](wj-h27-linux/task2-outlet-campaign/cp1-bfl/outlet-Bfl-audit.json) · [Native and analysis receipts](wj-h27-linux/task2-outlet-campaign/cp1-bfl/operations.json)

All 690 Bfl−Bf rows and the fixed experimental region pass. The stress and shear maxima on the upper r=.06 outlet owner are −.000551678 and +.000580718 tau0. That endpoint has zero RMS weight and is retained in MAX. The Bf radial owner is [.06,.064375] and the Bfl owner [.06,.0621875], with the same final x column. The whole-domain pressure maximum is .000160731 tau0. The largest gradient change is +.00293127 in Wr at x=−9.93967077, r=.06438252. The fixed-region stress, pressure, shear and gradient maxima are .00003640, .00013836, .00001836 and .00250036. [All signed rows and both physical owners](wj-h27-linux/task2-outlet-campaign/cp1-bfl/outlet-qual-Bfl-Bf.json)

| r=.06 at x=15 | Ur term / tau0 | Wx term / tau0 | Viscosity term / tau0 | Total signed shear / tau0 |
| --- | ---: | ---: | ---: | ---: |
| Lower owner | −.0002643830 | +.0005458250 | +8.089e−10 | +.0002814428 |
| Upper owner | +.0000348923 | +.0005458250 | +8.085e−10 | +.0005807181 |

Wx supplies about 94% of the upper owner's change. The independently reconstructed gradients agree with the raw accepted gradients to 1.77e−13 on Bf and 2.63e−13 on Bfl over the tracked and neighboring samples. At r=.06875 the outlet Ux changes by +.00139313. At the archived inlet probe x=−9.93967077, r=.05996564, Wr changes by −.00153039 in the unchanged lower owner. These results support the lower-band resolution hypothesis. The change in effective viscosity contributes less than 1e−9 tau0 at either face owner. [Independent coefficient and shear proof](wj-h27-linux/task2-outlet-campaign/cp1-bfl/outlet-qual-Bfl-owner-probe.json)

## Next direct challenge

Bf is the lowest-cost candidate supported by this direct repeat. Br1 challenges Bf with one bisection of every radial cell from r=.05 through .11541061001735, holding all x edges and H-adjacent owners. This 111×128 witness tests the whole blend, both flanks and the coupled response of the already challenged subintervals. The prospective expectation is that the confirmed lower-band peaks remain within bounds, with the upper transition and modified pressure providing the independent challenge. The 3,718 s complete forecast reserves 16.63 GiB against a named 17 GiB cap. The planned storage and live resource conditions fit. H convergence and independent x, quadrature and algebraic controls remain open. [Frozen wider support](wj-h27-linux/task2-outlet-campaign/cp1-bfl/outlet-campaign-br1-design.json) · [Measured full-cost forecast](wj-h27-linux/task2-outlet-campaign/cp1-bfl/outlet-qual-br1-forecast.json)

**Campaign OPEN. TASK2_COMPLETE=NO. TASK3_READY=NO.**
