# Flat outlet B collar replay and B3

## Replayed B2−B1c evidence

The independent replay exactly reproduces all 690 archived signed rows and the fixed experimental region. The failed rows contain 128 failing component maxima. Six stress-component misses sit behind a larger component maximum in the same row. The component counts are 60 Pi, 56 Rrr, four Rxx, four Rxr and four Vxr. Ninety-two of the failing component maxima lie in x∈[−2.5,8]. The failure map gives both physical owner intervals, signed differences and the raw state at every failed component maximum. The first replay compared in-memory owner tuples with JSON arrays. The corrected check normalizes the replay through JSON and requires exact equality. Both attempts are charged. [Signed replay and owner map](wj-h27-linux/task2-outlet-qual/cp1/outlet-qual-cp4-replay.json)

At the upper r=.08125 outlet owner, the signed B2−B1c shear change is −.002299986 tau0. Evaluating νΔ(Ur+Wx)−ΔRxr on those same owners gives −.000112160 tau0 from molecular shear and −.002187826 tau0 from the Reynolds shear term. This puts about 95% of the change in the Reynolds term. Both states have χ=.32570 and Wr=0 there. The raw D values are .00384743 and .00598639. The peak therefore concerns a lower collar shear change. The r=.10 trace diagnostic concerns Ryy at a different location. [Raw shear identity and state](wj-h27-linux/task2-outlet-qual/cp1/outlet-qual-cp4-replay.json) · [Prior r=.10 sensitivity](wj-h27-linux/task2-outlet/cp4/outlet-collar-normal-stress-diagnostic.json)

## B3 solve and recovery

B3 uses the full 111×77 nested mesh from the contract. It bisects all eight original B2 collar descendants and preserves all 112 streamwise edges. The first solve takes one accepted Newton update and writes a typed zero-time iterate at residual 2.25438e−8. Its wall charge is 233.662 s and its peak process-tree RSS is 8.79 GiB.

External load forced two watchdog stops during continuation. They charged 778.701 and 98.397 s, including the bounded pauses. The sampled shared census reached 94 during the first stop. The recovery restarts from the saved B3 Newton iterate, with a one-update cap in its next short segment and a 45-second quiet host admission window after the second stop. Every attempt retains the original 6×8 quadrature, coupled tolerance and raw constitutive gradients. [Original admission](wj-h27-linux/task2-outlet-qual/cp1/outlet-qual-b3-admission.json) · [Measured readmission](wj-h27-linux/task2-outlet-qual/cp1/outlet-qual-b3-admission-1.json) · [Watchdog records](wj-h27-linux/task2-outlet-qual/cp1/operations.json)

The first launch after the quiet interval stopped in 3.687 s when OpenMPI placed the rank outside its requested affinity. The corrected command disables MPI binding and pins inside the rank wrapper. A preflight verified CPU masks 0, 2 and 14. Later native receipts record the exact rank mask. The following three segments reached residuals 4.55837e−9, 1.28408e−10 and 2.75901e−11. The final segment converged in 308.300 charged seconds. All B3 native attempts charged 2,336.504 s in total. The largest sampled tree RSS was 9.34 GiB under the 12 GiB admission. [Placement, checkpoints and native receipts](wj-h27-linux/task2-outlet-qual/cp1/operations.json) · [Final readmission](wj-h27-linux/task2-outlet-qual/cp1/outlet-qual-b3-admission-3.json)

The accepted B3 coefficients pass the independent current, weak accounting, physical, boundary, interface and inlet-projection checks. Momentum, continuity, k and outer epsilon residuals are 2.64e−10, 4.68e−13, 3.91e−10 and 3.62e−10. The minimum covariance-over-k and principal-over-viscosity margins are .2810 and .004068. The global strong divergence RMS is 1.190e−4, compared with 1.194e−4 on B2. At H, the shared mass and k-advection owner jumps are zero, and the observed raw traction and k-diffusion jumps are comparable to B2. A broader radial witness should also bisect the two H-adjacent intervals to challenge those one-sided fluxes. [Independent assembly, physical bounds and owner fluxes](wj-h27-linux/task2-outlet-qual/cp1/outlet-B3-audit.json) · [Inlet projection](wj-h27-linux/task2-outlet-qual/cp1/outlet-B3-audit-inlet.json)

## Spatial outcome

The three radial differences on the original collar support reduce the missed row counts from 180 to 122 to 8. B3−B2 meets every inherited fixed-region row, with stress and Pi maxima of .000485 and .000571. Four global/outlet stress rows and four global/outlet shear rows still miss. The full-domain Rxr and Vxr peaks are at x=14.9999979 and r=.08135285. On the x=15 station their magnitudes reach .001363 and .001427 at the same radius. Both peaks sample the upper side of r=.08125 in the final x column. The compared owner radial bounds are [.08125,.0859375] and [.08125,.08359375]. Each norm contracts in 680 of the 690 rows. [All 690 rows, peak signs and both owners](wj-h27-linux/task2-outlet-qual/cp1/outlet-qual-B3-B2.json) · [Both sides of the blend and H](wj-h27-linux/task2-outlet-qual/cp1/outlet-qual-b3-blend-probe.json)

| Full-domain field | B2−B1c RMS / MAX | B3−B2 RMS / MAX | RMS and MAX limit |
| --- | ---: | ---: | ---: |
| R / tau0 | 1.62e−5 / .002188 | 5.38e−6 / .001356 | .001 |
| Pi / tau0 | 2.20e−5 / .001353 | 6.23e−6 / .000890 | .001 |
| Vxr / tau0 | 1.35e−6 / .002300 | 3.19e−7 / .001420 | .001 |
| Velocity gradient | 2.76e−5 / .003235 | 1.18e−5 / .002450 | .01 |

B2 fails its finer collar witness. B3 is the next candidate for a direct finer comparison. Scaling the 9.34 GiB measured B3 peak by the full 93/77 radial ratio and a 10% margin yields 12.42 GiB. The proposed B4 alternative bisects the eight B3 children in [.08125,.10] to give 111×85, with a forecast near 11.35 GiB. The lower-only change starts a new refinement sequence. Preserve the independent blend/flank and streamwise tests for candidate qualification.

The accepted B3 state SHA-256 is `4cea90ce0861e8fcd37532051e2f7ae44d6b5ac8135106292e57553a9d2ea5d6`. The mesh SHA-256 is `ba5a6a2af50c4071e8741fe3791c77a2390c9df0f917b0a15af59a54e3f8e55d`. The operation receipts record every executable command and parameter, and the ext4 archive includes the typed seeds, exact session XML, option files, outputs and drivers. CP1 settles at 52,004.174 cumulative charged seconds with 9,195.826 nonprotected seconds available. [Archive digest, zero-time commands and settled balance](wj-h27-linux/task2-outlet-qual/cp1/operations.json)

**TASK2_COMPLETE=NO. TASK3_READY=NO.**
