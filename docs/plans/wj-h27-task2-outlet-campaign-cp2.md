# CP2 qualification of the steady flat Bf recipe

Bf is the least expensive demonstrated mesh with passing direct spatial witnesses and both numerical controls. It uses 111×92 rectangular cells, Q3/Q2, native 6×8 quadrature and 1e−12 linear and 1e−10 coupled stops. The explicit flat outlet policy is `ZeroTransverseV1`, with the retained `WJNativeCoupledV2`, `ExactInnerFiniteNuV1` and `ScalarCdV2` choices. The successful native chain cost 1,059.481 s from its typed seed and peaked at 11.36 GiB. Its independent audit cost 207.663 s. The archived coarser Bn2 recipe missed 12 signed rows against Bf. [Selected recipe and exact witness bindings](wj-h27-linux/task2-outlet-campaign/cp2-selection/outlet-campaign-selected-recipe.json) · [Measured Bf operations](wj-h27-linux/task2-outlet-qual/cp2/operations.json) · [Bn2–Bf signed results](wj-h27-linux/task2-outlet-qual/cp2/outlet-qual-Bf-Bn2.json)

| Gate | Passing evidence against the same Bf state |
| --- | --- |
| Radial collar and neighbors | Bfl, Br3 and Br4 cover 36 original radial owners over [.05,.11541061001735]. [Direct radial coverage](wj-h27-task2-outlet-campaign-br4.md) |
| Streamwise | Bx2 and Bx3 cover the 25 owners touching the 23 experimental stations, the first inlet owner and the final three owners. Independent momentum is 7.965e−10 and 7.597e−10. [Staged streamwise coverage and signed overlap](wj-h27-task2-outlet-campaign-bx3.md) |
| H flux | Bh1 and Bh2 use the same original H support. Finest maximum jumps are 5.006e−8 in k diffusion and 1.753e−5 in total traction. Both shared flux checks pass. [Frozen H proof and signed direct rows](wj-h27-task2-outlet-campaign-bh2.md) |
| Quadrature | Bq1 at 8×12 and base stops. Largest original-inventory inner-k RMS and MAX are 9.366e−10 and 1.417e−9. [Quadrature control](wj-h27-task2-outlet-campaign-bq1.md) |
| Algebraic | Bt2 at 6×8 and 3e−13 / 3e−12. Largest original-inventory inner-k RMS and MAX are 1.064e−9 and 1.856e−9. [Algebraic control](wj-h27-task2-outlet-campaign-bt2.md) |
| Signed comparisons | Every listed direct witness passes all 690 rows on the original 29 stations and the fixed interval [−2.5,8]. [Witness inventory with exact comparisons](wj-h27-linux/task2-outlet-campaign/cp2-selection/outlet-campaign-selected-recipe.json) |
| Current and physical | Bf and all qualifying witnesses pass independent weak balances, current residuals, boundary work and reactions, positivity, root, covariance, principal, current-K and inlet projection checks. [Bf audit](wj-h27-linux/task2-outlet-qual/cp2/outlet-Bf-audit.json) |

Bf independent momentum, continuity, k and outer epsilon are 2.636e−10, 5.703e−12, 3.919e−10 and 3.616e−10. The full-domain physical RMS/MAX relative to the inlet are .000262/.00301 for U, .0000359/.000186 for W, .00176/.0403 for k/Kin and .00136/.0129 for local outer epsilon. Ground shear departure from tau0 has RMS .01209 and MAX .02198. Its separate departure from Tref, scaled by tau0, has RMS .01021 and MAX .01389. [Bf physical and inlet audits](wj-h27-linux/task2-outlet-qual/cp2/outlet-Bf-audit.json)

The first selection verifier applied the ordinary k limit to the whole-domain inner row because its key is `inner/k`. The measured row already used the required 1e−4 limit. The corrected verifier checks the final region and field in every row and passes the same archived evidence. Both validation receipts remain charged. [Exact commands and receipts](wj-h27-linux/task2-outlet-campaign/cp2-selection/operations.json)

The CP3 handoff starts from Bf state SHA-256 `23ee2191573f48a4312352ea5549fb1b37b76a44347eb102f17593b000e85782` and mesh SHA-256 `3a14dd102a8f1c829b7a5331bc7d5bedc35664c3a3d736caa6c0bc3bcd0b0ff4`. The tested binary SHA-256 is `3e8db421575122676808adc6c9783f0c897136e9d9d283b9d0cb200c5ec2df45`. The exact-inner table SHA-256 is `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`. CP3 next demonstrates x=15 parity and measures the matched x=17.5 extension. Physical time advanced is zero.

**CP2 MET. Campaign OPEN for CP3. TASK2_COMPLETE=NO. TASK3_READY=NO.**
