# Task 2 outlet and collar diagnostic contract

The baseline is `cabca8e91a76ad3fa365692612dc16c9312fcaf4` on `codex/wj-earsm-awf`. Linux work uses `jay@100.88.70.60` and the established ext4 run root. A0 and A1 are the accepted zero-time S2 states on the exact 110×62 M0 and strictly nested 111×63 M1. Retain Q3 velocity and transported scalars, Q2 pressure, native 6×8 panel quadrature, linear tolerance 1e−12 and coupled tolerance 1e−10. Use the exact-inner table, ScalarCdV2, finite molecular viscosity, complete WJ stress and production, positive implicit stress split with its remainder, sampled AWF force and current sampled-K cutoff. The modeled wall elevation is 3.5 mm. The matching height is .15 code units or 6.0 mm above that wall. [R2 and R3 contract](wj-h27-task2r-contract.md) · [S2 findings](wj-h27-task2rs-s2.md) · [Linux baseline](wj-h27-linux/task2-outlet/baseline.json)

| Case | Mesh | Outlet | Admission |
| --- | --- | --- | --- |
| A0 | M0 | A with frozen full traction | Reuse subject to archive parity |
| A1 | M1 | A with frozen full traction | Reuse subject to archive parity |
| B0 | M0 | B with W=0 and original normal traction | Required |
| B1 | M1 | B with W=0 and original normal traction | Required |
| A1c | M1c | A with frozen full traction | Required |
| B1c | M1c | B with W=0 and original normal traction | Candidate viability and measured budget |

M1c bisects exactly the two M1 radial intervals adjoining r=.10 and preserves every x edge, both terminal corner widths and all other radial edges. Use the archived states as typed steady guesses where the native inventory and provenance checks allow. Keep A as the default. Both laws leave axial velocity free on the outlet away from the shared no-slip corner and set the pressure level through normal traction. B removes the applied transverse load and reports the outward reaction from the full constrained momentum rows. Report raw traction separately. Group the global shared corner rows explicitly; each is counted once.

Use all 690 inherited rows and all 29 stations with their signed one-sided owners. Preserve the original scales and numerical, physical, current-equation and inlet-projection gates. Record all 16 open A rows and the four A corner maxima. Compare A1−A0 and B1−B0 as resolution sensitivity and B0−A0 and B1−A1 as boundary effects.

The fixed interior comparison region is x in [−2.5,8], with full height and separate inner and outer support. It includes the virtual crest x=0 of the Fig. 3 flat reference and the complete x=3.75–8 downstream measurement station series. Retain original measurement elevations including 7 mm from the original floor, equivalent to r=.0875. Report the other stations separately, including x=10, 12.5 and 15. At every station and in the interior region, compare componentwise RMS and MAX of each boundary change against the sum of the corresponding A and B within-policy spatial differences. Label changes at or below the sum unresolved at these two levels. A resolved effect on this screen has the same interpretation on both meshes and is assessed separately against the inherited physical limits. These indicators describe numerical sensitivity at the observed levels. [Experimental location, units and norm definitions](wj-h27-boundary-audit.md#5-prospective-physical-development-envelope-distinct-from-numerics)

Use weak pressure-tested continuity and raw strong D at interior quadrature points, on wall and outlet owners, and on fixed physical annuli about (15,0). Use .0025–.006 code units with radii .0025, .003, .0035, .004, .0045, .005, .0055, .006 and angles 12°, 18°, 24°, 30°, 36°, 42°, 48°, plus the broader inherited .0125, .025, .05, .075 radii at 15°, 30°, 45°, 60°, 75° to cover multiple owners. Judge B off the corner as well as at the four raw extrema. For B, Ux=0 along the ground and Wr=0 along the outlet give D=0 at their vertex. Retain the original gradients in all constitutive and acceptance evaluations. One B refinement difference screens sensitivity. A trend claim requires another independent difference or evidence of comparable strength.

At r=.10 on M1, use both original one-sided owner regions. Include all stations upstream and the previous collar extrema. Report D, all R components, k, dependent E, transported epsilon where supported, sampled K, gradients and production. Retain all 690 comparison rows. Any trace-free constitutive evaluation is a diagnostic beside the raw-gradient rows. A1c measures A's radial sensitivity. A candidate-specific conclusion requires B1c.

Use the higher-quadrature and tighter-algebraic controls, including each inner-k RMS/MAX limit of 1e−4, for any final-recipe claim. An extension, if admitted, must have its length, retained mesh, prior exit observations and physical tolerances frozen first. Retain the Task 3 launch hold. [Inherited gates](wj-h27-task2rs-contract.md) · [Staged operation forecast](wj-h27-linux/task2-outlet/forecast.json) · [One-time amendment receipt](wj-h27-linux/task2-outlet/outlet-cp0-amendment.json)

**TASK2_COMPLETE=NO. TASK3_READY=NO.**
