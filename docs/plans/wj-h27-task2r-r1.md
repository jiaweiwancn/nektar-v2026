# Task 2R R1 reference qualification

**R1=PASS.** The user supplied the recalibration ZIP in Downloads after R0. The full archive matched SHA-256 `1d4fff1e0efe09820bf3189d2c670ac4b55303536f7af45664760fad4421a619` and all 25 manifest entries verified on the Mac and on `jay@100.88.70.60`. The reduced shooting, independent collocation and full-tensor programs reran on the new host. The stationary analytic lower budget uses $T_\mathrm{ref}^2/a$ at finite viscosity. The collocation run solves full split column equations with roughness and cutoff as unknowns.

| Kin reference | Code value |
| --- | ---: |
| Calibrated roughness, reduced shooting | 0.11541061001726277 |
| Calibrated roughness, collocation | 0.11541061001728480 |
| Fixed roughness used across K | 0.11541061001735275 |
| Represented cutoff | 0.02439603217825280 |
| Collocation cutoff | 0.02439603217824522 |
| Roughness intercept | −0.002926631422545874 |
| Ground k | 0.012784290318074194 |
| U at H | 0.43724424118668737 |
| G at H | 0.8902437734714551 |
| Reference shear | 0.002658181904845143 |
| Integrated P and E | 0.0009781519636903385 / 0.0009781519636903675 |

The code matching height of 0.15 maps to 6.0 mm above the modeled wall and 9.5 mm above the tunnel floor. The fixed roughness and Kin cutoff map to 4.6164244 mm and 0.9758413 mm. Modified pressure and tractions use $U_0^2$. No physical velocity or time scale was inferred from the numerical normalization.

A degree 12 Chebyshev representation of `log(c/H)` used 18 training values and 47 independent held-out columns, including both endpoints and interleaved points near the edges. The second calculation integrates the upper ODE in log k with RK45, solves the shear in q coordinates, independently quadratures the lower production, and integrates the full lower state with the represented cutoff. The maximum c error is 1.61e−14 relative. The two fourth-order finite-difference sensitivity scales give errors up to 1.60e−11 in `Kin*dlog(c)/dK` and 1.76e−11 in `K*dlog(c)/dK`. At Kin, `c'=-0.6764507344720717` and `Kin*c'/c=-0.2442554888278606`. The maximum fixed-radius epsilon K derivative error is 2.01e−11 in `Kin/E` units. The represented inverse residual is 6.94e−18 code length and the independent K recovery error is 6.34e−14 in Kin units. The Chebyshev coefficient bounds including roundoff give a strictly negative log-cutoff slope in [−32.256469, −23.680550]. With positive endpoint radii below `.4H`, this certifies the represented interpolant over the interval. Held-out endpoint, integrated budget and total shear defects are at most 5.73e−14 in their fixed scales.

The independent general 3D stress reconstruction checks covariance and the minimum over all supported planar directions at each evaluated radius. Kin minima are 0.2817664872 and 0.0255601826. The lowest sampled minima across the family are 0.2817664872 and 0.0245836963. The radial and K samples provide resolved numerical evidence and retain the runtime guards.

The frozen inlet panels preserve all material breaks. Integrating fits of U', k' and Q' yields U and k value/derivative errors of 9.15e−14/2.61e−10 and 4.26e−13/2.70e−10 relative to independent collocation, in scales $1,1/H,K_\mathrm{in},K_\mathrm{in}/H$. The highest fixed-scale defect among the eight profile channels is 2.96e−9. The independent continuum matched shear at H is 0.002658181904845143. The top value is 0.0026014307488123924, with outward epsilon flux 2.613996809487407e−7. The outlet streamwise traction is `-Kin/Ne-2*Pi` below H and `-Kin/Ne` above it. The outward normal is $e_x$ at the outlet and $e_r$ at the top. Reference inner flux and pressure data are frozen with the boundary artifacts. Native assembled actions and on-mesh projection remain R2 gates.

The independently integrated outer finite-nu RHS per unit x width is approximately −5.6751156e−5 for streamwise momentum, zero for k, and 9.9147797e−7 for epsilon. These are the prospective diagnostic integrals for the unchanged outer logarithm and dissipation profiles.

The qualified table SHA-256 is `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`. The ext4 evidence archive SHA-256 is `4b6800ced6fbea6e13e2eb40f14d8bad5edb464fd50d9ba6ba16a9b36a852dcb`. Accounting preserves the previous 22,438.643 seconds. The settled charge is 22,496.738 seconds, with 17,103.262 nonprotected seconds and 5,786.552 seconds in the initial tranche. The first archive preflight failed before launching its command and was charged. The launcher corrected argument classification and the rerun verified both the research archive and the retained mesh. No reservations remained at R1 settlement.

**R2=NOT_RUN. R3=NOT_RUN. R4=NOT_RUN. RECALIBRATED_NATIVE_REFERENCE=NOT_RUN. SHORT_SCREEN=NOT_RUN. TASK2_COMPLETE=NO. TASK3_READY=NO.**

[Qualified table](wj-h27-linux/task2r-r1/tables-exact-inner-v1.json) · [Independent results](wj-h27-linux/task2r-r1/qualification.json) · [Archive manifest](wj-h27-linux/task2r-r1/r1-input-MANIFEST.sha256) · [Settled account](wj-h27-linux/task2r-r1/accounting.json)
