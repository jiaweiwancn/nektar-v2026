# Task 2R R2 native qualification

**R2=PASS.** The registered `WJNativeCoupledV2` application selects `ExactInnerFiniteNuV1` explicitly. The default historical selection passes under the original artifact. The unchanged tensor and `ScalarCdV2` use a single new roughness for auxiliary viscosity and the complete sampled-force support. The existing native rules sort the interior material breaks, including the new .11541061001735275 force height above .75H. The new pressure seed comes from the frozen profile. The top tangent and outer outlet shear use the local finite-nu extension. No reference shear or k flux is imposed on the evolving state at H.

The constructor authenticates the distinct table and its calibration, cutoff, inlet, loads, geometry, units and source bindings. The native history signature binds that reference, the reference policy, transport, source, actual geometry, fields and solver configuration. Legacy-to-new and new-to-legacy accepted histories both failed at the scientific identity check. A table with a modified force-support value and a recomputed component checksum also failed. The retained physical geometry has `H_code=0.15`, or 6.0 mm from the modeled wall and 9.5 mm from the tunnel floor. The original outer logarithm and epsilon are verified against the historical native implementation.

| Check | Witness |
| --- | --- |
| Native R1 profile, gradients, pressure seed, load and current derivative comparison | 171 inlet/load and 153 current rows. Maximum fixed-scale difference 8.87e−15. Eight out-of-domain or binding negatives rejected. |
| Historical preset parity | Independent legacy reference check passes. The new and already installed Linux executables produce identical U, Pi, k and epsilon coefficients with the legacy preset. |
| Nonparallel sampled-force and interface checks | General-tensor continuous MMS and independently assembled current equations pass. Both old and new support edges are distinguishable by the support witness. Separate prescribed body and sampled-H actions agree within 1e−10. |
| A→B→A moving cutoff and cache | Crossing travels from x=.57044948 to .44007954 and returns to the same original owners. Native/independent fixed-current residual difference ≤4.64e−14. Independent advancing momentum ≤3.20e−12 and continuity ≤2.23e−17, k ≤1.09e−13 and epsilon ≤2.22e−13 on the two return witnesses. |
| History and rejection | Recomputed current actions detect the stale scalar/cutoff control. Native history rejects a checksummed stale explicit action and a separately modified dependent-epsilon trace. Injected retry, cold reconstruction and A return agree within 1.22e−13 in fixed scales. The failed step preserves the accepted history bytes. |
| Evolving temporal MMS | Matched t=.4, dt=.08/.04/.02 on the inherited 4×8 mesh, with 8×16, quadrature and algebraic controls. Orders for U/k/epsilon are 2.126/2.015/2.019 and 2.058/2.011/2.013. Independent current, history and owner-gradient checks pass. |
| Independent dense nonlinear solution | The two-cell dt=.005 experiment passes separate fixed-state/current and dense nonlinear checks. Maximum fixed-scale RMS difference is 6.85e−11, in Pi. |

The independently integrated unforced continuum weak actions agree with the exact-inner and finite-nu outer balances to 1.23e−14 in fixed scales. The outer strong RHS integrals per unit x width are −5.6751156033e−5 for axial momentum, zero for k and 9.9147796497e−7 for epsilon. The steady weak residual uses the opposite sign and includes the outward loads. The outer outlet tangent meets the interior value at H, and the top/outlet vertical and horizontal shear components have a common corner magnitude. These diagnostic residuals are excluded from application forcing.

| Native FE inlet | RMS | MAX | Limit |
| --- | ---: | ---: | ---: |
| 28×112, U | 6.00e−8 | 3.48e−7 | 1e−5 |
| 28×112, k/Kin | 4.18e−8 | 4.67e−6 | 1e−4 |
| 28×112, radial velocity gradient in Uref/L0 | 8.42e−6 | 7.84e−4 | 1e−3 |

The retained 28×100 mesh missed the gradient MAX at .0146868556. Twelve radial bisections selected by an independent endpoint-constrained Q3 calculation fix that gate while retaining x, the matching plane, all old radial edges and the original mixed orders. The eight-cell crossing mesh was too coarse for the temporal spatial-control gate, although its temporal order and other controls passed. The inherited 4×8/8×16 comparison reduced the subordinate U RMS to 3.55e−10 against 1.17e−9. A bounded dense comparison on the moving eight-cell case reached its cap. Its independently assembled current equations passed, and the dedicated two-cell dense comparison passed.

Source identity `061c48fe6805b1b1bb66a8effffc3bfc28c5f86777c5e669f52afd9f8c0d861f`, application SHA-256 `7a94420aaa31f41d9d8b6754d20c16f11d1398e894f243eaa6beaa3af7924819`, and Linux evidence archive SHA-256 `2dc1a9986a95844bdcd216184d7a40d604c2a21c8bf0c19dd6bc7c1b8f81f3c3` bind the run on `jay@100.88.70.60`. The Makefiles build uses four workers, takes 260.123 wall seconds and peaks at 4,612,952,064 bytes process-tree RSS. The ledger records a measured shared-worker peak of 52. All numerical runs use one process and six single-thread limits. R2 settles at 23,614.203 total charged seconds with 15,985.797 nonprotected seconds, 4,669.087 seconds in the initial tranche and no remaining reservation. The initial Ninja configuration, history fixture and analysis retries are retained in the ledger.

**R3=NOT_RUN. R4=NOT_RUN. RECALIBRATED_NATIVE_REFERENCE=NOT_RUN. SHORT_SCREEN=NOT_RUN. TASK2_COMPLETE=NO. TASK3_READY=NO.**

[Native reference actions](wj-h27-linux/task2r-r2/reference-new-check.json) · [Unforced continuum](wj-h27-linux/task2r-r2/continuum.json) · [Projection](wj-h27-linux/task2r-r2/projection-28x112.json) · [A/B/A replay](wj-h27-linux/task2r-r2/cache-check.json) · [Temporal study](wj-h27-linux/task2r-r2/temporal-refined.json) · [Live-ledger settlement](wj-h27-linux/task2r-r2/accounting.json)
