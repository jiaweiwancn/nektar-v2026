# Bt2 exact-mesh algebraic control

Bt2 retains the exact 111×92 Bf geometry and native 6×8 quadrature. It tightens the linear and coupled stops to 3e−13 and 3e−12. The first two updates in 791.715 s save a typed iterate at 3.78039e−12. Two updates after reloading that packet take 524.175 s and converge to 2.30396e−12. The native cost is 1,315.889 s and the peak is 11.39 GiB within 13 GiB. Independent audit and comparison take 207.957 and 104.080 s. [Measured native operations and exact recovery](wj-h27-linux/task2-outlet-campaign/cp2-bt2/operations.json)

The independent audit passes momentum at 2.63494e−10, continuity at 1.74048e−14, k at 3.90629e−10 and outer epsilon at 3.61790e−10. Boundary, weak accounting, shared flux, physical, inlet projection and native guards all pass. [Independent Bt2 audit](wj-h27-linux/task2-outlet-campaign/cp2-bt2/outlet-Bt2-audit.json)

All 690 signed Bt2−Bf rows and the fixed-region rows pass. The largest inner-k RMS and MAX over the original inventory are 1.06371e−9 and 1.85606e−9 in Kin units. The MAX is negative at x=−2.5, r=.01416145, on the common owner [−3.125,−2.5] × [.0134375,.0153125]. In the fixed window, inner-k RMS and signed MAX are 6.11495e−10 and −2.08067e−9. Both comparisons pass the 1e−4 requirement. The measured quadrature and tolerance controls support retaining the nominal Bf recipe. [Full signed and fixed-window results](wj-h27-linux/task2-outlet-campaign/cp2-bt2/outlet-qual-Bt2-Bf.json) · [Independent quadrature control](wj-h27-task2-outlet-campaign-bq1.md)

**Campaign OPEN pending canonical binding and the measured outlet-location test. TASK2_COMPLETE=NO. TASK3_READY=NO.**
