# F1 fixed-choice reference columns

F0 separates the four choices on a single historical profile. F1 keeps
those same choices while solving the nonlinear column for each row. In
all four columns prescribe $k(H)=K$, $Q(0)=0$, constant total shear $t=K/q_e$
and $U(0)=0$. The endpoint flux $Q(H)$ and $U(H)$ are outputs. The full-v2
collar production depends on the column solution. For row D, calibrate the
cutoff against historical auxiliary production $p_{A,0}$, with no extra
source in the collar. No column flux supplies an actual 2D interface load.

Use the exact v2 cutoff table for A and B and the analytic cutoff for C and
D. Because both cutoffs are below .4H, integrate the increasing-shear
branch in the shared collar and express the lower-layer contribution to
the shooting residual analytically. On $Q(H)/A\in[-.4,.4]$, scan equally
spaced .05 increments, retain all invalid trials, isolate adjacent
admissible sign changes and independently evaluate their local sensitivity.
Record the dimensionless derivative of ground flux with respect to the
endpoint flux. A finite scan and sensitivity checks supply evidence on the
sampled branch. Verify each resulting whole-domain IVP with fixed cutoff,
ground flux and independent split budget error at most 1e-10 in units of A.
Scan K/Kin=.90,.95,1,1.05,1.10 and held-out .925,.975,1.025,1.075. Check
root uniqueness within observed brackets and look for additional sign
changes. Confirm that row A recovers its v2 near-zero endpoint flux.

Use DOP853 at rtol 2e-11 and maximum normalized step .025. At Kin repeat the
complete four rows at rtol 5e-13 and maximum step .0125. Integrate physical
budgets with split Gauss48 and repeat at Gauss80 for the refined Kin rows.
Include every cutoff, .4H, h, .75H, ground and H in comparative samples, and
use 257 uniform samples plus the budget nodes, refined to 1025 plus the
fine budget nodes at Kin. Compare to the exact retained Kin profile. The
Kin refinement requires the common-profile k/U differences and RMS changes
at most 2e-9, and change in sampled k/U MAX at most 2e-5. These are
subordinate to the original numerical and physical budgets.

Confirm the historical cutoff and diffusion derivatives with independent
50-digit differentiation on the admitted interval and reject invalid K.
For covariance and the supported planar principal response, use the
isotropic auxiliary block plus full WJ stress. Evaluate both common pure
shear and off-shear gradient changes. Use two centered fourth-order
step sizes and stationary directions of the corresponding planar
rank-one quartic. Record the finite radius coverage, positive shear,
covariance/k above 1e-10 and principal/mu above 1e-8, and a derivative
refinement agreement below 1e-8. Retain the exact discrete action limit
1e-10 for the ground and integrated budgets, the shear and the independent
coefficient actions. A positive mechanism decision for D requires its
Kin inner k RMS/MAX within .05/.10 and velocity RMS/MAX within .01/.02,
with the finite-interval guard and bracket evidence in this contract.

After F0 passes, admit 500 seconds of analysis including all scans and the
Kin refinement, and 35 seconds of F1 evidence. Read the live H27 account
before execution and carry the actual charges forward under the same
6,000-second fixture cap. This phase performs no native solve or new build.
[Fixed-profile contract](wj-h27-inner-k-fixture-contract.md)
