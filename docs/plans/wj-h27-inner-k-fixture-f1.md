# F1 fixed-choice nonlinear column results

**COLUMN_REMEDY=SUPPORTED.** The proposed row D solves the declared full-v2
shear column with auxiliary inner diffusion and historical analytic cutoff.
All four rows use the same full stress, finite molecular viscosity, no-slip
ground and constant total shear. Each solves for its own $Q(H)$ and $U(H)$
with $k(H)=K$ and $Q(0)=0$.

| Row at Kin | k RMS / sampled MAX in Kin | U RMS / sampled MAX | $U(H)$ | $Q(H)/A$ |
| --- | ---: | ---: | ---: | ---: |
| A, scalar / v2 cutoff | .09938390 / .32642486 | .00676160 / .01573733 | .42150691 | −1.0e-14 |
| B, auxiliary / same v2 cutoff | .18717354 / .42219409 | .00679690 / .01578352 | .42146072 | +.00031460 |
| C, scalar / historical cutoff | .11653008 / .48628496 | .00229809 / .00627167 | .43097257 | +.04451675 |
| D, auxiliary / historical cutoff | .01893856 / .02146724 | .00231731 / .00632234 | .43092191 | +.04487878 |

Row B uses the unchanged v2 cutoff, giving a distinct result from the
earlier auxiliary-diffusion column that recalibrated to $Q(H)=0$. The
cutoff in row D matches the integrated historical auxiliary production
$p_{A,0}$, while its full-v2 shear column evaluates $P_{\rm ref}$ on the
actual solved profile. This leaves positive $Q(H)$, which corresponds to
negative upward diffusive flux. The coupled 2D fixture solves its own
interface flux. The reference column by itself cannot join a constant-k
outer column with common flux. [F1 admission](wj-h27-inner-k-column-admission.md)

All 36 rows across nine values of K/Kin from .9 to 1.1 isolate one adjacent
admissible sign change. The sampled endpoint sensitivities at the roots
are 2.1846–2.9368 and the smallest scanned admissible slope is 1.2504.
The archive retains rejected collar trials and independent endpoint-slope
differences. The largest fixed budget and ground errors in the baseline
grid are 7.09e-12 and 6.99e-12 in A units. Row A recovers zero endpoint
flux within 5.81e-14 across the interval. At Kin the independently refined
IVP changes common-point k by at most 1.87e-11 in Kin units. Row D's
sampled MAX change is below 1.0e-11.

The second roughness branch stays distinct throughout the interval.
The historical cutoff and diffusion sensitivities agree with independent
50-digit differentiation. The planar guards use the amended isotropic
lower block and full WJ collar. Across the tested radii the minimum
covariance/k is .2952703 and the minimum principal/mu is .0465017 with
stationary-direction search on each sampled rank-one quartic. The radial
and K coverage is finite. The increasing shear branch, cutoffs, bracket
trials, derivatives, guard samples and exact numerical limits are retained.
[Column report and input hashes](wj-h27-linux/inner-k-f1/inventory.json)

The F0 and F1 package charge is 230.694418 seconds, leaving 5,769.305582
seconds within the fixture cap. The same H27 account is OPEN at
20,003.070967 of 43,200 aggregate seconds with 19,596.929033 seconds
available before the 3,600-second protected reserve. Both phases used one
rank, all six library thread limits one, five actual workers at peak and
253,149,184 bytes peak RSS. [Settled account](wj-h27-linux/inner-k-f1/settled.json)
