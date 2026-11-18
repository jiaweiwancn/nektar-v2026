# H27 CP2 source19 physical reference decision

**Classification C — the developed reference remains spatially unresolved
within the admitted v2 matrix.** All converged fields exceed the physical
inner-k envelope. The independent spacings still disagree beyond the frozen
numerical budgets, so the gate stays `NOT_MET` and physical long pilots remain
on hold. `TASK3_READY=NO`.

The active application is `WJNativeCoupledV2` with
`FiniteNuIsotropicWallFullWJCollarV2`. All candidates retain the original
same-inlet, unforced finite-viscosity problem, including the open total
traction, top, ground, roughness, full WJ collar and dependent epsilon
trace. Source19 is
`56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf`.
The ELF is `37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b`
and the table is
`b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d`.
The compiled inputs remain at their verified CP1 identities.
[Numerical foundation](wj-h27-continuation-cp1.md)

## Current inner-layer balance and inlet mechanism

The exact source19 28×200 reference has inner-k RMS .1327946821 in units of
Kin. The independent volume orders 20 and 32 agree in global RMS to
1.2e-16 and in normalized individual column balance terms to 2.5e-14.
The original audit reports sampled MAX .4416832429. A diagnostic rule
including endpoints on every original owner locates a slightly larger
sampled deficit .4418922347 at x=-8.75, r=0, from k=.0133453937 and the
retained k=.0172380291. First and second x columns contribute 5.39% of
the global squared departure together. The departure persists downstream,
where many column RMS values are about .14. The band from the historical
cutoff .02642 to .4H contains 48.1% of the squared departure. The band
below the historical cutoff contains another 24.3%.

| Inner-area term | Code units on the fixed 28×200 state |
| --- | ---: |
| Integrated current production | .0243968807591 |
| Integrated dependent destruction | .0241267797682 |
| Production minus destruction | .000270100990935 |
| Net exterior x advection | .0000542429372980 |
| Net exterior x diffusion | −1.83786878e-8 |
| Upward advective H exchange | −4.83272996e-6 |
| Upward diffusive H exchange, inner and outer owners | .000221353073 and .000221334652 |
| Upward ground flux | −8.95481883e-10 |

The table uses $F=Uk-D_k\nabla k$. The inner and outer sides are evaluated
with each element's own gradient. The integrated H flux jump is
1.8421e-8, and the summed absolute radial interior diffusion jump below H
is 1.9711e-5. A signed column surface defect reaches .1821 of $A$ times
column width near x=7.5. The original weak balance and its essential
reactions satisfy the discrete identity. The individual physical side
defects and gradients supply spatial approximation diagnostics alongside
that identity. [Localization protocol](wj-h27-cp2-localization-admission-v1.md)

At fixed historical parallel inlet data the historical cutoff is
.02641989831 and the active v2 cutoff is .02954486569 at Kin. Ground v2
$D_k$ is 19.736 times the historical auxiliary diffusivity applied to that
profile. The integrated current $P/A$ and $E_w/A$ are 2.753836094 and
2.624725315. A locally parallel steady column on those fixed profiles
would require $Q(H)/A=\int_0^H(E_w-P)\,dr/A=-.1291107790$ with zero
ground $Q$. Applying v2 diffusivity to the historical gradient gives
$Q(H)/A\simeq0$. This demonstrates a mismatch of the inherited parallel
profile with the v2 local budget. The actual two-dimensional solution
can adjust both production and transport, and exchanges flux across H.
The native coefficients and independent samples remain positive. In
particular, the 28×200 state has certified native minimum k .0087600377.

## Resolution and physical envelopes

The source19 experiment uses Q3 velocity/k/epsilon, Q2 pressure, current
moving cutoff panels, `FullGradientV1` steady inverse and zero elapsed
physical time. All new fields must satisfy fresh steady acceptance and
independent order32 equations. The common quadrature keeps both original
gradient owners, the open boundaries, the inlet/interface corner and every
frozen physical station. The source11 fields were initial guesses for
steady iterations. [Frozen refinement inputs](wj-h27-cp2-refinement-admission-v1.md)

| Source19 steady field | Independent momentum | Continuity | k | Epsilon | Physical inner-k RMS / MAX |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Base 28×100 | 1.590e-12 | 1.541e-15 | 1.584e-10 | 5.166e-14 | .13279450 / .44142666 |
| Base 12×16 quadrature | 6.935e-14 | 1.467e-15 | 6.628e-12 | 5.056e-13 | .13279450 / .44142666 |
| Base tighter solve | 1.590e-12 | 1.541e-15 | 1.584e-10 | 5.166e-14 | .13279450 / .44142666 |
| Radial 28×200 | 1.622e-12 | 1.622e-15 | 5.816e-10 | 5.240e-13 | .13279468 / .44168324 |
| Streamwise 56×100 | 1.532e-14 | 2.972e-15 | 2.937e-12 | 1.114e-13 | .13264342 / .42365273 |
| Targeted 44×132 | 6.663e-14 | 2.652e-15 | 5.177e-11 | 2.461e-13 | .13267988 / .40626250 |
| Limits | 1e-9 | 1e-9 | 1e-8 | 1e-8 | .05 / .10 |

The `targeted` solve took five Newton updates and 118 Krylov iterations,
with 8.208e-13 native true residual and no physical time advanced.
Minimum native k was .0088026839, epsilon 1.7439425e-5, covariance/k
.27967747 and supported principal/mu .00426732. The boundary and weak
accounts passed. The independent physical inner U RMS/MAX were
.00531435/.01002361, W .00011749/.00061990, outer epsilon
.00922143/.06832151 and ground shear .03808521/.05394833. Those whole-layer
envelopes pass. The k gate still fails for the inner domain and at repeated
frozen downstream stations. The owner-gradient k budget normalized by its
fixed whole-domain scale is 3.415e-5 for this field, separate from the
1.717e-15 weak accounting identity error.

The pure quadrature control passes every fixed comparison norm. Its
inner-k RMS/MAX changes are 3.56e-11/8.68e-10 in Kin units. The tighter
control returns identical field coefficients from the base initialization.
Both controls pass the independent steady equations. Increasing radial
resolution alone passes the inner-k and U field budgets. Increasing x
resolution changes inner k and U by considerably more and exposes a
larger one-sided gradient difference. The targeted mesh reduces the
physical peak departure to .40626 in Kin units and its physical RMS
stays above .1326. Its common-owner differences still miss several budgets.

| Common owner comparison (order20) | Inner-k RMS / MAX | Inner-U RMS / MAX | Inner velocity-gradient RMS / MAX | Failed rows of 690 |
| --- | ---: | ---: | ---: | ---: |
| Base → higher quadrature | 3.56e-11 / 8.68e-10 | 1.74e-12 / 2.46e-11 | 1.03e-10 / 3.90e-9 | 0 |
| Base → tighter solve | 0 / 0 | 0 / 0 | 0 / 0 | 0 |
| Base → radial | 6.62e-6 / 5.66e-5 | 8.18e-7 / 4.18e-6 | .001800 / .01392 | 240 |
| Base → streamwise | .008029 / .13088 | .0005324 / .004579 | .01661 / .60986 | 549 |
| Streamwise → targeted | .004775 / .20418 | .0003125 / .003731 | .01019 / .46901 | 382 |
| Radial → targeted | .009417 / .26041 | .0006824 / .006825 | .02175 / .84093 | 548 |
| Limits | .001 / .001 | 1e-4 / 1e-4 | .01 / .01 | 0 |

The independent exact-polynomial integration of common rectangles confirms
where these field differences concentrate. The ground inner-k sample
near x=-9.9681 differs from the streamwise and radial fields by about .2031
and .2604 respectively. The streamwise-to-targeted differences are also
significant over x=4.75 to 12.5. Their inner-k RMS in the two constituent
bands is .005952 and .008566. They are .000313 and .000237 in the central
bands x=-5 to 0 and x=0 to 4.75. Radial refinement reduces RMS owner k
diffusion jumps at H from 1.8674e-5 to 2.1784e-6 in the fixed oracle scale.
The targeted H traction RMS jump falls to .0002234 with sampled MAX .01550.
This separates persistent physical k departure from numerically sensitive
pointwise features and separates improved weak/mean behavior from unresolved
owner maxima. The 2D continuum envelope cannot yet be classified as B.

## Remedy and scientific decision

Further v2 resolution work should address streamwise changes in the early
ground layer, inlet/H corner and downstream width transitions. The actual
thin-layer and owner stress/flux measurements are part of that check.
The targeted 5,808-cell solve reached 6.773 GB peak tree RSS. A uniform x
doubling of that recipe would exceed the current forecast for the 8 GiB
cap. A further v2 admission needs a newly justified mesh allocation or a
measured smaller-memory inverse, with current equations and guards preserved.

The two-change v4 proposal remains a separately reviewable model experiment.
It would use inner $D_k=\nu+\nu_A$ together with cutoff
$c(K)=H\exp(1-J(K,H))$ from the historical auxiliary-column budget. Its
independent Kin-column result gives inner-k RMS/MAX .018939/.021467 and
counterfactual $Q_{\rm ref}(H)/A=.04487878$. The counterfactual root solves
for this endpoint flux with $Q_{\rm ref}(0)=0$ and $k_{\rm ref}(H)=K$.
The two-dimensional equations would still solve common physical flux at
H, permitting distinct one-sided gradients. The computed counterfactual
flux belongs to reference/outlet construction. The proposed change needs
explicit scientific approval, a newly qualified bounded endpoint-root and
outlet table, new model/history identities, independently differentiated MMS
and a new 2D reference and pilot assessment. The Kin result supplies
a mechanism test for that decision. [Inactive v4 equations, roots and tests](wj-h27-historical-inner-transport-proposal-v4.md)

## Provenance and resources

Executed on `jay-Super-Server`, one rank and six numerical-library thread
limits of one. CP2 charged 1,976.829730 s including the localization and
publication. The live account was OPEN at 19,093.301722 of 43,200 s after
publication, with 20,506.698278 s available before the 3,600 s evidence
reserve and zero active reservations. Peak workers were 13 of 72. Peak
process-tree RSS was 6,773,403,648 bytes. Science and source/build/runtime
used 3,245,778,998 and 3,007,184,492 bytes. Available ext4 space was
68,894,314,496 bytes. CP2 needed no build.

Linux archive `evidence/continuation-cp2-source19.tar.gz` under the live H27
root contains 92 checked members plus the inventory. It is 1,103,326 bytes,
SHA256 `c8662eeb760ebd6479d27a4499cbf4cb7b266615dac09e5cdf2051a7011afad6`.
The full independent owner, profile, coefficient and field files stay in
the Linux run directories. The verified inventory binds their paths and
hashes. The archived sessions and `inputs/cp2-refinement-v19/run.py` invoke
the same immutable runtime via the finite campaign runner. Reproduction
on the live account requires new output prefixes and operation labels.
[Compact assessment](wj-h27-linux/continuation-cp2-source19/assessment-summary.json)
[Evidence inventory](wj-h27-linux/continuation-cp2-source19/inventory.json)
[Settled receipt and account](wj-h27-linux/continuation-cp2-source19/post-publication.json)
[Launch hold](wj-h27-task3-launch-hold-v2.json)
