# F0 signed mechanism on the historical inlet profile

**MECHANISM_DIAGNOSIS=PASS.** The inherited ground profile has
$k_A=.01723802909157$, $k_A'=0$, and $k_A''=-111.5520289651$.
The historical and v2 cutoffs are .02641989831 and .02954486569.
At the ground $P_A=.02596488930$, $E_A=.01314647786$ and
$E_{v2}=.01175597181$. The v2 scalar diffusivity amplifies the negative
inherited curvature by a factor of 19.73620851 relative to $D_A$.

| Fixed-profile row | Ground $D$ | Ground $(D k_A')'$ | Ground $P-E$ | Ground total $L$ | Required $Q(H)/A$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| A, scalar / v2 cutoff | .00226788202 | −.252986841 | +.0142089175 | −.238777923 | −.1291107790 |
| B, auxiliary / v2 cutoff | .00011490971 | −.0128184114 | +.0142089175 | +.0013905061 | −.1291107790 |
| C, scalar / historical cutoff | .00203480141 | −.226986226 | +.0128184114 | −.214167815 | −.0173182502 |
| D, auxiliary / historical cutoff | .00011490971 | −.0128184114 | +.0128184114 | 0 | −.0173182502 |

Changing diffusion at fixed v2 cutoff gives a ground tendency increment
of +.2401684295 (B−A). Changing only the sink when diffusion is auxiliary
gives −.0013905061 (D−B). At scalar diffusion, the cutoff change alters
both destruction and diffusion, giving +.0246101087 (C−A). The four-term
interaction D−C−B+A is −.0260006147. The lower-layer decomposition agrees
with direct differentiation to 3.60e-14 in the fixed action scale.

The historical reference has integrated $P_A/A=E_A/A=2.736517844$.
Full-v2 stress on the fixed profiles gives $P_{v2}/A=2.753836094$ for every
row, including the signed collar contribution .01731825025. Active v2
destruction integrates to 2.624725315 in those units. The inherited
$D k_A'$ gives $Q(H)/A$ at roundoff. This leaves a positive cumulative
flux defect of .1291107790 for A and .01731825025 for D. Row D restores
the pointwise lower-layer balance to 7.49e-16 in the action scale. Its
remaining tendency equals the collar production change, including the
full-WJ region above the blend. These are fixed-profile interior tendencies.
The imposed inlet trace remains fixed. [Signed profiles and plot](wj-h27-linux/inner-k-f0/inventory.json)

The exact CP2 inner RMS is .1327946821 and order48 sampled MAX .4416832429.
Reintegration of its original-owner squared terms agrees to 1.39e-16.
Its separate owner-endpoint sample reaches .4418922347 at the ground
near x=−8.75. The first two x columns account for 5.39% of squared error,
and the downstream columns remain near .14 RMS. The localized ground
diffusion signature explains a strong early tendency. Transport, collar
production, and interface exchange also enter the developed solution.
CP2's independent streamwise/refinement differences leave the maxima
unresolved. [CP2](wj-h27-continuation-cp2.md)

Split orders 32/64 agree to 8.88e-16 in the budget scales. The original
one-sided ground curvature refinement missed 1e-10. Its preserved failure
motivated two further step halvings at the same precision and stencil.
The largest independent multiprecision derivative discrepancy on the
refined one-sided set is 9.10e-12 in the prescribed action scales. F0
charged 22.129 aggregate seconds including evidence. The compact archive
binds the failed and accepted attempts, raw signed rows, local literature
hashes, exact inputs and resource receipts.
[Assessment](wj-h27-linux/inner-k-f0/mechanism.json)
[Account](wj-h27-linux/inner-k-f0/settled.json)
[Contract](wj-h27-inner-k-fixture-contract.md)
