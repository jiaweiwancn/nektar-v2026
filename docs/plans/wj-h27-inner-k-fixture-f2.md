# Inner-k fixture F2 result

The changed native operator gate passed. The isolated Release executable
WJH27InnerKFixture reused the native coupled operators and steady solver.
Its default-v2 result reproduced all four source19 field and advection
history hashes on the CP2 28×100 base mesh. The native k and outer epsilon
residuals were 1.883e-13 and 4.771e-14. The exact source19 independent
current-state audit remains bound to the reproduced state. The ordinary
production factory remains on ScalarCdV2. [Prospective checks](wj-h27-inner-k-operator-admission.md)

| Identity | SHA-256 |
| --- | --- |
| Diagnostic executable | 148e33dc937beb6998509de194387f5714f958df867c91ec1baa33bc171d3706 |
| Diagnostic source | f6000e775564bbdb665d8cce292511eeb39324ea68509309f75216de145ac90b |
| Shared kernel source | 9dc5aec5f2e5493b37ed4b026398cb2e5fc19b1fe3c035919775f646c79eb6b4 |
| AuxiliaryHistoricalInnerKFixtureV1 reference | cdbabf20a178cf8d1ea89e9160517b8e2276050f1efd5dc0072e8df6a61cb2ae |
| Retained v2 reference table and physical boundary identity | b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d |

The 25 coefficient, inverse and sensitivity checks had a maximum scaled
error of 1.58e-15. Both nonstationary positive states passed complete
signed scalar-action checks against independently split order20 and order32
rules. The largest signed-action error was below 7e-15. Using stale
scalar-Cd inner diffusion changed the independent action by .0143 in the
same scale. The history byte reader rejected both cross-policy directions.
The independent checker also rejected a candidate field presented under
the v2 reference. The maximum accounting error was 8.14e-15 across
both independent orders and the minimum tested current native
principal response divided by mu was .00182.

At x=.45, the lower and upper H owners had diffusivities .0030365814
and .0029859034 and normal k gradients .001875856 and .001907694.
Both upward fluxes agreed with the nonzero target 5.6961892e-6 within
6.7e-13 in units of A. The moving K state changed the dependent Gamma by
.00263687 of its nominal value. Both owner rules split the shared
cutoff face at x=.30000000000015, within 3.0e-15 of the independent
position. The current dependent epsilon projections passed the retained
coefficient gate. These action states exercise the interface without
prescribing any H flux in the actual unforced steady solve.

The first input attempt detected an existing historical cutoff face, so
the rerun reused it. The first executable's signed-source export held
an Eigen expression referencing a temporary. The full residual comparison
and the signed-source failure are both preserved. The corrected exporter
materializes the vector before serialization. The fresh v21 build retested
parity and all changed actions using identical shared operator code.

F2 consumed 352.373 seconds, including two builds, rejected attempts,
independent verification and publication. The v21 four-worker build ran
72.128 s and peaked at 4.835 GB RSS. Its 700-cell action run took
43.797 s at .598 GB, and the 2,800-cell default-v2 parity took 31.582 s
at 1.365 GB. All science used one rank and six numerical thread settings
of one. F0–F2 account for 583.068 of the 6,000 package seconds.
The existing live account has 19,244.556 non-protected seconds and
3,600 evidence-protected seconds. The fixture has 5,416.932 seconds left.
On the F2 costs and the prior 5,808-cell peak, a first matched 2,800-cell
candidate run with a conservative 4.5 GB forecast and explicit wall caps
fits. [F2 validation](wj-h27-linux/inner-k-f2/assessment.json)
[Raw evidence inventory](wj-h27-linux/inner-k-f2/inventory.json)
[Settled account](wj-h27-linux/inner-k-f2/settled.json)
