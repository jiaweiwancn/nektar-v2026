# Current-data inverse preparation and offset cold restart

This checkpoint follows `c88e26ac7f9d55c18b4c5139949de628015cc182`.
The active approved model remains `FiniteNuIsotropicWallFullWJCollarV2`:
physical nu in both regions, isotropic lower-wall normal stress plus auxiliary
shear, the same C2 collar into full WJ, raw actual production, and inner
`Dk=nu+Cd*k*k/(sigma_k*Ew)`. `Ew=A(K)/max(r,c(K))`, the v2 reference cutoff,
the retained inlet, finite-nu boundaries, conservative transport and current
admissibility/residual gates are unchanged. V4 remains inactive.
`TASK3_READY=NO`; full H27-C/D and the complete production-cost gate remain open.

## Implementation and identities

The new opt-in `CurrentDataDriftV1` compares a dimensionless current-data tuple
on the fixed original quadrature owners with its prepared tuple. Its exact
entries/positive reference units are in the
[prospective admission](wj-h27-current-data-cache-admission-v1.md). It includes
current k/E/aux/full planar G and all inverse coefficients. The relative L2
norm has no denominator floor. This heuristic has its own diagnostic name,
`ScaledCurrentQuadratureDataL2V1`; it is not an operator-error bound.

On reused bases, skip frozen velocity-tangent and scalar inverse-matrix
preparation. Keep the actual-alpha/.02/eight-base/cold refresh rules and all
three factors, matching saved epsilon full operator, current Gamma derivative,
typed same-dt cold retry and single step deadline. Current residuals, FD Jv,
forcing and accepted-state guards never read the cache. Defaults and the
previous actual-matrix `OperatorDriftV1` control remain available.

Zero-NumSteps terminal serialization is moved under final export, repairing
source14's unqualified zero-step timing path. Advancing output behavior is
unchanged; derived caches remain absent from restart files.

Linux namespace: `/home/jay/build/wj-h27-coupling-r1-20261117`.
Source15 scientific identity:
`47fb3597c4e56348361cfab6eadf9a34ba8d89027e57ed01153ca4d9c5ab10d1`.
Runtime ELF SHA256:
`be3dd8fa87236e3c25a5e6b890077090587dd6f85fa83cc309a31a46415b1a05`.
Reference-table SHA256 remains
`b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d`.
All1717 locally present C++/header/CMake files in the executed capsule match
this Mac checkpoint. Build: normal production main, Release/GCC11.4,
`-ffp-contract=off -fno-fast-math`, four compiler jobs,199.824672s, at most15
actual shared-host workers. All science uses one serial process/rank and
all six numerical-library thread counts1. Exact commands/inputs/receipts and
runtime dependencies are preserved in both linked evidence inventories.

## Independent numerical and history results

PolynomialScalarsV1/body1 four steps at dt=.02 and Historical MMS eight
steps at dt=.0025 both pass fixed<=1e-8 fresh-versus-current-data coefficient,
owner field, complete history, essential lift, Gamma/Ew/cutoff and explicit
mass-dual comparisons. Both current-equation oracles pass. Threshold-zero
control refreshes on each nonzero measured data change and its independent
current equations pass. No mathematical acceptance threshold is changed.

On the original14x50=700-cell Q3/Q2 physical case,12 startup steps at dt=.005,
nonlinear3e-12/linear1e-12/no Anderson, both OperatorDriftV1 and CurrentDataDriftV1
pass the fixed<=1e-8 comparison against fresh factors. At the current-data
endpoint scaled coefficient MAX is U9.29726e-15, Pi8.17128e-10,
k5.19885e-14 and epsilon2.38363e-15. Pressure is not shifted.

The new-policy control has7 triple refreshes and41 reused Newton bases.
Cold, actual-alpha and age refreshes are observed. Its physical .02 drift
trigger is not observed; the separate zero-threshold control exercises
nonzero-data-change refresh. All per-base reasons/ages remain in the logs.

Prefix3+cold9, actual retained recovery9+cold3, recovery12+zero steps and
one retry at step4 pass the independent physical history foundation. Recovery
writes3/6/9/12 and actual deletion3/6 retain usable9/12. All four numbered
scientific snapshots remain, with the descriptor bound separately. Prefix3
and recovery9 continuations happen to coincide with the uninterrupted
inverse refresh schedule and match all six available state/history hashes
byte-for-byte. The zero-step final binary also matches its input exactly.

An additional prospectively admitted prefix4+cold8 test addresses that
alignment: uninterrupted step5 reuses age4 factors; cold step5 starts at age0
and performs a real triple rebuild. It passes all complete state/explicit/
older-history, owner field, essential and dependent-state limits<=1e-8.
Its endpoint is not byte-identical: current coefficient MAX is U1.49004e-14,
Pi1.37051e-9, k5.78898e-15 and epsilon1.72901e-15. The pressure level is retained.

| Advancing endpoint | Independent momentum | Continuity | k | epsilon |
| --- | ---: | ---: | ---: | ---: |
| Current-data control | 1.08642e-14 | 6.62851e-16 | 1.32075e-13 | 1.17847e-14 |
| Cold3→12 | 1.08642e-14 | 6.62851e-16 | 1.32075e-13 | 1.17847e-14 |
| Recovery9→12 | 1.08642e-14 | 6.62851e-16 | 1.32075e-13 | 1.17847e-14 |
| Retry4→12 | 8.66787e-15 | 6.48460e-16 | 1.30096e-13 | 1.17910e-14 |
| Offset cold4→12 | 1.00145e-14 | 6.45930e-16 | 1.21710e-13 | 1.16512e-14 |

All supplied current-equation/history-boundary/quadrature/cache/accounting
subgates are `MET`. Overall physical-oracle reports remain `NOT_MET` because
reference, discretization and long-time gates remain unresolved.

## Measured short complete costs

These blocks include actual scheduled logs/recovery/scientific output.
Each final physical step has four Newton updates, seven Krylov iterations
and30 current residual evaluations. Preserve the full per-step variation;
the three post-startup blocks are still a short unforced diagnostic.

| Steps | Fresh factors, s/step | Operator drift, s/step | Current data, s/step |
| --- | ---: | ---: | ---: |
| Startup1–3 | 11.22374720 | 8.73468855 | 6.60053037 |
| 4–6 | 11.04010699 | 7.09556100 | 6.33374149 |
| 7–9 | 10.99479771 | 7.46398634 | 5.83401900 |
| 10–12 | 10.98425307 | 7.10885258 | 6.30173802 |
| All12 | 11.06072625 | 7.60077212 | 6.26750722 |

Whole operations are137.544478s fresh,95.958922s operator drift and80.075016s
current data. The mean advancing-step ratios are1.76477x fresh/current-data
and1.21273x operator/current-data. The two drift metrics differ and do not
promise equal refresh frequencies; actual events are retained.

Prefix3 takes25.329635s; cold9 continuation61.878692s; recovery9's cold3
continuation26.024490s; zero-step recovery12 import/export6.486796s; the12-step
retry trajectory86.390377s. The extra offset prefix4 takes32.108565s and
cold8 takes55.699432s. These include full process construction/init/export.
They must not be added to a forecast that already counts their accepted
steps again. No complete long-trajectory or campaign forecast is qualified.

Current-data native advance categories attribute99.9673% of elapsed native
advancement. Guards cost22.243407s over12 steps. The raw assessment's
`known_equation_work_fraction_lower_bound` value .652531 is a diagnostic
category ratio: constitutive, wall, assembly, solves, factors, residuals and
nonlinear orchestration over the whole short native operation. Assembly
includes current-data metric comparisons; it is not a separately instrumented
formal genuine-equation-only fraction. Initialization equation work is omitted
from that numerator and qualification audits are outside both numerator and
denominator. The production50%/preferred70% objective remains unevaluated.

Zero-step terminal output now measures .079654866s frame serialization plus
.000058602s index within .112786761s final export, with zero advancement and
no inverse preparations. The sole frame is explicitly imported history;
default zero residual/guard fields are not accepted-step evidence. The one
final frame is not a complete terminal-window assessment.

## Evidence and corrected outcomes

Main [evidence inventory](wj-h27-linux/native-current-data-cache-v2g/inventory.json):
16,953,921byte immutable Linux archive, SHA256
`59497819621482b72a0e2f55a00026828aadd83e41a5aa646e71dccd387f0111`.
The [offset addon](wj-h27-linux/native-current-data-offset-v2g2/inventory.json)
is3,971,072bytes, SHA256
`aa2d5820e9ca81ec8bee1479b2103f78a75994fc9787ba793a2f5cc26152e6cf`.
Full raw physical coefficients/owners/history/profile data remain on Linux.

The first MMS current-oracle invocation fails its nested time-budget guard
before equation assessment (1.597125s); the corrected explicit25s total/20s
direct command passes. The first history assessment fails after27.677273s
because it counts `.science-index.json` as an extra numbered frame. Its
test-only correction binds that descriptor separately and still requires
exactly3/6/9/12, with no unknown scientific files. The corrected complete
history assessment passes in35.381772s. Both failed commands, reports and
receipts remain immutable and charged. No native result or gate was edited.

The main and addon publications take9.637170s and3.233172s. The campaign
account remains OPEN at11873.077442873851 aggregate operation seconds,
including all carried-forward work and failures. The43200s ceiling and3600s
protected evidence reserve are unchanged. This entire source15 matrix and
addon fit the original3060s phase admission.

## Open closure gates

The v2 same-inlet reference still misses the inner-k envelope (last5600-cell
RMS .1326434/MAX .4236560 versus .05/.10), and its numerical refinement budget
is not met (inner-k .00802872/.13088236 versus .001/.001 on the latest x-only
refinement). These are unresolved finite-resolution results, not proof of a
mesh-resolved limiting barrier. The separately proposed v4 transport/cutoff
candidate remains pending its scientific decision.

Required evolved duration>=26.13350407707765, baseline/startup/both-sign finite
perturbations, selected h/p/q/dt, monitoring escalation, scientific retention/
storage reduction and complete per-case/campaign cost remain open. Raw full
scientific output is still retained. The next operational work addresses that
storage gap. The original Task3 horizon505.02466521278404 and terminal interval
[447.8483832973273,505.02466521278404] remain unchanged. No qualified full-Task3
launch command exists yet.

The subsequent [scientific-retention checkpoint](wj-h27-scientific-retention-checkpoint-v2.md)
qualifies persisted short-window statistics, actual owned raw/statistics
rotation and both retained cold imports on this same executable. It preserves
the unresolved physical, selected-recipe, monitoring and complete-cost gates.
