# Explicit moving-cutoff x panels: numerical repair admission

The source18 radial16x12 state converges in420.124789s with6,036,103,168bytes
peak tree RSS. Independent32/48 current k residuals8.466472223361823e-9 and
7.147672466181481e-9 pass1e-8, and the fixed common-owner field comparisons
pass. Physical inner-k development still fails the unchanged.05/.10 limits.
The conditional native32 branch is not triggered and has not been executed.

The separately admitted first cutoff-crossing witness locates
`r=.02936741102104816, x=-9.43122807524855`, on adjacent original owners896/924.
Independent full conservative k volume-row errors in Ki*cell-area units are
approximately1.00e-6 (unsplit Gauss8),4.11e-8 (16),5.61e-9 (32) and3.91e-9
(48), all above the fixed1e-10 small-action target. Splitting x at the crossing
and using Gauss8 gives<=1.35e-14; independently refined split64/128 references
agree within1.75e-15. Ew-load results independently exhibit the same issue.
Keep the two different gates separate: a passing1e-8 whole equation residual
does not erase the failed local1e-10 operator test.

Implement numerical integration repair in a new immutable source19 capsule.
The approved-v2 physical tensor, transport, dependent Ew, table/cutoff,
inlet/geometry/boundaries, forcing policy and all tolerances remain fixed.
No v4 science or full Task3 run is admitted.

* Certify the frozen log-cutoff derivative is negative throughout its K
  interval using its Chebyshev coefficient enclosure. Invert this unchanged
  monotone cutoff only to precompute K levels at native radial faces.
* For each current Q3 H trace, isolate cubic K-level crossings using its
  derivative's monotone subintervals and guarded bisection. Include stationary
  points and handle constant/degenerate traces. All roots retain their native
  signed owner. No state/cutoff is clamped or averaged.
* Integrate positive x panels at the crossings, then the existing positive
  r panels. Consider every cell overlapping the qualified cutoff range so a
  narrow crossing cannot be missed by the original x Gauss nodes. Keep the
  static rules elsewhere. Rebuild these dependent panels with each current
  state and finite-difference action; cached inverses remain approximations.
* Use the exact quadratic manufactured H trace for its separately derived
  continuous forcing. Do not substitute the numerical cutoff surface into
  verification loads. Keep all application forcing unchanged.
* Upgrade the independent physical oracle with a separately implemented
  direct-cutoff/SciPy-root x split. Preserve source18 oracles and all their
  failed or passing outcomes. Independent split-order refinement, native
  signed volume-row agreement and full current equations remain required.

Before native radial reuse, qualify cubic/multiple/tangent/constant/no-crossing
and invalid-input controls; native two-owner crossing rows against the frozen
independent split witness at<=1e-10; two continuous MMS families; the source18
versus19 physical/history comparison at<=1e-8 plus independent current
equations. Require all positive scalar/root/covariance/principal guards.
Then solve the same radial5600 steady equations from the source18 candidate
only as a fresh initial guess, using native8x12 split panels and full independent
32/48 audits. Measure complete costs and preserve every failure. This cannot
by itself close the engineering, x-resolution, long-pilot or production gates.
Then cold-reload the new radial binary at zero additional physical steps.
Require all four histories to match exactly, no advancement/factorizations,
and fresh reconstruction of the nonempty face-level and panel data.

Admit materialization30s/build360s, standalone cubic compile40s/run10s,
preparation30s/native crossing40s/independent operator30s, two MMS natives60s
and audits70s, physical12 native120s/current audit150s/cross-source comparison90s,
radial native900s/current audits400+900s/same-grid comparison120s,
cold radial zero-step reload60s and evidence120s:
<=3530 additional aggregate operation seconds. This replaces the untriggered
conditional32 branch of3220s and adds310s to the prior6620s quadrature plan;
the cumulative radial/trace/quadrature/panel phase is bounded by11910s.
Carry forward every actual charge under the same43200s campaign ceiling and
3600s evidence protection. No reset. Use four compiler workers, one science
rank, all six library thread limits1, actual shared-host workers<=72,
science RSS8GiB/build16GiB, source/build8GiB/science8GiB/free15GiB caps.
