# CT4 observation-coordinate correction — CT5, 2026-11-12

This notice supplements CT4. It does not replace its sealed observations,
fits, mathematical-reference results, native trajectories or failed accuracy
decisions. The corrected version is `ClosestGraphNormalV2` /
`CT5ObservationAuditV1` in the adjacent CT5 evidence.

CT4's closest-normal calculation searched only strict sign changes and omitted
stationary roots exactly on its search nodes. At `(x,z)/Hh=(0,1.25)` it selected
the foot `x/Hh=-2.5` and returned `2.7950849718747373 Hh`. The closest foot is
the crest, and its distance is exactly `.25 Hh`. CT5 includes critical and
endpoint roots and compares every valid graph/flat candidate. It does not
clamp a distance to the vertical height.

Forty previously exported nominal crest distances change. Four nominal smooth
observations formerly labelled outer are inner at each of 5.5, 6 and 6.5 mm.
The full point IDs, old/new roots, distances, classifications and extraction
envelopes are retained in `observation-correction.json` and its hashed full
Linux receipt. The 36 visible Fig.7 observations, previously absent from this
audit, now include both RMS-ratio and shear-ratio series. Their `yprime/delta`
coordinates and extraction bounds use `delta/Hh=6.25`; original quantities,
roles, values, weights and units are preserved. The rough Fig.7 lower bounds
are approximately 9.335–9.375 mm. The rough Fig.12 7 mm traverse remains the
limiting known turbulence height. Fig.11 physical heights remain unavailable
where `delta_7.5` is unspecified.

For the compact cosine-squared graph, analytic critical points partition the
normal-stationarity equation into monotone intervals. CT5 explicitly tests
sample-node, endpoint and tangential roots, competing candidates, crest and
flat controls, foot branches, reflection symmetry, crest continuity, extraction
envelopes, and the vertical-distance upper bound. Independent distance
minimization is a separate check. Binary64 coalescent roots within the stated
arithmetic tolerance are not separately certified for multiplicity.

The CT4 label `nu_over_nut_H` means the **inner auxiliary** ratio
`nu/nu_t,A(H_A)`. CT5 retains it under that explicit name and separately
reports `nu/(kappa_m*u_tau*H_A)` for the outer scale. Outer molecular diffusion
remains zero in the selected high-Re equations. These are scale indicators,
not permission to change the equations.

No experimental/probe-position uncertainty is invented. Extraction envelopes,
the nominal boundary-layer depth and the paper's assumed displacement retain
their separate meanings. No inlet refit, wall-height change, native curved
evolution or hill prediction follows from this correction.
