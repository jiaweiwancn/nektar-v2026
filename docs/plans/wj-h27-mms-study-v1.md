# H27 independent MMS endpoint study, prospective version 1

This specification precedes inspection of the study endpoint results. It adds
analysis code only; it changes no native equation, source, reference table,
boundary condition, integration rule, or convergence tolerance. Numerical
execution is Linux-only and is admitted/launched by the coordinator. The script
does not launch a solver. Its output is create-only JSON, including failures.

## Inputs and separate evidence scopes

`solvers/IncNavierStokesSolver/Tests/RANS/h27/mms_study.py` reads completed
`WJH27RegisteredCoupledRunV2` endpoints with coefficient/history export. Every
run must have the same model source identity, reference-table identity, MMS
variant, prescribed body amplitude, and physical terminal time. Stored state,
prior-state, and advection identities are checked. Steady artifacts are rejected.
Each endpoint has at least two completed steps, so its final equation is
BDF2/EX2, with actual retained BE startup. No missing endpoint is interpolated,
extrapolated, or treated as completed.

The formal temporal study uses `PolynomialScalarsV1`, terminal time 0.4, and
steps 0.08, 0.04, 0.02 (5, 10, 20 completed steps); further exact halvings are
permitted. Its mesh, quadrature, and solver settings are identical. This variant
has velocity Q3×Q2/Q2×Q3, pressure Q2, k Q2×Q2, and outer epsilon Q3×Q2.
Exact spatial representability does not make nonlinear constitutive integration
or projection of the current nonlinear dependent epsilon trace exact.

A separate matched-time `Historical` study may have no temporal-role rows. Its
rational outer epsilon and quartic radial k factor provide genuine spatial
approximation content. The formal temporal gate then reads `NOT_EVALUATED`, and
the overall temporal-study pass is null. Observed mesh rates are total field
error rates until temporal contamination is independently bounded. Optional
`kind: time` comparisons halve dt at that short time on an unchanged mesh and
report a 10% subordinate temporal-difference check. Such a pair is not a
three-level temporal-order certificate.

Neither scope qualifies physical inlet development, long-horizon behavior,
principal coercivity between samples, the full H27 gate, or Task 3.

## Fixed norms, resolution and gates

All field errors are recomputed using the actual signed native local maps and
the analytic functions in the independent continuous MMS checker. No target
projection, native matrix-times-target, or fitted endpoint is used. Pressure is
compared in its prescribed gauge; no post hoc mean alignment is made. Physical
gradients are evaluated from one element's signed basis. On shared boundaries,
both owners are sampled; gradients are never averaged.

For a field q with reference qref, RMS is
`sqrt(integral |(q−qexact)/qref|² / domain_area)`. Vector velocity uses Euclidean
length and its gradient uses Frobenius length. The fixed scales are Uref=1,
Lref=1, Tref=1, kref=Kin, epsilonref(r)=Ain/r, and Piref=tau0. Gradient scales
are the corresponding field scale divided by Lref. Scalar-gradient errors mean
`|grad(q−qexact)|/qref`, not the derivative of a spatially normalized field.
The script also reports unscaled physical RMS and physical sampled MAX for U,
k, epsilon, Pi, and all their planar gradients. Epsilon and its gradient use
only the outer domain; other fields use the full domain.

Independent Gauss 20×24 and 24×32 rules split radial integration at both current
and exact reference cutoffs, collar endpoints 0.4H/0.75H, h, H, and cell edges.
Vertical boundary integration is split by the same rules. Each MAX samples all
Gauss nodes and an endpoint-inclusive owner grid, refined from 17×17 to 33×33.
These maxima are sampled diagnostics, not certified global suprema. Their
change is reported against a 1% criterion, without imposing MAX temporal order.
RMS integration is resolved when refinement changes the error by at most 0.1%
of its refined value or the fixed reference-scaled band 1e−12.

The inherited temporal gate remains RMS order at least 1.8 for U, k and epsilon
on both successive dt pairs. Errors at/below 1e−12 do not generate an infinite
order or a fabricated pass. Pressure, gradients, and sampled MAX orders are
diagnostics. No MAX-order condition is added to the inherited gate.

Quadrature, mesh, and tolerance controls use the finest temporal dt, changing
only their declared setting. Their common-point RMS differences must be at
most 10% of the finest temporal analytic RMS error for U/k/epsilon, with a
fixed 1e−12 band. All other field and sampled-MAX differences are quantified.
Mesh controls require nested rectangular edges, more cells, and a smaller
maximum cell diameter. Quadrature controls raise at least one native order;
tolerance controls tighten at least one convergence tolerance. A temporal-study
pass requires all three control kinds. Common-point comparisons use a union
mesh, unique element owners, and both independent rule/grid levels.

Current-equation verification retains momentum, continuity, and gauge limits
1e−9 and k/epsilon limits 1e−8, with independent quadrature differences no more
than 20% of those limits. Actual EX2 histories and current/prior/older essential
traces retain 1e−10 limits. Scalar dual residuals use the true weighted Riesz
matrices (Kin² or (Ain/r)²); no alpha scaling or `max(1,RHS)` normalization is
introduced. Stored native admissibility minima are identified as native sample
coverage, not as an independent global certificate.

## Physical balances and algebraic reactions

Write `delta_t q=(3q_n/2−2q_(n−1)+q_(n−2)/2)/dt` and
`EX2 F=2F_(n−1)−F_(n−2)`. For k, physical balance is

```
B_k = integral delta_t k
    + integral_boundary EX2(U k)·n
    + integral_boundary (−Dk grad k·n)
    − integral P + integral E − integral f_k.
```

For outer epsilon, replace k by epsilon, Dk by De, production by
`1.44 P epsilon/k`, and destruction by `1.82 epsilon²/k`. Actual current k and
destruction E determine Dk/De; inner E is `A(K)/max(r,c(K))`, never a dummy
outer-epsilon value. The H epsilon flux uses the outer owner's gradient and
outward normal −ez even though H is essential. Polynomial k has
`k_r(0)=−K eta/R`, so its nonzero ground flux is explicitly integrated. All
exterior advective fluxes are retained.

The prescribed MMS scalar source integral is obtained from independent
continuous weak forcing plus its independently integrated prescribed natural
flux. This is an integration-by-parts identity of analytic continuous fields,
not cancellation against a native discretized target. For k the callback load
includes all exterior faces, including rows later eliminated at the inlet. For
epsilon it excludes H; the independent analytic H term remains in the volume
forcing identity. Each face reports physical owner flux and prescribed flux
separately, including essential faces.

Let `V=nu(grad U+grad U^T)−R+2kI/3` and let `Pi` be modified pressure.
Mechanical balance, tested by current U, is

```
B_m = integral U·delta_t U
    + integral U·EX2[(U·grad)U]
    + integral grad U:V − integral Pi div U
    − integral_boundary U·(V−Pi I)n
    − integral U·f_M − integral U·b − W_AWF.
W_AWF = integral_(r<h) u_n { EX2[2u_H (u_x)_H]
                           + (Pi_x)_H,n − (2/3)(k_x)_H,n − b_x,H }.
```

The convective, pressure, known-k, and sampled-body work terms are individually
reported with their physical signs. Bulk body work is separate. `f_M` is the
continuous verification source with prescribed bulk and sampled-body actions
removed from the full RHS accounting; application forcing is not replaced.
All current stresses and boundary tractions use owner gradients and full v2
stress. The BDF2 kinetic identity is checked independently:

```
integral U_n·delta_t U = [ ||U_n||² + ||2U_n−U_(n−1)||²
                       − ||U_(n−1)||² − ||2U_(n−1)−U_(n−2)||²
                       + ||U_n−2U_(n−1)+U_(n−2)||² ] / (4dt).
```

Balance rate scales are fixed: Kin×full_area/Tref for k,
`Ain log(R/H)/Tref` for outer epsilon (unit x length), and
`Uref²×full_area/Tref` for mechanical work. Reported quadrature uncertainty is
the absolute change between the two independent integration rules. On h
refinement the absolute normalized owner balance must decrease beyond both
quadrature uncertainties, or both values plus uncertainties must lie within
the prospective fixed 1e−10 balance band. This is a numerical study criterion,
not measurement uncertainty or a change to the equation residual limits.

The constant-test scalar essential reaction, free-row defect, and assembled
diffusion cancellation are reported separately. Mechanical essential reaction
power and free velocity-test defect are likewise separate. Their algebraic
identities retain a 1e−10 fixed-scale check. In particular,
`B_scalar−one^T residual = owner_flux−prescribed_flux` and
`B_m−U^T residual = −owner_complete_traction_power`. These identities explain
why an algebraic reaction balance is not a physical owner-flux proof. The
script additionally reports signed and absolute two-owner interior diffusive
flux jumps and complete-traction power jumps. Shared CG test functions do not
make these pointwise physical fluxes equal. The current nonlinear H trace's
pointwise target discrepancy is also reported independently of its projected
essential-coefficient constraint.

## Manifest and execution interface

Paths are absolute or relative to the manifest. The dense independent analysis
is prospectively bounded to at most 128 cells and 4096 coefficients in each
native field. The coordinator must admit memory and a finite wall-clock cap;
single-thread numerical-library environment variables are set by the script.

```json
{
  "schema": "WJH27MMSStudyManifestV1",
  "terminal_time": 0.4,
  "runs": [
    {"label":"dt08", "native":"dt08.result.json", "role":"temporal"},
    {"label":"dt04", "native":"dt04.result.json", "role":"temporal"},
    {"label":"dt02", "native":"dt02.result.json", "role":"temporal"},
    {"label":"qfine", "native":"qfine.result.json", "role":"quadrature"},
    {"label":"hfine", "native":"hfine.result.json", "role":"spatial"},
    {"label":"tolfine", "native":"tolfine.result.json", "role":"tolerance"}
  ],
  "comparisons": [
    {"label":"q", "coarse":"dt02", "fine":"qfine", "kind":"quadrature"},
    {"label":"h", "coarse":"dt02", "fine":"hfine", "kind":"mesh"},
    {"label":"tol", "coarse":"dt02", "fine":"tolfine", "kind":"tolerance"}
  ]
}
```

For a spatial-only Historical study, use its common short terminal time, omit
temporal roles, and include mesh comparisons. An optional comparison
`{"label":"short_time", "coarse":"hfine", "fine":"hfine_dt_half", "kind":"time"}`
bounds remaining time error by the common-point field difference; no temporal
order gate is evaluated. A successful endpoint analysis is not a spatial-order
certificate unless the accompanying spatial-rate and contamination evidence
also support that claim.

The coordinator runs, under an admitted Linux analysis entry:

```sh
python3 -B h27/mms_study.py --manifest study.json --table reference.json \
  --seconds ADMITTED_SECONDS --out study.result.json
```

The output binds the manifest, every endpoint, reference, script, independent
continuous MMS, basis, constitutive oracle, and coupled oracle by hashes. It
always retains `full_coupled_gate_pass=false`. A study result cannot override
an earlier failed current-equation, direct-reference, or physical gate.
