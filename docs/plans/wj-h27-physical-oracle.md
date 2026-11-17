# Independent H27 physical coefficient oracle

`h27/physical_oracle.py` reads the registered normal-main physical result or its
exported coefficient state. It rejects MMS states, a pressure gauge, the wrong
model/table identity, incomplete signed maps, and non-outer epsilon ownership.
It contains no time mass, historical EX2 advection, manufactured source, native
matrix or native constitutive call. Global Riesz matrices use SciPy sparse
assembly/factorization, with an 8 GiB virtual-memory ceiling and explicit
cooperative elapsed allowance. The parent resource runner remains authoritative.

The only reused numerical basis helper is the independently defined modified
Jacobi basis in `h26/fixtures.py`; its source hash is recorded. Full three
dimensional WJ uses batched cubic companion eigenvalues, selects the unique real
root at least Cp, and checks residual and positive implicit derivative. The
three dimensional deviatoric strain and raw `P=-R:G` retain weak divergence.
Lower-wall isotropization and the fixed collar polynomial are evaluated directly.

For each current coefficient state the momentum residual is

```
integral v.(U.grad U) + grad(v):V - Pi div(v)
 - integral_{r<h} v_x [2 U_H U_H,x + Pi_H,x - 2 K_H,x/3]
 - integral_{natural boundary} v.t_prescribed.
```

The scalar residual is conservative:

```
- integral grad(v).U q + integral grad(v).D_q grad(q)
 - integral v.source_q + integral_boundary v U_n q
 + integral_natural_boundary v j_prescribed,
j=-D grad(q).
```

The k source is current `P-E`; the outer epsilon source is
`(1.44 P-1.82 epsilon) epsilon/k`. Epsilon has no inner field or PDE. The
dependent H essential lift uses current `K_H`, with endpoint-constrained
Gauss16 projection and independently evaluated basis. The physical top epsilon
load is the fixed finite-nu flux. Outlet traction is the pinned table's full
total traction. Top normal velocity is essential; there is no added top
normal traction or pressure row.

Volume and vertical-face integration split at the current cutoff, .4H, .75H,
h and H. Scalar and momentum Riesz scales remain Kin, Ain/r, and unit
acceleration/velocity; there is no alpha or right-hand-side normalization.
Declared residual limits are momentum/continuity 1e-9 and k/epsilon 1e-8.

Scalar and mechanical budgets have two distinct reports. The weak accounting
identity uses only the essential reaction coefficients and must agree with the
remaining free-test action to 1e-10 in fixed physical scales. This is an
assembly identity, not evidence that the actual P-E budget vanishes. Separate
budgets integrate actual owner gradients on every physical boundary, and mixed
budgets combine physical essential gradients with prescribed natural loads.
These physical balances and interior one-sided flux/traction jumps require
resolution assessment. Inlet/H epsilon reactions are reported together because
their shared-corner modal reactions are not uniquely separable. Their physical
one-sided H fluxes remain separately observable.

The mechanical account retains convection, molecular work, full stress work,
raw turbulent production, modified-pressure work, the weak-divergence isotropic
k term, sampled roughness work, and boundary traction work. Its sum with the
actual k account cancels P algebraically; it never sets P=E.

The fixed station union and engineering envelopes come from the H27 contract
and boundary audit. Profiles retain both owners, H limits, all native face and
material/cutoff breaks, Fig.3 physical heights and the physical 7 mm point.
Their comparator is the retained imposed inlet. Pi is reported without an
invented inlet-pressure target. MAX includes split Gauss points, break
endpoints, panel midpoints, and relative 1e-8 inward offsets; it is a finite
sampled maximum. Actual native solution-grid points are not exported by the
current artifact and are not claimed to be included.

Optional finer comparison integrates common physical element intersections,
preserving both original gradient owners at every common endpoint. It also
reassembles the finer steady residual rather than trusting its native status.
Coarse/fine physical budget changes are reported separately. A field comparison
cannot alone close a gradient-flux budget or admissibility gate.

`physical_oracle_test.py` constructs signed conforming Q3/Q2 coefficients for a
nonzero analytic test field over the physical rectangle. This fixture does not
satisfy the physical essential data and is explicitly not a physical solution.
Exact moments test convection, the factor-two sampled convection, pressure and
cubic-known-k signs, every scalar conservative/advective moment, the top epsilon
flux and outlet traction signs. Missing-term mutations must be detectable, and
the fixture must fail the steady residual gate. General 3D tensor values are
cross-checked against the independent scalar polynomial-root implementation.

Source preparation and Python AST checks were performed on the Mac. All
numerical outcomes are pending authoritative Linux execution. Example commands:

```
python3 h27/physical_oracle_test.py --table TABLE --out NEW-controls.json --seconds 60
python3 h27/physical_oracle.py --native RUN.result.json --table TABLE \
    --out NEW-physical.json --seconds 300
python3 h27/physical_oracle.py --native RUN.result.json --table TABLE \
    --compare-finer FINER.result.json --out NEW-comparison.json --seconds 600
```

The overall report remains `NOT_MET` while separate full-state admissibility,
physical gradient-budget convergence, temporal and finite-perturbation gates
remain unresolved. Core equations, boundary representation, physical envelope,
and field-refinement evidence each have explicit separate statuses. A completed
measurement may exit successfully while reporting scientific `NOT_MET`;
exceptions exit 2 and retain partial evidence in the immutable output.
