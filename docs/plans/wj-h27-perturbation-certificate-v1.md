# H27 prospective perturbation certificate and native projection probe

This is verification infrastructure only. It changes no existing solver,
application data, perturbation helper, reference, CMake target, or active state.
No numerical execution was performed on the Mac. The coordinator admits and
launches Linux work and records all resource use. Qualification of these shapes
does not establish an evolved physical reference, admissibility of a seeded
state, consistent new histories, successful pilots, or the full H27 gate.

## Authoritative layout and immutable normalization

The native layout export is authoritative. Its initial input is the retained
700-cell physical fixture supplied by the coordinator:
`/home/jay/build/wj-h27-coupling-r1-20261117/science/coupled-fixtures-v9/physical-aa.xml`.
The table is bound to SHA256
`b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d`.
The probe requires an existing `.xml.opt` file and reads it without alteration.

The probe constructs native flow and full/outer scalar fields, exports the
actual signed native maps and geometry, and proves that every element is one
cell in the common tensor inventory. Merely taking a union of unrelated edge
coordinates is insufficient. Layout mode requires exactly 700 elements.
`Perturbation::SelectLayout` is called for `PrimaryQuartile` and
`SecondaryLowerSkew`; no radial coordinates or nx value are guessed. The
independent Python certificate repeats selection from this native export.

Each mathematical certificate freezes the complete layout, base mesh SHA,
source/reference identities, and unit-velocity global peak enclosure. Its file
SHA is computed after creation. Project mode requires that SHA and reads the
same layout from the certificate. `SetNormalization` receives that exact file
SHA; a second load is tested and must be rejected. Finer meshes must contain
the original physical knots and use the original certificate. They never call
`SelectLayout` again. Velocity/k/independent epsilon signs are separately ±1.

## Global polynomial proof and rounding

On each of the 16 spline knot rectangles,
`|UnitVelocity|²=Nx² Nr'²+Nx'² Nr²` is a polynomial of bi-degree at most (6,6).
The certificate reconstructs the Cox-de Boor pieces using exact rational
arithmetic from the exported binary64 knot values. It converts the power
coefficients to Bernstein coefficients using the exact binomial formula.

The initial exact Bernstein coefficients are converted to binary64 intervals
with outward `nextafter` reserve. Every midpoint addition and multiplication
in subsequent de Casteljau subdivision is separately rounded outward. A box's
largest upper coefficient bounds its polynomial everywhere; lower-enclosed
values at its corners/center are attained-point lower bounds. A priority queue
subdivides boxes which can still improve the global maximum. No sampled or
local-optimizer maximum is used as an upper bound. The arithmetic used to
choose a subdivision direction affects efficiency only, never enclosure.

The squared-peak relative enclosure target is 5e−9, leaving margin for the
native requirement that the velocity peak enclosure width divided by its
estimate be at most 1e−8. Square roots use a 180-bit integer enclosure, then
outward binary64 serialization. Exact rational construction has no polynomial
roundoff; the interval and serialization reserves are explicit. This is a
certificate for the mathematical shape, not a claim that native floating-point
evaluation/projection is exact.

The certificate first checks its enclosure algorithm against independent
closed-form examples: an interior maximum at non-grid rational coordinates,
a corner maximum, a half-level area exactly 7/20, and a polynomial energy
exactly 37/300. These controls execute only in the admitted Linux analysis.

The projection checker independently constructs the actual native field from
the signed Q3 coefficients. In local t∈[0,1], its Modified-A basis is
`1−t, t, t−t², −2t+6t²−4t³`. All coefficient operations are exact rational
operations on the exported binary64 coefficients. It subdivides at every
native and prescribed shape break, including tiny geometry-roundoff slivers.
The projected velocity/k minus target is polynomial. The normalized epsilon
`r delta_e/Ain` is also polynomial, of degrees at most (3,4), and its difference
from the prescribed normalized cubic bump is polynomial despite the unscaled
target's 1/r factor. Bernstein bounds therefore give a global representation
error bound, not merely a nodal discrepancy.

The actual projected fields receive their own global peak certificates. For k
the squared-peak subdivision tolerance is tightened to 1e−11 so that enclosure
width cannot masquerade as failure of the arithmetic peak reserve. The native
velocity bound includes the actual loaded streamfunction coefficient and all
projection errors. Helper evaluations are independently checked at interior
points; those point comparisons supplement, and do not replace, global bounds.

## Fixed amplitudes, support, norms and areas

The following conditions are prospective and precede probe results:

- Uref=1, velocity Euclidean peak .001; actual peak must differ by at most
  1e−8 relative. Global projected velocity representation error is at most
  1e−12 Uref.
- k/Kin peak .01, with global representation/peak arithmetic reserve 1e−12
  in reference units.
- Independent `delta_e/(Ain/r)` peak .01. Its actual projected peak may differ
  by at most 1e−5 reference units, i.e. 0.1% of the prescribed increment. Its
  global representation-error bound and exact RMS error are reported
  separately. This tolerance was explicitly selected before outputs; it is
  not permission to retune or shrink the bump.
- Domain RMS minima are 1e−5 for velocity and 1e−4 for each normalized scalar;
  support RMS minima are .0001 and .001, respectively.
- Main support lies within x∈[−2.5,2.5], r∈[.075,.45], has area at least 1.5,
  and spans at least four native intervals in x and r. The independent outer
  epsilon bump also spans at least four outer intervals. Its own support is
  x∈[−2.5,2.5], r∈[H,b], and its area is reported separately.
- Each field's half-absolute-peak area must have a certified lower bound at
  least .05. The independent epsilon area is outer-domain area only.

Squared field integrals are exact rational polynomial integrals over all
pieces; RMS values are outward-rounded square-root enclosures with the actual
full or outer domain area. No `max(1,RHS)` or support-dependent normalization is
used for domain RMS. The checker also reports exact polynomial gradient
energies and global gradient upper bounds. The epsilon gradient reported here
is explicitly `grad(r delta_e/Ain)`, not `grad(delta_e)/(Ain/r)`.

For half-peak areas, a box is certainly inside if its Bernstein lower bound
exceeds one quarter of the squared-peak upper bound, and certainly outside if
its upper bound is below one quarter of the squared-peak lower bound. All other
box area remains unresolved. Subdivision uses exact rational box measures.
Report nested area enclosures at unresolved-area fractions .01, .003, and .001
of support area, with subdivision counts and remaining boxes. Acceptance uses
the lower area, not an unrefined indicator-quadrature estimate. Contraction of
these enclosures is numerical area-resolution evidence, not a physical h-order
claim. The same procedure is applied to actual native projected fields.

Outside-support projected leakage is bounded globally and must be at most
1e−12 in the relevant reference units. It must also be below half the certified
peak before restricting the half-peak set integration to the declared support.
The distinction between exact compact target support and possible arithmetic
leakage of the projected field stays visible.

## Divergence, boundaries and dependent Gamma

The exact streamfunction identity gives `u_x+w_r=0`. The checker bounds the
actual projected polynomial divergence globally and independently assembles
every native Q2 pressure moment using Gauss 8×8 and the signed maps. Its own
sparse pressure mass matrix supplies the moment residual dual norm. Both the
strong divergence bound and weak dual norm must be at most 1e−13, i.e. 1e−10
of the prescribed velocity peak with Lref=1. Native D action is compared with
the independent action. Essential velocity/k/independent-epsilon increments
must be zero within the 1e−12 reference reserve; pressure remains exactly zero
in this increment probe.

For a deliberately synthetic constant Kbase=Kin only, the probe projects the
1% k seed, forms Knew from the actual native k coefficients, and applies the
sole endpoint-constrained essential lift `Gamma=A(Knew)/H`. The independent
epsilon bump still has zero H trace. The probe exports the independent bump,
old/new Gamma lifts, Gamma increment, and their total separately. Incompatible
left/H corner data must be rejected.

The checker independently reconstructs Gamma using fixed endpoint values and
L2 projection of only the two interior edge modes, comparing Gauss 24 and 40.
Its native coefficient error limit is 1e−10 in Ain/H units; projection-rule
change is at most 1e−12. It reports the pointwise target discrepancy as a sampled
diagnostic, distinct from coefficient compatibility. It verifies the 1% K
increment at x=0,H, reports the induced cutoff change, and requires the induced
relative Gamma target change at most 2%. Scalar composition and the synthetic
constant base are checked independently. The Gamma extension is not masked to
fake compact total-epsilon support. Actual evolved-reference K, admissibility,
current-data refresh, and fresh coupled histories remain later engine gates.

## Interfaces and evidence sequence

New files only:

1. `h27/perturbation_probe.cpp`
2. `h27/perturbation_certificate.py`
3. `h27/perturbation_check.py`

Native CLI:

```sh
WJH27PerturbationProbe case.xml reference.json request.json new-native.json
```

Layout request:

```json
{"mode":"layout","source_identity":"SOURCE_SHA256","reference_sha256":"TABLE_SHA256"}
```

For each shape, run the independently admitted mathematical certificate:

```sh
python3 -B h27/perturbation_certificate.py --mesh-native layout-native.json \
  --shape PrimaryQuartile --seconds ADMITTED_CAP --out primary-certificate.json
```

Project request (certificate path is explicit, and its complete-file hash is
required):

```json
{"mode":"project","source_identity":"SOURCE_SHA256","reference_sha256":"TABLE_SHA256",
 "certificate":"/absolute/path/primary-certificate.json","certificate_sha256":"CERT_SHA256",
 "velocity_sign":1,"k_sign":1,"epsilon_sign":1}
```

Then independently check the actual native projection:

```sh
python3 -B h27/perturbation_check.py --native primary-native.json \
  --certificate primary-certificate.json --table reference.json \
  --seconds ADMITTED_CAP --out primary-check.json
```

Repeat for secondary and required sign/finer-grid cases, retaining each frozen
layout and certificate. Every output is create-only. Native output uses
exclusive POSIX creation and fsync; Python writes with mode `x`. Source, mesh,
reference, request, certificate, and verification-code identities remain in
the evidence. All records keep `full_coupled_gate_pass=false` and
`developed_reference_qualified=false`.
