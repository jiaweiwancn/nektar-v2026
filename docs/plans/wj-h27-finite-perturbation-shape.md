# H27 finite perturbation shape V1 — prospective construction

This specifies source shapes and a normalization/evaluation procedure before
finite-pilot outcomes. It does not report normalization, support norms, native
admissibility, or successful perturbation evolution. Actual native 700-cell
radial knots have not been invented or evaluated here. No Mac numerical work
was performed. The helper changes no engine, application data, model, histories
or pressure field.

## Frozen native-edge layout

Use actual common tensor-product x/r edge layers, with `r=z-d0`. The caller
must establish their correspondence to every actual element in the support;
the union of unrelated local edge coordinates is insufficient. The application
domain remains `x∈[-10,15]`, `r∈[0,19.9125]`, with H=.15 and d0=.0875.

The fixed streamwise spline knots are `[-2.5,-1.25,0,1.25,2.5]`. Every one must
already be a native edge; finer meshes may add edges. No nx=20 or 700-cell
topology is inferred from a filename or assumed here.

For the initial selection, take the smallest/largest actual native r edges
inside `[.075,.45]` as a,b. Reject unless `a<H<b`, `b-a≥.3`, and at least five
distinct radial edges are available. Choose three interior knots minimizing
the sum of squared physical fractional distances to `[.25,.5,.75]`; ties retain
the lexicographically first native triple. This is `PrimaryQuartile`.
`SecondaryLowerSkew` instead targets `[.2,.4,.65]` and excludes the primary
triple. Reject if no distinct triple exists. These choices use mesh geometry,
not pilot stability or cost results.

Select the independent epsilon peak at the actual outer native edge nearest
`H+.5(b-H)` for primary or `H+.65(b-H)` for secondary, with an interior-edge
requirement and lower-coordinate tie break. H, a, b and this peak are native
edges. The constructor requires at least four native intervals in x, in the
cross-H normal support, and in the outer epsilon normal support. The last check
makes the same prospective resolution minimum explicit for that scalar bump.
Unsupported layouts require mesh/support review; they are not repaired by
reducing amplitudes or inventing edge coordinates.

Serialize and freeze the selected complete layout. A finer mesh must contain
these SAME physical knots; do not rerun selection and compare different
perturbations. Main support area is exactly `5(b-a)≥1.5` in code units. The
independent epsilon support area is `5(b-H)`. The support is strictly interior
to inlet/ground/top/outlet, so all exterior essential traces remain unchanged.

## Exactly representable streamfunction velocity

Let Nx, Nr be the single degree-three Cox–de Boor B-splines on the respective
five distinct knots. Their support is the interval from first to last knot;
they and their first two derivatives vanish at both endpoints. Define

```
psi = C Nx(x) Nr(r)
deltaU = ( C Nx Nr', 0, -C Nx' Nr )
grad(deltaU)_xx = C Nx' Nr'
grad(deltaU)_xz = C Nx Nr''
grad(deltaU)_zx = -C Nx'' Nr
grad(deltaU)_zz = -C Nx' Nr'.
```

Consequently divergence is identically zero, including the exact cancellation
of the two diagonal derivative expressions. The velocity is C1 and compact;
its nonzero pieces have degrees `(3,2)` and `(2,3)`, both within tensor Q3 on
every native cell. Skipping intermediate native edges does not harm this
property, since polynomial breaks still occur only on selected native edges.
No coordinate-gradient averaging or numerical Helmholtz projection is needed.

For every native Q2 pressure test, integration by parts gives
`∫q div(deltaU)=0`; internal traces cancel and the exterior velocity is zero.
The actual signed native coefficient projection and assembled D action must
independently verify this identity to arithmetic/representation error.

`UnitVelocity` means C=1 with positive sign. Its global Euclidean vector peak
must be determined on Linux before constructing histories. A normalization
record contains the exact frozen layout, positive lower/estimate/upper peak
bounds, and an immutable evidence SHA256. Relative bracket width must be
`≤1e-8`; this is a prospective numerical construction tolerance, not a measured
physical uncertainty. The runner must verify the record identity and evidence;
the C++ container only checks its binding and numerical invariants.

Load this record once, then set `C=sign_velocity*1e-3/peak_estimate`.
The helper refuses normalized velocity evaluation without it and rejects any
second normalization load. The unit streamfunction coefficient, knots, signs,
peak bracket and final C are reviewable through the API. No scalar or velocity
amplitude is tuned after viewing pilot outcomes.

A suitable independent global peak procedure works on all 4×4 spline knot
rectangles. There `|UnitVelocity|²=Nx² Nr'²+Nx'² Nr²` is a polynomial of maximum
bi-degree (6,6). Search interior stationary roots and each edge polynomial,
then bound the global supremum with Bernstein/interval subdivision, including
all endpoints and one-sided cell owners, until the prescribed bracket closes.
A sampled local optimizer alone must not be described as a global bound.
Native projected-field peak error must be added separately to the bracket.

## Scalars and the dependent interface

Let `fx=1.5 Nx`; the fixed uniform x knots give `max fx=fx(0)=1`. For any native
edges `l<p<u`, define a C1 Hermite bump which is zero outside [l,u], rises as
`3s²−2s³`, `s=(r-l)/(p-l)`, and falls with the same polynomial at
`s=(u-r)/(u-p)`. Its peak is one and its derivatives vanish at l,p,u.

Set `fr=Hermite(a,H,b)` and

```
delta k = sign_k * .01 Kin * fx fr.
```

This is piecewise Q3, has zero value/normal derivative at support endpoints,
and `deltaK(0,H)=sign_k*.01 Kin`. For a developed base with Kbase≠Kin, the
relative interface change is `.01 Kin/Kbase`, not automatically one percent.

Set `fe=Hermite(H,epsilon_peak,b)` and independently prescribe

```
delta epsilon_ind = sign_e * .01 (Ain/r) fx fe, r>H,
```

with zero extension otherwise. Its value and first derivative vanish at H
and b. The support closure touches H but the independent trace is exactly
zero. This target is rational in r; unlike the velocity and k seed, it is
NOT claimed to be exactly Q3. Its actual native projected normalized peak,
support, RMS and interface trace must be measured. Failure is a construction
or resolution failure, not permission for amplitude reduction.

After perturbing k, the engine recomputes current K, c, auxiliary viscosity
and Ew, then applies the sole dependent epsilon-H lift. Use the exact ratio
`Ew_new(H)/Ew_base(H)=(Knew/Kbase)^(3/2)`. Do not replace that lift by the
independent 1% epsilon bump or force its amplitude to 1%. Require induced
relative interface change ≤2%, and record total epsilon separately. A state
leaving the qualified K interval is rejected, not clipped. Dependent column
changes can extend below a even though the prescribed velocity/k support is
compact; they must not be masked to make reported norms smaller.

The exact ratio describes the physical dependent trace target before native
edge projection. A general developed Q3 K trace raised to the power 3/2 is
not itself Q3. Preserve the engine's declared endpoint-constrained current
Gamma projection, verify its coefficients and corner values, and report its
pointwise target error separately. Do not silently demand an incompatible
second epsilon trace or hide the projection error in the independent bump.

Signs for velocity, k and independent epsilon are separately ±1. An opposite
signed seed changes no support or normalization. The secondary shape changes
the radial velocity spline and potentially the outer epsilon peak, while
retaining the same prescribed k/interface seed. It needs its own bound record.

## Prospective qualification measurements

Measure the actual native increments at common physical coordinates after
projection and dependent trace reconstruction. Use code-area integrals and
fixed reference scales, with the true full-domain/outer-domain measures.

| Increment | Peak | Whole relevant-domain RMS minimum | Support RMS minimum |
|---|---:|---:|---:|
| Euclidean delta velocity / Uref | .001 | 1e-5 | .0001 |
| delta k / Kin | .01 | 1e-4 | .001 |
| independent delta epsilon / (Ain/r) | .01 | 1e-4 (outer domain) | .001 |

Use the prospective peak and independent normalization/representation error,
never a pilot-dependent tolerance. For all three quantities require area of
the half-absolute-peak superlevel set ≥.05; epsilon's area is measured outside
H. For velocity, incremental kinetic energy per domain area in Uref units is
at least `0.5*(1e-5)^2=5e-11`. Report component norms and all gradient changes.

The unprojected velocity energy integral is separable:

```
integral |deltaU|² = C² [integral Nx² * integral Nr'²
                      + integral Nx'² * integral Nr²].
```

Four-point Gauss on every spline knot interval integrates these polynomial
moments exactly; the same order suffices for unprojected k and normalized
epsilon shape moments split at their Hermite breaks. For native projected
epsilon, `(r delta_e/Ain)²` can have degree eight in r, so use at least
five-point physical Gauss per native interval and independently refine it.
Half-peak set area needs level-set/interval refinement; nonzero nodal count
or one unrefined indicator quadrature is insufficient. Include both sides of
all shape and native element breaks in gradient/MAX observations.

Require actual velocity gradients, k, independent epsilon, interior K and
cutoff changes to exceed ten times independently measured representation noise.
Admissibility and all boundary/current-data checks apply to the final seeded
state everywhere, including the dependent wall columns. Keep pressure
unchanged by this helper. Construct fresh, consistent coupled time levels and
explicit histories in the engine only after seed admission; preserve the
unperturbed developed state and immutable qualification record.

These support, amplitude, RMS and resolution criteria are prospective
engineering choices inherited from the H27 audit. They are not measured
disturbance spectra or experimental sensor uncertainty. No perturbation gate
is closed by this source artifact alone.
