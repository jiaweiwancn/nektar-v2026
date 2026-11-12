# Cao–Tamura budget AWF H2b

2026-11-12, Asia/Shanghai. The composite solver path is implemented and its
native operator checks pass. **The proposed hill initialization fails the
retained covariance gate at T=0. No physical hill step was accepted.** The
failure is present before the two-dimensional projection, and is not explained
by a side-inconsistent audit or k negativity. The saved candidate is retained.

The bounded post-startup flat triplet completed 70 steps. Its differences show
approximately first-order scaling, but no tested larger timestep passes the
predeclared combined absolute/response-relative screen. No long campaign was
launched.

## Configuration and provenance

The starting local/fetched head was `5cbaed815ff54d46312af657663aa4f346a5fe40`.
Source/Git work stayed on `codex/cls-epawf-pimple-one-sweep` in the existing Mac
worktree. Linux supplied all compiled, numerical and timing results. Exact
implementation revision, executable/library closures and file hashes are in
[source revisions](cls-cao-tamura-budget-h2b-evidence/source-revisions.json).
The implementation commit is `be5e856b5ba9d4708540d455b19b4a04c3e83ee7`;
all 109 relevant native source/header files match the Linux export. Final
solver SHA256 is `2a4fc727c6273548d1b77c05f8fb95aa49233344b29662ea174cdac1e8f5aac9`;
probe SHA256 is `92e891aea510e37112f9524a91f73de405262578afa7e1b849a4b230375125ac`.

The H1 finite completion is unchanged: profile SHA256
`9375601441e2c70fb87cf8f31de305afd4b67fd2f82a750f9274f6bb68eed697`.
The frozen Cao–Tamura dataset was not edited or redigitized. Fig.3d and every
hill observation remain held out; no speed-up or reattachment score is reported.

The selected problem retains Regularized CLS/Linear inner stress,
beta=.4, Gamma=1, gamma0=0, Hw=z_K=z_m=.015 m, endpoint K=k(Gamma_H),
h_eff=.000878311859293 m, nu=3.80952380952e-6 m²/s and
u_star=.0742679379242 m/s. Hh=.04 m, L=.1 m, cube height b=.005 m,
aerodynamic z0=.0002 m and displacement d0=.0035 m remain distinct quantities.
The actual inner source is P_actual-E_C. The cutoff uses the counterfactual
flat parallel p_A0(K); neither the finite-reference production nor terrain
production replaces it. The floor, reserve, source signs, pressure spaces,
bound controls and one-sweep schedule are unchanged.

The data volume was remounted after a server reboot. Its verified UUID is
`5AB8A669B8A642FD`; the actual mount is `/media/jay/新加卷1`, while historical
reports retain `/media/jay/新加卷`. No checkpoint bytes, signatures or histories
were edited. Heavy evidence is under:

```text
/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2b-20261112
```

[Artifact index](cls-cao-tamura-budget-h2b-evidence/artifact-index.json) records
payload paths and hashes. Historical H1/H2a and neutral runtimes were preserved.
The index covers 3,653 files, about 5.58 GB, including frozen dependency copies.
The isolated Release build used four compiler workers, rank one and one
computational-library thread per native process; at most four native jobs ran
concurrently, within the aggregate 72-worker limit. Flags retain
`-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS`.

## Transactions and finite diagnostic IVPs

Rollback restores coupled coefficients, native histories, stage time/alpha,
physical samples and accepted diagnostic vectors before fallible packet output.
The pending writer is transferred and cleared before invocation. The original
rejection and packet-output exception are recorded separately. A failed writer
cannot cause repeat flushing. The rollback receipt also retains the actual
attempted callback history and its time/alpha, separately from restored history.

Focused native tests passed for an unselected hard rejection plus writer
failure, and for an explicit-callback rejection **after candidate assignment**
plus writer failure. The latter changes the assigned coupled-state hash before
rollback. Both restore all tested accepted fields/history/diagnostics exactly.
Normal selected and unselected packet paths also pass.

`FiniteDevelopingIMEX1V1` is a separate affine diagnostic IVP contract. It
validates the original parent timestep/signature/checksum, preserves the complete
coupled seed, absolute pressures, solved epsilon, fixed reference/defect and
physical starting time, and forms each integrator's own IMEX1 histories without
Lift or reprojection. Same-dt imported seed/reference/history and first-step
governing vectors are bitwise equivalent to ordinary restart. Ordinary
changed-dt restart still rejects. Packet cadence now uses the explicit native
attempt count; a nonzero-time epoch gives selections `[true,true,false]` at
interval two, independently of time*alpha.

See [verification](cls-cao-tamura-budget-h2b-evidence/verification.json).
The focused Linux `RANSKBoundQPChecks` and `RANSCaoTamuraBudgetH1Checks`
CTest targets both pass, with test parallelism one.

## Composite geometry and native operators

`RANSWallGeometry=CompositeHillV1` selects the constant normal strip below Hw
and vertical extrusion above it. Native outer zHat is converted through
eta=(zHat-Hw)/(Z-Hw); its Jacobian column is divided by Z-Hw and inverse-gradient
row multiplied by Z-Hw. The normal map is not extended to the roof.

Physical maps, native-element derivatives, quadrature/masses, face cofactors,
source rays, endpoint functionals, pressure/roughness action, residuals and
readbacks use the composite interpretation. Outer CLS/Yap uses certified
nearest-floor distance; immutable FiniteInflowV1 uses Cartesian incoming height.
Neither query uses outer eta as wall distance. Native region/element and
one-sided hill-foot provenance select derivatives and constitutive law.

Here Hw*max|curvature|=.296088. Positive x_H', positive Z-z_H and nested normal
level graphs establish injectivity and a common interface without gaps/overlap.
The relative AWF normal Jacobian and base surface measure each enter once.
No extra curvature force is added to Cartesian weak stress divergence.

The mode is restricted to the finite/open Regularized/Linear accepted-budget
one-sweep pilot. Unrelated legacy mapped guards remain. Geometry/chart,
boundary policy and physical parameters participate in strict restart identity.
Native GetCoords overshot a fitted endpoint by at most 1.11e-16 m; a bounded
16-epsilon geometric endpoint normalization handles that roundoff. It changes
no transported value and rejects genuinely outside coordinates.

The 6x13 pilot has **2918 coupled slots**: cubic velocity/scalars and linear
regional pressures. Tangential knots are [-.2,-.1,-.05,0,.05,.1,.6] m, with
explicit feet and crest. It preserves the H2a vertical partition, including
.25/.4 completion panels. This is an operator/startup mesh, not experimental
resolution. New manifests report actual topology/dimensions; earlier inherited
template labels are preserved as historical metadata, with XML/native dimensions
authoritative.

The open policy keeps the composite inlet, no-slip ground/natural k, impermeable
full-stress slip top/natural scalars, floor-only Yap distance and conservative
upwind/backflow outlet. The absolute **total traction** datum fixes pressure;
there is no closed-box mean gauge or top net-flux adjustment. With modified
pi=p_mechanical+2k/3, sigma=-R+(2k/3-pi)I+nu*S (nu inside only), the outlet
combines un+*q+un-*q_seed and Pout*n. Its existing backflow contribution is not
added twice. The inlet/interface epsilon corner uses E_C(Hw), while gamma's
immutable reference table retains its different raw epsilon value.

| Native check | Observed result |
|---|---:|
| Curved minimum sampled physical Jacobian | .7039133322 |
| Map/inverse maximum error | 8.94e-15 m |
| Independent physical-gradient FD error | 3.34e-9 |
| Volume, GL16/32 | .636 m², agreement about 2e-15 |
| Interface position / opposing area-vector errors | 0 / 2.22e-16 |
| Full pressure action: matrix / native FD errors | 5.38e-17 / 2.53e-12 |
| Native residual/component-sum error | 2.82e-16 |
| All-mode frozen selected-stress/scalar flux errors | at most 1.91e-16 |
| Nonzero outlet datum Pout=.013 | traction integral .0104 |
| Outflow/backflow independent flux errors | at most 4.45e-16 |

GL8/16/32 checks refine both directions. Independent source integration converts
physical finite-profile knots back to each regional vertical chart. Flat
polynomial exactness is not assumed for curved metrics.

The shared inventory includes mapped native/source/reporting rows, actual and
dependent epsilon separately, endpoint K, both radial sides, one-sided foot
profiles, moving source points and continuous cubic ground extrema. A negative
outlet-ground k fixture repairs with one KKT solve using the physical consistent
mass; full inventory and ground feasibility pass, stationarity/dual defects
are zero at reported precision, and complementarity is about 9.5e-31. Its
raw-centered objective, 2.57826e-15, agrees with independent physical GL32
integration. Non-k and essential coefficients remain bitwise unchanged. A
negative actual H-plus epsilon is detected despite positive dependent E_C.
This protects the declared inventory/continuous ground, not every interior point.

In the composite flat limit, mass and pressure actions agree with affine values;
full spatial-load error is at most 3.83e-14. Two physical steps reach T=1e-5 s.
Maximum coefficient difference is 7.64e-15; physical pressure L2 difference is
6.57e-17 and traction difference zero. A 1+1 split is bitwise equivalent;
post-assignment rejection restores state/history and incompatible geometry
restart rejects. These flat-limit steps do not count as nonzero-hill stepping.
Each accepted flat step retains one U, pressure, outer-epsilon and k solve and
one logical bound call. Solved epsilon coefficients remain authoritative after k.

[Geometry evidence](cls-cao-tamura-budget-h2b-evidence/geometry.json).

## The actual hill gate

The seed uses the native radial completion, a physical streamfunction lift
matching the strip, and a Hermite continuation with no terrain lift in the
constant freestream tail. Physical consistent-mass projection enforces inlet,
ground, matching and tail constraints; a setup Schur projection enforces native
weak continuity. Setup algebra is separately counted: three radial and two
field projections, six projection mass solves, four fixed-constraint mass
solves, one continuity mass factor and one Schur solve (91 independent rows).
It performs zero governing solves and no physical advance.

The resulting seed hash is
`7adda6c42bd9d6e8fae1c5259055dbd60290888593c0fbc626d027b39b2203e4`.
Its native weak continuity maximum is 1.73e-16; strong divergence L2 is .1113
inside and .2882 outside. The native solver and diagnostic fixture give the
same complete coupled vector. Scalars and cutoff budgets pass, including a
budget ratio above 3.4, but the hard covariance audit rejects:

| Hill-side witness at x=-.1 m | Before 2D projection | Native seed |
|---|---:|---:|
| Inner, r=.015, geometry19 | -.00235109 | -.00545013 |
| Outer, physical z=.09, geometry37 | -.05149109 | -.06532722 |

The pre-2D column uses the **native 1D radial representation**, not analytical
gradients substituted into a discrete residual. Step-halved differentiation
errors are below 1.6e-8 in these gradients. Independent coefficient contraction
and tensor eigendecomposition reproduce the native witnesses. Flat-side
neighbors at the same physical points remain positive; no gradients are averaged.

An additional analytical inner check removes the radial projection error. The
exact incoming AWF radial gradient is 9.38948/s; following the wall through the
strip adds dW/dx=21.53198/s at the foot. With the unchanged local k and wall
viscosity, the Linear covariance minimum is still **-.00137022**. Thus this
chosen tangential-following seed is already incompatible with the selected
covariance domain. Representation worsens it. At the native outer witness,
regularization contributes -.0709279 to the critical eigenvalue, versus
+ .0088384 from isotropy and the smaller Original-CLS terms.

This is a limitation of the proposed initialization under the selected law,
not proof that every hill state is inadmissible or that ellipticity is lost.
No floor, gamma, source, covariance tolerance or physical height was changed.
The initial rejection, actual coefficients/reference/histories and witnesses
remain archived. The requested 8/16-step hill pair was **not executed**.

[Seed evidence](cls-cao-tamura-budget-h2b-evidence/seed.json).

## Representation, timestep and cost

The zero-step flat control adds the positive viscosity-intercept knot
3.82864e-5 m and removes .5 m from the initially constant tail, retaining .25/.4.
The rapid molecular/linear-eddy scale is 9.76599e-5 m; neither is h_eff or z0.

| Initial representation | Retained mesh | Targeted split |
|---|---:|---:|
| Slots | 2000 | 1961 |
| Sampled traction | .00625786 (+13.45%) | .00566078 (+2.63%) |
| Velocity L2 error | 2.14601e-4 | 2.14485e-4 |
| Velocity-gradient L2 error | .530546 | .393778 |
| Inferred weighted ground reaction | .00557170 | .00558825 |
| Integrated actual inner production | .001266889 | .001266545 |

The reaction uses the same Q2(s)*max(0,1-r/.00025) lifting, zero at side corners,
with free-row work explicitly retained. Its behavior does not justify selecting
a mesh from sampled traction alone. Native scalar margins remain positive;
GL16/32 source integrals agree closely. No changed-mesh trajectory was launched.
[Representation evidence](cls-cao-tamura-budget-h2b-evidence/representation.json).

The local flat IVPs start from the actual H2a T=.0008 checkpoint SHA256
`7485f1b78243f6cb09d7406f3c685caa8779cf868c5826f9d8f9aca025dd6605`.
dt=4e-5/2e-5/1e-5 s complete elapsed .0002/.0004, with zero correction. Absolute
pressure and the complete Hw-to-.8 outer-epsilon domain are retained. At the
last observation, observed orders are .978 (velocity), .980 (k), 1.005 (epsilon),
1.022 (pressure) and .959 (mean traction). Fine-step Richardson estimates are
about 4.04e-8, 5.13e-9, 5.90e-9 and 7.89e-8 in dimensional field norms, and
1.13e-8 for traction. All absolute normalized targets pass. The separate 1%
response-relative traction screen fails even at the finest level (about 12.6%
at the earlier observation and 4.7% at the later one). No larger dt is adopted.
This local experiment does not recover accumulated T0 error or qualify a
developed precursor. Actual inlet bulk speed is 1.38360 m/s: domain and
inlet-to-observation flow-through times are .28910 and .14455 s. The local
elapsed horizon is only .001384 domain flow-through times; physical t* moves
from .02 to .03. Hw/u_star=.20197 s.
[Temporal evidence](cls-cao-tamura-budget-h2b-evidence/temporal.json).

Three alternating eight-step H2a cache pairs preserve governing stages,
coupled fields and native histories bitwise. Median end-to-end speedup is
1.097: approximately 36.1–36.3 s uncached versus 33.0–33.3 s with scoped reuse.
The 32-step block takes 111.4 s, peak RSS about 885 MiB. Exclusive CPU totals
include 46.1 s scalar inventories, 18.9 s constitutive observations, 6.23 s
residuals and 6.16 s physical audits. Factorizations/applications remain small.

Fixed anchor unions and immutable surface views are shared without dropping
observations. Exact scope-local constitutive reuse distinguishes geometry,
element, side and actual/dependent epsilon. Static nearest-wall certificates
are retained; additional queries have bounded storage. A certified chord
envelope reduces representative distance searches from 12k–60k boxes to 23–26,
at the unchanged 1e-10 m tolerance and 200k interval limit. Independent stationary
roots and multiple-footpoint tests pass. The old combined curved diagnostic
timed out at 1800 s; this is an infrastructure limit, not a numerical rejection.

Removing eager error-string formatting on successful mapped queries reduces a
matched complete curved diagnostic from 66.2 to 44.5 s and composite-flat
two-step execution from 67.0 to 35.7 s. These are single matched replays;
numerical actions/audits/correction results and accepted stages/histories are
bitwise unchanged. Curved diagnostic peak RSS is about 1802 MiB. No actual
curved physical-step cost is claimed. At the measured flat 3.5 s/step, another
640 steps at retained dt=5e-6 would cost roughly 37 minutes to reach t*=.1;
this forecast did not authorize or launch a continuation.
[Performance evidence](cls-cao-tamura-budget-h2b-evidence/performance.json).

## Reproduction and next gate

Use unused output roots. Activate the existing Linux MPI environment and set
OMP/OPENBLAS/MKL/BLIS/VECLIB/NUMEXPR thread counts to one. The following paths
refer to the isolated run root above:

```bash
python3 source/solvers/IncNavierStokesSolver/Tests/RANS/run_cao_tamura_budget_h2b.py \
  --phase curved-operators --solver curved-lazy-final/runtime/IncNavierStokesSolver \
  --probe curved-lazy-final/probe-runtime/RANSWallResolutionProbe \
  --h1 ../cls-cao-tamura-budget-h1-20261111/compatibility-resolved --output NEW_OPERATORS

python3 source/solvers/IncNavierStokesSolver/Tests/RANS/run_cao_tamura_budget_h2b_flat.py \
  --validation-only --solver rollback-context-final/runtime/IncNavierStokesSolver \
  --probe curved-lazy-final/probe-runtime/RANSWallResolutionProbe \
  --h1 ../cls-cao-tamura-budget-h1-20261111/compatibility-resolved \
  --parent ../cls-cao-tamura-budget-h2a-20261111/startup-v2/primary-step160/checkpoint.16 \
  --output NEW_TRANSACTIONS
```

The artifact index retains every actual command, session, runtime and native
payload. Source collection scripts reproduce the compact reductions. An operator
probe's zero exit code means its diagnostic protocol completed; its separate
`accepted_awf_hard_admissible=0` still blocks hill evolution.

The independent P2 campaign was read once: stopped at step7232/T=4.52 after
three confirmed complete-field windows ending at step7200/T=4.5. Its checkpoint
SHA256 is `b9d7e139840b757030996028564031d951cc538637fd65d4a0119dc0938f243e`.
Runtime, cap and stopping policy were unchanged. No H2b campaign is left evolving.

The next implementation priority is a **covariance-admissible mapped seed with
interface velocity allowed to adjust**, checked before two-dimensional
projection and then under native continuity and all hard audits. Keep the
selected model/completion fixed. Only after that gate passes should the bounded
8/16-step hill pair execute. Curved restart/rollback and finite-step sensitivity
remain unqualified until an admissible nonzero-hill leg runs. Developed flat-flow
calibration, affordable larger assembly/storage, spatial/temporal resolution and
held-out hill validation remain separate gates; raising the 3000-slot limit is
not a scaling solution.
