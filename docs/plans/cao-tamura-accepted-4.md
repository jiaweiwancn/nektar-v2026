# CT-Accepted-4: continuous neutral inlet and rough mean-profile fit

This package changes the declared approach-flow idealization. It uses the
original Natural A93, A95 and E94E95 closures with the Accepted budget AWF,
constant outer equilibrium k, and transported nonconstant inner k. CT1–CT3
finite-depth/Hermite results, failures and resource receipts retain their
original identities. They are regression evidence, not results for this inlet.

The implementation worktree started clean at the reviewed and fetched remote
head `50bd04c315283e996b53c360d8fa7a67877f02f3`. Source/Git work used the
existing Mac branch worktree. All builds and numerical results used the isolated
Linux experiment at
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-4-20261112`.
`findmnt` confirmed the active `/dev/sda2` mount at `新加卷1`; the shared
installation, other worktrees and running campaigns were preserved.

The [original contract](cao-tamura-accepted-4-evidence/contract.json),
[pre-native numerical supplement](cao-tamura-accepted-4-evidence/native-contract.json)
and [pre-evolution projection contract](cao-tamura-accepted-4-evidence/projection-contract.json)
record the decisions in order. The supplements preserve earlier failed
preparations and the original accuracy limits. The
[continuation supplement](cao-tamura-accepted-4-evidence/continuation-contract.json)
applies the requested conditional T=.02 gate to every admitted 6 mm pair,
including the P4 control; it was declared before inspecting that endpoint.

**Decision:** implementation and the nine central flat mathematical references
pass. All 20 timestep members finish T=.008 without a native domain rejection,
and their complete budgets, field, mass and timestep gates pass. All ten
model/height/space pairs fail the unchanged ground-traction accuracy limit.
No pair is admitted to T=.02. Selected degree sensitivity and future terrain
admission remain unresolved. This is bounded verification with an explicit
accuracy failure, not neutral-preservation or experimental validation.

## Physical input and calibration

Coordinates are physical height above the numerical flat wall. The numerical
wall remains at zero; neither the outer logarithm nor epsilon is displaced.
Hh=40 mm, rough Re_H=10500, and reported u_tau/U_inf=.051 are retained.
With length Hh and the original interpolated U(Hh) as scales,
nu=1/10500, U_inf=1.4562340769492705, and u_tau=.07426793792441279.
This Reynolds-matched scaling does not assign an unreported dimensional air
viscosity or velocity to the experiment.

| Model | q | Natural kappa | sigma_epsilon | K=q u_tau² | Fitted B=predicted U(Hh) | Unshifted z0 (mm) |
|---|---:|---:|---:|---:|---:|---:|
| A93 | 3.401265417397 | .428323628708 | 1.3 | .0187604501485 | .939137432047 | .177748710 |
| A95 | 3.822000426077 | .404060995299 | 1.3 | .0210811094289 | .947684430747 | .230598353 |
| E94E95 | 3.324506456314 | .412087471625 | 1.2 | .0183370687047 | .944745541279 | .211571661 |

q is evaluated from the specified roots and kappa from
kappa² q=sigma_epsilon Delta. The compatible wall constants are q^-2 and
kappa q^(3/2). No common .41, measured-RMS fit of K, finite-Re E-model
denominator or outer molecular diffusion is introduced.

The new CT4 mean-fit policy gives equal weight to frozen `f3a_R` point IDs
0–3. Their full-precision physical heights are 10.240861264, 15.311355771,
20.232034934 and 30.254565303 mm. The fit changes only the intercept B of
U=B+u_tau/kappa log(z/Hh). For each height a separate effective h_s then solves
the original finite-nu wall-resistance matching equation, with bracket
0<=h_s<H_A and an explicit root/branch check. These are model-conditioned
equivalent inputs, not newly measured surface properties.

| Model | h_s(5.5 mm) (mm) | h_s(6 mm) (mm) | h_s(6.5 mm) (mm) |
|---|---:|---:|---:|
| A93 | 3.509953 | 3.498719 | 3.489255 |
| A95 | 4.457119 | 4.438703 | 4.423218 |
| E94E95 | 4.219487 | 4.203862 | 4.190714 |

These reproduce the preliminary values from full-precision source coordinates.
The fitted U(Hh) remains below the measured normalization datum; it is not
renormalized to one. For IDs0–3 the signed velocity residuals, in U(Hh) units,
are A93 `[.042461,-.0000945,-.010027,-.032340]`, A95
`[.036822,-.001546,-.008576,-.026700]`, and E94E95
`[.038761,-.001047,-.009075,-.028639]`. All unused mean-point residuals remain
in the evidence. A fitted intercept is not an experimental-validation pass.

Joint extraction sensitivity uses the frozen ±.5 PDF-point envelope on both
axes of IDs0–5: all 4096 joint vertices, 256 deterministic Sobol interior
realizations and the central realization. Membership and equal fit weights
remain fixed. Each realization is reused in the fit and 40 mm interpolation,
including the shared ID3. The interpolation bracket changes from [3,4] to
[4,5] in some realizations. The evaluated h_s ranges are:

| Model | 5.5 mm | 6 mm | 6.5 mm |
|---|---:|---:|---:|
| A93 | 3.05694–3.97432 | 3.04761–3.96072 | 3.03974–3.94928 |
| A95 | 3.90212–5.03711 | 3.88716–5.01456 | 3.87457–4.99561 |
| E94E95 | 3.70266–4.75647 | 3.68989–4.73740 | 3.67914–4.72137 |

All evaluated roughness roots exist. These are evaluated deterministic
extraction envelopes, not statistical intervals or a proof of extrema over
the continuous input box. Experimental/probe-position uncertainty, approximate
Reynolds number, assumed displacement and roughness homogenization uncertainty
remain separate and unquantified. No hill observation enters this fit.
Every candidate interface is below the first visible rough mean observation;
the matching speed at H_A is therefore inferred from the stipulated logarithm,
not directly measured at 5.5–6.5 mm.

## Wall height, observation coordinates and mathematical references

For H_A=5.5/6/6.5 mm, H_A+ is 107.2243/116.9720/126.7197; using the paper's
assumed d_obs=3.5 mm gives (H_A-d_obs)+=38.9907/48.7383/58.4860.
The nominal physical cube is 5 mm high, leaving only .5/1/1.5 mm of
cube-top clearance. h_s, cube height, aerodynamic z0, assumed displacement and
numerical wall origin are distinct quantities. The interface's roughness-sublayer
adequacy is unestablished.

The compact cosine-squared hill has maximum |curvature|H_A of
.1085656/.1184353/.1283049. All fail the retained historical .1 flag and satisfy
the separately reported .25 flag. The latter is not a replacement acceptance
criterion or an error bound. Positive local normal Jacobians and the sufficient
global normal-map reach check do not establish curved equilibrium or stability.
No native curved execution is enabled.

The lowest frozen flat stress extraction bounds are about 3.849 mm for S and
9.121 mm for R, with the Fig.3 axis/displacement conventions retained. Rough
flat stress checks therefore lie above all three candidate interfaces. Fig.12
streamwise RMS traverses retain physical heights 2 mm (S) and 7 mm (R).
The observation audit transforms every applicable hill-profile location to a
closest actual normal distance, retaining lower/upper coordinate envelopes,
foot limits and multiple-root information. Vertical height is not substituted
for normal distance. Fig.11 has no inferred physical height because its
downstream boundary-layer depth is unspecified.

The independent finite-nu oracle splits integration at viscosity and
destruction breakpoints, finds candidate k extrema from Q0=0, and compares
48/96-order integration and the original native Accepted kernel. It uses
counterfactual zero-forcing production for the cutoff. Its central results are:

| Model | H_A (mm) | Lambda | z_c/H_A | minimum k_w/K | ground k_w/K | minimum flat inner covariance/K |
|---|---:|---:|---:|---:|---:|---:|
| A93 | 5.5 | 2.801239 | .165094 | 1 | 1.756531 | .378853 |
| A93 | 6.0 | 2.884734 | .151870 | 1 | 1.772754 | .378353 |
| A93 | 6.5 | 2.961793 | .140606 | 1 | 1.786970 | .377927 |
| A95 | 5.5 | 2.666509 | .188905 | 1 | 1.745033 | .410809 |
| A95 | 6.0 | 2.749829 | .173804 | 1 | 1.763400 | .410347 |
| A95 | 6.5 | 2.826736 | .160938 | 1 | 1.779522 | .409952 |
| E94E95 | 5.5 | 2.714000 | .180144 | 1 | 1.830974 | .372416 |
| E94E95 | 6.0 | 2.797379 | .165733 | 1 | 1.850473 | .371890 |
| E94E95 | 6.5 | 2.874338 | .153456 | 1 | 1.867578 | .371443 |

All nine central rough references are mathematically admissible. Minima occur
at the endpoint; near-endpoint numerical Q0 roots are retained in the receipts.
The integrated production/destruction and finite-nu mechanical-work identities
are verified separately. nu/nu_t(H) is .01825–.02261; sqrt(nu epsilon(H))/K is
.03656–.04525. These nonzero finite-Re indicators are reported without changing
the stipulated high-Re outer equations.

The inner parallel reference has equal normal variances 2k_w(z)/3. Outer
predicted RMS/u_tau in (streamwise, spanwise, vertical) order is
A93 `(1.618905,1.459005,1.432822)`, A95
`(1.893175,1.596246,1.229590)`, and E94E95
`(1.488737,1.488737,1.488737)`. The unused-stress receipt retains individual
experimental discrepancies; neither K nor h_s is refitted to conceal them.
Using all retained rough flat stress heights and their unchanged weights,
longitudinal/vertical RMS discrepancies are respectively .84351/.79869 (A93),
1.02644/.62013 (A95), and .77576/.84977 (E94E95), in u_tau units. These include
the outer layer where a constant-K idealization cannot reproduce the measured
decay. The separate Fig.3e duplicate mean check retains its original coordinate
conversion and normalization ambiguity; its RMS discrepancy is 1.47428,
1.18719 and 1.28568 in U/u_tau for A93, A95 and E94E95. It is not a fit datum.
These whole-profile checks evaluate the declared continuous idealization;
heights above the native 80 mm top are not native observations.

## Continuous provider and native representation

`NeutralCompositeV1` supplies initialization, live inlet/top data, absolute
outlet traction and probe exports from one immutable continuous calculation.
It uses `ReferenceK`, `ReferenceKGradient`, `ReferenceKFlux`,
`ReferenceKDiffusivity` and the original finite-nu velocity reconstruction.
The outer equations use constant K and epsilon=u_tau³/(kappa z); the inner
k is the integrated nonconstant budget reference. The current evolving
endpoint K still controls the actual wall kernel and its frozen prescribed
sink. The reference u_tau is not imposed as evolving ground traction.

Pi_w=0 and Pi_o=-Qzz. At +x the imposed total traction is
`(Qzz-Qxx,u_tau²)` outside and `(0,u_tau²)` inside. For the 6 mm references,
outer Pi is .001183289539 (A93), .005714899175 (A95), and zero (E94E95).
The outer normal outlet traction is -.003132237014, -.011429798349 and zero,
respectively; the shear component is .005515726604 for every model.
Independent closure evaluations verify these nonzero pressure offsets and
tractions. Scalar outlet advection and zero streamwise diffusive data remain;
there is no added pressure gauge or outlet pressure condition.

The provider, projection, model/policy, physical parameters, coordinate convention
and absolute-pressure convention are bound into the existing open identity.
Known outer FwdTrans endpoint ulps permit analytic continuation at the original
positive coordinate only, with explicit logging. No physical coordinate,
generating volume owner, quadrature weight or native side semantics changes.

The domain is x=[0,.125,.25], z=[0,2] in Hh units: two columns, 10 mm length
and an 80 mm compatible top truncation. It is not a tunnel-roof model, a
free-stream cap or a developed precursor. Inner knots include the predeclared
H fractions, h_s and cutoff; the pre-native supplement fixes outer grading
factor1.25. P3/P1 is primary; A95 6 mm also has P4/P2 on the same partition.
Basis points are 8x8. Native normal/tangential/mass/volume/interface orders are
32/8/10/10/10. Frozen budget controls are native, q48f16, q64f16 and q64f32,
with all four equations and both reconstruction records mandatory.

The first unconstrained FunctionsL2V1 preparations failed the declared traction
accuracy despite small field-value errors. A95 6 mm P3 gave ground error
.0114397 and interface RMS .00200281 against the unchanged .001 limits.
These inputs/histories were not evolved. A separately versioned
ConstrainedRadialV1 preparation uses the existing split radial mass projection
with initial value/one-sided gradient conditions. Its live inlet/top rows are
the same immutable projected reference trace in both data objects. Existing
dependent epsilon lifting precedes freezing the initial corner map. This
operation never restores an interior evolved state or repairs epsilon after k.
All new histories are fresh; the old projection history is rejected.
An earlier coordinate/side exception and secondary collection error are also
archived in the preparation receipts. The endpoint handling and failed-report
collection were corrected before constructing the admitted fresh histories.

## Native outcomes and limits

All nine model/height combinations have 1607 stored native slots on P3/P1;
the A95 6 mm P4/P2 control has 2839. Every count is measured, below the unchanged
3000-slot guard. Every same-space .002/.001 pair has identical actual initialized
coefficients. Driver/probe effective initial states and explicit/implicit live
boundary maps agree. All ten fresh initializations pass domain, source, field,
wall/interface traction and continuous-provider checks.

All nine independent continuous references also pass the native-test-space
weak balances for momentum, shared k, model-specific outer epsilon and
continuity at order64. For representative A95, scaled maximum free-row defects
are 2.05e-14 (momentum x), 2.72e-16 (momentum z), 1.79e-15 (shared k),
3.13e-16 (epsilon), and zero (continuity), against 2e-9. Interface rows are
counted once. The accepted inner higher source moments cancel with diffusion,
not separately. Wrong equilibrium data/provider, omitted nonlinear pressure
offset and constant-inner-K substitutions have discriminating negative tests.

These analytical results do not make the constrained polynomial reference a
discrete equilibrium. For A95 the initialized order64 native alpha-zero
mass-dual norms in momentum x/z, shared k and epsilon are
`[.122792,.00142068,.0202246,.000760459]` on P3 and
`[.0236343,.000340766,.00549357,.0000206931]` on P4. Raw unforced residuals,
initial value/gradient errors and actual evolved residuals remain separate.
The existing packet also exposes diagnostic-only `*_cancelled` variants; they
are not used in this qualification or in governing evolution. In particular,
they must not be read as a qualified open-boundary reference-cancellation route.

Both timestep members complete T=.008 in every row below. There are no rejected
constitutive steps; 120 accepted steps across these 20 trajectories are retained.
The table reports dt=.001; complete bodies include both members and all common
times. Errors are normalized by u_tau²; the limit for each is .001.

| Model / H_A / space | Mean wall-traction error | Full interface-traction RMS | T=.008 numerical gate |
|---|---:|---:|---|
| A93 / 5.5 / P3 | .0122033 | .00337508 | fail traction |
| A93 / 6 / P3 | .0143930 | .00536158 | fail traction |
| A93 / 6.5 / P3 | .0166160 | .00761121 | fail traction |
| A95 / 5.5 / P3 | .00688143 | .000729192 | fail ground traction |
| A95 / 6 / P3 | .00827326 | .00131265 | fail traction |
| A95 / 6.5 / P3 | .00972243 | .00227222 | fail traction |
| E94E95 / 5.5 / P3 | .00836245 | .00100347 | fail traction |
| E94E95 / 6 / P3 | .00999054 | .00188629 | fail traction |
| E94E95 / 6.5 / P3 | .0116720 | .00317234 | fail traction |
| A95 / 6 / P4 | .00146759 | .0000818372 | fail ground traction |

The solver's hard-domain checks run at the existing production stages. Accuracy
is the declared paired-milestone gate, not an added predictor rejection rule.
Read-only accepted-state reconstruction retains earlier accuracy crossings:
for A95 P3 dt=.001 the ground error is .00346912 at T=.001; for P4 it is
.000809643 at .001 and .00113757 at .002. The subsequent .008 states are failed
accuracy evidence, not qualification past these crossings. No continuation,
reference reset, changed sink, epsilon repair or bound correction is used.

The complete four-control/four-equation physical-budget classifier passes for
every member. Both reconstruction records, actual state/configuration/boundary
identities, finite positive scales and endpoint coverage are required. Actual
physical face gradients and frozen stage histories are used; prescribed-row
reactions do not replace fluxes. Sample A95 P3 dt=.001 q64f32 signed physical
defects are approximately -3.8e-7 (inner k), 2.4e-9 (outer k), 7.1e-8
(outer epsilon), and -4.3e-7 (whole k), on the recorded physical scales.
Quadrature/face spreads pass separately. Endpoint mass defects remain below 9e-16 on
the declared scale, against 1e-8. Native and independent physical inventories
retain positive k/epsilon, covariance and active-cutoff margins. A95 P3's
audited minimum radicand is about .691469; no domain-margin clipping is applied.
Weak continuity does not imply pointwise incompressibility: A95 dt=.001 has
sampled signed divergence ranges [-.00652990,.00544207] on P3 and
[-.00109372,.000987476] on P4. The largest reported outer isotropic production
work magnitudes are 2.18107e-6 and 1.60945e-7. The raw gradients and signed
divergence work remain in the actual stress and production; no gradient
projection is applied.

At T=.008 the primary 6 mm scaled reference errors are:

| Model / space, dt=.001 | U | w | k | outer epsilon | absolute Pi |
|---|---:|---:|---:|---:|---:|
| A93 P3 | 9.384e-6 | 3.758e-7 | 3.226e-5 | 3.617e-6 | 1.390e-6 |
| A95 P3 | 6.526e-6 | 1.963e-7 | 3.274e-5 | 3.824e-6 | 8.320e-7 |
| E94E95 P3 | 6.921e-6 | 2.330e-7 | 3.541e-5 | 3.803e-6 | 1.029e-6 |
| A95 P4 | 8.076e-7 | 2.270e-8 | 5.363e-6 | 2.683e-7 | 3.955e-8 |

These are physical weighted RMS, with fixed U(Hh), K, epsilon(H) and pressure
scales, not local relative errors. The limits are 1e-3, with w additionally
limited to 1e-4. Initial projection error and drift from the actual initial
numerical field are exported separately. All pairwise timestep comparisons
pass the1e-3 field gate; no temporal-order claim is made.

P4 materially improves the selected representation and wall traction but does
not clear its accuracy limit. A95's U-gradient physical RMS error drops from
.00851939 to .00111311 initially and from .00485832 to .000684273 at .008,
for P3 versus P4. Its dt=.001 evolved alpha-zero norms are
`[.00604764,.000352222,.00216118,.000302658]` and
`[.000553330,.0000386188,.000284671,.00000830092]`, respectively. These are
full native norms on different test spaces, not a common-test convergence rate.
The common physical-value contrast has resolved timestep uncertainty for U,
w, k and epsilon; Pi remains unresolved: maximum timestep difference
2.22203e-7 exceeds max(1e-7,.1*spatial difference)=1e-7. This bounded contrast
does not establish spatial convergence or admit a longer run.

Same-method uninterrupted/restart states and histories agree; complete rollback
and changed physical/provider/projection identity rejection pass. The old
open-boundary pressure-offset control gives coefficient defect4.86e-15 for
offset .125; internal reversed flow, outlet restriction, ordinary CLS/cache
and existing native diagnostic controls pass. The open route still rejects
its unsupported cache mode. Existing kernels and specifications are reused
only through their unchanged source hashes.

The evolution controller exits1 for the required accuracy failures. Analytical
audit and compact collection exit0 for completed diagnostics; these do not
override the failed native decision. The explicit conditional-selection receipt
rejects all four 6 mm continuations. No T=.02, T=.2, developed precursor or hill
case is claimed.

## Cost, resources and provenance

The numerical run uses the frozen `build-v3/runtime` compiled from
`123f31c974e0b75e551e43a5203749a89378fade`. Fresh evolution/frozen-budget
orchestration used `c1578b47372884a0a88736e92c3f94f27f693c80`; the read-only
all-reference/collection finish used `c42825aa6ef0ebaa269788dbdd6215e1342cda3a`.
The additional unused duplicate-mean audit used
`4164c2ccdacf56f5d3fd31d2c20bcc7516aa3e8c`; it changes no native source.
The final native C++ bytes are unchanged from the tested build source.

| Executable | SHA256 |
|---|---|
| IncNavierStokesSolver | `9f78fa07ed1eab1f1463bc3a7b546a039fb27d8195c757bd6577590c41bd6496` |
| RANSPreservationProbe | `9291343c45ea984eaa04132e5c9df52bc8b9662e772127a2768bf4bace6c4fc0` |
| RANSNewModelsProbe | `34dab79530e4e2f51b510d7d5709dbd727bd77b86da5c7a3f728cb6da0c6a669` |

Compiler flags, resolved and actually mapped library paths/hashes, Python and
NumPy/SciPy versions, actual mount, configurations, boundary maps, coefficient
vectors, checkpoint/native-history hashes and snapshot links are bound in the
compact source scope and full per-member receipts. Builds use two jobs;
native execution uses one serial rank of the MPI-enabled build and one-thread
numerical libraries. All new families use the owner-aware launcher, 16 project
including helpers and 72 aggregate. Fresh measured peaks are 14 project and 15
aggregate; all final resource receipts pass. Collector/observation/archive
finalizers are supplied separately, closing the collector's pending-own-finalizer
note. No suspension/allocation wait was required. Historical CT1–CT3 resource
failures remain unchanged.

Representative ordinary solver process costs, including startup/output but
excluding independent Python/probe audits, are:

| A95 6mm space / dt | accepted steps | process CPU (s) | elapsed (s) | peak RSS (KiB) | exclusive instrumented CPU/step (s) |
|---|---:|---:|---:|---:|---:|
| P3 / .002 | 4 | 19.918 | 20.484 | 227792 | 3.693 |
| P3 / .001 | 8 | 35.203 | 35.786 | 236040 | 3.758 |
| P4 / .002 | 4 | 25.239 | 25.788 | 333988 | 4.846 |
| P4 / .001 | 8 | 44.871 | 45.452 | 341960 | 4.883 |

The full evolution/analysis family costs 3323.665 process CPU seconds and 2590.274
elapsed seconds; final read-only weak checks/collection cost 361.199/294.777.
CPU includes owned helper processes and is distinct from elapsed time. Exclusive
timings are summed once, not nested inclusive totals. For A95 P3 dt=.001,
state/domain/constitutive work accounts for 21.322 of 30.062 exclusive CPU seconds;
factorization/solve validation accounts for .216 seconds. This remains a small
dense reference case, not a backend-throughput qualification.

## Smooth limitation and next gate

No smooth native case is admitted. At .5/1 mm the independent smooth columns
fail cutoff existence for all models. At 1.5 mm A93/A95 have negative reference
k while E94E95 still fails the cutoff. At 2 mm and 8 mm all have negative
reference k; at 8 mm minimum k/K is -.323602 (A93), -.057157 (A95), and
-.388142 (E94E95), despite Lambda>1. The unshifted mean-matching defects remain
visible. This diagnoses these references, not all possible finite-depth or
nonequilibrium smooth flows. Below the 2 mm RMS observations there is no
demonstrated comfortable high-Re smooth operating interval. No clipping,
redistribution, EP fallback or other repair is introduced.

The strict logarithmic/constant-outer-K profile is an explicit approach-flow
idealization. It is not the complete measured boundary layer. Curved evolution,
EP, developed precursors, hill predictions and backend expansion remain outside
this package. CT1–CT3 and the eight historical unresolved epsilon-order
observations are unchanged.

## Reproduction and evidence map

The [compact decision](cao-tamura-accepted-4-evidence/compact-v1/decision.json)
and complete per-case bodies separate mathematical, native accuracy,
spatial uncertainty, diagnostic completion and future-terrain decisions.
[Seals](cao-tamura-accepted-4-evidence/compact-v1/completed-artifacts.json)
hash the exported files; the
[export map](cao-tamura-accepted-4-evidence/compact-v1/export-provenance.json)
links them to full Linux originals. Full matrices, states, initializations,
failed preparations, histories and profile packets remain outside Git.
Compact stage inventories retain extrema with times, coordinates, native sides
and geometry, plus full-array hashes. The observation supplement is separately
sealed in `observation-checks-v2`; it adds Fig.3e checks without altering the
older compact collection or its numerical decisions.

For a fresh isolated reproduction, synchronize this committed source to a new
Linux `CT4_ROOT/source`, retaining the recorded CT2 dependency/runtime libraries.
Use the same active mount and one numerical family at a time. The scientific
commands are below; use new destinations, not the sealed experiment folders.
The builder records the new executable and loaded-library identities. A new
build is a reproduction, not reuse of the pinned runtime's receipts.

```bash
export CT4_ROOT=/media/jay/新加卷1/nektar-v2026-runs/ct4-reproduce-NEW
export CT2_ROOT=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-2-20261112
export S="$CT4_ROOT/source/solvers/IncNavierStokesSolver/Tests/RANS"
export D="$CT4_ROOT/source/docs/plans/cao-tamura-validation-data"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$CT4_ROOT/tmp"
mkdir -p "$CT4_ROOT/tmp"
launch() {
  CT4_LABEL="$1"
  shift
  python3 -B "$S/ct_accepted3_command.py" --control "$CT4_ROOT/launch-$CT4_LABEL" -- "$@"
}
launch build python3 -B "$S/ct_accepted4_build.py" --root "$CT4_ROOT" --ct2 "$CT2_ROOT" --output "$CT4_ROOT/build-v3"
export RT="$CT4_ROOT/build-v3/runtime" LD_LIBRARY_PATH="$CT4_ROOT/build-v3/runtime/lib"
launch reference python3 -B "$S/ct_accepted4_reference.py" --data "$D" --contract "$CT4_ROOT/source/docs/plans/cao-tamura-accepted-4-evidence/contract.json" --native-probe "$RT/RANSNewModelsProbe" --output "$CT4_ROOT/reference-v1"
launch initial python3 -B "$S/ct_accepted4_jobs.py" --phase initial --root "$CT4_ROOT" --runtime "$RT" --initial-root "$CT4_ROOT/initial-v3" --output "$CT4_ROOT/initial-controller-v1"
launch controls python3 -B "$S/ct_accepted4_controls.py" --root "$CT4_ROOT" --ct2 "$CT2_ROOT" --runtime "$RT" --initial-root "$CT4_ROOT/initial-v3" --output "$CT4_ROOT/controls-v1"
launch evolve python3 -B "$S/ct_accepted4_jobs.py" --phase evolve --root "$CT4_ROOT" --runtime "$RT" --initial-root "$CT4_ROOT/initial-v3" --evolution-root "$CT4_ROOT/evolution-v1" --output "$CT4_ROOT/evolution-controller-v1"
launch weak python3 -B "$S/ct_accepted4_all_weak.py" --runtime "$RT" --initial-root "$CT4_ROOT/initial-v3" --output "$CT4_ROOT/all-reference-weak-v1"
launch analysis python3 -B "$S/ct_accepted4_analysis.py" --root "$CT4_ROOT" --runtime "$RT" --initial-root "$CT4_ROOT/initial-v3" --evolution-root "$CT4_ROOT/evolution-v1" --require-evolved --output "$CT4_ROOT/evolved-analysis-v1"
launch observations python3 -B "$S/ct_accepted4_observations.py" --references "$CT4_ROOT/reference-v1/results.json" --data "$D" --output "$CT4_ROOT/observation-checks-v2"
```

The evolve command is expected to return 1 for the retained accuracy failures;
run the read-only audits separately and preserve that exit. The immutable
requests/resource execution receipts record the exact original build,
incremental-build, controls, collector and conditional-continuation commands,
including the preceding failed preparation versions. Reproduction of the old
projection-history negative uses its archived CT4 `initial-v2` artifact; it is
not synthesized by rewriting a checkpoint signature.

The next bounded numerical gate is the ground-layer representation: hold these
continuous functions and fits fixed, localize the assembled momentum load, and
test one near-ground radial refinement of the P4 reference under the 3000-slot
guard. The degree contrast supports that hypothesis but does not yet prove
the cause or qualify a repaired trajectory. Do not extend duration to compensate
for the failed traction gate.

A subsequent curved Accepted package can independently implement and verify
the metric-weighted column/provider and conservative Cartesian momentum/scalar
coupling on a small manufactured curved patch, including one-sided full traction,
shared k flux, generating ownership and absolute pressure. Keep production
curved evolution guarded until that work, the flat traction/representation gate,
the retained .1 thin-layer decision and physical roughness-sublayer/interface
assessment are resolved. S/R preparation and measured complete approach-flow
adequacy remain separate prerequisites for a developed precursor or hill
prediction. No experimental validation is admitted here.
