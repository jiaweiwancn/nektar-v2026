# CT-Accepted-8: displaced A95 reference

Work starts at `6f545f53299915657a18368ac18879a7fda935e5`, clean and equal
to the fetched remote. CT7 and all earlier sealed evidence remain immutable.
The [contract](cao-tamura-accepted-8-evidence/contract.json) precedes native
execution. A95 physical height7.5mm is explicitly authorized for this package;
the earlier below7mm constraint remains part of its historical contract.

`SamePhysicalWallDisplacedAcceptedV1` is an explicit extension of the Accepted
appendix. The physical wall remains at r=0. It changes the prescribed sink and
its dependent trace consistently with aerodynamic displacement, while retaining
the original finite-viscosity momentum resistance, transported inner k and
one-sweep freezing. Its native admission is restricted to original Natural A95,
serial flat2D. Metric columns are read-only kernel checks, not terrain admission.

The central displaced reference is mathematically admissible and the native
implementation passes the independent kernel/operator and history controls.
Both fresh physical histories complete T=.008. Domain, traction, field-drift,
mass and complete physical scalar budgets pass, but the absolute-pressure
timestep-uncertainty gate fails. **T=.02 is not admitted or executed.** This
is a bounded completion with failed preservation uncertainty, not a qualified
precursor. No flat result replaces the CT6 joined-map temporal/degree failures.

## Input and explicit specification extension

The [original contract](cao-tamura-accepted-8-evidence/contract.json) fixes the
three physical-height-selected Fig.3a_R IDs0,1,2, equal weights, d=3.5mm,
Hh=40mm, u_tau/U_inf=.051 and nu*=1/10500. Membership is selected at physical
z<=25mm before displacement is subtracted. Point3 and the four-point CT7 fit
remain held-out/historical data. The measured40mm normalization is retained;
the fitted U(40mm) is not reset to one. Neither k nor constitutive coefficients
are fitted to RMS observations.

The full-precision [physical reference manifest](cao-tamura-accepted-8-evidence/physical-reference-manifest.json)
copies the authoritative fitted inputs and their dataset/specification hashes.
The [sealed compact evidence](cao-tamura-accepted-8-evidence/README.md) and
[closed publication scope](cao-tamura-accepted-8-evidence/publication-v1/scope.json)
separate completion from admission. The inherited independent oracle retains
its historical `...CandidateV1` label; native signatures/configuration use
the implemented `SamePhysicalWallDisplacedAcceptedV1`, not that oracle label.

For zeta=z-d, the intercept is fitted in consistent dimensionless coordinates:
`C = mean(U_j - u_tau/kappa log(zeta_j/Hh))`,
`z0/Hh = exp(-kappa C/u_tau)`. The inner finite-nu resistance equation is then
solved once for the central fixed surface parameter h_s, with 0<=h_s<H_A.

| Central A95 quantity | Independently calculated value |
|---|---:|
| Natural q | 3.8220004260766363 |
| Natural kappa | .404060995298521 |
| u_tau in retained normalization | .07426793792441279 |
| K0=q u_tau² | .021081109428869052 |
| Fitted intercept C | .9919449895413243 |
| z0 | .18124943074282322mm |
| h_s | 6.844974060108776mm |
| Fitted U(40mm)/measured U(40mm) | .9751145933293871 |
| Physical H_A | 7.5mm |

All64 joint coordinate/velocity extraction corners and the centre are evaluated
with fixed membership, u_tau, d, nu and normalization. Eight corners have no
admissible resistance root; they are retained failures. The evaluated z0 range
is .15511893426191492--.2086512192524378mm. This deterministic extraction
sensitivity is neither a confidence interval nor experimental/probe-position
uncertainty. The central candidate is not claimed robust to every corner.

The new policy is an explicit extension of the Accepted appendix, with the
physical wall remaining at r=0. It keeps the original finite-nu resistance,
roughness correlation and counterfactual p0. For actual sampled K it uses
`t0=K/q`, `A=(K/q)^(3/2)/kappa` and `E=A/(max(r,c)-d)`. The derivative of the
flat normalized budget is `-c/(c-d)^2<0`; its value at c=H is H/(H-d) and it
diverges as c approaches d from above. Thus the declared budget inequality
gives a unique cutoff, including c=H at equality. For positive J the analogous
derivative is `-integral_0^c J dr/(c-d)^2<0`. The native solve uses this
monotonicity, a representably positive bracket, and explicit rejection. It
does not clip an inadmissible budget. At d=0 the existing arithmetic path is
retained and checked against the original policy.

The shifted sink is finite at r=d because c>d. Its trace is A/(H-d); the
nonconstant reference k follows Q0 and D_A, not a constant-K substitution.
K0 fixes initialization/inlet data; the actual evolving endpoint K continues
to enter wall viscosity, cutoff, sink and dependent epsilon. h_s is fixed.

| Central column check | Result |
|---|---:|
| Cutoff c | 6.125936352851381mm |
| P0/A | approximately2.753714718 |
| Minimum admissible P0/A | 1.875 |
| Reference k(0)/K0 | 3.5544801988435704 |
| Minimum reference k/K0 | 1, at H |
| Epsilon(H) | .01013811394236224 |
| Epsilon plateau | .01544304595403272 |
| Integrated sink/production defect | 4.34e-19 |
| Finite-nu mechanical-work defect | -4.34e-19 |
| Upper k flux | 1.30e-18 |
| Minimum flat inner covariance/K0 | .40925154516634393 |
| Outer covariance minimum | .006111572556576726 |

Split independent integration includes roughness/viscosity/cutoff breaks and
stationary extrema. Signed hill-curvature metric-column controls pass, as do
budget-boundary and invalid displacement/height/roughness/metric/budget/model
controls. These are frozen-column results, not native curved admission.

## Native implementation and boundary identity

`EquilibriumAWFProfiles` owns the named shifted policy, metric budget,
source moments, Q/reference-k integrals and dependent dissipation. The patch
propagates it to actual sampled-K profiles and fingerprints the sink policy
and displacement. Native physical admission is deliberately restricted to
Natural A95, serial flat2D, Accepted one-sweep, cacheOff and the qualified
Richards--Norris boundary. Other model/geometry/action requests remain guarded.

`RANSWallEvolutionOpenBoundary` uses the distinct immutable
`DisplacedNeutralCompositeV1` for initialization, inlet, gradients, top and
full outlet traction. Its identity also binds the fit policy, physical floor,
zeta convention, all scales, H, h_s and displacement. The original
`NeutralCompositeV1` still rejects nonzero displacement.

The stored pressure is Pi_inner=0 and Pi_outer=-Qzz; centrally
Qxx=.005714899174565924, Qzz=-Qxx, Qyy=Qxz=0. Neutral Qxz=0 is separate from
the complete Rxz=-u_tau². The outlet uses `(Qzz-Qxx,u_tau²)` outside and
`(0,u_tau²)` inside. The RN roof has essential w only, natural tangential
traction +u_tau², g_k=0 and prescribed physical
`g_epsilon=-u_tau^4/(sigma_epsilon*(Z-d))`. This datum is independent of
evolving diffusivity. No top pressure/normal-traction condition or mean gauge
is added. Explicit and implicit live boundaries agree.

Complete one-sided interface values, fluxes and full normal/shear traction
are checked. The high-Re outer molecular-diffusion policy, finite inner nu,
raw-gradient signed production, full3D covariance, stage order and frozen
prescribed sink are unchanged. Corrections and extra governing solves are zero.

## Native static qualification and retained failed preparations

Resolved539-slot P3/P1 and945-slot P4/P2 fixtures exercise actual displaced
production rows, independent Nektar basis/derivative/face/mass checks, all-mode
top loads, full pressure action/rank, interface terms and frozen matrices.
Maximum raw all-row disagreements are2.78e-17 and1.82e-17; the largest frozen
matrix disagreement is6.94e-17. The continuous native-test balance maximum is
3.04e-14. Pressure singular ratios .00567/.00138 exceed the1e-10 criterion;
absolute traction fixes the pressure level without an extra null mode.

The [static export amendment](cao-tamura-accepted-8-evidence/static-verification-amendment-v1.json)
raises three read-only matrix/tensor packet limits from500 to1000. It does not
change any governing matrix, arithmetic, equation or the3,000-slot evolution
guard. The amendment's initially predicted539/945 counts are subsequently
confirmed natively. Failed preflights are retained: old diagnostic export
guards; an unaligned cutoff's independent face integration defect; and a
compressed-cell projection with negative covariance. The latter is a genuine
finite-element representation failure, despite the admissible continuous
column. The final resolved fixtures pass without physical modification.

Same-method uninterrupted/restart equivalence, complete rejected-step rollback,
stale fit-history rejection, wrong model/provider/cache rejection, pressure
offset invariance and state-independent prescribed roof flux pass. Constant
inner-k and omitted pressure-offset controls fail for their intended reasons.

## Fixed full-height short strip

The physical strip retains x=[0,.125,.25]Hh and roof20Hh (800mm). Isotropic
P4/P2 Modified_A fields and the fixed breakpoint-resolved radial mesh have
**2,957 stored slots**: u765,w765,k765,outer-epsilon441,inner-pressure95,
outer-pressure125 and one excluded unused gauge slot. Physical geometry,
degree/points and live boundaries are recorded in the preflight and cases.
The two fresh zero-step native histories have identical actual coefficients.
The shared continuous input is also compared read-only with P3/P1 on the same
mesh; its different projection is reported separately, not as evolved spatial
convergence or a measured degree order.

Initial P4 scaled reference RMS errors are U1.77302e-5, k4.06928e-7,
outer-epsilon1.05667e-5 and Pi2.85613e-16; w is zero. P3--P4 static physical
differences are U1.54033e-4, k5.24964e-6, outer-epsilon6.34212e-5 and Pi/w
zero. These are representation errors/differences with fixed physical input,
not subsequent drift. Initial physical ground/interface traction errors are
6.29e-16/6.53e-15, while the separately reconstructed outlet mismatch is
9.76e-5. Natural outlet data are enforced weakly, not by pointwise resetting.

The reference volume-average speed is1.3539582080070294; the short-strip
bulk flow-through time is .18464380844368136. T=.008 covers .04332666265622494
of that transit (4/8 accepted steps). It is more than a two-step smoke but is
neither a full transit nor a developed precursor.

Both physical histories complete with all accepted-state accuracy checks,
native/independent domain checks, signed production, covariance, cutoff,
mass, frozen matrix/reconstruction and complete scalar-budget checks passing.
Budgets retain exactly native,q48f16,q64f16,q64f32 and inner_k,outer_k,
outer_epsilon,whole_k, with k/epsilon reconstruction and independent physical
face fluxes. Prescribed-row reactions are not substituted for physical fluxes.

| T=.008 endpoint, scaled physical quantity | dt=.002 | dt=.001 |
|---|---:|---:|
| U drift | 2.63947e-6 | 2.74485e-6 |
| w drift | 5.56040e-7 | 5.72809e-7 |
| k drift | 2.94000e-7 | 3.01355e-7 |
| Outer epsilon drift | 1.43159e-6 | 1.47106e-6 |
| Absolute Pi drift | 8.15427e-6 | 8.41188e-6 |
| Mean ground-traction error | 9.09376e-5 | 9.03650e-5 |
| Full interface-traction RMS error | 7.46358e-4 | 7.73743e-4 |
| Net mass defect | 9.31269e-15 | 1.11655e-14 |

Absolute endpoint U reference errors are1.66747e-5/1.66432e-5, k errors
5.02014e-7/5.06357e-7 and outer-epsilon errors1.06692e-5/1.06751e-5.
Pi reference error agrees with drift within roundoff. Generic solver-console
L2 entries are not used for these independent provider/reference comparisons.
All complete-budget
controls pass the unchanged .02 criterion; scaled signed-defect quadrature
spreads are at most1.06e-14. Actual total defects and reconstructed physical
faces remain visible, rather than only reporting this integration agreement.
Native/independent constitutive inventories remain admissible. The fine final
inventory still reports raw divergence up to9.81884e-4 and isotropic production
work up to9.03524e-6: weak continuity is not asserted to imply pointwise zero
divergence. No gradient projection is applied.

The unchanged timestep classifier recomputes physical RMS primitives:

| Field | Scaled timestep spread | Allowed spread | Decision |
|---|---:|---:|---|
| U | 1.1901146269884219e-7 | 2.7448525483885866e-7 | pass |
| w | 2.073090188406101e-8 | 1e-7 | pass |
| k | 8.606596415770304e-9 | 1e-7 | pass |
| outer epsilon | 4.089372266715992e-8 | 1.47106301050867e-7 | pass |
| absolute Pi | 4.804457393280823e-6 | 8.411880056796802e-7 | **fail** |

Traction-specific temporal spreads also pass: ground profile1.26385e-6,
signed mean5.72587e-7 and full interface3.93728e-5, against respective
limits9.19539e-6,9.09376e-6 and4.16455e-5. The interface uncertainty margin
is relatively narrow. Admission uses absolute pressure; removing a mean,
resetting the field or relaxing the relative criterion is not permitted.
The failed pair is sealed and the driver exits nonzero. NoT=.02 continuation
or additional timestep/mesh search is performed.

Read-only pressure offsets, variation, region means and actual stored flow/
epsilon/k stages are exported separately. These quantities diagnose the
failure; they do not change the absolute-pressure metric. There is no evidence
here that a closure change or epsilon retuning would be justified.

The final pressure difference has physical RMS2.65000735e-8, mean3.34349e-9
and variation RMS2.62883e-8. It is therefore not a removable common gauge
offset. Inner/outer RMS differences are2.54597e-7/9.77072e-9; their mean
offset difference is2.40150e-7 in magnitude. The full-height norm weights the
thin inner region lightly, so these regional quantities are retained.

The reference pressure projection is already at approximately1.3e-18 RMS.
The first corrected-flow pressures depart by8.77845e-8/1.01768e-7 RMS at
their respective first step times, then fall to4.49767e-8/4.63976e-8 atT=.008.
Momentum predictor pressure is unchanged. Epsilon and k stages also leave
the corrected pressure coefficients unchanged while changing the constitutive
traction slightly. First-step full-G/splitting/linear identity defects are
of order1e-21--1e-23, and boundary defects are zero. Frozen momentum defects
2.57e-8/1.90e-8 are the recorded one-sweep splitting defects, not failed linear
solves. Refreshed epsilon-trace lags decrease from1.24e-8/8.61e-9 at the first
step to3.23e-9/9.55e-10 at the endpoint; they are measured, not repaired.
These observations localize a transient spatial/pressure-coupling uncertainty,
but do not uniquely separate projection, quadrature and splitting causes.
The reference is steady, so retiming its pressure cannot remove this error.

## Observations, approximation and retained temporal evidence

The corrected sealed CT5 closest-normal audit supplies all796 observations
without changing source data, extraction envelopes or roles. At7.5mm every
Fig12_RR7mm row is an inner-AWF/coupled comparison; the lowest Fig5b_RR
envelopes at x/Hh=+/-1.25 overlap the interface. Actual transported-state
flat-strip stresses at7mm are exported. They are not predictions at experimental
hill stations outside this strip. Fig11 physical heights remain unavailable
where delta7.5 is unspecified.

The physical cube height, assumed displacement, fitted z0, effective h_s and
numerical wall origin remain distinct. H_A clears the5mm cube tops by2.5mm;
this alone does not establish roughness-sublayer adequacy. High-Re/metric
indicators and the separate inner-auxiliary/outer molecular-viscosity ratios
are recorded. The historical curvature-height .1 screen remains exceeded;
positive metric does not prove curved-model accuracy. Smooth preparation is
still unresolved and no other-model displaced physical route is enabled.

Numerically H_A+=146.2150, (H_A-d)+=77.9813,
nu/nut_outer(H)=.0317367 and nu/nut_inner_auxiliary(H)=.0164247.
The two viscosity indicators are different quantities. The maximum hill
curvature-height product is .1480441 and its interface stretch bound is
.8519559; the retained .1 screen fails despite positive stretch.

Retained CT6/7 mass-norm difference vectors are evaluated without evolving or
altering them. U/w ratios2.38977/2.35335 accompany alignments .97113/.92859
and relative norms of a-2b .27372/.37936. Their absolute manufactured errors
decrease, but those sizable second differences and rotations leave a concrete
finite-step/nonasymptotic risk. They do not demonstrate a violated discrete
identity or establish asymptotic first order. k/epsilon ratios are near two;
Pi ratio is1.85629. The old upper-band failures remain unchanged. CT8's flat
pair does not close the joined-map spatial/temporal question.

## Capacity and next implementation

The new paired terrain/flat geometry inventory fixes H_A/Hh=.1875,
x/Hh=[-25,30] and roof20Hh with an explicit radial grid. It checks both16/32
point geometry inventories, positive collar/join-to-horizontal-roof geometry
and actual flat native counts **33,267(P3/P1)/60,907(P4/P2)**. These are geometry/count-only sessions with
the displaced/RN boundary identity recorded separately from production
admission. They exceed the unchanged dense limit; no full-domain evolution
is launched. The experimental hill's pure normal map cannot reach this roof.
The joined map still needs the corresponding coupled numerical qualification.

The analytic `C1CollarRoofV1` certificate has collar0.3Hh, transition end1Hh,
minimum interval Jacobian .259323 and minimum x_q .571062. Native polynomial
geometry is independently sampled at16/32 orders; its extrema/volume remain
separate from the analytic certificate. Terrain field counts are not measured
by the current affine inventory (its region classifier assumes physical y=H).
The paired logical topology is exported, but a mapped region-aware field
inventory belongs with coupled map qualification. The measured flat counts
alone already block the required full-domain evolution; no curved count or
coupled operator qualification is inferred from them.

The [CT9 design](cao-tamura-accepted-8-evidence/ct9-backend-design.md) audits
startup masses/projections, frozen volume/face assembly, consistent-free-mass
inverses, predictor/scalar solves, pressure-Schur action, mandatory feasibility
and diagnostic mass operations, plus snapshots/replay. Installed PETSc is
**OFF**; source adapters are not evidence of available runtime support. CT9
must retain `D[(alpha*M_free)^-1(Gp)]-Cp`, absolute pressure and every physical
check, with no hidden dense conversion or mass lumping.

The next package may implement and qualify that separately versioned serial
backend on the small dense references. Before physical full-domain preservation,
it must also close the retained absolute-pressure timestep uncertainty under
a prospectively declared bounded study. Backend equivalence alone cannot
admit preservation. Terrain production, smooth preparation, physical approach-
layer fidelity, roughness-sublayer accuracy and curved thin-layer accuracy
remain separate gates. EP and full variable-height evolution stay deferred.

## Provenance and reproduction

Authoritative work is isolated under
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-8-20261114/` on Linux.
The mount, source/configuration hashes, executable and actual loaded-library
identities, complete commands, resource receipts and heavy archive locations
are recorded with the compact evidence. Shared installations and all historical
archives are unchanged. Mac source/Git work uses the dedicated implementation
worktree. Tailscale SSH was used after the USB address timed out.

The governing implementation begins at4889a8f6; the final runtime is built
at e6bbe66a, whose later changes only bound read-only verification exports.
The physical controller is ba86ac0c. Subsequent diagnostic/reporting commits
do not change governing C++ or the sealed physical histories. Exact full
revisions and per-file hashes distinguish these identities.

The runtime is Release GCC11.4.0 against the installed MPI-enabled SDK,
executed with one serial rank and one-thread numerical libraries. Actual
mapped-library hashes are captured during execution; a Git hash alone is
not used as binary provenance. Solver SHA256 is
`a95fc6037946ce088185337fbf92d3d143a1e12b0dfa2d62ed591bf2f9fcdfa1`.

The two physical solver executions take31.379/54.163 elapsed seconds and
28.683/50.479 process-CPU seconds, with peak RSS399,240/437,044KiB. Those
figures include their startup/output, not the separate quadrature/probe work.
Exclusive phase timings and all family CPU/elapsed/allocation-wait values are
preserved. The fresh old curved-history regression costs721.55 elapsed seconds;
it is compatibility work, not part of the physical strip cost. No nested
inclusive timing totals are added.

Closed numerical/build families pass resources with high waters13 project
and20 aggregate workers, including the two control helpers. Collector and
publication verification are separately closed and recorded. Builds use two
compiler jobs; all numerical families run sequentially with library threads1.
Historical resource failures are not replaced by these fresh compliant receipts.
The five focused CTests (independent model oracle plus CT6/7/8 checks) pass.
Fresh CT5 A95 states and history payload, unchanged curved restart and complete
rejected-step rollback are bitwise compatible.

Run each command through `ct_accepted3_command.py --control NEW -- COMMAND`
with all numerical libraries restricted to one thread, the isolated runtime's
lib directory and TMPDIR, as in the recorded launch receipts. The reproducible
sequence is `ct_accepted8_reference.py`, `ct_accepted8_native.py --mode preflight`,
`ct_accepted8_controls.py`, then `ct_accepted8_native.py --mode advance --target
.008` with its immutable reference/contract/preflight/compatibility inputs.
The latter must exit1 for this sealed pressure-uncertainty failure. Independent
completion uses `ct_accepted8_pressure.py`, `ct_accepted8_assess.py` and the
focused CT6/CT5 regression. `ct_accepted8_collect.py` recomputes primitive
classifications and preserves separate mathematical, numerical, resource and
physical-applicability outcomes. Missing members/equations, foreign identity,
incomplete endpoint, false matrix/temporal labels and incompatible histories
are negative controls, not alternative admissible physical cases.

For a fresh Linux reproduction with the pinned runtime, run
`bash docs/plans/cao-tamura-accepted-8-evidence/reproduce-strip.sh NEW_ABSOLUTE_DIRECTORY`.
It rejects an existing destination and stops on the expected failed pair gate;
it does not attempt .02 or change a history's timestep. The SDK dependency
audit separately records Release/MPI=ON/PETSc=OFF with its cache hash and closed
worker receipt. Python3.10.12, NumPy1.21.5 and SciPy1.8.0 are recorded in the
publication environment. Final source changes after the diagnostic receipt
are documentation and collector wording only, with no numerical reclassification.
