# H2h: fixed-velocity standard-model scalar feasibility

13 November 2026, Asia/Shanghai. **The fixed standard-B velocity cannot be made
covariance-admissible in the declared scalar space and bounds.** A checked LP
dual bound establishes a sampled zero-margin obstruction. Independently, a
negative H-plus tensor has identically zero response to every permitted scalar
direction. The terminal native and fresh GL24 audits fail. No certificate was
issued; physical hill attempts and accepted steps are both **zero**.

This completes the bounded H2h experiment. The dedicated calibrated-standard
versus regularized-CLS comparison remains open for admissible initialization,
common flat spatial qualification, developed flow, temporal/spatial sensitivity
and held-out experimental assessment.

## Preserved problem and exact state

The branch/worktree is `codex/cls-epawf-pimple-one-sweep` at
`.worktrees/pimple-one-sweep`. Fetch matched reviewed head
`f30c309d6944dc7a117f14d6d7e2b929bbcf0d54`; the H2g implementation is
`07b680966481d331d1d1d727d118e425829b3f51`. Other worktrees, the shared Linux
checkout's local changes, historical jobs and all frozen runtimes were preserved.

Validated implementation through: `58c67776299bb2dfa9e5f97771dee43e90902715`. The
[commit binding](cls-cao-tamura-budget-h2h-evidence/commit-binding.json) confirms
that its recorded source blobs match the Linux test manifest.

The [plan](cls-cao-tamura-budget-h2h-plan.md),
[predeclared policy](cls-cao-tamura-budget-h2h-evidence/policy.json),
[comparison manifest](cls-cao-tamura-budget-h2h-evidence/comparison-manifest.json),
[context verification](cls-cao-tamura-budget-h2h-evidence/context.json),
[optimizer result](cls-cao-tamura-budget-h2h-evidence/optimizer.json),
[active constraints and independent tensor replay](cls-cao-tamura-budget-h2h-evidence/active-constraints.json),
[native result](cls-cao-tamura-budget-h2h-evidence/native.json),
[focused verification](cls-cao-tamura-budget-h2h-evidence/verification.json) and
[source/runtime identities](cls-cao-tamura-budget-h2h-evidence/source-revisions.json)
are the compact record. [The artifact index](cls-cao-tamura-budget-h2h-evidence/artifact-index.json)
locates heavy data in:

```
/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2h-20261113
```

UUID `5AB8A669B8A642FD` was verified at that mount before work and at final
collection. Native compilation and conclusions are from Linux, with four build
workers, one native rank and computational-library threads fixed to one. The
aggregate remained below 72; other project work was preserved. USB SSH was
unavailable; the established Tailscale fallback was used.

| Preserved quantity | Deployed CLS | Calibrated standard |
|---|---:|---:|
| H_A=z_K=z_m, m | .0065 | .0065 |
| u_tau, m/s; nu, m²/s | .07426793792419734; 3.8095238095238094e-6 | same |
| kappa; fitted z0_eff, m | .4405539051866587; .00012222264395414748 | same |
| q; K0, m²/s² | 3.1906750816242875; .017598891430879275 | 3.3333333333333335; .01838575534504143 |
| sigmaEpsilon | 1.3 | 1.3478315512167742 |
| Physical nonlinear / regularization / Yap | retained; beta=.4, Gamma=1, gamma0=0 | all absent; constant Cmu=.09 |
| Deployed h_s, m | .0022472751006336437 | .0022035111852050995 |
| Wall cMu; cLength | .09822789532821126; 2.5108624752197968 | .09; 2.6811256851969465 |
| Solved flat cutoff, m | .0008312858196541334 | .0008312858196541342 |
| Mesh / native slots | exact retained common mesh / 2804 | same |

Ce1=1.44, Ce2=1.92 and sigmaK=1 remain common. The full-precision H2g calibration
manifest is retained by hash; this package performs no refit. The common neutral
mean/epsilon, unshifted physical incoming-height convention, inlet/top/outlet
roles, absolute pressure convention, Linear inner stress, inner molecular
viscosity and budget-cutoff source are unchanged. Hill observations remain held
out. The physical roughness-sublayer and upper-profile qualifications remain open.

Verified immutable state hashes:

- Standard cold: `038fc076cebbea41bcb7cc349147eacc2b938a9051827bbb386baf7c3cb3716c`.
- Archived standard B: `c83da50b65d2e9980325d365779a9d53394fbf10e6431d1b87f79728366ae13d`.
- Terminal phase-I state: `f7ad3a6403f45398ae5da899f812d9a1848d9f9374e0ab211241dc8e3e7b3e25`.

Every trial starts as an exact B copy. Only permitted outer-k/epsilon rows are
updated using the saved eight-plus-eight directions and original conditioning.
Velocity and both pressures are **bitwise unchanged**, including signed zeros.
All inner k, endpoint K, prescribed/interface scalars, dependent epsilon and
reference/defect data retain their identities. Scalar envelopes refer to cold,
not B. These are initializer restrictions; production evolution has no new clamp.

## Native scalar context and implementation

`RANSWallScalarFeasibilityContext.h` adds a probe-only `FixedBScalarContextV1`.
It evaluates actual B gradients directly, rather than cold plus reconstructed
velocity increments. Binary64 arrays retain all 450871 native memberships and
31200 GL20 volume points, with element, geometry, radial/tangential sides,
actual/dependent epsilon, family and physical coordinates. Named blocks identify
B/cold scalars, the sixteen responses and conservative centring-error estimates.
The reader checks the exact native receipt, model, array and saved-basis hashes.
GL24 is absent from this context and from optimization.

Three legal nonzero retractions toward the old scalar origin give maximum
relative affine/native error **1.90e-15** and zero gradient change. The independent
standard minimum-eigenvalue formula agrees with native tensors to **1.13e-17**.
The complete protected-profile hash is unchanged:
`5009b2e0ec55d25b567713b94cd8f0cf62d4c25716ce8f2a69d55c213eaf6bc7`.
Native contracts also check frozen rows, endpoint functionals, dependent traces
and observation membership. No evolving reference or covariance repair is used.

For fixed G, S=G+Gᵀ and s=lambda_max(S), the selected standard law is
R=(2/3)k I−.09 k² S/epsilon. Thus zero margin is equivalent to
(2/3)epsilon−.09 s k>=0. Positive margin delta=1e-4 K0 uses the convex function
.09 s k²/epsilon−(2/3)k+delta. Supporting cuts use its exact first derivatives;
nonpositive-s rows are checked separately using positive scalar lower bounds.
No extra factor of two or deviatoric projection is introduced.

The HiGHS dual-simplex LP uses signed variables and hard original [.5,2] scalar
bounds. Scalar rows are normalized by cold values; zero-margin covariance by
cold epsilon; positive-margin cuts by K0. Row scaling is recorded in the sparse
dual replay. The conditioned coefficient box +/- .740806494509108 is **implied**
by independent physical scalar rows, not an added amplitude restriction. Their
condition number is 3.106 and inverse residual is 8.39e-16.

Phase I minimizes covariance slack alone. Distance penalties cannot influence
its obstruction diagnosis. Positive-margin cut generation and unbuffered
positive-margin obstruction checks are implemented and tested, but unnecessary
on the hill data because even zero margin fails. No secondary objective ran.

## Checked obstruction and its physical location

One LP with 720147 rows and 17 variables completes in 8.06 s. It has 653828 hard
scalar rows and 66319 positive-s covariance memberships. The checked results are:

| Phase-I quantity | Value |
|---|---:|
| Primal slack | .16835786057497626 |
| Dual objective | .16835786057497631 |
| Primal–dual lower-bound gap | 3.67e-14 |
| Maximum scaled primal violation | 5.20e-17 |
| Dual stationarity infinity norm | 1.15e-15 |
| Conservative native/data rounding allowance | 6.12e-12 |
| Verified lower bound after allowance | **.168357860568823** |

Slack is dimensionless: a common allowance on
(.09 s k−(2/3)epsilon)/epsilon_cold. It is not a covariance eigenvalue. Seventeen
nonzero dual rows independently reproduce the lower bound. They combine k
lower/upper bounds, epsilon upper bounds and covariance rows near the slope,
upstream blend, crest/outer transition, foot and roof. The GL20 constraint set
contributes active scalar bounds; the search is not roof-only. The sparse rows,
multipliers, physical locations and their cancellation are saved for replay.

There is a simpler independently checked obstruction. At physical
(x,z)=(.23549580558569316,.0065), outer element/geometry41, H-plus:

- k=K0=.01838575534504143 and epsilon=.1430512592910451 are protected.
- Every allowed scalar response is exactly zero.
- dU/dz=57.92565906673285/s, nu_t=.00021267369553732568 m²/s.
- Covariance eigenvalues are **[-.00007280610118924, .0122571702300, .0245900616023]**.

The same physical point on inner element35, H-minus, has dU/dz=57.48728001135/s,
auxiliary nu_t=.0002129398143720893 m²/s and minimum **+1.83857535726e-6**.
Tangential derivatives agree to roundoff; normal derivatives differ. These are
separate native evaluations, with no averaging or side substitution. Independent
isotropic/linear tensor contractions, eigenvectors and residuals below 1e-18 are
retained. Actual and dependent outer epsilon agree there and are recorded
separately. There are 18 negative fixed-response memberships, eight distinct
native sided/context evaluations. Enlarging scalar freedom while retaining
these fixed interface values cannot repair this particular B velocity.

If unrelated scalar values could be chosen independently within [.5,2] at each
point, the optimistic zero-margin test is positive. That relaxation ignores
interface freezing and shared basis coupling; it does not rescue this contract.
The result proves neither global physical infeasibility nor loss of ellipticity.

## Terminal audits and comparison

The phase-I candidate minimizes the declared normalized slack, not the worst
stress eigenvalue; its roof minimum is worse than B. B also had small GL20
scalar-envelope violations, now explicitly visible. No rejected candidate is
called an improved or accepted physical seed.

| State / set | Inner minimum | Outer minimum | Outer minimum/K0 |
|---|---:|---:|---:|
| CLS retained H2f B, native | +1.759889e-6 | −.0298040470 | −1.6935 (CLS K0) |
| Standard cold, native | +.002126938 | −.0358364390 | −1.9491 |
| Standard B, native | +1.838575e-6 | −.00152399208 | −.08288983 |
| Standard B, GL20 | +3.890556e-5 | −.00153609961 | −.08354835 |
| H2h phase I, native | +1.838575e-6 | **−.01004178649** | −.54617209 |
| H2h phase I, GL20 constraints | +3.890556e-5 | **−.00952358826** | −.51799 |
| H2h phase I, fresh GL24 | +9.474957e-7 | **−.00965347314** | −.52505 |

Native H2h roof witness: (x,z)=(.0799863864322,.8), k=.033966763923,
epsilon=.00066189887514. GL20 and GL24 outer witnesses remain near this roof
region; their full tensors and independent eigenvalue checks are saved.

| Outer scalar range relative to cold | k | epsilon |
|---|---:|---:|
| B, native | [.5,2] within rounding | [.802296,2] within rounding |
| B, GL20 | [.499640449,2.000195718] | [.803680879,2.000979443] |
| Phase I, native | [.5,2] | [.569478,2] |
| Phase I, GL20 constraints | [.5,1.999749871] | [.570556,2] |
| Phase I, fresh GL24 | **[.499548706,1.999573826]** | **[.570345737,2.000223756]** |

All fresh sampled scalars are positive and finite. GL24 nevertheless fails the
unchanged envelopes, and its inner minimum/K0=5.1534e-5 misses the positive
startup target 1e-4. GL24 was used once after terminal selection; no candidate was
subsequently changed or tuned against it. GL20 is a constraint set, not withheld
validation. None of these finite samplings proves positivity everywhere.

TKE integral is .0115867499829 versus cold .0117020062643 (−.98493%) and B
.0119310431966. These are initializer changes, not correction energy supplies.
GL20 time-scale ratios span [.25148,3.23866], length-scale ratios
[.17795,4.40248]. Full dimensional ranges and witness eddy viscosities are saved.
The unchanged velocity has strong-divergence L2=.0603080741529 against
.0603080994220; GL20 and GL24 integrals agree to about 3e-16. Native weak
continuity/constraint norms remain below 3.6e-14. Pressure bits and all sampled
velocity response caps retain their B values.

Selected alpha=0 scaled U/W/k/epsilon residual norms are
[1.2441133,2.4817564,2.2870106,.08852130]; these are not equilibrium residuals near
zero. The unchanged interface/source budget passes, with maximum relative
budget defect 3.17e-16. Full-vector traction RMS is .00296236, normal-traction RMS
7.39249e-5, maximum normal jump .000318415 and k-flux RMS 2.39992e-10.

No new counterfactual law was evaluated on standard fields. The valid H2g
six-law frozen-CLS-field table is retained; standard's Gamma=0 was not passed
silently into a regularized-CLS counterfactual. Historical CLS and standard
optimizer costs are not a matched performance benchmark. Standard neutral
isotropic normal stresses remain a model property, not an RMS-data calibration.

## Verification, cost and stopping gates

The final focused tests cover signed variables, known feasible and independently
bounded infeasible problems, convex supporting cuts, zero versus positive-margin
obstruction, nonpositive strain, numerical failure versus budget exhaustion,
and missing/failed/wrong-state/model GL24 issuance/launch receipts. The unchanged
H2g model/AWF calibration test passes. The updated probe reproduces every old
complete affine-audit and GL20 field bitwise, including state/reference/defect
hashes. Production solver code and runtime are unchanged; completed H2g native
restart, nonzero V2 import and rollback evidence is reused.

H2h issuance and its runner require the bound fresh receipt and positive margin,
in addition to native/GL20, scalar, budget, algebraic, divergence and exact-runtime
checks. A JSON attachment alone cannot enable the H2h launcher. The failed receipt
prevents issuance before any integrator history or physical step is created.

| Work | Wall time | Peak RSS |
|---|---:|---:|
| Native context export | 36.00 s native; 44.30 s driver | 1047676 KiB |
| LP/scalar separation | 9.18 s; 9.93 s whole driver | 3051892 KiB |
| Two read-only dual/localization passes | 1.54 s total | included separately |
| Complete native + GL20 + fresh GL24 | 32.73 s native driver; 34.37 s process | 1359492 KiB |
| H2h physical hill step | unexecuted | unmeasured |

The separate [I/O check](cls-cao-tamura-budget-h2h-evidence/io-cost.json) measures
one warm-cache context read/hash at 0.1777 s and a durable candidate JSON write
at 0.1003 s. This is not a cold-disk benchmark or an exclusive decomposition
of native export time; export I/O remains included in its reported wall cost.

Binary context arrays occupy 246.42 MiB. The scalar basis occupies about351 KiB.
There is no repeated velocity search or full-context rebuild per trial. Actual
optimizer work, including both diagnostic replays, is **one LP round, six complete
scans and 10.72 s** against 24/80/300. One of three permitted native candidate
checks was used. Independent final audit cost is excluded from optimizer time.
SciPy 1.8 uses serial dual simplex; OpenMP and computational-library limits are
one. Its ignored `threads` option warning remains in the original log; serial
strategy and environment controls are recorded rather than assuming that option
was recognized.

| Gate | Result |
|---|---|
| Exact B, basis, calibration/runtime identities | passed |
| Affine responses, protected rows/profile, complete membership | passed |
| Zero-margin scalar feasibility | **verified limited obstruction** |
| Positive-margin hill optimization | unnecessary after zero-margin obstruction |
| Complete native covariance | failed; scalar/budget/algebraic gates passed |
| GL20 envelopes | passed for terminal candidate |
| Fresh GL24 covariance target / scalar envelopes | failed / failed |
| V2 certificate issuance | blocked as required; none issued |
| H2h hill 8/16, 4+4 restart, curved rollback | unexecuted; zero physical attempts |
| Common flat U/k spatial target, developed flow, temporal/spatial accuracy, experimental validation | open |

## Reproduction and next package

The source/runtime manifest records the new frozen probe SHA256
`9cfdc803d847ef63185a84310010c36b622a9f96fac8bdca16229c6b27e2dbe0`.
The unchanged selected solver is H2g `runtime-v2`, SHA256
`844871bb5c70224658c88cca3f8088d7738e6f50a08ae521c74c7d02f4aec0a0`.
The [runtime binding receipt](cls-cao-tamura-budget-h2h-evidence/runtime-bindings.json)
verifies both selected closures against their frozen manifests under the actual
library path. Final handoff review fixed comparison of ASLR-dependent `ldd` text
and ambient library lookup: issuance/launch now compare stable executable/path/
library hashes. Original failed-receipt ambient metadata is preserved; this
sidecar supplies the corrected provenance. The native invocations already used
the frozen probe libraries, so their fields and audit results did not change.
Repeated runtime identity and changed-library rejection checks pass.

`LP-used-source/` preserves the exact initial optimizer source. Subsequent reader
hardening and report-label changes leave the optimization classes/functions
unchanged; the sparse dual and actual tensors were independently replayed.

```sh
R=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2h-20261113
G=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2g-20261113
H1=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111/compatibility-resolved
cd "$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python3 test_cao_tamura_h2h_scalar.py --output /tmp/h2h-check.json
# For a fresh unused result root, preparation receives the saved H2g scalar basis:
python3 run_cao_tamura_budget_h2h.py prepare --parent "$G" \
  --probe "$R/probe-runtime/RANSWallResolutionProbe" --h1 "$H1" \
  --preflight "$G/comparison.json" --output /path/to/unused-context
python3 run_cao_tamura_budget_h2h.py search --parent "$G" \
  --context "$R/scalar-context/context.json" --preflight "$G/comparison.json" \
  --output /path/to/unused-search
```

Do not repeat the withheld audit as another optimization opportunity. Its exact
command, inputs, state hash and failed receipt are indexed under `terminal-audit`.
No H2h continuation is running, and no long run was launched.

**Next recommendation:** one bounded, native-continuity-preserving outer-velocity
repair that explicitly controls both H-minus/H-plus covariance and the roof/slope
witnesses, retaining the scalar envelope and endpoint data. The fixed H-plus
witness rules out scalar-only enrichment as a complete remedy for this B state.
Keep the present scalar LP as a conditional subproblem; do not release K, widen
bounds or retune physics automatically. Affordable common flat spatial
qualification remains necessary before developed-flow or experimental claims.
