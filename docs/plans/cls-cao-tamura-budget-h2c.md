# Cao–Tamura budget AWF H2c

2026-11-13, Asia/Shanghai (work executed November 12–13). The analytical neutral provider, frozen-data height
inventory, per-height roughness calibration, native flat controls and moving-foot
coverage repair are implemented and tested on Linux. **No hill physical step
was accepted.** All three tangential-following seeds fail outer covariance.
One explicit free-interface initializer improves the primary candidate but also
fails before physical advancement. Its analytical, pre-projection counterpart
is already inadmissible. The conditional 8/16-step pair was not run.

## Identity and retained formulation

Work stayed in the existing `codex/cls-epawf-pimple-one-sweep` Mac worktree,
starting from local/fetched `85ce59e95a2e44a5217216ae1842e248d8a24ec9`.
The implementation commit is `f4a886071485dedb48168dcfe6ad2baf17dfcaa0`.
Read-only witness side labels and evidence reduction were finalized in
`c44c4bdda418abb16567c4f4eba26022cedf90ff`; no native arithmetic changed.
All 20 changed implementation/test files match the Linux export. Final solver
SHA256 is `3b1ba4f2b2b1b48fbe6e980d82482f382f68917d541e4e68918310f2dc9ab2d8`;
final probe SHA256 is `59035e0ea945732c344c3819e62a7142c96b74cf3960f44527040589005c8889`.
[Source/runtime identities](cls-cao-tamura-budget-h2c-evidence/source-revisions.json)
and the [artifact index](cls-cao-tamura-budget-h2c-evidence/artifact-index.json)
identify every executed variant. Source hashes were compared with the isolated
Linux export; builds used four workers, native pilots one rank and one
computational-library thread. Other project jobs and the shared Linux checkout
were preserved. The verified data-volume UUID is `5AB8A669B8A642FD`; its actual
mount remains `/media/jay/新加卷1`.

Heavy artifacts are under:

```text
/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2c-20261112
```

The equations remain Regularized CLS/Linear inner stress, beta=.4, Gamma=1,
gamma0=0, calibrated budget-cutoff AWF, physical inner molecular diffusion,
and high-Re outer equations. The sink is E_C; the actual source is P_actual-E_C.
The counterfactual flat parallel p_A0(K) determines the metric-weighted cutoff.
There is no EP blending, E_C/k reaction, stress clipping or increased floor.
Evolution uses native IMEX1, one U/p/epsilon/k solve and one logical bound call.
Solved epsilon remains authoritative after k.

H1/H2a/H2b finite-profile cases remain separate, unchanged histories. The old
H2b seed reproduces its exact coupled hash
`7adda6c42bd9d6e8fae1c5259055dbd60290888593c0fbc626d027b39b2203e4`
and its negative covariance. The new foot coverage does not explain it away.

## Analytical reference and data contract

`RANSWallRegularizedReference=AnalyticalNeutralInflowV1`, paired with the same
named native initializer, selects an immutable analytical provider. It calls
`RANSNeutralReference::AtHeight()` with the selected CLS parameter object:

- q_star=3.1906750816242875, kappa_CLS=.4405539051866587;
- K0=q_star*u_tau², epsilon_N=u_tau³/(kappa_CLS*z);
- U_N=(u_tau/kappa_CLS)*log(z/z0_eff), nu_t,N=kappa_CLS*u_tau*z;
- S_ref=U_N'(z)(e_x tensor e_z+e_z tensor e_x), gamma=Gamma*nu_t,N.

The coordinate is physical incoming-frame height, with a globally fixed
horizontal frame and zero normal origin. It is distinct from inner r, outer
zHat/eta and nearest-floor distance. No table interpolation or constant tail is
used. The exact log/1-z outer law continues to the roof. The inner velocity uses
the native AWF, inner k uses its **nonconstant ReferenceK**, and the interface
uses actual E_C(H_A). During evolution K remains the transported endpoint trace.

Closure, implicit viscosity, cancellation-stable residual components and
constitutive diagnostics share the selected provider. Legacy R3 Upstream actions
reject it. Legacy analytical weak-load substitution is explicitly affine-only.
The identity includes mode/version, all selected model/reference parameters,
origin/frame, fit-content hash, H_A, h_s, mesh and boundary interpretation.
Ordinary incompatible fit/inventory restarts reject; unsupported fresh normal
origin shifts reject. Default NeutralLogV1 and FiniteInflowV1 paths retain
bitwise governing vectors, states and histories against the frozen H2b runtime;
the neutral split restart also agrees bitwise.

The frozen `cao-tamura-2006-v1.1` snapshot passes its auditor: 98 series, 1,151
points and unchanged definitions. Only rough-flat Fig.3a enters fitting. Fig.3d
and all hill observations remain held out. The physical convention retains
u_tau=.0742679379242 m/s, nu=3.80952380952e-6 m²/s and nominal U(Hhill)=1 m/s.
These are the existing explicit dimensional equivalent, not independently
measured absolute air properties. Neither u_tau nor nu was fitted.

## Common fit, height and roughness qualification

With frozen weights, fit only the intercept of U=(u_tau/kappa_CLS)log(z/z0_eff).
The predeclared 10–51 mm window contains seven points, including 50.27 mm.
The result is **z0_eff=.122222644 mm**, weighted RMSE **.025492133 Uinf**.
The six-point 15–51 mm sensitivity is retained with its own residuals and h_s
roots. No window was selected using hill or stress outcomes.

[Calibration evidence](cls-cao-tamura-budget-h2c-evidence/calibration.json) records
all points, weights, residuals, root brackets and residuals, calibrated wall
coefficients, signed-curvature/K screens and GL16/32/64 checks. Each h_s solves
u_tau² J0(H_A;K0,h_s,nu)=U_fit(H_A), using the actual roughness correlation and
cMu_wall=q_star^-2, cLength_wall=kappa_CLS*q_star^(3/2).

| H_A (mm) | h_s (mm) | cutoff (mm) | H_A+ | nu/D_A(H_A) | clearance below 7 mm |
|---:|---:|---:|---:|---:|---:|
| 5.5 | 2.259266 | .825999 | 107.224 | 2.070% | 1.5 mm |
| 6.0 | 2.252768 | .828856 | 116.972 | 1.901% | 1.0 mm |
| 6.5 | 2.247275 | .831286 | 126.720 | 1.758% | .5 mm |

The sensitivity-window h_s values are 2.151537, 2.145138 and 2.139727 mm.
U_fit(H_A) is an extrapolation below measured means. All three outer analytical
laws agree on their common domain, so this common mean fit cannot rank heights.
The fitted law implies U(Hhill)=.97620304 m/s and Re_H=10250.13, 2.38% below
the nominal 10500. The existing rounded experimental Reynolds-number
inconsistency remains; no simultaneous calibration was imposed.

The [measurement inventory](cls-cao-tamura-budget-h2c-evidence/measurements.json)
retains 283 rough/RR stress/RMS points, their original conventions, visibility,
weights, extraction bounds and physical classifications. 271 have recoverable
physical heights. The flat stress/RMS minimum is 10.089 mm; Fig.12's RR wake
traverse at 7 mm is the known minimum. Fig.11's 12 points remain unclassified in
physical height because local delta_7.5 is unavailable. Incoming delta=.25 m
was not substituted. Frozen displacement conversion was not applied twice.
Flat-reference points are measured from the flat floor, not from a fictitious
hill at their x=0 convention. Hill distances use independent nearest-floor
minimization; vertical experimental heights are not treated as normal r.

H_A+>=100 and molecular fraction<=5% are **provisional project screens**.
They all pass. The planes are only 1.1–1.3 physical cube heights; this does not
establish roughness-sublayer clearance. h_s, cube height b=5 mm, nominal
z0=.2 mm, fitted z0_eff and d=3.5 mm are distinct. Displaced wall units are
reported only as diagnostics; no z-d law was inserted into epsilon or Yap.
A physically validated overlap between measurement coverage and wall matching
is therefore unresolved. These are exploratory cases, not validated interfaces.

## Native representation and flat verification

The final six-streamwise-panel mesh has **2804 coupled slots**, below the
unchanged 3000 limit, cubic velocity/scalars and linear regional pressures.
Feet and crest remain explicit. Radial panels resolve the positive rapid
viscosity scale, cutoff, h_s and H_A; no negative intercept is made into a knot.
The first knot is one eighth of nu/slope-z_v. The fitted z_v values are
-8.57 to -8.13 micrometres and are retained without clipping. Positive actual
viscosity, ReferenceK and cutoff budgets pass. Curvature ratios are .1086–.1283;
positive normal reach, monotone x_H and positive outer Jacobians certify the
same composite domain construction. Cutoffs are recomputed for each curved
profile/current K; they are not frozen to the flat values.
The molecular fraction accounts for the small inner/outer **turbulent** shear
difference. Analytically the inner molecular-plus-turbulent traction equals
the outer high-Re turbulent traction, u_tau². Reference k fluxes at ground and
matching are zero to budget/root accuracy; compatible regional modified-pressure
offsets are retained.

The affine control uses reference-compatible prescribed streamwise/top traces,
including the nonzero neutral top shear/epsilon-flux context. It is distinct
from the hill's impermeable free-slip top and absolute total-traction outlet.
No artificial volume cancellation source is supplied. Projected boundary traces
are part of representation error. Direct analytical values **and** gradients
under native tests are reported separately from coefficient-derived residuals.

At H_A=6.5 mm, the native h0 errors are U L2=6.518e-5, grad-U L2=.259916,
k L2=1.005e-6 and complete outer-epsilon L2=1.118e-5 in their dimensional norms.
Initial sampled traction is .005522619 versus u_tau²=.005515727 (+.125%).
Analytical inner P/E is one; native integrated P/E is 1.000633. Full weak source
modes, diffusion, transport, storage and constrained-row reactions remain
separate in [flat evidence](cls-cao-tamura-budget-h2c-evidence/flat.json).
Near-unity integrated production/destruction does not establish balance of
higher modes.

At h0, direct analytical values and analytical gradients under the same native
homogeneous tests give scaled free-load norms
[1.83e-13,5.75e-15,9.90e-15,6.73e-16] for U/W/k/epsilon. Essential rows are
zeroed for this diagnostic; it is not the residual of the represented field.
This separates analytical reference balance from its nonzero native projection
and differentiation error. The prescribed boundary projection is included in
the representation qualification, not hidden by a volume cancellation term.

Eight steps at dt=2e-6 s reach **T=1.6e-5 s**, with zero bound activation/supply.
The 4+4 restart agrees bitwise; deliberate post-assignment rejection restores
coupled state and all histories and leaves the parent intact. Each accepted
step has exactly one U/p/epsilon/k solve and one logical bound call. Epsilon
retention is bitwise confirmed. This short exact-reference fixture is not a
transient temporal-order test or an open experimental precursor.

Normalized h0-to-step8 changes are U/W combined 1.29e-6, k 4.84e-6,
epsilon 4.81e-7 and common-gauge pressure 4.54e-5. They pass the predeclared
1e-4 absolute field screen; no response-relative drift screen was used.
Mean traction changes by .1282% of u_tau². Native selected alpha=0 block norms
at h0 are [.13430,.004721,.047615,.001790]; after eight steps they are
[.083337,.004058,.021680,.001788]. **The 1e-3 spatial equilibrium target fails.**
Small short-time drift is not proof of an exact discrete neutral equilibrium.

## Moving-foot coverage and projection

`AcceptedBudgetV2` / `NativeTangentialSideV3` are explicit new identities.
Moving foot-minus and foot-plus rays select their intended native element
independently of radial interface side. Observations retain tangential side,
geometry, exact point/mean rows and actual/dependent epsilon context. Covariance
uses the same element's physical derivative. Exact audit reuse keys already
include this provenance; no coordinate-only deduplication was introduced.

The synthetic moving point is outside the fixed reconstruction grid. Its two
sides have different gradients, with negative covariance only on the selected
opposite side at that point. Both native inventories retain distinct geometries.
Cached/uncached primary states, scalar inventories, complete audits, witness
results and this fixture are bitwise identical.

The new inventory's active projection fixture passes continuous-ground and full
inventory feasibility, exact non-k/essential preservation and independent GL32
physical-mass objective verification. It uses one KKT solve: objective
1.5996180e-16, stationarity mass-dual 1.22e-24, primal/dual defects zero and
complementarity 5.92e-32. This is a projection-only test, not a physical step
or acceptance of a covariance-inadmissible seed. The guarantee remains the
declared finite inventory plus continuous ground, not global volume positivity.

## Hill initialization: completed gates and blocking witnesses

Each height was screened at zero native steps. With the tangential-following
streamfunction initializer, inner covariance stays positive, but outer minima
are -.0766133, -.0735663 and -.0711611 for 5.5/6.0/6.5 mm. At the primary
outer witness (x,z)=(.1,.1) m, the independent **analytical radial** construction
has minimum -.0901318. Native radial/projection error does not create this
failure; the 2D projected minimum is less negative in this case.

One bounded alternative, `PhysicalHeightInterfaceV1`, lets interface speed
follow the incoming profile at its physical height. A complete streamfunction
supplies its normal velocity; physical-mass projection enforces inlet/ground/top
and native weak continuity. Scalar/interface data and the immutable reference
are retained. The upper target remains the neutral law above the declared
z=.2 blend height, without a constant tail. Both initializer versions have
explicit identities. No stress/floor/source parameter was adjusted.

For this alternative at 6.5 mm:

- Inner projected covariance minimum is +.000335983.
- Outer projected minimum is **-.053444902** at approximately
  (x,z)=(-.0249201,.1306473) m.
- The independent analytical radial target at that witness is already
  **-.041237863**; native representation/projection worsens it.
- Independent tensor/eigenvalue replay agrees with native tensors to about
  1e-17; gradient step-halving differences are about 5e-8/s.

[Witness evidence](cls-cao-tamura-budget-h2c-evidence/witnesses.json) separates
analytical radial fields, native 1D targets at retained scalar values, and the
actual projected 2D fields, with native sides, gradients, tensors, eigenvectors,
term contractions, epsilon contexts, divergence and scalar/budget margins.
Regularization is evaluated with physical incoming height throughout.
At the adjusted projected outer witness, critical-direction contractions are
isotropic +.0117326, linear -.00356713, quadratic -.00308783, cubic -.000361314
and regularization -.0581612, summing to the reported negative eigenvalue.
Native weak continuity is 2.60e-16 after setup, independently of its nonzero
strong divergence and local constitutive strain.
These findings concern the chosen seeds. They do not prove nonexistence of an
admissible hill solution or loss of ellipticity.

The actual native solver reproduces the alternative's rejected coupled hash at
T=0, with zero implicit callbacks, zero governing solves and zero accepted
steps. Its initialized native histories/reference/defect and full failure are
retained. **The 8/16-step hill pair, curved physical split restart and hill
finite-step sensitivity remain unexecuted because the initial gate fails.**
No unchanged numerical rejection was retried as a physical campaign.

## Cost, reproduction and next decision

The flat h0 takes 14.50 s (peak 966 MiB). The eight-step block takes 70.41 s
wall/67.11 s user CPU, including restart, reporting and checkpoints, peak
1598 MiB. Its inclusive cost is 8.80 s per accepted step, not an isolated
factorization timing. Native mapped initialization/rejection takes 14.52 s,
peak 1001 MiB. The mapped projection/operator diagnostic uses about 27.3 CPU s
and 1566 MiB. **No physical hill-step cost was measured.**

Initial software-only stops are archived separately: the inherited constant-tail
freeze was incompatible with native continuity for a log profile; an isolated
synthetic fixture needed a clean velocity background/moving point; the old
reaction lifting assumed a .25 mm interface absent from this mesh. These were
fixed. Completed accepted trajectories were reused through strict readbacks,
not rerun or relabeled as numerical failures.

Use unused roots, the recorded frozen runtimes, rank one and all library
threads one. Representative commands, from the run root:

```bash
build/RANSCaoTamuraBudgetH2cChecks \
  source/docs/plans/cao-tamura-validation-data/profiles.csv NEW_PREFLIGHT.json
python3 source/solvers/IncNavierStokesSolver/Tests/RANS/run_cao_tamura_budget_h2c.py \
  --solver build/IncNavierStokesSolver --probe build/RANSWallResolutionProbe \
  --preflight preflight-final.json \
  --h1 ../cls-cao-tamura-budget-h1-20261111/compatibility-resolved \
  --phase screen --output NEW_SCREEN
python3 source/solvers/IncNavierStokesSolver/Tests/RANS/run_cao_tamura_budget_h2c.py \
  --solver build/IncNavierStokesSolver --probe build/RANSWallResolutionProbe \
  --preflight preflight-final.json \
  --h1 ../cls-cao-tamura-budget-h1-20261111/compatibility-resolved \
  --phase primary --seed PhysicalHeightInterfaceV1 --output NEW_PRIMARY
```

All H2c experiments are bounded and finished; no H2c long campaign was launched.
The three focused kernel/physical-mass test targets pass. The final executable
also reproduces the retained four-step neutral state, governing stages and
native histories bitwise. The artifact index covers 1,911 files, approximately
4.05 GB; large fields, histories, matrices and frozen libraries remain on Linux.
Unrelated live experiments and historical neutral/finite campaigns were left
untouched. Quantitative experimental validity, roughness-sublayer matching,
spatial accuracy, developed flow and hill statistics remain open.

The next package should implement **a bounded covariance-feasible initializer
with free interface and outer streamfunction degrees of freedom**, preserving
inlet/ground/top constraints, native continuity, scalar completion and the
selected equations. Check its analytical target before native projection.
Only an admissible seed should unlock the deferred 8/16-step pilot. Another
height or model-parameter sweep is not justified by this initialization failure.
