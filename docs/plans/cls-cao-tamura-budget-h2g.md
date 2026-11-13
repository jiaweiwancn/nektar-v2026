# H2g: calibrated standard k–epsilon and regularized CLS

13 November 2026, Asia/Shanghai. **The calibrated standard path is implemented
and passes its native flat/model/identity tests. Its bounded hill initialization
still fails outer covariance. No physical hill step was attempted.** The deployed
CLS base and short native trajectory reproduce their retained bits.

## Definition and tested implementation

Work started at `f13c9240a1f21d3c9550509c2ca317ef05c543be` on
`codex/cls-epawf-pimple-one-sweep`, in `.worktrees/pimple-one-sweep`. Fetch matched
that reviewed head. The main checkout, other worktrees and shared Linux source
were preserved. The latter had unrelated local documentation work.

Validated implementation: `07b680966481d331d1d1d727d118e425829b3f51`.
The [Git/source binding](cls-cao-tamura-budget-h2g-evidence/commit-binding.json)
links the commit to the unchanged native source/runtime manifest.

The [plan](cls-cao-tamura-budget-h2g-plan.md),
[predeclared policy](cls-cao-tamura-budget-h2g-evidence/policy.json),
[comparison manifest](cls-cao-tamura-budget-h2g-evidence/comparison-manifest.json),
[native findings](cls-cao-tamura-budget-h2g-evidence/native.json),
[frozen comparison](cls-cao-tamura-budget-h2g-evidence/frozen.json),
[flat results](cls-cao-tamura-budget-h2g-evidence/flat.json),
[verification](cls-cao-tamura-budget-h2g-evidence/verification.json) and
[source/runtime manifest](cls-cao-tamura-budget-h2g-evidence/source-revisions.json)
are the compact record. [The artifact index](cls-cao-tamura-budget-h2g-evidence/artifact-index.json)
locates native fields, contexts, packets, matrices and rejected states under:

```
/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2g-20261113
```

UUID `5AB8A669B8A642FD` was verified. Linux Release builds use four workers,
`-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS`, one native rank and one
computational-library thread. Other project jobs remained under their existing
16-worker guard; the aggregate remained below 72. USB SSH timed out and the
established Tailscale fallback was used. No historical runtime was overwritten.

The final production/test closure is `runtime-v2` / `probe-runtime-v2`:

- Solver SHA256 `844871bb5c70224658c88cca3f8088d7738e6f50a08ae521c74c7d02f4aec0a0`.
- Probe SHA256 `240e15cf219fd64085dcac21bd5bb471e747728e70bed0fead3727ccc51df7ec`.

The first frozen closure is retained. It produced the unchanged CLS readbacks
and frozen-field data, then exposed a CLS-only forensic counterfactual invoked
on standard. That diagnostic dispatch defect was fixed before standard
certification or physical stepping. Its rejected setup receipt remains archived.
The V2 closure was frozen before certificate issuance and stayed immutable.

## Two complete calibrations on one mesh

Common inputs are H_A=z_K=z_m=.0065 m, u_tau=.07426793792419734 m/s,
nu=3.8095238095238094e-6 m²/s, kappa=.4405539051866587 and
z0_eff=.00012222264395414748 m. The frozen seven-point mean fit, weights and all
Cao–Tamura data hashes are unchanged. No new fit or hill observation was used.
Both outer inlets are the same unshifted physical-height logarithm and 1/z epsilon
through the roof. The upper-profile and roughness-sublayer qualifications remain
open; H_A is only 1.3 cube heights. The existing open hill boundary roles and
pressure convention are unchanged. The flat verification fixture uses its
separate compatible prescribed traces.

| Quantity | Deployed CLS | Calibrated standard |
|---|---:|---:|
| Native model | CLSHighRe | StandardKEpsilon / StandardNeutralCalibratedV1 |
| Cmu law | CLSCmu(max(qS,qOmega)) | constant .09 |
| Cepsilon1 / Cepsilon2 / sigmaK | 1.44 / 1.92 / 1 | same |
| sigmaEpsilon | 1.3 | 1.3478315512167742 |
| Nonlinear stress / Yap | quadratic + beta=.4 cubic / existing Yap | absent / zero |
| Physical Gamma / gamma0 | 1 / 0 | 0 / 0 |
| q | 3.1906750816242875 | 3.3333333333333335 |
| K0, m²/s² | .017598891430879275 | .01838575534504143 |
| Wall cMu | .09822789532821126 | .09 |
| Wall cLength | 2.5108624752197968 | 2.6811256851969465 |
| Deployed effective h_s, m | .0022472751006336437 | .0022035111852050995 |
| Viscosity intercept, m | about −8.133457321054e-6 | same |
| Independently solved flat cutoff, m | .0008312858196541334 | .0008312858196541342 |
| J1/J0, m | .00045714516865748154 | .0004685594287524002 |
| Actual native slots | 2804 | 2804 |

The native standard sigma-epsilon follows
kappa²/[(Cepsilon2−Cepsilon1)sqrt(Cmu)]. The native roughness root uses a bracket
strictly within (0,H_A), molecular viscosity and the actual AWF correlation.
Its residual is zero at binary64 reporting precision. Re-solving the CLS root
gives .0022472751006336117 m, within 3.2e-17 m of its deployed value; **the CLS
case keeps the deployed value**. This is numerical root arithmetic, not a refit.

U_N(H_A)=.6698827910508315 m/s. Independent molecular-retaining budget evaluation
P=u_tau^4(J0−nu J2) gives .002842113026557174 m³/s³. A_epsilon is
.0009298331853917931 m³/s³ and E_C(H_A)=.1430512592910451 m²/s³.
Both constructors solve their own cutoff; it is not assigned from the other
model. Curved checks at signed curvatures ±19.7392088/m and K/K0=.8,1,1.2 retain
positive reference k and explicit local budget roots.

On the flat neutral reference the maximum paired differences in U, nu_A,
epsilon, inner k gradient and additive k offset are respectively about
1.1e-16, 2.7e-20, 1.1e-15, 2.4e-14 and 3.2e-17 in their physical units.
Delta K0=.000786863914162155 m²/s². The inner profile remains nonconstant.
Reference production/dissipation agree to roundoff. The J1/J0 difference is
about 2.50%; calibrated neutral equality does not make the off-reference
forcing/pressure-response operator identical. Off-reference and curved budgets
continue to use each model's current transported endpoint K.

The common mesh is the exact retained six-streamwise by thirteen-radial-panel
H2f allocation, with cubic velocity/scalars, regional linear pressure and 16×16
production quadrature. Standard roughness/cutoff roots do not move its knots.
Each profile still inserts its own required integration breakpoints and builds
its own inventory. The 3000-slot production guard remains unchanged.

## Native dispatch, identity and flat qualification

`RANSWallModelSpec.h` is the shared evolution/probe loader. It passes effective
coefficients into the actual native standard constructor, validates the
calibration-manifest content, and checks expression-level Qstar, Kappa and K0.
Both solver and probe emit native coefficients and active-term readback.
`StandardCalibratedBudgetAWF` has a distinct identity; the historical CLS
formulation and signatures retain their meaning. Standard uses the `Original`
variant selector to mean its unmodified standard stress; this is not an
Original-CLS evolution campaign.

The analytical inlet is shared independently of physical CLS regularization.
Standard selects no quadratic/cubic stress, Yap or regularizing stress. The
positive implicit viscosity, complete stress-work production, scalar sources,
normal pressure reference and audits use its selected model. The standard-safe
forensic diagnostic reports its CLS-only Regularized counterfactual as
inapplicable. The legacy R3 action remains explicitly excluded for prescribed
inflow; it is not silently evaluated with a logarithmic surrogate.

The standard context is linked to its exact native export receipt, array hashes,
model and immutable base. The independent evaluator dispatches by that native
specification. Standard fixed-scalar Rayleigh derivatives are analytic; sampled
checks differ from finite differences by at most 2.6e-13. CLS evaluation still
calls the retained formula. Reference, model, AWF, mesh, inventory, cache and
checkpoint/certificate identities remain distinct.

Focused tests reject wrong model manifests, wrong effective coefficients,
Qstar mismatch, solver/probe model mismatch, a relabelled CLS context/evaluator,
cross-model checkpoints/certificates and ordinary changed-dt restart. The touched
CLS h0 and two-step native histories are bitwise identical to retained evidence;
no H2f optimizer campaign was repeated. Existing consistent-mass QP checks pass.

The standard nonzero V2 fixture changes only permitted outer-k/epsilon rows.
Complete native and GL20 gates pass before issuance. Exact import preserves the
candidate/reference distinction, defect, endpoint K, prescribed data and
protected epsilon. It performs zero Lift/projection. Two direct steps agree
bitwise with 1+1 split stepping; post-assignment rejection restores the coupled
state, native histories, time/alpha, samples and parent checkpoint bytes.
Every accepted step has one U/p/outer-epsilon/k solve and one logical bound call;
solved epsilon remains authoritative after k. No correction activates or supplies
energy. Refreshed matching is reported separately, not used to overwrite epsilon.
In the first unperturbed standard step, the cutoff range moves from its
reference value to [.00083128581965,.00083128657596] m as endpoint K changes.
The inlet/reference K0 is not an evolution clamp.

An additional unperturbed standard control and the short retained-CLS regression
each accept two steps at dt=2e-6 s through T=4e-6 s. These are compatibility/startup
results, not developed flow or steady convergence.

| Flat metric | CLS | Standard |
|---|---:|---:|
| Initial sampled traction, m²/s² | .00552261884447 | .00552261884447 |
| Initial velocity-gradient L2 error | .259916409246 | .259916409244 |
| Free scaled U spatial residual | .1342977314 | .1341183090 |
| Free scaled W spatial residual | .00472120461 | 5.03e-16 |
| Free scaled k spatial residual | .04761493032 | .04557832032 |
| Free scaled epsilon spatial residual | .00178977690 | .00289100605 |
| T=4e-6 normalized U drift | 5.0396e-6 | 5.0299e-6 |
| Normalized k drift | 1.5476e-6 | 1.4815e-6 |
| Normalized full-outer-epsilon drift | 9.6346e-8 | 1.4082e-7 |
| Normalized outer-pressure drift | .00827250 | .00593566 |

Dimensional norms, domains, native field scales and both regional pressures are
in `flat.json`. One common pressure gauge preserves their relative offset.
Pressure has an immediate algebraic response; the maximum drift must not be
reported as only the small velocity drift. Exact analytical weak loads under
the homogeneous tests are near roundoff; the final constrained-row entries in
the raw norm arrays are reactions, not unforced evolution equations.
Both represented U/k residuals still miss the existing 1e-3 spatial target.
Removing the nonlinear flat normal-stress error does not remove the common
radial representation error or qualify the standard model's physical accuracy.

## Frozen-field constitutive experiment

These rows use the **same native CLS base/A/B tuples**, not model-consistent
standard initializations. All gradients, local k, actual/dependent epsilon,
coordinates and native sides are retained. Six native laws and independent
algebra agree on every minimum and negative-point membership count. Native
versus independent tensor/production error is about 1.3e-15; combined velocity/
scalar directional errors are below 4.2e-11, including formal strain/rotation
branch limits.

| Frozen law: minimum covariance, m²/s² | CLS base | H2f A | H2f B |
|---|---:|---:|---:|
| Deployed regularized CLS | −.05289055 | −.03426784 | −.02980405 |
| beta=.4 CLS, no regularization | −.00201757 | −.00076126 | −.01303724 |
| Original CLS, beta=1, no regularization | −.00690906 | −.00433812 | −.01303724 |
| Local-CLS-Cmu linear | +.00160236 | +1.75989e-6 | +1.75989e-6 |
| Local-CLS-Cmu linear plus regularization | −.05047703 | −.03765474 | −.03314650 |
| Calibrated constant-Cmu standard | −.03233253 | −.03206052 | −.07550303 |

GL20 confirms the same signs. Complete unique/membership counts, every law's
own witness/eigenvector, common-witness tensor contributions, production and
invariants are retained. The physical GL20 domain volume is .636 m² per unit
span. Negative-volume measures are quadrature estimates: for the frozen base,
regularized CLS and standard give .0331434/.0264682 m²; on H2f B they give
.188419/.173774 m². An improved worst eigenvalue does not imply a smaller
negative-volume measure.

Covariance and tangent diffusion are distinct. On the exact CLS neutral tuple,
the 64-angle frozen-scalar symmetric incompressible-symbol scan gives Original
CLS about −.003156 m²/s while beta=.4/no-regularization is positive in the tested
branches. The standard scan is positive even at its negative-covariance
witnesses. The variable-Cmu linear law's positive covariance on these fields
is not a qualification for evolution. These finite directional scans do not
prove global ellipticity or nonlinear stability. Physical incoming height and
nearest-floor distance remain separate in the source/tangent diagnostic.

## Fresh standard hill initialization and bounded searches

The standard base is regenerated from its own reference, inner completion,
pressure balance and defect. Its hash is
`038fc076cebbea41bcb7cc349147eacc2b938a9051827bbb386baf7c3cb3716c`;
the CLS base remains `83b1a4fd…045e79b5`. Standard cold velocity coefficients
differ from CLS by only 3.48e-13 m/s. The velocity direction pool and physical
mass Gram are bitwise unchanged; divergence/strain forms differ by at most
2.8e-13. Both searches use the **same numerical divergence limit**
.06030809942199109 and common initializer strain weight kappa*u_tau*z, despite
standard's zero physical Gamma.

| Model-consistent hill state: outer minimum | Native | GL20 |
|---|---:|---:|
| CLS cold | −.0528905535 | −.0526123788 |
| Standard cold | −.0358364390 | −.0357093840 |
| Retained H2f CLS B | −.0298040470 | −.0300134540 |
| Standard B | −.0015239921 | −.0015360996 |

These rows compare feasibility, not matched optimizer work or evolved flow.

The 41-direction pool contains the 32 original directions and nine local
H2f directions selected using CLS witnesses. No enrichment was rerun or added.
The fresh standard context has 187395 numerical observation records, with all
memberships preserved. It is not interchangeable with the CLS inventory.
The scalar basis remains eight outer-k and eight outer-epsilon directions;
inner k, inlet/interface values, endpoint K, dependent epsilon and pressure
are frozen bitwise only during initialization.

| Standard state | Native inner / outer | GL20 inner / outer |
|---|---:|---:|
| Cold base | +.002126938 / −.035836439 | +.002144972 / −.035709384 |
| A: velocity only | +1.838575e-6 / −.014825005 | +4.859429e-5 / −.014922766 |
| B: outer scalars free | +1.838575e-6 / −.001523992 | +3.890556e-5 / −.001536100 |

A stops after 19 majors at its predeclared relative-progress threshold; B uses
all 40 majors and is still improving at the end. They consume 140/267 complete
scan work units and 76/320 active derivative batches, within the respective
180/320 and 3000/6000 limits. No search was extended. Inner cuts have no outer
slack; every saved/terminal state passes the complete inner scan. At fixed
standard scalars the tensors are affine and the declared finite-inventory
constraints are convex. No global feasibility/infeasibility certificate was
computed; B remains nonlinear. The inherited raw search label is retained with
this explicit mathematical interpretation. The positive
1e-4 K0 search margin and unchanged hard audit remain separate.

B's velocity L2 departure is .05060 (limit .12), maximum component increment
.46729 m/s (limit .5), and divergence L2 .0603080742. The inner margin and
strong-divergence restriction are close; scalar ratios reach k=[.5,2] and
epsilon=[.8023,2]. Response/scalar rows enter the local working sets, but the
accepted velocity-distance and component caps are not at their limits.
Initial integrated TKE rises by about 1.96%; this is seed construction, not
positivity-correction supply. Time-/length-scale ratios and exact native
frozen-row checks are in the compact evidence.

The final native outer witness is at physical (x,z)=(.0799864,.8), on the roof,
with k=.0366434, epsilon=.000970184 and nu_t=.124560 m²/s. Its critical isotropic
and linear contributions are +.02442894 and −.02595293; quadratic/cubic and
physical regularization are absent. A forensic reconstruction leaves only an
ulp-scale remainder (~2e-18), not an active regularization term. Independent
tensor/eigenvalue replay agrees to roundoff. At this frozen gradient, modest
local k reduction or epsilon increase would help that tensor; global basis,
other witnesses and active bounds still restrict any such field change.

All complete scalar/budget, frozen-data, boundary, continuity and divergence
gates pass. Outer covariance fails. **No hill certificate was issued; the
8/16-step pair, curved 4+4 restart and hill rollback were not executed.**
This is a bounded local search failure within this CLS-selected span and
specified bounds, not a global physical incompatibility certificate. The
retained H2f CLS searches have different starting histories and do not provide
a matched optimizer-performance comparison.

## Costs, gates and handoff

Measured cold native readbacks cost about 21–28 s each. Standard context
construction costs 142.54 s; A/B searches cost 33.39/123.89 s with peaks about
1.20/1.19 GiB. The standard nonzero V2 two-step driver costs 42.22 s including
setup/audits/output, versus roughly 4 s per continued native physical step.
The unperturbed two-step block costs 29.25 s and peaks at 1518644 KiB. These
are measured job costs, not an alternating performance benchmark or a claimed
model speedup. The largest changed-candidate fixture peak is about 1.50 GiB.
There is no measured physical hill-step cost because no hill step qualified.

| Gate | Result |
|---|---|
| Native model coefficients, analytical neutral balance and AWF roots | passed |
| Shared mesh and common initializer metrics | passed |
| Complete frozen-field and branch-aware independent algebra | passed diagnostic |
| CLS h0/two-step state/history regression | bitwise passed |
| Standard flat native stepping, nonzero V2, split and rollback | passed; zero correction |
| Model/manifest/context/expression/restart incompatibility tests | passed |
| Represented flat U/k spatial-equilibrium target | remains failed |
| Cold/enriched standard hill covariance | failed after bounded search |
| Conditional hill pair and quantitative experimental validation | unexecuted / unqualified |

The [final observed status](cls-cao-tamura-budget-h2g-evidence/final-status.json)
at 19:36 Asia/Shanghai confirms all H2g jobs have finished; no long continuation
was launched.
Standard remains a numerical and mean-flow control: its parallel neutral
normal stresses are isotropic and are not a fit to the measured unequal RMS
components. All held-out hill observations remain held out.

**Next recommendation:** one bounded native active-constraint/accessibility
analysis at the standard B roof witness, including divergence and scalar-basis
coupling, before choosing additional optimization work or one local response.
The last accepted search step still improved the objective; this evidence does
not justify widening the envelope, releasing endpoint K, changing the closure
or declaring that another direction is necessarily required.

## Reproduction

Use a fresh root and the indexed source overlays; never repurpose a historical
checkpoint. Configure the external RANS test build against the established
Linux installation with `RANS_BUILD_SOLVER=ON`, the recorded flags and immutable
prebuilt SpatialDomains library, then build with `--parallel 4` these targets:

```
RANSCaoTamuraBudgetH2gChecks RANSH2gConstitutiveProbe
RANSWallResolutionProbe IncNavierStokesSolver RANSKBoundQPChecks
```

`RANSCaoTamuraBudgetH2gChecks retained-H2c-preflight.json calibration.json`
performs native calibration. `run_cao_tamura_budget_h2g.py prepare` combines it
with the retained CLS preflight and exact H2f mesh definition. Its `cold`,
`frozen`, `context`, `flat`, `unperturbed`, `certify` and conditional `pilot`
operations are explicit. Freeze solver/probe closures before certification.

The bounded searches use `search_cao_tamura_budget_h2f.py --context
standard-context/context.json --preflight comparison.json --start ... --stage
A|B --divergence-limit .06030809942199109 --output ...`; the native model
specification selects standard evaluation. A starts from the new cold state,
B from its completely audited A candidate. No modified CLS checkpoint is used.

`check_cao_tamura_budget_h2g_context.py`,
`verify_cao_tamura_budget_h2g.py`,
`diagnose_cao_tamura_budget_h2g_flat.py` and
`analyze_cao_tamura_budget_h2g_frozen.py` execute the focused checks/readbacks.
`collect_cao_tamura_budget_h2g.py` verifies Mac/Linux source hashes and builds
the compact evidence/index. Exact input packets, commands, native histories,
usage logs, runtime/library hashes and terminal receipts are indexed on Linux.

## H2h follow-up, 13 November 2026

The [executed H2h package](cls-cao-tamura-budget-h2h.md) freezes the exact standard-B
velocity and pressure and uses the existing sixteen scalar directions. It finds
a checked sampled zero-margin obstruction and an independently negative H-plus
tensor whose scalar responses are identically zero. The terminal native/GL24
audits fail, so no hill certificate or physical step follows. This advances the
continuing comparison without changing either calibration or repeating H2g's
searches. Common flat spatial qualification and later developed-flow, accuracy
and held-out experimental gates remain open.
