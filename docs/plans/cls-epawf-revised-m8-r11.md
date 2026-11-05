# M8-R11 objectives at the native trial scale

R11 completes the bounded four-policy comparison at the unchanged R9
initial state, the physical R8 iterate-10 endpoint. The retrospective
check confirms that R10's anchor ranking reverses at f=.125.
N1 improves the full-scale model objective but fails practical promotion.
Lf and the first Nf tangent exhaust the retained cut-exchange limit.
The usefulness of the trial-scale objectives remains unresolved at
direction construction. L1 is retained and no continuation is run.
M8 remains open.

The reviewed head is a44639b99f9edd7f0bb4f8f2134a8065ee03b132.
There were no newer remote commits before implementation. The producing
Linux revision is 136b09ef38b26fd0817b9dc427fbf03d0482fcec.
The producing checkout and build remain preserved in
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r11/epoch-001.
The separate segregated-Picard worktree, its untracked work, and the
shared Linux source checkout were preserved.

The [decision ledger](cls-epawf-revised-evidence/m8-r11/decision-ledger.json)
contains exact executable, runtime, XML, source compatibility, history,
state and result bindings. Full arrays, cut-exchange records, matrices,
field vectors and detailed logs remain on Linux.

## Fixed base and numerical contract

| Bound quantity | SHA256 |
|---|---|
| Physical R8 iterate-10 state | eb063e15f16b2610f1c5a832280d87894f567e8dab321560202e6dd9af95194c |
| Current physical J, historical R9 solve1.J0.bin | 2194965dee9c5960ba4e1c1604f76e8e85372cabef2d5ced7c8a394f593c8dc6 |
| Original initialized h0 | b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71 |

The replay reproduces E=2.911641510023811, continuous ground
minimum 6.733016805109465e-6 and scaled actual matching discrepancy
4.2429689461870394e-4. The original box and scalar inventory pass.

The matrix is bound to the initial state. The saved R10 linear-control
trial has a different state. R11 imports the verified D(q) and J_D(q)
through the historical ledger and rechecks the current source value
against its bound bytes. The residual, Norms, source value, nonpressure
actions and physical transition/source implementations pass the
preserved compatibility checks.

The problem retains 2,165 outer6 p3/p1 unknowns and the original
quadrature, alpha 400, target time .0025, eta zero, original h0 and zero
forcing. SaturatedR2 retains rho .8, beta .4, GammaScale 1,
GammaUniform zero and the original outer C(q). Viscosities, scalar
diffusion, transition, samplers, profiles and projected trace are
unchanged. Capped destruction uses Linear smoothing, theta .2 and
delta .1. Radius .125, 5% reserves, the disabled positive-value
inward-only rule, the original box and the 3,000-unknown guard apply.
Root, action and original-row gates remain 2e-12, 1e-7 and 5e-12.

RANSWallR11Experiment=On is the isolated numerical selector with
JointRestorationR5, RefinedR6, ScalarFeasibleR3 and VerifiedR3Action.
The fixed objective fraction is f_star=.125.

## Retrospective prediction and anchor ranking

The native source extractor reconstructs
$$
 M_q(fd)=R(q)+fJ(q)d+D(q+fd)-D(q)-fJ_D(q)d.
$$
Eight source-only values, including the base, cover the nine replay
records through cache reuse. This work precedes all new native search
trials and takes 3.89 s.

| Saved R10 trial | Actual E | Model E | Absolute E error |
|---|---:|---:|---:|
| Linear, .125 | 2.910253517943 | 2.910252782498 | 7.3545e-7 |
| Source, .125 | 2.911671604342 | 2.911671199969 | 4.0437e-7 |
| Source, .0625 | 2.911629893491 | 2.911629792332 | 1.0116e-7 |
| Source, .04541015625 | 2.911628012134 | 2.911627958733 | 5.3401e-8 |

The predeclared replay test compares governing E error with
max(2e-6, .01*abs(E(q)-E_actual)). All four pass. This finite-step
prediction test is separate from the 1e-7 directional-action gate
and the 2e-12 native root criterion.

The blockwise errors remain visible. At the linear .125 state the
streamwise momentum, normal momentum, k, epsilon and constraint
errors are 3.8022e-5, 3.6116e-7, 1.2320e-5, 9.4169e-8 and
5.9842e-6. At the source .04541015625 state they are 2.3891e-6,
2.9888e-8, 8.7733e-7, 9.0041e-10 and 5.0056e-7.
The verified model retains only the destruction nonlinearity.
Momentum and projected-trace nonlinear remainders remain in the
native evaluation. The governing k norm is predicted more accurately
than these smaller-block vector errors alone would suggest.

| R10 anchor | Retained Phi_1 | Reconstructed Phi_.125 |
|---:|---:|---:|
| 0 | 4.3210357075 | 4.2660586417 |
| 1 | 4.3190921048 | 4.2678888543 |
| 2 | 4.3171253688 | 4.2694377441 |
| 3 | 4.3154992235 | 4.2700450347 |
| 4 | 4.2992543536 | 4.2701884264 |

Every retained full-scale improvement in R10 worsens the objective at
the usable fraction. The source ranges and ordered limiter, cap,
profile and panel changes for every replay are retained in the ledger.
The calculation includes sampled K, moving profile quadrature and the
complete normalization. Lambda one does not freeze m=b.

## Four policies and the mixed tangent

L1 minimizes the full affine differential objective. N1 minimizes the
full nonlinear-destruction objective with the native-slope row from its
first tangent. Lf minimizes the affine differential objective at .125d.
Nf minimizes the nonlinear-destruction objective there, initially using
Lf when a verified Lf anchor is available.

For objective scale s, the internal objective is
Psi_s=Phi_s/s^2. On free differential rows, the nonlinear tangent uses
$$
 A_D=[J-J_D(q)+J_D(q+sa)]_D,\qquad
 r_D=M_q(sa)_D/s-(Aa)_D.
$$
All algebraic and overwritten rows retain A_C=J_C and r_C=R_C.
Lf uses the same J factorization with r_D=R_D/s and r_C=R_C.
The physical prediction is M_q(sd), including its fractional
algebraic residual.

All model anchors retain full total-displacement constraints relative
to q. The full-displacement source domain and continuous wall, the
fractional source sample, original-row scans and physical-base slope
are checked. Armijo backtracking follows a+eta(d_star-a), preserving
the verified incumbent and every accepted intermediate anchor.
The slope bound is g_E(d)<=-2e-4 E(q). L1 computes
g_E=-.0149347535214. The retained N1 direction lies on the safeguard
at -.000582328302005. Nf reuses the L1 slope with its seed.

The frozen mixed-row regression is also checked on the unfinished
Lf and Nf candidates without another solve or native trial. Their full
algebraic norm maxima are 1.45e-18 and 1.49e-18. At fraction .125,
their predicted trace block is .000371259782791366, matching
(1-.125) C(q). The fractional C identity errors are below 1.87e-19.
The incorrectly scaled C(q)/s formulation would instead leave a
.00297007826233093 defect at full d. These candidates remain unqualified
because their inequality construction is unfinished.

## Bounded results and construction limits

The predeclared limits are one affine solve per affine policy,
three source-model tangents per nonlinear policy, ten Armijo fractions
per tangent, 50 new source values per nonlinear policy, nine selector
trials per policy and a 600 s cooperative comparison cap.

| Policy | Direction construction | Actual chosen fraction | Actual E | Gain |
|---|---|---:|---:|---:|
| L1 | Established affine control reproduced | .125 | 2.910253517943 | .001387992081 |
| N1 | Three verified monotone tangents | .03125 | 2.911631752682 | .000009757342 |
| Lf | 96-cut exchange limit | — | — | — |
| Nf | L1 seed fallback, then first tangent hits 96-cut limit | .125, cached L1 state | 2.910253517943 | .001387992081 |

N1 reduces Phi_1 from 4.3210357075 to 4.2950038123,
4.2936400055 and 4.2848645362. Its moved source-action checks pass
at 8.37e-9, 1.02e-8 and 8.69e-9. It stops at the three-tangent
budget with continuing model progress. Its selected nonlinear-model E
is 2.911631709333, with governing error 4.33e-8 and k-vector error
6.59e-7. Its actual trace block is .000411499233423.

Lf's last original-row violation is .1275846026 in a retained-epsilon
box row at (1.5,0,2), element 34, outer side. Nf's is .4946742192
in the streamwise wall-normal gradient row at (-3,0,2.13178346368),
element 31. At these incoming cuts, the physical, coefficient and
transformed gaps agree within 1.4e-15. The archive preserves the active
rows, entering cut, candidate, QR and exchange data. These are unfinished
direction calculations with substantial remaining inequality violations.
The small arithmetic feasibility barrier from R5 is a separate
historical finding.

Nf uses the verified L1 total displacement after the Lf construction
failure. Its first fractional tangent then stops before a usable new
direction is constructed. Its moved-direction action qualification
therefore remains open. The selector's two evaluable seed trials are
exact cache hits. Their accepted state supplies no new nonlinear
direction result.

| Policy | Native trial queries / new states | New source values / actions | Policy cost |
|---|---:|---:|---:|
| L1 | 2 / 2 | 0 / 0 | 26.58 s |
| Lf | 0 / 0 | 0 / 0 | 29.97 s, unresolved |
| N1 | 4 / 4 | 21 / 3 | 125.10 s |
| Nf | 2 / 0 | 1 / 1 | 90.00 s, fallback |

Each nonlinear policy is charged its seed construction and its own
model and trial work. Nf's cost includes the failed Lf construction,
the required L1 seed and original evaluation cost for the cached
native states. Shared work is reused in the study. The per-policy
counterfactual costs consequently do not sum to the study wall time.

The comparison phase takes 196.52 s. Retrospective replay takes
3.89 s separately. The native process exits normally after 236.49 s,
including initialization and final reporting. There are six distinct
new native states and two exact reused trial states. The native R8
selector retains its fitting, held-out, proposal and fallback rules.
No full f=1 native trial is added.

## Selection, promotion and audits

The selector explicitly publishes its chosen fraction, state hash and
residual hash for each policy. Collection and accounting resolve that
exact tuple. Lf has no chosen native state. L1 and Nf each have
multiple eligible seed trials, and N1's selected result is the .03125
fallback after the .03662109375 proposal fails sufficient decrease.

N1 fails the gain, absolute-gain and efficiency gates against L1.
It passes the ground and trace gates. Nf retains the same gain as L1
while costing more and fails promotion. Lf has no eligible native
trial, so the matched Nf/Lf value comparison remains unresolved.
No policy passes its declared promotion tests. No continuation is run.

The retained L1 state and residual are byte-identical to the independently
audited R10 retained state with the same h0. That complete audit is
reused with verified bindings. The Nf seed result reuses the same audit.
The newly selected N1 proposal receives an independent audit tied to
its selected tuple. Its independently reconstructed E is
2.91163175268201. The finite scalar inventory and fixed box pass.
The independently reconstructed cubic ground minimum is
6.5331327285e-6 at x=1.5638412564. Source and same-state energy defects
are 1.44e-15 and -3.73e-14. This successful read-only audit takes
242.21 s.

The first read-only N1 audit stopped with KeyError('phase') when the
range check re-entered domain reconstruction with a minimal observation
table. The corrected adapter reuses its already verified scalar-domain
result within that audit. Its failed log and script are preserved.
This analysis failure creates no additional native search trials.
Successful audit timings and this separate failed attempt are recorded;
a precise combined wall time for the failed audit is unavailable.

The isolated authoritative Linux build uses 16 compiler workers.
Native runs use one rank and one thread each for OpenMP, OpenBLAS, MKL,
BLIS, vecLib and NumExpr. The active separate Picard work was checked and
left running within the total 72-worker ceiling. The three focused
native CTests, nine predictor-binding tests and five selected-candidate
tests pass.
Compile-time checks cover mixed-row scaling, fractional C prediction,
directional scaling and incumbent retention. The native runtime
records verify the full and fractional tangent identities and
original-unit slope scans.

The retained state has E=2.910253517943, above 2e-12, and matching
trace discrepancy .000377244015022. Its finite scalar inventory and
continuous cubic ground trace pass. The ground minimum is
5.9334752793e-6 at x=1.5638416679. The fixed box passes with
gradient maximum .0009999849430 and retained/projected epsilon minima
.9800000328669 and .9800000329396. Measured covariance minima are
.0001190994061 and .1768770013736; minimum lambda is one.
Original startup slots and the native accepted-state table are verified.
Zero physical time steps are accepted.

## Recorded commands and reproduction

The producing executable SHA256 is
98a7e65d801baf270e4c8f17e10f4a768d6634e9df11fe39820d71d4faec8afd.
The ledger pins the build receipt, runtime fingerprint, XML and
supervisor receipt with complete argv, environment and logs. The
following abbreviated invocation uses the archived epoch; a rerun
requires a fresh external epoch and the same bound historical payloads.
It must not overwrite this evidence.

```bash
r11_epoch='/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r11/epoch-001'
r11_evidence="$r11_epoch/source/docs/plans/cls-epawf-revised-evidence"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python3 -B "$r11_evidence/m8-r3/build-epoch.py" \
  --source "$r11_epoch/source" --build "$r11_epoch/build" \
  --prefix /home/jay/.local/opt/nektar-v2026-mpi \
  --output "$r11_epoch/build-run" --jobs 16 \
  --expected-source 136b09ef38b26fd0817b9dc427fbf03d0482fcec \
  --target IncNavierStokesSolver --target RANSNonpressureActionProbe \
  --target RANSEquilibriumTransitionProbe --target RANSEquilibriumTimeProbe \
  --target RANSEquilibriumEvolutionProbe
python3 -B "$r11_evidence/m8-r11/experiment.py" prepare \
  --epoch "$r11_epoch" --out "$r11_epoch/local-bound"
python3 -B "$r11_evidence/m8-r3/run-stage.py" \
  --epoch "$r11_epoch" --session "$r11_epoch/local-bound/case.xml" \
  --manifest "$r11_epoch/local-bound/predictor.json" \
  --output "$r11_epoch/local-run" --phase r5-local --seconds 650 \
  --expected-source 136b09ef38b26fd0817b9dc427fbf03d0482fcec \
  --cpu 0 --metadata-cpu 50
taskset -c 51 ctest --test-dir "$r11_epoch/build" -j 1 --output-on-failure \
  -R '^(RANSEquilibriumTimeNative|RANSEquilibriumEvolutionNative|RANSEquilibriumWallStageSmoke)$'
python3 -B "$r11_evidence/m8-r3/test-predictor.py"
python3 -B "$r11_evidence/m8-r11/test-selection.py"
```

The runtime fingerprint and isolated geometry-library setup follow
the preserved R3 epoch workflow before preparation and native execution.
Postprocessing uses the final `analysis/experiment.py summary`,
`analysis/audit-selected.py`, `analysis/check-fractional-rows.py` and
`analysis/review.py` copies. They leave the producing source untouched.
The committed R11 drivers reproduce their collection and selection
logic; their evidence pins identify the exact copies used. No new
search trials or matrices are needed to rebuild the review tables.

## Separate conclusions and next decision

| Question | Bounded result |
|---|---|
| A. Slope row from the first tangent | Better full-scale model reduction within three tangents, but no competitive native update |
| B. Affine objective at f_star | Unresolved because Lf did not finish direction construction |
| C. Added nonlinear value at that scale | Unresolved because Nf's first tangent did not produce a new usable direction |
| D. Coupled root or physical step | Neither obtained |

The scale mismatch in anchor ranking is measured. The practical
trial-scale test is limited by constrained-direction construction.
Further iterations or a larger exchange allowance have not been
assumed to resolve it, and no stationary or infeasibility conclusion
is drawn.

The next task should diagnose the retained Lf active-set subproblem
and obtain one independently verified Lf direction before spending more
nonlinear-model or native trajectory work. Preserve the original
constraints and tolerances while identifying why the 96-cut exchange
bound is exhausted.

M8 remains open. Equilibrium preservation, constitutive applicability,
nonlinear convergence, full-element positivity and long-time physical
robustness retain separate status. Curved-wall MMS, active-cap and
variable-height coupling, fully discrete positivity, refinement,
long-time robustness, scalability and physical floor/hill validation
remain outstanding.
