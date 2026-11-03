# M8 compact proof and diagnostic evidence

## Material Passport

- Origin: approved M8 model proposal and academic-research-suite experiment-agent workflow.
- Date: 2026-11-03.
- Mathematical/implementation verification: **completed for the receipts listed here**.
- Broad robust-model qualification: **FAILED for the screened variants on the declared state box**.
- Physical terrain accuracy: **not qualified by these experiments**.
- Authoritative numerical environment: Linux through Tailscale, `jay@100.105.247.110`.

Passing a verification receipt means that its identities, numerical comparisons,
or explicit failure witnesses were reproduced. It does not turn a failed model
property into a passing qualification. The milestone report owns the complete
native/regression campaign and final integration status; this index describes
the compact mathematical, reduced-mode and state-coverage evidence.

## Files and interpretation

| Evidence | What it establishes |
|---|---|
| [Fresh kernel qualification](proof-kernel-fresh.json) and [runner](proof-kernel-fresh-runner.json) | 37,310 independent checks, including 27,907 native probe checks; explicit common-cubic scaling, consistent modified production, capped-source identities and one-sided actions, required invalid-input rejection. |
| [Kernel replay](proof-kernel-replay.json) and [native provenance](proof-native-replay-provenance.json) | 37,309 checks with the later interval-proof code and the same stored native report. One native-command exit check is absent because this is a replay. Original native execution `07363199` is distinguished from the AA compiled campaign and `11aa94c0` proof source. |
| [Transverse qualification](transverse-validation.json) and [frozen plan](transverse-plan.json) | 305,254 checks of the reduced Fourier/native-CG three-component transverse subsystem, with independent stress, weak matrices, native histories, energy, covariance and refinement oracles. Local decay is qualified; coupled/terrain stability is not. |
| [Joint local proof, AA](proof-local-joint-aa.json) | Original local momentum/covariance certificate used by the recorded native state-coverage audits. |
| [Joint local proof, 11aa](proof-local-joint-11aa.json) and [runner](proof-correlated-runner.json) | Sharper directed-rounding coefficient enclosures, the same declared local boxes, fixed-gamma covariance proofs, and an explicit outside-envelope control. |
| [Conservative wide bound](proof-conservative-wide.json) | Sufficient gamma values from conservative interval bounds. These values are **not necessary or minimal**; a failure of an overdiffusive sufficient choice alone would not rule out every gamma. |
| [Fixed-gamma feasibility proof](proof-fixed-gamma-feasibility.json) and [plan](proof-fixed-gamma-plan.json) | Separate exact rational witnesses impose contradictory **necessary** lower/upper gamma constraints. This stronger result rules out the specified fixed scalar gamma law on the declared broad box. |
| [Invariant-response analytic controls](proof-invariant-controls.json), [runner](proof-invariant-controls-runner.json), [independent screen](proof-invariant-screen.json) and [plan](proof-invariant-plan.json) | Explicit directed-interval constitutive/covariance failures of the two optional C(q) research variants, with all C′ and β′ terms included. No V3 production model was adopted. |
| [p3 reference coverage](coverage-reference-p3.json), [perturbed coverage](coverage-perturbation.json) and [runner](coverage-runner.json) | Independent reconstruction at the actual declared outer source, mass, reconstruction, interface and transfer evaluation sets. These p3 cases leave the local certificate box. |
| [p5 coverage](coverage-reference-p5.json), [p7 coverage](coverage-reference-p7.json), [refinement runner](reference-refinement-runner.json), [p5 energy](reference-p5-energy.json), [p7 energy](reference-p7-energy.json) | Actual zero-step reference projections with independent energy checks. p7 covers the finite outer certificate sets and has positive sampled inner covariance; this is an initial-state result. |
| [Runtime epoch summary](runtime-epoch-summary.json) and [build runner](runtime-build-runner.json) | Compact source/runtime binding for the producing build epoch. The complete fingerprint remains on Linux; this summary does not replace final source/runtime equivalence checks. |

## Reduced transverse experiment

The figure selects the actual `p7-n8-xi4`, IMEX2, amplitude `1e-4` rows at
`dt=0.00015625 s`, `T=0.05 s`. The disturbance is a reduced Fourier/native-CG
experiment supported above the wall layer, not a nonlinear terrain trajectory.

| Model in the reduced experiment | Final energy ratio |
|---|---:|
| Original CLS, reused M7 control | 8.62184899621 |
| Common cubic β=0.4, Γ=0 | 0.589289214077 |
| Common cubic β=0.45, Γ=0 | 0.734610537997 |
| Common cubic β=0.4, Γ=0.25 | 0.333442685761 |

Here Γ multiplies the prescribed neutral eddy viscosity in the persistent
reference-preserving stress increment. Γ=0.25 is a reduced-mode sensitivity
case. It is distinct from the Γ=1 native reference-projection cases and from
the fixed-gamma feasibility problem. The original control is the preserved
[M7 receipt](../m7/transverse-validation.json), not a new M8 solver execution.

## Local certificates and broad fixed-gamma failure

Use the upstream orthonormal frame and normalization
`A*=A/G0`, `k*=Qstar*k/K0`, `epsilon*=epsilon/epsilon0`, where
`G0=u_star/(kappa*r)`, `K0=Qstar*u_star^2` and
`epsilon0=u_star^3/(kappa*r)`. The tangent has physical units νt,0 and
the Reynolds covariance has physical scale `u_star^2`.

Both limiting max-invariant branches are covered by the interval certificates.
The following 11aa bounds use `k/K0, epsilon/epsilon0` in `[.99,1.01]` and
entrywise gradient distance from `e_x tensor e_z`:

| Fixed Γ | Gradient radius | β=0.4 lower principal bound / νt,0 | β=0.45 lower principal bound / νt,0 | Covariance |
|---|---:|---:|---:|---|
| 1 | .001 | .471272810698 | .359364129975 | Strictly positive throughout the declared box by interval LDL |
| 3 | .01 | .547381238192 | .252563404664 | Strictly positive throughout the declared box by interval LDL |

These are local constitutive properties with frozen k and epsilon. They do not
prove state-box invariance, outer turbulence positivity, a high-order discrete
maximum principle, long-time coupled stability, or terrain accuracy. Covariance
positivity for one fixed Γ does not automatically extend to larger Γ.

The wider feasibility box has entrywise gradient radius `.25`, exactly
trace-free witness gradients, and `k/K0, epsilon/epsilon0` in `[.5,1.5]`.
The selected exact rational witnesses imply these necessary constraints:

| β | Momentum requires Γ at least approximately | Covariance requires Γ at most approximately |
|---|---:|---:|
| .4 | 2.44640920758313 | .177525831861587 |
| .45 | 2.27794509320359 | .156117116286741 |

The intersections are empty. The actual directed interval endpoints, rational
states and directions, branch checks and gamma-shift identities are retained
in `proof-fixed-gamma-feasibility.json`. The plotted intervals are necessary
conditions, not asserted sufficient or safe regions. This result concerns
the specified fixed scalar stress increment on this box, not every possible
constitutive regularization or global model impossibility.

## Optional invariant-response alternatives

Both research alternatives retain the original C(q) up to `qc=3.3` and use
the C1 power tail `C(qc)*(q/qc)^(-alpha)` above it, with
`alpha=-qc*Cprime(qc)/C(qc)`. Candidate A keeps common β constant; B uses
`beta0*[qc*C(qc)/(q*Ctail(q))]^2` above the join. B is the explicitly
specified research completion of the adaptive-coefficient direction. Neither
alternative was adopted in production.

Candidate A has an explicit adverse transverse coefficient at q=8. Candidate B
repairs that selected direction, but the independent strict-active screen
finds other negative incompressible principal forms for both beta0 values.
These forms fail the declared constitutive coercivity gate; they are not
new native coupled PDE trajectories or nonlinear instability theorems.

There is also a simple common covariance witness for both alternatives:
`A02=5/4`, other gradient entries zero, `k=19/4`, `epsilon=1/2`,
`q=95/8`, and covariance direction `(1,0,5)/sqrt(26)`. It lies within the
declared broad box and gives variance `-0.2832467961122523`. The screen's
boundary control with `k=1.5*Qstar` is a different state and gives approximately
`-0.310052142763282`; the two results are consistent.

## Initial projection and finite coverage

These results use β=.4, Γ=1 and the same physical model inputs:

| Degree | Maximum normalized gradient-entry deviation | Minimum inner covariance at source quadrature | Finite outer joint-certificate coverage |
|---|---:|---:|---|
| p3 | .0642009880560146 | -.177666349051788 | Outside |
| p5 | .00264835174984812 | -.0198538557236063 | Outside |
| p7 | .000102433249649092 | +.0806939841621756 | Covered at the declared finite sets |

The p3 figure row selects t=0 from the stored anchored control. Its subsequent
control coefficients are exactly identical, so the stored aggregate coverage
extrema equal its initial projection. The p5/p7 sessions have zero steps.
This separates resolution of the projected reference from the broad
constitutive failures. A positive p7 initial sample does not qualify an evolved
trajectory or guarantee the minimum over every point of the polynomial.

## Figures, curation and reproduction

- [Three-panel PNG](m8-qualification.png) and [vector PDF](m8-qualification.pdf).
- [Exact plotted values and input hashes](figure-data.json).
- [Figure receipt](figure-receipt.json) records the plot-source hash, producing
  data commits, versions, one-worker environment, output hashes and actual command.
- [Plot/curation source](plot-m8.py) reads actual immutable JSONs. It executes no solver.
- [Curation index](curation-index.json) records each full original path, byte
  count and SHA256 alongside the curated-file hash and declared transformation.

The final scientific figures are generated on Linux under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/figures-002`.
The initial layout attempt remains under `figures-001`. The final compact
archive is `curated-proof-002`. All rendering and mathematical proof jobs use
one worker with numerical library threads one; the figure run is CPU 48.

To render the exact extracted data again on Linux, use a fresh output path:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  taskset -c 48 python3 -B plot-m8.py --data figure-data.json --output FRESH_FIGURE_DIRECTORY
```

Large native fields, complete histories, solver logs and full 9x9 interval
tangent arrays remain in the Linux run directories. The three local/wide
proof JSONs remove only the recursively named `tangent_intervals` arrays and
add explicit curation metadata; their exact scalar bounds, null-form weights,
interval LDL pivots, states and covariance witnesses are retained. Other
listed proof/coverage receipts are copied byte for byte. The runtime summary
declares its omitted full manifests and the authoritative original hash.
