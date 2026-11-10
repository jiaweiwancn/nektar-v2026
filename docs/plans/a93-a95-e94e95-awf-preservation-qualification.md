# Native audit and neutral-reference preservation qualification

Starting revision: `93b6c6b31a58a976ce4462526114243646e6552d`, clean local branch
`codex/a93-a95-e94e95-awf`; remote fetch agrees. Previous primary native binary:
`1950dedf4e04a903a3e919b64bc3ce7a25c25358`. The numerical source hashes at the
reviewed head agree with that binary's manifest. Its immutable private runtime
is reusable for before/after controls, not relabelled as a new build.

Mac source worktree remains `.worktrees/a93-a95-e94e95-awf`. Linux work is under
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-preservation-20261110`, with
separate source, build, runtime and experiment subdirectories. The older A93
and active CLS directories are read-only dependencies/evidence. Initial Linux
conservative process census: 18 relevant threads. Build allocation: eight;
native jobs initially serial, MPI ranks 1, all numerical-library thread counts
1. Additional concurrency requires a fresh census and remains below 72.

## Contract and defects

The first-work contract and supplied TeX specifications govern unchanged.
Their SHA-256 remain respectively
`056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64` and
`17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40`.
No continuum coefficient, source, wall calibration, epsilon semantics or
one-sweep ordering is changed. No anchoring, cancellation forcing, field reset,
radicand clipping or new positivity device is permitted.

Audit corrections:

* Required cutoff checks distinguish accepted budget profiles, active EP
  corrections and inactive EP baselines. Actual scalar/constitutive/covariance
  checks remain independent. Activity, fallback and limiting reason are logged.
* New-model forward audits omit unsupported tangent stencils. The exact A95
  near-boundary physical witness goes through the production diagnostic helper;
  unavailable tangent evidence is explicitly not an ellipticity certificate.
* All covariance callers use one scale-normalized symmetric eigensolver,
  including an exactly zero branch. Hard covariance eligibility uses the same
  relative reference-k tolerance under coherent rescaling.

## Predeclared numerical gates

Each result distinguishes execution, expected rejection, analytical identity,
representation accuracy, quadrature uncertainty, spatial convergence and
temporal order. Missing coverage or a failed required check returns a nonzero
process status after available reports and reproducing XMLs are preserved.

Analytical weak balances use native test-function values/gradients and physical
integration points, with independently evaluated physical constitutive laws and
wall profiles. Scaled weak defects must be below 1e-9 (absolute unscaled floor
1e-12). Shared interface rows occur once. Accepted-wall source and diffusion
moments are retained separately; they need not separately vanish.

Frozen-state quadrature series: wall-normal 32/48/64, tangential and outer-volume
8/12/16, corresponding mass and interface 8/12/16. The coefficient state and
polynomial spaces stay identical. Radial reference projection uses its own
reported native Gauss32 rule and geometric subpanels. Require the last
quadrature change below 1% of an interpreted spatial-error change, with a
scaled absolute floor 1e-10; otherwise classify the spatial conclusion as
quadrature unresolved. Common evaluation-grid uncertainty is checked separately.

Independent h series: three nested radial levels at p=3, with declared physical
breakpoints (cutoff/intercept/roughness/H) aligned in the base mesh. Independent
p series: p=3,4,5 on one fixed mesh, with pressure degree p-2. Values, gradients,
alpha=0 spatial residual groups and pointwise interface quantities use fixed
physical scales. A resolved observable requires sustained reduction at both
increments (ratio <=0.85), beyond five times quadrature/evaluation uncertainty.
Near-zero observables are reported as floors, not convergence orders.

Temporal triplets use dt=1e-4,5e-5,2.5e-5 at T=.0016 with identical initial
coefficient states, fixed quadrature/spaces and consistent native histories.
Comparisons use matching physical grids and weighted RMS, plus one common
global pressure gauge. Differences above five times their uncertainty use the
initial first-order band 0.8<=p<=1.2; smaller signals are unresolved. A targeted
finer dt or bounded extension may diagnose a specific failed/unresolved order.
Perturbation response is measured against the unperturbed numerical trajectory
at the same time, separately from reference drift; no monotone-decay criterion
is imposed.

The six natural-kappa model/wall combinations receive temporal triplets. All
twelve model/wall/policy configurations receive selected spatial and short
native checks; prescribed-kappa temporal coverage is explicitly tabulated.
Capabilities not exercised are not inferred from successful endpoints.

## Results

In progress. Before/after witnesses, frozen-state geometry/basis packets,
uncertainty and convergence tables, runtime/source identities and final gates
will be recorded here and in the adjacent compact evidence directory.
