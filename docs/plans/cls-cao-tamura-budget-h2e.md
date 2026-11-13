# Cao–Tamura budget AWF H2e

13 November 2026, Asia/Shanghai. **Initialization remains blocked by covariance.**
The mass-conditioned velocity control and the conditional outer-scalar search
were executed and independently re-audited on Linux. No physical hill attempt
was made. Physical equations, calibrated roughness, floors, closure parameters,
physical boundaries and the 3000-slot guard are unchanged.

## Fixed problem and provenance

The Mac branch worktree is `.worktrees/pimple-one-sweep`, branch
`codex/cls-epawf-pimple-one-sweep`, starting at `29b6da5`. Tested implementation: `5cb7ce4d583dd864f443670e032ee9019dc48583`. The shared Mac checkout
and the unrelated Linux checkout/work were preserved. Linux storage is on UUID
`5AB8A669B8A642FD`, mounted at `/media/jay/新加卷1`. The isolated run root is:

```
/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2e-20261113
```

[Policy](cls-cao-tamura-budget-h2e-evidence/policy.json),
[native reduction](cls-cao-tamura-budget-h2e-evidence/native.json),
[verification](cls-cao-tamura-budget-h2e-evidence/verification.json),
[source/runtime identities](cls-cao-tamura-budget-h2e-evidence/source-revisions.json)
and [artifact index](cls-cao-tamura-budget-h2e-evidence/artifact-index.json)
provide the reproducible record. Heavy arrays, candidates and native output stay
on Linux. Source archives/overlays and build logs preserve the intermediate
probe identities as well as the final runtime.

H_A=z_m=z_K=.0065 m; h_s=.0022472751006336437 m; the seven-point mean fit,
u_tau/viscosity, AnalyticalNeutralInflowV1, Regularized CLS beta=.4/Gamma=1/
gamma0=0, Linear inner stress and calibrated budget cutoff are retained.
Endpoint K is transported in evolution. It is frozen only during this controlled
initializer experiment. The corrected PhysicalHeightInterfaceV2 base has hash
`83b1a4fd9a4fcda694f6e6f10bfb789700d1fc6fc83443b91f67bfcb045e79b5`.
The H2d best seed was not the starting point; its excessive divergence and all
historical rejected states remain archived.

## Native basis and physical constraints

The 32 existing essential/weak-continuity-nullspace directions were imported
from the bounded prefix of the H2d context. They were checked against the newly
constructed native physical mass and continuity constraints. No old search was
rerun. Symmetric eigendecomposition of G=Z^T M Z retains eigenvalues greater
than 1e-10 times the maximum. All 32 directions survive: eigenvalue range
9.23777e-6–7.00104, condition number 757871 before conditioning, and maximum
mass-orthogonality defect about 1.27e-11 afterward. There is no outcome-based
rank truncation. Old physical increments, gradients, stresses and mass distances
are reconstructed in the new coordinates; frozen rows remain exact.

The velocity trust metric is physical L2 distance in mass-orthonormal
coordinates. The unchanged limits are .12 in that distance and .5 m/s for every
sampled velocity component increment. The .12 scale keeps H2d's nominal
flat-volume/Uhill convention; actual curved-volume/RMS measures are reported
separately. Nearly active response rows enter the local subproblem; newly
violated rows are added, and every actual trial still receives a complete cap
scan. No observation is removed.

Strong divergence uses a physical-volume GL16 quadratic, with actual native
basis derivatives, rather than inventory counts. Its bound is
base_L2*(1+1e-6)+1e-9 = .0603080994. The base is .0603080381. Native weak
continuity remains an independent constraint. A gamma-weighted strain-distance
predictor and a soft roof streamwise-gradient penalty use the same physical
quadratures. No roof tangential Dirichlet condition or production grad-div term
was introduced.

The phase-I objective is the true maximum normalized covariance deficit plus
the predeclared small physical-distance/strain/roof penalties. The desired
interior margin lambda_min/K0>=1e-4 is distinct from the unchanged native hard
tolerance. Local eigenvalue linearizations plus slack are nonconvex-search
subproblems, not feasibility certificates. Every trial recomputes the true
selected tensors.

Two initial implementation checks exposed poorly scaled local optimizer
variables. Scaling constraint residuals alone did not fix it; solving in unit
trust coordinates did. Those two seven-major-iteration diagnostics are retained
as `stage-A-v1/v2`. Their work counts against the 20-major Stage-A allowance;
the final control therefore uses six additional major iterations, not another
20. Total Stage-A work remains below its 300-s/180-full-scan/3000-active-batch
limits. No physical limit or objective weight was changed to improve survival.
The final control's subproblems converge and its accepted changes are saved
immediately, including the last major iteration. True before/after objectives,
radii, active rows, complete scans, active derivative batches and termination
reasons are durable.

## Executed searches and independent certification

| Quantity | Corrected base | A: velocity only | B: velocity + outer k/epsilon |
|---|---:|---:|---:|
| Native inner covariance minimum | historical H2d base | +.000792568 | −.001006332 |
| Native outer covariance minimum | approximately −.0534449 | −.038148745 | −.029992657 |
| Independent GL20 minimum | — | −.038162358 | −.030095454 |
| Physical divergence L2 | .0603080381 | .0601611911 | .0603080602 |
| Velocity increment L2 | 0 | .03251855 | .05350842 |
| Maximum inventory component increment, m/s | 0 | .24508835 | .5 |
| Complete scalar/budget/boundary/weak-continuity checks | pass | pass | pass |
| Covariance certificate | fail | fail | fail |

The A control ends at its allocated iteration budget. Its complete native
rejection triggers B. B uses eight smooth outer-k and eight outer-epsilon
native directions, separately conditioned in relative-field physical mass.
Their inlet/interface coefficient rows are exactly zero. Every inner-k row,
endpoint functional, prescribed row and pressure coefficient stays fixed,
including signed-zero bits. Smoothness and near-reference penalties act on the
actual projected scalar directions; no positivity claim is based on the
unprojected analytic shape.

Each coupled trial consumes fresh native k/epsilon values, recomputing CLS
invariants, C_mu, viscosity and stress. gamma/S_ref remain immutable at physical
incoming height. The native validator independently confirms identical moving
inventory coordinates/cardinality, cutoff data and dependent epsilon, rather
than assuming frozen K was sufficient. Actual transported and dependent epsilon
remain distinct. The old H2d velocity-only validator is unchanged.

B completes 40 major iterations, 152 full tensor scans and 3880 active batches
in 92.97 s, below 900 s. Local subproblem iteration limits occur near the many
active scalar/response constraints; only actual fully checked improving trials
are retained. This is a work-limited local search, not an incompatibility
certificate. The terminal accepted state is saved and independently reconstructed
bitwise. The complete native scalar envelopes are approximately k/base
[.500091,2] and epsilon/base [.540098,1.99936]. They were not widened.

B's whole-domain initial TKE integral rises from .01120156 to .01237889
(+10.51%); the outer integral rises from .01109927 to .01227659. Independent
sampling gives k/epsilon time-scale ratios [.5193,2.0768] and
k^(3/2)/epsilon length-scale ratios [.3675,2.9274]. Velocity kinetic energy
changes from .5643663 to .5644814. The roof sampled delta-U maximum is about
.04695 m/s; the far-field and roof-gradient measures are in the native reduction.
Divergence RMS is .0756218/s and the GL20 local maximum is 15.6334/s. Small weak
continuity defects do not imply small pointwise divergence.

The native B inner failure is at (s,r)≈(.235496,.0065), downstream on the
matching plane. Its outer failure is near (s,zHat)≈(−.022625,.135920).
Independent GL20 sampling finds (physical x,z)≈(−.025588,.167338), with k=.0349818,
epsilon=.00821337 and lambda_min=−.03009545. The full physical gradient,
geometry/element, tensor, eigenvalues and independent critical eigenvector are
retained. The indexed tensor replay agrees to implementation precision.

At that independent B witness, the contributions along the critical eigenvector
are isotropic +.0233212, linear −.00613986, quadratic −.00684528, beta-weighted
cubic −.000646423 and regularization −.0397851. Their signed sum is the negative
eigenvalue. These are local tensor contributions, not additive percentages of
an integrated flow error. At fixed limiting velocity gradient, a 41-by-41
scalar-envelope diagnostic followed by a local bounded optimization still gives
about −.0299086, at k/base=2 and epsilon/base≈1.21247. This does not exclude a
different velocity field or scalar spatial distribution. The saved H2d roof
same-state equal k/epsilon scaling diagnostic reproduces 5.812321; that multiplier
was never applied to a field.

The evidence supports insufficient mean-strain improvement within this bounded
span and active response/scalar limits. It does not establish global physical
infeasibility, nonexistence of a hill solution, or loss of ellipticity.

## Exact-state, boundary and runtime binding

`OuterScalarInitializerV1` declares native permitted rows and rejects any change
to inner/interface/inlet/pressure/frozen signed-zero coefficients. Inlet and
interface functionals are checked independently using immutable-base boundary
data. Native fixtures reject inner k, trace epsilon, pressure, inlet,
signed-zero and old-frozen-scalar-contract violations. Zero scalar amplitude
reproduces the old validator's state, complete audit, scalar inventory, physical
samples, constitutive diagnostics and native budgets bitwise.

`CertifiedNeutralStartupV2` preserves the immutable base reference/defect and
boundary data separately from the candidate. It binds state/base/frozen-row
hashes, allowed rows, basis/context/envelope identities, geometry/fit/reference/
inventory signatures, source provenance, native report/receipt, validation probe,
solver executable and linked-library closure. A different runtime requires
exact-state revalidation and a new certificate. Import performs no optimizer,
Lift or projection and forms fresh IMEX1 histories. Ordinary incompatible
restart checks remain strict. V1 remains supported with its original signature.

Both rejected hill receipts are refused by the certificate issuer. A positive
affine verification fixture passes exact seed/reference/history import,
bitwise same-dt governing stages and solved epsilon, split restart, deliberate
post-assignment rollback, parent-checkpoint preservation, and rejection of
changed dt and corrupted runtime/receipt/boundary identities. Initialization
algebra is not counted as physical solves. These affine checks are not a curved
hill startup result.

## Cost, disposition and next package

Native context construction is 39.08 s, 46.78 s inclusive, peak about .765 GiB.
Binary64 arrays and membership records occupy about 330 MiB instead of the old
789-MB JSON. There are 181975 unique evaluations and 433207 retained inventory
memberships. Native element, geometry, both sides, coordinates and epsilon
context are part of the reuse identity. Fixed-K profiles/cutoffs are reused only
after native invariance checks. No full context is exported/reparsed per trial.

A full tensor scan is about .34 s. Search peaks are approximately .89 GiB,
versus the historical H2d 4.85 GiB. These are measured initializer costs, not a
matched production-solver speedup. Active derivative batches are reported
separately from full scans. Affine regression costs and exact runtime identities
are in the verification receipt; **there is no measured physical hill-step cost**.
Builds use four workers; pilots and numerical libraries use one rank/thread.
Other project jobs and frozen historical runtimes were left intact.

| Gate | Result |
|---|---|
| Mass conditioning, old-span reconstruction, fixed bits/weak continuity | passed |
| Physical divergence and active cap implementation | passed on saved candidates |
| Fresh outer scalars, fixed endpoint/cutoff/trace, native envelope checks | passed |
| Complete A/B covariance and independent extra sampling | failed |
| Exact-state/runtime/boundary provenance and affine transactions | passed |
| Conditional 8/16 hill steps and curved 4+4 restart | unexecuted: no certificate |
| Flat 1e-3 spatial target | remains unclosed; no new mesh trajectory |
| Developed flow or experimental validation | unqualified |

The [assembly plan](cls-cao-tamura-budget-h2e-evidence/assembly-plan.md) identifies
local CSR/block assembly and exact sparse physical-mass applications with a
small dense pressure Schur as the bounded route toward 4399–5339 slots. Estimates
are distinguished from measurements; no guard was raised and the rejected
single-knot transfer was not adopted.

**Next recommendation:** one bounded enrichment of the continuity-preserving
velocity span localized around the outer slope/transition witness and downstream
matching plane, with explicit protection of initially admissible inner tensors.
Keep the present divergence/response limits and scalar envelope for that test.
The local scalar diagnostic points first to the velocity strain, rather than a
reason to widen scalar bounds or retune Regularized CLS. This remains an
initialization experiment; affordable spatial/development qualification is still
required before quantitative hill claims.

## Reproduction commands

Use a fresh output root for any reproduction; never overwrite these receipts.
The indexed source archive starts at 29b6da5 and its overlays/hashed files identify
all tested implementations. The authoritative build used:

```sh
cmake -S "$H2E_ROOT/source/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$H2E_ROOT/build" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS' \
  '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake' \
  -DRANS_BUILD_SOLVER=ON \
  '-DRANS_PREBUILT_SPATIAL_LIBRARY=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2b-20261112/curved-lazy-final/runtime/lib/libSpatialDomains.so.5.9.0'
cmake --build "$H2E_ROOT/build" --parallel 4 \
  --target RANSWallResolutionProbe IncNavierStokesSolver RANSKBoundQPChecks
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
```

Here H2E_ROOT is the run root above; H2E_TESTS is its
`source/solvers/IncNavierStokesSolver/Tests/RANS`, H2E_PREFLIGHT is the retained
H2c `preflight-final.json`, and H2E_H1 is H1 `compatibility-resolved`.

```sh
python3 -B "$H2E_TESTS/run_cao_tamura_budget_h2e.py" context \
  --probe "$H2E_ROOT/build/RANSWallResolutionProbe" \
  --h1 "$H2E_H1" --preflight "$H2E_PREFLIGHT" \
  --old-context '/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2d-20261113/setup-v1/65-hill-probe/result.json' \
  --output "$H2E_ROOT/context-v1"
python3 -B "$H2E_TESTS/search_cao_tamura_budget_h2e.py" \
  --context "$H2E_ROOT/context-v1/context.json" --preflight "$H2E_PREFLIGHT" \
  --stage A --max-major 6 --output "$H2E_ROOT/stage-A-final"
python3 -B "$H2E_TESTS/run_cao_tamura_budget_h2e.py" certify \
  --probe "$H2E_ROOT/build/RANSWallResolutionProbe" \
  --context "$H2E_ROOT/context-v1/context.json" --preflight "$H2E_PREFLIGHT" \
  --candidate "$H2E_ROOT/stage-A-final/best-candidate.json" --output "$H2E_ROOT/certify-A"
# B was authorized only after the complete A covariance rejection:
python3 -B "$H2E_TESTS/search_cao_tamura_budget_h2e.py" \
  --context "$H2E_ROOT/context-v1/context.json" --preflight "$H2E_PREFLIGHT" \
  --stage B --start "$H2E_ROOT/stage-A-final/best-candidate.json" --output "$H2E_ROOT/stage-B"
```

The B native command is the same `certify` operation with B candidate/output.
`analyze_cao_tamura_budget_h2e.py --root ... --preflight ... --old-native
.../source/docs/plans/cls-cao-tamura-budget-h2d-evidence/native.json --output ...`
reconstructs each saved state exactly and independently contracts/eigendecomposes
the tensors. `verify_cao_tamura_budget_h2e.py` takes `--solver`, `--probe`,
`--context`, `--h1`, `--preflight`, `--output`; final output is
`verification-final`. The indexed native input packets, receipts and logs
record complete actual commands and runtime paths.

The final affine two-step exact-import block is 32.91 s inclusive with peak
1542140 KiB; its zero-step import is 18.44 s. These include startup, audits,
checkpoint/packet I/O and teardown, not exclusive equation-solve times.
`RANSKBoundQPChecks` and the two focused branch-aware Python tests pass.
Native tensor discrepancies are around 7e-18, whereas the B negative eigenvalue
is about −.03 (hard covariance tolerance about −1.8e-12). The rejection is far
larger than replay roundoff. Signed cancellation in local stress contributions
is preserved throughout.
