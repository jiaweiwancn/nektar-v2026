# CT14: short full-domain A95 flat-pilot qualification

Status: **IN PROGRESS; no selected-mesh evolution admitted**.

Fetched `codex/a93-a95-e94e95-awf` on 2026-11-15: actual remote and clean
dedicated Mac worktree both at `f121dd1e5dc14ed147714876ea182a913e6fcae5`.
The parent Mac worktree and Linux main checkout contain unrelated work and
will not be switched, reset, cleaned or built in place. USB SSH timed out;
the authorized Tailscale fallback works. The verified heavy-storage mount is
`/media/jay/新加卷1`. CT14 uses a new isolated experiment directory there.

The prospective [contract](cao-tamura-accepted-14-evidence/contract.json)
precedes costly execution. The frozen CT8 physical manifest remains SHA256
`004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296`.
CT8/11/12/13 results, failures, component identities and histories are retained.
The original target's failed pressure spread is 3.824614134304275e-6 against
5.71117429897371e-7. The selected 62,533-slot, 961-element P4/P2 mesh remains
at T=0. Its 96.27% initial shear-pressure-image reduction is diagnostic only.

| Gate | Fixed criterion | Intended evidence |
|---|---|---|
| Predictor capacity | 70k slots, 20M matrix nnz, 40M factor entries, 16GiB address space | Four hash-bound CT13 canonical A/RHS controls, then every evolved factor |
| Linear accuracy | infinity residual <=2e-12 infinity(rhs)+1e-25; at most three reused-factor refinements | Original-system extended-precision actions, actual/tiny signed/cancelling RHS |
| Component | scaled actions/reconstruction <=2e-10; state/absolute-pressure/traction equivalence <=1e-8 | Original and small controls, canonical equality, inverse permutation, genuine restart/full rollback/incompatible-history rejection |
| Static admission | All target native modes, alpha1000/2000, flow/scalar/consistent mass/full pressure; physical projection, gradients, signed loads, boundaries, domain/covariance/radicand/cutoff | Complete actual-state/stage-bound native and quadrature packets, CT13 evidence reused only with matching scope |
| Physical pair | Fresh dt=.001/.0005 to T=.002; all CT11 classifiers unchanged | Own-initial drift, reference error and common-time differences separately; absolute and regional pressure, current-K physical tractions, mass and four-control scalar budgets |
| Resources | <=24,000 process-CPU s, <=28,800 active-wall s, <=16GiB new output; <=16 project/72 aggregate workers, <=2 compiler jobs | Closed owned-family receipts including failures/controllers/audits/publication; one numerical family, one thread per library |
| Reproducibility | Four normal checkpoint commits/pushes when their prerequisites pass | Exact source/runtime/library/input/history hashes and remote containment |

The first numerical method is deterministic Scotch nested dissection of the
undirected structural union A+A^T with no diagonal graph edges, sorted unique
neighbors, seed zero and one fixed explicit strategy. Its precise strategy,
available headers/linkage and forward/inverse convention will be recorded
before the first factor trial. Only internal factor numbering changes; every
solution and residual returns to canonical IDs. Scalar/preconditioner COLAMD,
consistent-mass LDLT, coefficients, explicit-zero policy, pivots and scaling
stay fixed. No governing component is admitted by an ordering test alone.

Failed attempts will be preserved and diagnosed. A subsequent attempt requires
a prospective causal amendment and remaining allocation; there is no blind
ordering/mesh/timestep sweep. Static admission and a pushed short-pair contract
must pass before advancing fresh physical histories. No horizon beyond .002,
terrain, other-model displaced evolution, smooth production, EP, MPI/3D or
variable-height production is authorized by this package. Lower fill does not
establish long-run affordability or curved readiness.

## Qualified capacity and component

`CT14-VectorCOLAMDV1` qualifies the predictor backend with the unchanged
matrix/zero/mask/scale/pivot/residual/refinement contract. The opt-in
`RANSWallPredictorOrdering=VectorCOLAMDV1` orders the directed quotient of
the two velocity components sharing each native FE ID, then expands groups
in canonical within-group order. Both matrix indices are permuted internally;
every original coefficient is checked exactly after permutation and all
solutions return to canonical IDs. Scalar/preconditioner COLAMD and
consistent-mass LDLT remain unchanged. Default component histories retain
their old identity; the new component has its own fingerprint.

| Frozen system | Matrix nnz | Factor entries | Margin to 40M | Process CPU s | VmPeak KiB |
|---|---:|---:|---:|---:|---:|
| Original initial, alpha1000 | 2,081,794 | 35,438,135 | 4,561,865 | 2.787 | 1,118,888 |
| Original evolved, alpha2000 | 2,082,125 | 35,399,996 | 4,600,004 | 2.844 | 1,119,488 |
| Selected initial, alpha1000 | 2,152,435 | 37,266,312 | 2,733,688 | 2.571 | 1,176,632 |
| Selected initial, alpha2000 | 2,152,435 | 37,225,528 | 2,774,472 | 2.529 | 1,176,632 |

CPU values cover the complete ordering probe, including CSR input, permutation
checks, actual/tiny/cancelling RHS solves and output. They are not isolated LU
timings. Actual RHS infinity residuals are at most 2.885e-19; tiny residuals
at most 3.885e-31. Every case meets `2e-12*maxabs(rhs)+1e-25` with at most
three reused-factor refinements. The independent nonsymmetric 48-row vector
control and permutation inversion pass. Selected solutions are qualified
against canonical matrices/RHS and independent controls, without requiring
an old selected solve that violates its fill guard.

All four P3/P1/P4/P2 zero/displaced small complete-operator controls pass.
The largest independent action defect is **4.313e-15**, below 2e-10. Fresh
original full-domain .001 two-step evolution agrees with the qualified CT11
update at every common state: maximum scaled field difference **1.355e-14**,
including absolute pressure, and physical traction difference **4.071e-14**,
both below 1e-8. CT11's historical physical pressure failure stays unchanged.
Genuine same-component restart and full rejected-step rollback reproduce
all native histories bitwise. Old-component, new-to-old, stale-fit histories
and unknown ordering reject. No old history is relabelled.

The tested native source is `6477d9d9`. Governing executable SHA256 is
`9ae7d3c63492c82a4c9f979616f673d8e4909a5c41ee8069fc1592f2fb23da06`;
diagnostic probe SHA256 is
`95d6a0258926ab250b69ea0940976f43a993bf4bdccc3001a5819f4f36835824`.
These are private Linux runtimes; shared installations are preserved. Exact
loaded libraries, input/runtime/source bindings, complete numerical bodies,
history witnesses and closed resource receipts are in
[capacity-component-v1](cao-tamura-accepted-14-evidence/capacity-component-v1/).

Three rejected attempts remain archived: scalar Scotch ND gives 40,862,187
entries; native-vector Scotch ND gives 40,119,892; retaining its external order
gives exactly the same fill. The latter rules out Eigen postordering as the
fill cause. The directed COLAMD quotient was then selected using the measured
success of original-system COLAMD, rather than a seed/threshold sweep. Each
attempt has a prospective amendment and an immutable failed witness. All
attempts count toward the one CT14 allocation.

Capacity/component qualification admits the prospectively recorded
[static tests](cao-tamura-accepted-14-evidence/static-tests-v1.json).
**Complete selected static admission and the physical pair remain pending.**
