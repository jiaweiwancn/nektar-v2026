# CT9: serial sparse one-sweep and full-domain flat preservation

This is a prospective implementation contract, not a CT8 scalability result.
The CT8 7.5mm radial inventory supersedes the CT7 capacity *input* only; CT7
counts and blocked decisions remain immutable. No dense guard is raised here.

CT8's central displaced strip completes both histories through T=.008, but
its absolute-pressure timestep uncertainty fails: scaled spread
4.804457393280823e-6 exceeds 8.411880056796802e-7. T=.02 is withheld.
This is independent of capacity. Backend implementation/equivalence may proceed;
physical full-domain preservation must additionally close that numerical gate
under a prospectively declared, bounded pressure/time-level study. A sparse
solver is not presumed to remedy the failure. No additional CT8 timestep was
introduced to search for a pass.

The new native affine target counts are33,267(P3/P1) and60,907(P4/P2)
on930elements. A single hypothetical combined dense binary64 matrix would
occupy8,853,546,312/29,677,301,192bytes, respectively, as reported by the
native counter; no such matrix is allocated in CT8. This is a capacity
illustration, not a measured peak-RSS estimate for the present block algorithm.
The paired terrain's geometry is checked independently, but its physical-H
region classifier is not supported by the affine field-count path. A later
mapped region-aware inventory must not reuse global-y classification.

The inspected authoritative Release SDK has `NEKTAR_USE_MPI=ON` and
`NEKTAR_USE_PETSC=OFF`. `cmake/ThirdPartyPETSc.cmake` and the MultiRegions PETSc
adapters exist in source; this does not make PETSc available in the tested
runtime. CT9 must either verify and isolate a PETSc-enabled dependency build or
use the installed Nektar sparse/matrix-action facilities. It must not link an
unverified optional library or alter the shared installation.

## Reachable dense paths to replace or bound

| Phase | Present path | Required CT9 treatment |
|---|---|---|
| Native numbering | `RANSWallGlobalField`, signed native CG maps | Retain native IDs, orientation and exact coefficient restriction; sparse graph from these IDs. |
| Startup mass | `RANSWallGlobalField::Mass`, `PhysicalMass`, constructor `m_mass` | Assemble element contributions to COO/CSR, with the same physical quadrature and deterministic summation; no global long-double N by N staging. |
| Initial projection | `Project`, constrained `RadialProject`, boundary segment transforms | General field projection uses sparse consistent mass. Keep the small radial KKT and face transforms bounded by radial/trace size, never global size; verify that they remain small on the target grid. |
| Frozen flow/scalars | `AssemblePicardBlock` in `RANSWallEvolutionPicard.cpp` | Direct sparse element/face assembly and residual loads; never construct dense then convert. Preserve roughness-volume pressure coupling to actual outer H rows, dependent epsilon constraints, top and outlet loads. |
| Predictor | PIMPLE dense A factor | Nonsymmetric sparse solve with true unpreconditioned residual receipt; retain the same frozen matrix and RHS. |
| Mass inverse | `Factor(m)` on free velocity rows | Consistent SPD free-mass solve, not diagonal/lumped mass. Reuse immutable symbolic pattern only under identical free-ID/mass identity. |
| Pressure correction | Dense `bg`, full dense Schur, dense factor | Implement **action** `p -> D[(alpha M_free)^-1(Gp)] - Cp`; no storage of all M-inverse-G columns or dense pressure matrix. Preserve absolute outlet traction and exclude only the existing unused gauge slot. |
| k/epsilon | Frozen scalar matrices/factors | Sparse conservative advection/diffusion and existing signed source splitting; flow/epsilon/k ordering and frozen prescribed sink unchanged. |
| Mandatory checks | `ScalarDomain`, covariance, A95 radicand, cutoff, continuous ground trace, outlet backflow | Keep every consumed evaluation and its current full dependencies. Row geometry may be immutable; constitutive validity is state/gradient dependent. No reduced inventory to save cost. |
| Norms/audits | `Norms`, diagnostic restricted-mass factors, regional physical masses, resolution components | Sparse restricted consistent-mass solves or direct physical integration. Keep free IDs, reactions, full/common-test definitions and equation coverage. Disable only optional global dense matrix *exports* on large cases with an explicit message, not mandatory physical checks. |
| Snapshots/replay | Dense flow/mass JSON and full-row packets | Small fixtures keep exact dense exports for equivalence. Large cases use selected-stage sparse triplets/actions, vectors, hashes and bounded streaming inventories; native state/history/boundary identities remain strict. |
| Nonlinear/QP/research paths | Dense Jacobian, Gram/KKT and finite-difference routes | Remain guarded as unsupported under the new backend; original bounded dense implementation is retained for their historical scope. |

The pattern must include element-local terms, one-sided interface rows,
roughness sampling links to the outer H pressure/k rows and all natural faces.
Changing the top from essential to natural changes the free-ID set and thus
the mass/Schur space. Do not reuse closed-boundary masks or pressure gauges.

Nektar already provides `NekSparseMatrix`/COO/BCO/BSR storage,
`SparseUtils`, `NekLinSysIterCG`, `NekLinSysIterGMRES` and
MultiRegions iterative-full/static-condensation operators. They are useful
building blocks, not a ready adapter for this patch's combined regional pressure
and sampler couplings. Start with a serial CSR assembly/action adapter. Choose
and verify a sparse SPD mass factor or tightly controlled iterative inverse.
If its tolerance makes the pressure action variable, use a qualified flexible
outer Krylov solve; do not assume the installed GMRES is flexible. Record
inner/outer true residuals, iteration counts and action-linearity uncertainty.
An approximate preconditioner may change convergence cost, never the declared
governing M inverse or stored one-sweep equations.

## Required equivalence and admission sequence

1. Isolated dependency/backend availability audit; explicit memory accounting
   for every reachable allocation. No global quadratic allocation under the
   scalable policy. Keep the 3,000 guard for the old dense policy.
2. Small CT7 d0 and CT8 displaced A95 fixtures: compare all native volume,
   interface/top/outlet loads, free matrices/actions, true linear residuals,
   consistent-mass corrections and absolute pressure level against dense.
   Use all-row scaled2e-10 and state/history differences bounded by an explicitly
   declared solve-error budget. Include nonzero roughness pressure action and
   `p=x^2`, top epsilon flux sign, wrong orientation/metric, stale masks and
   wrong model/policy/boundary negative controls.
3. Same-method restart and rejected-step rollback; cacheOff reference and
   immutable-row reuse comparison. No reuse of state-dependent admissibility.
4. Three increasing *flat* meshes, ending at the actual CT8 target domain:
   record native unknowns/nonzeros, peak RSS, CPU/elapsed, assembly, mass,
   predictor and pressure costs, conditioning/iterations. Only admit the
   scalable route beyond3,000 after equivalence and measured resource gates.
5. Reconstruct initial full-height reference error and unforced spatial loads.
   Preserve the continuous CT8 physical fit, displacement, surface parameter,
   Richards--Norris top and absolute outlet traction. No residual restoration.
6. First close the retained small-strip absolute-pressure uncertainty with
   an independently specified diagnostic/qualification contract; retain
   the CT8 failed pair. Predeclare independent timestep histories and duration in target-domain
   transit times. Advance `.01,.1,1` transit stages only on passing field,
   traction, domain, mass, complete scalar-budget and temporal-uncertainty gates.
   The small strip is an equivalence case, not full-domain preservation.

Do not include terrain evolution in CT9. The collar-to-horizontal-roof map still
needs coupled Cartesian gradient/cofactor/transmission qualification. CT6's
old U/w temporal gate, smooth preparation, roughness-sublayer/interface adequacy,
thin-layer accuracy and full measured approach-layer fidelity stay open.

## Recommended CT9 prompt

Implement a separately versioned serial scalable backend for the exact CT8
Natural A95 displaced Accepted/Richards--Norris one-sweep route. Preserve the
native fields, physical functions, consistent free mass, pressure-Schur action,
absolute traction, stage order and all physical checks. Audit and remove every
reachable global dense allocation listed above; keep the old dense backend and
its guard unchanged. Verify installed dependencies before selecting a sparse
solver. Demonstrate dense/scalable operator, pressure-level, trajectory,
restart/rollback and complete-budget equivalence on the sealed small cases.
Then perform measured increasing-mesh preflight and, only when independently
admitted, including closure of CT8's small-strip pressure uncertainty, the
staged full-target-domain flat neutral-preservation test with
fresh paired histories and declared flow-through horizon. Preserve failed gates;
do not begin curved production, retune the reference, or claim hill validation.
