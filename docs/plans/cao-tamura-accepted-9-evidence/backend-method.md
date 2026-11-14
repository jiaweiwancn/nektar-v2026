# AcceptedSerialSparseV1 implementation and allocation audit

The governing implementation is pinned by the private build receipt. Its model,
reference, pressure/traction convention and one-sweep order are the CT8 central
displaced Natural A95 definitions. The backend identity is separate from the
immutable model descriptor. Checkpoints and preservation/replay signatures bind
it; cross-backend histories are rejected rather than relabelled.

The dependency audit finds Eigen sparse headers and the installed Nektar
matrix-action/GMRES interfaces. The installed SDK remains MPI ON, PETSc OFF.
`SparseMatrix.*`, `SparseUtils.*`, `NekLinSysIter*` and the MultiRegions sparse
adapters were inspected. A narrow sorted-row/CSR adapter handles the combined
regional pressure fields, native signed CG restrictions and roughness sampler
links directly. It does not route the coupled patch through a single-field
assembly map or convert a global dense matrix to sparse storage. Eigen supplies
private header-only SparseLU and sparse LDLT; installed Nektar supplies GMRES.
Header hashes, compiler flags, binaries and loaded libraries are in the sealed
runtime/experiment receipts. No PETSc build or shared installation is changed.

| Reachable phase | Implementation and bound |
|---|---|
| Native maps and rows | Existing signed field maps and generating-owner/side semantics; linear index vectors. |
| Startup mass | Direct sparse accumulation, retaining the dense path's long-double per-entry mass summation; immutable CSR. No global long-double square array. |
| General projection | Local native FwdTrans seed, then the existing three consistent-mass defect corrections with sparse positive LDLT. |
| Radial/face projection | Existing local systems; radial KKT dimension capped at 2048. Degree is uniformly P3/P1 or P4/P2, so face transforms stay small. |
| Frozen flow and scalars | Original element/interface/natural-face contribution order, direct sparse row insertion; epsilon dependent rows replaced explicitly. |
| Predictor and scalars | SparseLU; original frozen matrices and right-hand sides; true residual recomputed. |
| Consistent correction mass | Exact free-ID restriction of physical consistent mass; positive LDLT. No diagonal mass in a governing update. |
| Pressure | Nektar flexible right-preconditioned GMRES, restart 80, maximum 1000 iterations; vector actions only. No dense Schur or array of mass-inverse-G columns. |
| Preconditioner | Sparse LU of `D diag(M_free)^-1 G / alpha - C`. Its diagonal mass is confined to preconditioning. |
| Required norms/audits | Sparse restricted consistent masses and original full constitutive/domain inventories. State/gradient validity is never inferred from a k/epsilon-only cache. |
| Optional exports | Dense matrix exports and per-element global-vector localization limited to 1000-slot fixtures; aggregate vectors, CSR entries/actions and hashes above it. |
| Native snapshots | Sparse flow, mass, G/D/C and k matrices plus actual stage/history vectors. Readers validate encoding, dimensions and hashes. |
| Retained histories/records | IMEX1 and at most 64 steps per segment; at most 65 accepted state vectors and bounded attempted-stage records. Longer durations require genuine same-method restart segments. |
| Legacy nonlinear/action/bound routes | Explicitly unsupported. The old dense backend and its 3000-slot limit remain unchanged. |

The sparse route additionally caps stored slots at 70000, each CSR matrix at
20 million nonzeros. The original 30-million factor cap is retained in the
failed capacity evidence. The prospective measured-fill amendment permits
40 million reported entries, after the unchanged P4 momentum block measured
38165601. Direct/factor
storage and local dense work are accompanied by a 16 GiB process address-space
ceiling in the qualification drivers. LDLT fill is the library's stored lower
factor count; its diagonal and permutations add O(N) vectors. Actual peak RSS
and factor/nonzero counts are measured rather than inferred from these caps.

For each pressure vector the governing operation is
`D[(alpha M_free)^-1(Gp)] - Cp`. Direct consistent-mass factor application makes
this a fixed action apart from roundoff. Every direct solve recomputes its
unpreconditioned residual. Pressure RHS normalization avoids the installed
GMRES small-RHS magnitude floor; the returned **unscaled** action residual
must satisfy `max|r| <= 2e-12 max|rhs| + 1e-25`. The internal Krylov tolerance is
stricter. Frozen packets also verify the complete independently assembled
native G action, consistent-mass correction residual and true continuity.

Small-case complete-matrix/action comparisons use the independent dense oracle
and the 2e-10 row limit. Fresh physical trajectories use the stricter 1e-8
scaled field/traction equivalence budget. This qualifies numerical equivalence,
not physical pressure uncertainty. The CT8 failed classification must remain
failed under the new backend. Final gate outcomes and actual capacity metrics
are supplied by the compact qualified evidence, not by this implementation note.

The final frozen C++ revision is `1e73bf53353490e7117f36799250ef79208e8df6`.
The qualified runtime is `build-v8/runtime`; its complete binary/header/loader
identities are in `qualified-v2/runtime.json` and the loaded-file catalog.
Direct residual refinement reuses a factor at most three times under the same
true-residual bound. Its actual use is recorded separately from governing
equation stages. The intermediate preconditioner needs one correction; the
other final capacity cases need none. The original failures are not replaced.

The ordinary solver build defaults this feature OFF. The production CMake
switch is `NEKTAR_RANS_SERIAL_SPARSE`; the isolated test/runtime build uses
`RANS_ENABLE_SERIAL_SPARSE`. Both require Eigen3 and preserve the installed SDK.
The private build driver records the exact tested configuration. A qualified
case explicitly selects:

```xml
<I PROPERTY="RANSWallLinearBackend" VALUE="AcceptedSerialSparseV1" />
<!-- In PARAMETERS: -->
<P> FlexibleGMRES = 1 </P>
```

All other admitted provider, wall, boundary, method and basis guards still
apply. Old dense or earlier sparse-component checkpoints are not portable to
this identity; create fresh histories. Within a component, genuine unchanged
method/timestep restart is verified. Optional dense exports do not remove any
mandatory physical-state or constitutive check.

This is a serial flat prototype, not a distributed backend or an inf-sup,
positivity, realizability, long-time or terrain theorem. The intended terrain
still needs mapped region-aware numbering, joined-map Cartesian gradients,
cofactors/transmission and the outstanding temporal/spatial qualification.
