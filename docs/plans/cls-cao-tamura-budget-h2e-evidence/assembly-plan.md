# Bounded representation/assembly handoff

This is a design assessment using the retained H2d zero-step results, not a new
mesh experiment or a performance measurement. The 2804-slot production pilot
and 3000-slot guard remain unchanged. The rejected single-knot transfer is not
adopted. The 4399–5339-slot allocations remain targets for a subsequent package;
the analytical weak loads near roundoff and improved diagnostic radial errors
do not close the native 1e-3 U/k spatial gate.

## Actual dense paths

- `RANSWallEvolutionPicard.cpp::AssemblePicardBlock` builds a dense selected-row
  system even though its volume and boundary contributions are element-local.
  Regional interface/essential elimination is a graph operation, not an
  inherently dense physical interaction.
- `RANSWallEvolutionPimple.cpp` copies the flow system into A, M, G and D,
  factors A, retains a second mass-key matrix and a mass factor, applies
  M_free^-1 to every pressure column, and forms/factors a dense pressure Schur
  complement. Preserve B=(alpha M_free)^-1 and the full roughness/open-boundary
  pressure action. The momentum inverse is a different operator.
- `RANSWallEvolutionComposite.cpp` retains physical velocity mass, a factor,
  independent continuity Q rows and M^-1 Q^T for initialization. The Schur
  solve is genuinely global, but only has 105 independent rows on this pilot.
- The scalar transport and consistent-k correction paths copy dense regional
  masses. The projection's mass whitening/factor applications are global
  algebra; changing their storage must preserve the same quadratic objective.
- Geometry, source/reconstruction rows and complete audit inventories remain
  material memory/work costs. Sparse factorization will not remove them.
  H2e's roughly 0.8-GiB native setup and 0.9-GiB search peaks include these
  costs; they are not isolated dense-factor measurements.

One binary64 n-by-n array requires 62.90 MB at n=2804, 154.81 MB at n=4399 and
228.04 MB at n=5339 (decimal MB). These are arithmetic storage estimates, not
measured resident memory. Multiple matrices, copies, factors, evaluation rows,
allocator overhead and output buffers coexist. A cubic dense-factor scaling
estimate grows by approximately 3.9–6.9 times; it is not a timing prediction.

## Smallest useful replacement

1. Assemble the unchanged physical element and face contributions directly
   into block CSR graphs. Start with regional/velocity consistent masses and
   the selected flow/scalar matrices. Eliminate essential values using the
   same global coefficient maps and retain signed CG sharing.
2. Use an exact sparse direct mass factor and exact sparse numerical factors
   for the changing momentum/scalar matrices. Keep the small pressure Schur
   matrix dense initially, applying the sparse mass inverse to its existing
   complete G columns. Do not reuse a changing numerical factor or replace
   the mass inverse by a lumped/approximate operator.
3. Reuse the sparse mass factor for native initialization and mass-dual norms.
   Adapt projection mass applications separately, retaining raw centering,
   exact free-row elimination, retained cuts and the same KKT/feasibility
   checks. Do not label a sparse solver behind the current dense assembly a
   storage solution.

For a bounded Q3/P1 graph estimate of 100–200 stored entries per coupled row,
5339 rows need roughly 6.4–12.8 MB of binary64/int32 CSR entries, plus row
pointers. Element/face work arrays and factors are additional. Pressure Schur
storage remains quadratic in pressure rank. Actual sparsity, fill, pivoting,
peak memory and audit costs must be measured on symbolic graphs before any
production-size limit is reconsidered.

## Required equivalence gates

Use the existing affine/composite flat limit and curved fixtures: complete
mass/inverse actions, volume and boundary loads, nonconstant pressure and
roughness coupling, absolute-pressure outlet, matching and epsilon mass
columns. Include actual/counterfactual production, all free and constrained
rows, complete residual/action comparisons, one-correction splitting,
consistent-mass bound objectives, active cuts and mapped sided witnesses.
Then require matched stage vectors/histories/rollback and strict checkpoint
identities. Deterministic summation can preserve bits; any changed assembly
order needs predeclared action/residual tolerances and a measured error bound.

Only after those gates should the retained 4399-slot radial allocation receive
a fresh initialization and zero-step spatial comparison on the original six
streamwise panels. A subsequent 5339-slot confirmation remains conditional.
No evolving knot sweep, new physical model or raised guard is part of H2e.
