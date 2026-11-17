# H27 native coupled EquationSystem and outer epsilon ownership

Read-only engineering audit of the current branch, after the user-approved
isotropic-lower-wall/full-WJ-collar amendment. No solver, CMake or model files
were changed by this audit. No numerical jobs were run. This note recommends
the smallest genuine normal executable integration, not a renamed probe.

## Recommended registration and class boundary

Introduce `EquationSystems/WJNativeCoupled.{h,cpp}`, deriving **directly from
`SolverUtils::EquationSystem`**. Give it a distinct factory name such as
`WJNativeCoupledV2`; keep all old factory names and defaults intact. The class
owns the current mixed flow kernel, full-domain k space, outer-only epsilon
space, current constitutive/reference state, coupled nonlinear iteration,
accepted BDF2/EX2 histories, restart and output. It does not inherit stock VCS
or the historical RANS evolution algorithms.

Evidence:

- `IncNavierStokesSolver.cpp:54-69` already creates SessionReader, MeshGraph and
  the selected Driver, then calls Execute. No new main is necessary.
- `library/SolverUtils/Driver.cpp:86-121` chooses the equation factory from
  SolverType (falling back to EqType), and uses the registered creator.
- `library/SolverUtils/DriverStandard.cpp:75-84` calls PrintSummary,
  DoInitialise, DoSolve and Output. Its lines99-107 subsequently request exact
  solutions and errors for every variable.
- `RANSEquilibriumWallSystem.h:14-33` and `.cpp:14-17,41-47` demonstrate the
  direct EquationSystem pattern, including `EquationSystem::v_InitObject(false)`.
- `EquationSystem.cpp:294-304,709-716` still sizes `m_fields` and loads Time,
  TimeStep and NumSteps when automatic field creation is disabled.

Implement the following overrides:

| Override | Required behavior |
|---|---|
| `v_InitObject(bool)` | Call base with false; require serial physical 2D flat conforming quads and supported orders; load immutable model/BC/integrator settings; create independent u,v,p,k fields and outer epsilon topology; build current MixedFlow/scalar actions and owner maps. Reject unsupported geometry/partitioning explicitly. |
| `v_DoInitialise(bool)` | Read fresh functions or the complete compatible coupled checkpoint; establish all scalar/flow histories and current dependent E/cutoff/trace; apply the correct startup policy; store native coefficient/physical views. A checkpoint must not be disguised as InitialConditions. |
| `v_DoSolve()` | Advance the actual U, pressure, k and outer epsilon equations with current coefficients and repeated coupled iterations; refresh all dependent data before accepting the full residual; preserve accepted histories on rejection. Update `m_time` only for accepted physical progress. |
| `v_Output()` | Write fields with their own element inventories; epsilon has only outer elements. Export inner dependent E separately as a diagnostic, not an epsilon unknown. Own the complete scientific/recovery output policy. |
| `v_EvaluateExactSolution(...)` | Resize output to `m_fields[field]->GetTotPoints()` and evaluate that field's own geometry; this is essential for outer epsilon because DriverStandard initially allocates the first field's point count. |
| `v_GenerateSummary(...)` | Report the explicit factory/model version, molecular viscosity, scalar convention, spaces, physical boundaries, coupled integration and qualified serial scope. |

Initially require `NumQuadPointsError=0` or implement the error overrides for
the restricted epsilon inventory. Generic extra-point routines can reconstruct
a full-domain expansion. The ordinary L2/Linf routines work with correctly
sized exact arrays on each field (`EquationSystem.cpp:800-849,876-906`).

Factory setup follows existing constructors:

```cpp
static EquationSystemSharedPtr create(const SessionReaderSharedPtr &s,
                                     const MeshGraphSharedPtr &g)
{
    auto eq = MemoryManager<WJNativeCoupled>::AllocateSharedPtr(s, g);
    eq->InitObject();
    return eq;
}
```

Register both `GetEquationSystemFactory().RegisterCreatorFunction(...)` and
`SessionReader::RegisterEnumValue("SolverType",...)`. Append a new SolverType
enumerator in `EquationSystems/IncNavierStokes.h:54-68`; do not renumber earlier
values. Registration is mandatory: `SessionReader.cpp:2609-2624` rejects an
unregistered SolverType even when the equation factory knows its string.
Use existing `EqType=UnsteadyNavierStokes`, `EvolutionOperator=Nonlinear`,
`Driver=Standard`, `Projection=Continuous`, and the new SolverType. Do not
register a new default SolverType or alter the default Driver.

The session variable order can remain `u,v,k,epsilon,p`. Construct MixedFlow
with an explicit vector `{u,v,p}`; do not pass `m_fields` wholesale. Store the
five native field views in the corresponding session order. Read-only native
output views are not the accepted history vectors themselves.

## Outer-only epsilon: smallest native construction

The existing code already supplies the necessary topology operation:

```cpp
std::vector<unsigned int> outerLocalExpansionIds = /* validated selection */;
auto subset = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(
    *kField, outerLocalExpansionIds, true);
auto epsilonField = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
    session, *subset);
```

Use the full-domain Q3 k field as the prototype when epsilon also uses Q3.
This avoids creating or advancing any inner epsilon unknown. A plain ExpList
prototype for a separate epsilon basis is also possible if a later qualified
recipe needs different orders. The subset argument is **local expansion
indices**, not geometry IDs.

- `MultiRegions/ExpList.cpp:131-149` copies only the selected local expansion
  pointers, preserving their geometry objects/global IDs, and allocates fresh
  subset coefficient and physical arrays.
- `RANSWallFieldLayout.cpp:102-107` uses precisely this construction.
- `MultiRegions/ContField.cpp:212-220` builds its native AssemblyMapCG from the
  subset. This is a topology-only constructor: despite its old comment, it
  does **not** generate XML boundary expansions or a native Dirichlet mask.
- `MultiRegions/DisContField.cpp:795-797` only copies the ExpList in this path;
  `AssemblyMapCG.h:77-89` explicitly supplies empty boundary data by default.

Therefore the new scalar kernel must own an explicit epsilon essential-DOF
set/lift and sparse free-DOF map. Do not assume the first N global coefficients
are epsilon Dirichlet values. Do not call generic native Helmholtz/FwdTrans
boundary enforcement on this topology and expect H to have appeared as a
Dirichlet boundary. Use native local basis/derivative/quadrature and signed CG
assembly, with the same explicit sparse elimination/lift pattern as MixedFlow.

This recommendation reuses the native **construction**, not the historical
`RANSWallFieldLayout` mixed system. That layout splits pressure and creates an
extra gauge (`RANSWallFieldLayout.cpp:327-335`), while H27 has the current
global Q2 pressure and boundary-dependent gauge. Its `RANSWallNativeField::Mass`
and `Project` use dense matrices (`:199-205,243-253`), unsuitable for the
production scalar integration. Its field wrapper has useful local sampling
and trace examples, but should not become the production scalar linear solver.

The topology-only approach needs no second SessionReader, second full mesh,
temporary native XML or duplicated initial-condition interpretation. Read
external epsilon inlet/top/outlet data from the original session on selected
outer edges. There is no physical epsilon boundary on the numerical ground
because the epsilon field does not reach it. Do not fabricate one simply to
satisfy an unused full-domain epsilon ContField constructor.

## Ownership, edges and dependent Dirichlet trace

During construction, enumerate geometry IDs, local expansion indices, physical
bounds, edge geometry IDs and signed local-to-global maps once per field.
Require that H is an existing complete mesh edge layer: reject any element
straddling `physical_z=d0+H`. Select outer elements by physical model distance
`r=physical_z-d0`, never by treating physical elevation as r. Keep the separate
full-field and epsilon local indices in a geometry-ID join.

The subset shares the original expansion objects, so its outer geometry IDs,
vertex IDs, edge IDs, basis and metric data are retained exactly. Global CG
coefficient numbering can change completely. Never equate full k and epsilon
global IDs, or assume their element-array positions match. The old layout's
checks at `RANSWallFieldLayout.cpp:319-325` use unshifted r coordinates; copying
them directly into the physical mesh with d0=.0875 would classify incorrectly.

Construct an edge adjacency table from `(geometry ID, local edge)` and edge
global ID. It must identify:

- every full-domain internal edge and its two owners;
- each H edge with one inner and one outer owner;
- each outer epsilon H edge with exactly one outer epsilon owner;
- the actual inlet, top and outlet exterior edges.

Use edge IDs and orientation, not coordinate-only hashing. At coincident points,
gradient evaluations carry their element owner. A single physical coordinate
can have two valid C0 gradients; it must not merge those rows or erase one-sided
traction/scalar flux data. K at H uses the common k trace; the prescribed outer
epsilon trace is assembled from the outer owner.

For each outer H edge and each current nonlinear iterate:

1. Sample current Q3 k on that exact edge through its full-k signed map.
2. Recompute the current wall/reference state and actual dependent E(H). For
   the unamended healthy source this is A(K)/H; use the declared amended E if
   a positivity policy actually changes it.
3. Build the epsilon edge lift with the declared trace projection. For Q3
   Modified_A, a boundary-constrained segment transform retains exact endpoint
   values and projects its interior edge modes.
4. Scatter edge coefficients through trace-to-element orientation signs and
   the epsilon CG map. Check duplicate corner assignments for consistency.
5. Rebuild the scalar RHS lift, current stresses, diffusion and all dependent
   residuals after this trace changes. Do not add an independent epsilon flux
   constraint at H.

`RANSWallNativeField::Boundary`, `RANSWallFieldLayout.cpp:257-290`, demonstrates
boundary-constrained segment transforms and corner checks. Its hard-coded
Modified_A mode indices assume canonical axis-aligned quads. Either retain
and validate those exact restrictions or use `GetTraceToElementMap` and
`GetTraceOrient` (`StdRegions/StdExpansion.h:775-781`,
`LocalRegions/Expansion.h:181-183`). Do not merely copy edge orientation from
the full-k global map: epsilon has its own map.

A(K(x))/H is generally not a Q3 polynomial. Define and test the projection as
part of the discrete trace contract; distinguish coefficient equality from
pointwise nonlinear trace error. Evaluate the latter independently under
quadrature/refinement. At the inlet/H corner, use the common exact endpoint
K_in and verify epsilon_in(H)=E(H;K_in); no boundary assignment may silently
override the other. Boundary-constrained projection is useful here because it
preserves endpoint values even when interior edge projection is approximate.

For k, one full-domain native CG space already enforces the common trace. Its
weak operator must use common flux accounting; no artificial H Dirichlet k
condition is added. For epsilon, H is a genuine essential boundary of its
restricted PDE. No epsilon mass, advection, source, diffusion row or history is
assembled below H.

## Alternative: a separate native epsilon boundary session

Use this only if the scalar implementation actually needs stock native
boundary-condition/Helmholtz machinery. Generate an outer-domain session mesh
with only the outer domain composite and correctly restricted exterior
boundary composites, plus H as a real boundary. Preserve original element,
edge and vertex IDs; do not renumber them. Build
`ContField(epsilonSession,epsilonGraph,"epsilon",true,false)`, rather than the
topology-only constructor. Validate the graph's geometry/metrics/basis against
the full mesh by original geometry IDs.

The H boundary has a declared Dirichlet type so the native map marks those
DOFs. Its changing coefficients must be overwritten from current E(H), or
lifted explicitly after native data evaluation. Do not let a later
EvaluateBoundaryConditions call restore an XML placeholder. The dependent
data still belong to the coupled solver, not a fixed session expression.
This route requires another geometry/session identity and exposes extra
boundary filtering/ordering risks. The subset-plus-explicit-lift route is
smaller for the existing Eigen sparse implementation.

## Normal build/run integration and pinned object reuse

Production source files to register in
`solvers/IncNavierStokesSolver/CMakeLists.txt` are:

- existing `EquationSystems/WJNativeFlow.cpp`;
- new `EquationSystems/WJNativeMixedFlow.cpp` and `WJNativeWallV2.cpp`;
- the new scalar/reference/coupled implementation source files and
  `EquationSystems/WJNativeCoupled.cpp`.

Add public headers to `IncNavierStokesSolverHeaders` when solver libraries are
enabled. The current top-level solver list does not contain these native
files. Its production macro creates `IncNavierStokesSolver-obj` and inserts
all its objects into the actual executable (`solvers/CMakeLists.txt:8-24`),
which reliably retains factory static initializers. Ensure Eigen include/link
usage reaches that object target as well as the final executable. The current
Eigen dependency is only declared in the external h26/h27 projects, not the
normal IncNavierStokesSolver CMake file.

A compile-time option for the new native implementation may default off if
preserving existing build dependencies is desired; when on, require Eigen3
explicitly and include the new files. The runtime factory stays opt-in either
way. Do not link a static archive containing only factory registration and
assume the linker extracts it: include the creator translation unit directly
or as an object library. Its methods then reference the MixedFlow library.

There is an existing isolated **normal-main** precedent:
`Tests/RANS/CMakeLists.txt:310-323` reads the production LIBRARY_SOURCES list,
adds the real `IncNavierStokesSolver.cpp`, and builds an isolated executable.
Unlike a bespoke driver, it exercises the standard SessionReader/Driver/factory
lifecycle. If production sources are added through conditional CMake variables
or target_sources, update this extractor: its regex only sees literal
`./path.cpp` entries inside the original LIBRARY_SOURCES block.

The H26 pinned stock-object option is narrower:
`Tests/RANS/h26/CMakeLists.txt:21-32` excludes RANS and the real main and builds
`stock.cpp`. That target is not the requested production entry point. Reuse
its object inventory only with a new target that compiles the **real main**
and new factory/core sources, and retains whatever historical registrations
are claimed to remain supported. Explicitly include the new factory object.
Do not relabel WJH27MixedProbe or WJNativeFlowDriver as IncNavierStokesSolver.

If unchanged stock objects are reused to contain build cost, pin their
compiler flags, ABI, source/header dependency hashes and link identities.
Exclude every object recompiled from a new source to avoid duplicate
definitions/registrations. An appended enum is ABI neutral but still changes
a header; record/review its dependency treatment rather than claiming the old
headers were byte-identical. Prefer a final normal incremental target build
once the integrated path is stable. Keep the historical source capsules and
installed default executable unchanged.

The current `h27/build.py:12` only admits three probe/driver target names; a
real new target requires a corresponding runner/build-manifest extension.
Retain immutable H27 source/runtime identities and the measured startup cost.

The reproducible invocation should ultimately be ordinary:

```sh
IncNavierStokesSolver case.xml
```

where the XML explicitly selects `SolverType=WJNativeCoupledV2`, carries the
physical u/v/k/epsilon data, mesh expansions, timestep and output/restart
contract. Debug request JSONs can remain verification tools, but the
application solver must not require a Python advancing wrapper or a probe
request that bypasses DoInitialise/DoSolve.

## Integration gates that catch the likely mistakes

- Factory smoke through the real main, first with NumSteps=0 and then a
  genuinely coupled accepted step. Summary must identify the new path.
- Enumerate epsilon geometry IDs, coefficient count and histories: no inner
  epsilon IDs/rows; full k and velocity still own every intended element.
- Verify different u/v essential masks, signed v permutation and top tangential
  traction. `WJNativeMixedFlow.cpp:49-95` explicitly handles this; do not
  bypass it by copying u's map for v or by prescribing top u.
- Verify outer H lift coefficients, nonlinear pointwise trace discrepancies,
  inlet/H corner consistency and owner-specific gradient evaluation.
- Verify complete conservative scalar boundary loads, H common k flux,
  epsilon H Dirichlet-only treatment, known-k force from the actual k space,
  and current full coupled residuals after every dependency refresh.
- Run final native output and the DriverStandard error loop with different
  point counts for u/p/epsilon, then cold restart all fields and histories.
- Ensure scientific checkpoints contain global accepted vectors and identities,
  not only whichever native physical view was most recently sampled.

These are prerequisites to the package's coupled numerics, physical pilots and
cost gates. Passing construction/factory smoke does not itself qualify them.
