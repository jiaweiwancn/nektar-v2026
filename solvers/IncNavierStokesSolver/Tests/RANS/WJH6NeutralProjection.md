# H6: one prescribed normal bisection, no governing steps

The continuum family remains WJNeutralBudgetReferenceV1. Use
`wj_h5_cases.case_xml(folder,c,bisect_once=True)` for the sole H6 mesh;
the default False retains the H5 generator. Every old normal interval is
bisected in both regions, with all endpoints and the H interface preserved.
This gives16 inner/34 outer cells and2328 native states on Nx1/Q3/Q1.

`study_wj_h6.py --root R --mode generate` writes exactly the three fixed-height
NumSteps=0 inputs and verifies their Functions against the imported H5 XML.
There is no timestep/height search or restart option. `run_wj_h6.py` enforces
the new1200s allocation and three unique-height CP0 reservations with no retry.
The H5 resource journal and its three original startup decisions are unchanged.

`wj_h6_projection.py` assembles the physical order16 mass/load equations from
the exported CG maps/ids/signs, solves them offline, reconstructs the native
GLL8 endpoint-constrained exterior transforms and then the actual-K dependent
epsilon trace. These are reconstructed stages, not native intermediate exports.
Two integral-defect refinements reduce hierarchical-basis cancellation while
leaving the mass equations and quadrature unchanged. No governing state is
modified. Rejected H5 candidates are read from genuine report/event coefficients,
never promoted to accepted checkpoints or passed to the strict checkpoint reader.

The `replay` mode qualifies final coefficients/fields/gradients first. The
`observe` mode reuses the qualified native reader and corrected projected-test
weak oracle; it performs no Project/Lift or governing solve. `compare` uses a
frozen physical union with all legitimate owners. Native admissibility and
the additional half-reference-margin/error-reduction gates remain separate.

Registered focused tests are WJH5NeutralInlet (nine methods) and
WJH6NeutralProjection (six methods). Scientific results, command identities,
the two different reading intervals and all limits are in the H6 contract and
evidence under docs/plans. H6 does not change the full-height control into a
finite-depth inlet or resolve the7mm observation/height conflict.
