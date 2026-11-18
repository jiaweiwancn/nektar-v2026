# F3 matched steady inner-k fixture

The first coupled row D run uses the exact source19 CP2 base 28×100 Q3/Q2
mesh and the retained v2 physical boundary table. Clone the base session.
Change only the steady wall-clock cap, output paths and the fresh steady
initial-guess declaration. Keep the original historical inlet, top and
ground, prescribed outlet traction and all scalar boundary loads. Keep
native q8×12, default linear tolerance 1e-12, steady tolerance 1e-10 and
the FullGradientV1 steady inverse. Source19 coefficients are only a
starting iterate at step and time zero. Use the F2-checked isolated ELF
and check each accepted state's reference, boundary and source identity.

Require independently assembled order32 current momentum and continuity
at 1e-9 and k and outer epsilon at 1e-8. Recompute the dependent epsilon
lift and report the one-sided representation defects separately. Retain
the current independent cubic-root and covariance evaluation, the native
positivity and supported planar principal guard, and all independent
accounting checks. The physical boundary budgets retain both essential
reactions and separately calculated owner-gradient fluxes. No prescribed
interface Q is added. Maintain the unchanged stress and roughness forcing.

Report the physical k departures normalized by constant Kin, with
inner RMS at most .05 and sampled MAX at most .10, over both the bulk
and every frozen station. Keep every prior observation rule, exact end
points on both owners, the inlet/H corner, and the cutoff, ground and
collar. Use the source19 base independent audit and a common-rectangle
physical difference using order20. Verify matching geometry, retained
essential velocity and k lifts, and the input table. Report U, W, k,
epsilon, stress, pressure, gradient and traction changes separately from
candidate discretization error. Use fresh signed production-destruction
and advective/diffusive owner flux integrals by original x column and
material panel to localize any failure.

Initial serial tranche caps are 35 seconds for inputs, 860 for the native
solve including startup, 235 for current analysis with localization and the
complete station inventory, and 180 for the same-mesh physical comparison.
Reserve 60 for publication. The CP2 2,800-cell resource data and F2 default
parity provide a 4.5 GB peak tree RSS forecast, with the existing 8 GiB
science stop. Admission reads the live allowance before launch. Each
operation uses one science process and all six numerical thread limits at
one. Make the outcome and spent allowance prerequisites for the next
tranche. Compare a promising state by targeted directional refinement,
with separately resolved quadrature and algebraic error. Use the CP2
common-owner numerical limits for all coupled quantities, with inner-k
RMS/MAX at most .001. The numerical-only controls must each contribute at
most .0001 in inner-k RMS/MAX. A promising state passes the independent
equations and substantially reduces inner k departure on the matched mesh.
Use its signed transport and owner error indicators to freeze the
directional mesh before running it.

The independently split postprocessing decomposes the candidate's actual
P−E into the historical parallel auxiliary budget, the change in
auxiliary shear production from the current velocity field, current
isotropic work, and the full-WJ production increment. Integrate at orders
20 and 32 by original owner and below, inside and above the blend.
Use independently integrated common-owner polynomial k differences at
Gauss8, which integrates the squared cellwise Q3 difference exactly.
Refine the MAX probes from Gauss20 to 40 with exact endpoints on both
owners. Keep the existing full common-owner comparison and the separate
quadrature and tight-solve controls. The fixed fields and boundary
recipes are unchanged by these checks.
