# H27 CP2 selected-source reference refinement

The fixed-field audit on source19 at 28×200 cells reproduced the inner-k RMS
departure of .1327946821 in units of Kin. Both independent rules agree to
1.2e-16 on that norm. The largest sampled deficit is .4418922347 at the
ground and x=-8.75. About 94.6% of the integrated squared departure lies
outside the first two columns. Most downstream column RMS values approach
.14. The column flux defects oscillate downstream, reaching .1821 when
normalized by A times column width near x=7.5. Refine both locations.
[Localization protocol](wj-h27-cp2-localization-admission-v1.md)

All steady cases retain Q3/Q2, the original flat domain, all prescribed and
natural boundaries, the retained inlet and the approved v2 law. Use source19
and its existing binary, table and cutoff panels. Reuse the two archived
source11 fields solely as explicit initial guesses for the steady iteration.
Fresh current equations, dependent traces and coefficients are required.
Solve at zero physical time, with FullGradientV1, Krylov dimension 80,
nonlinear tolerance 1e-10 and linear tolerance 1e-12 unless specified below.

| Case | Cells | Native x/r quadrature | Isolated input |
| --- | ---: | --- | --- |
| `base` | 28×100 | 8/12 | Current-source common h baseline |
| `quadrature` | 28×100 | 12/16 | Native quadrature relative to `base` |
| `tight` | 28×100 | 8/12 | Tolerances 3e-12 and 3e-13 relative to `base` |
| `streamwise` | 56×100 | 8/12 | All original x intervals bisected |
| `radial` | 28×200 | 8/12 | Retained exact source19 radial candidate and audits |
| `targeted` | 44×132 | 8/12 | Refine the inner layer and selected x intervals |

For `targeted`, retain every base edge and every `radial` edge below H.
Insert x=-9.921875, -9.6875 and -9.375. Bisect original x interval indices
1–3, 15–23 and 27, zero-based. This keeps the inlet corners, covers the site
of the early deficit and resolves the changing width and flux defect near
x=7.5. Above H retain the base radial edges. None of the three refined meshes
is uniformly finer than the others.

Independently audit each new converged field at order32 using the original
current equations, weak accounts, owner interfaces, stress/admissibility
evidence and frozen physical stations and norms. Compare `base` with both
controls and with `streamwise` and `radial`, and compare `targeted` with
`streamwise` and `radial`. Use the prior exact source19 radial audits at
order32 and order48. The common-cell comparisons use the unchanged oracle
at order20, retaining both actual one-sided gradients and all frozen
samples. Record budget changes on their fixed physical scales as well as
weak reaction identities. Any proposed positive classification needs
confirmation of the decisive comparison at order32 and all existing
resolution and admissibility gates. Report C if these remain unresolved.

The current 5600-cell radial case took 103.73 s and 5.470 GB peak tree
RSS. The older streamwise solve used 6.622 GB on 5600 cells. Forecast
6–8 GB for the new targeted case given 5808 cells and its intermediate
aspect ratio. The binding limit remains 8 GiB RSS. Run serially and let
the existing watchdog stop the case at the cap. No full joint doubling is
admitted. Use the live account, storage and worker census before admission.

| Operation | Aggregate seconds cap |
| --- | ---: |
| Fixtures | 40 |
| `base`, `quadrature`, `tight` native | 450 + 520 + 420 |
| `streamwise`, `targeted` native | 1100 + 1300 |
| Independent equations and boundary audits | 3 × 210 + 2 × 355 |
| Six order20 common-owner comparisons | 2 × 130 + 2 × 180 + 2 × 270 |
| Compact assessment and publication | 80 + 90 |
| Total | 6500 |

At the post-localization reading the same account was OPEN at
17166.227363280836 of 43200 seconds. The full cap leaves 15933.772636719164
non-protected seconds with 3600 for evidence. Keep four build workers only
if later needed, single rank and all six numerical-library thread limits
one, actual host workers <=72, science/build RSS 8/16 GiB, storage 8 GiB
each for science and source/build, and free space >=15 GiB. This tranche
has no build. Forecast less than .9 GiB of additional science storage.
Preserve all failures, skip their dependent comparisons and continue the
independent cases within the cap. An order32 comparison confirmation or
any further refinement requires an explicit bounded follow-up tranche on
the same account. No long physical pilot or Task 3 horizon is admitted.
