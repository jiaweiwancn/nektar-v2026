# WJ-H13: full-height neutral inlet and paired terrain-domain preparation

Frozen before the offline scientific batch, 2026-11-13. Starting and fetched
remote HEAD is `36115115e7d4e0ce2e98d9ed83c4c28bc2e9ec89`; the WJ worktree is
clean. No intervening commit, H13 source/artifact, existing H13 Linux root/job
or subsequent selected domain manifest was found. The unrelated main checkout
and its edits/untracked literature remain untouched.

This records the user's explicit change of direction: the **full-height neutral
WJ reference** replaces every earlier H13 finite-depth/Richards–Norris completion
design. Richards–Norris now supports the selected reference's top shear/flux,
not a finite-depth supplement. H12 audit PASS and its closed129.168471056/600 s
allocation remain immutable; H11 temporal qualification stays UNRESOLVED.

## Frozen case and artifacts

Use WJNeutralBudgetReferenceV1 with H_h=40 mm, U_infinity=1, u_tau=.051,
nu=.051*.005/3.97 (from reported .2 mm, not fitted z0), d0=.0875,
top=20, H_w=.15, and the exact stored primary H5 fit/intercept/roughness bytes.
The 6 mm primary has h_s=4.341801417014184 mm. Preserve the5/5.5 mm records in
metadata only. No fit, model, wall origin or matching-height search occurs.
Retain original full3D WJ, DevStrainRawWorkV1, raw full-stress production,
independent Cd=572/6561, qe=81/(2sqrt143), Ne=81/20 and the original constants.
Molecular viscosity remains in all AWF operations and absent only from the
declared high-Re outer equations.

Outer U=(u_tau/kappa)log((z-d0)/z0_e), W=0, k=K=qe*u_tau^2,
epsilon=A/(z-d0), A=u_tau^3/kappa, over the entire outer height. Reuse the
finite-nu inner U and nonconstant budget k, current-K cutoff and dependent
epsilon. The full-height logarithm is the selected idealized inlet, including
U(800 mm)/U_infinity about1.09003; no taper, rescaling or enforced unit top speed.

Shared chart x knots are
[-10,-7.5,-5,-2.5,-1.25,0,1.25,2.5,4,5.5,6.5,7.5,10,12.5,15].
Reuse `wj_h5_cases.mesh_knots(primary6mm,bisect_once=True)` without rounded
reconstruction:16 inner+34 outer cells, unique H, original breakpoints and
observation faces. Retain Q3 differential/Q1 regional pressure,8-point native
GLL representation/lifts, order8 canonical mass/projection and order16 physical
mass/projection; operator tangential12/normal16 are separate settings.

Generate one shared topology/space manifest, paired chart/physical geometry
artifacts and a guarded **NumSteps=0 NOT_RUN flat session**. Flat selects the
horizontal Plane/CompositeTunnelV1 and FlatInletTraceV1; hill selects amplitude1
CosineHill with half-width2.5, the normally displaced effective wall and the
certified composite outer continuation to physical top20. Only the geometry
specialization differs; do not fabricate a runnable curved OpenTunnel session.
Preserve element/edge IDs, boundary roles, region partitions and interface
owners. The initial domain is [-400,600] mm with physical top800 mm; it is not
a boundary-independence or resolved-hill claim.

## Boundary and geometry checks

Use physical Functions (`y` is physical z and `v` is W), not outer chart rho.
Inlet U/W/k and outer epsilon have ground no-slip and dependent-H precedence;
no inlet pressure Dirichlet datum. Ground k flux is zero; no inner epsilon PDE.
Reference pi_outer=0, pi_inner=-K/Ne is a seed, not an evolving jump constraint.
ReferenceTractionV1 top: U_n=0, t_x=u_tau^2, outward j_k=0,
j_epsilon=Cd*K^2/[sigma_epsilon*(top-d0)]>0. Normal traction is a reaction.
Complete outlet (t_x,t_z): inner(+K/Ne,+u_tau^2), outer(-K/Ne,+u_tau^2).
Zero scalar diffusive outflow and existing backflow stop; no extra gauge.
Distinguish fixed prescribed flux from gradient/Robin/Dirichlet laws.

Check the unspecialized coefficient B=Cd*K^2/sigma_epsilon, its equality to
(Ce2-Ce1)*A^2/K for the selected coefficients, derivative B/r^2 and signed
integrated epsilon flux. Reuse H5 stationary-wall identities, independently
sample common traces, total finite-nu shear, Q endpoints and pressure signs.
Expression accuracy remains2e-11. Do not demand zero projected spatial residual.

Geometry checks use H1's whole-height analytic normal-reach/composite-Jacobian,
injectivity and nonoverlap certificate, plus every element's tensor5 GLL
nodes, corners/faces and interface cofactors. This is25 points per700 elements
per partner (35,000 owned map samples), with both one-sided foot curvatures.
Shared faces coincide; H cofactors oppose; upstream/downstream boundaries are
flat and the physical top is exact. No whole-height normal rays or foot smoothing.
Retain numerical metric differences between partners, not a false equality.

## Capacity and future preservation protocol

Source-aware lightweight topology enumeration must confirm state formula
3*(3Nx+1)*(3(Ni+No)+1)+(3Nx+1)*(3No+1)+(Nx+1)*(Ni+No+2),
or1720*Nx+608 here:2328/4048/24688 at Nx1/2/14. Count individual fields,
essential/dependent trace DOFs, free velocity/pressure/scalar blocks and dense
storage lower bounds using integers only. Never allocate a global dense matrix.
Canonical topology IDs are not claimed to be native-exported numbering.

The actual reviewed predicates retain default3000 and exclude this mapped
CompositeTunnel route from the4331 opt-in. Native construction already creates
dense mass/projection structures before evolution admission. The external
preflight checks manifest/XML/topology/source capacity before any executable;
NumSteps=0 is not admission. H13 blocks every solver/probe/build/restart launch.
No guard increase, mesh coarsening, stretched single-column substitute or map
identity change is allowed. Produce a concrete mapped block/operator prerequisite
with source call graph and dense-reference equivalence tests, not a backend rewrite.

Machine-readable future stages: A fresh target CP0 and <=4-step admission pilot
(future independent authorization/cost admission); B domain-duration preservation;
C separately justified sensitivity. No H11 checkpoint conversion. Actual q_h^0
after projection/lifts/initializer differs from analytical q_0 and possibly the
frozen internal reference. Separate representation, drift, reference and fetch
errors. Use physical regional Gauss8 norms, owner-preserving finite maxima,
outer-only epsilon and height-relative epsilon errors; preserve traction-fixed
pressure. Freeze stations including inlet/outlet, feet, crest and wake through7.5.
Keep Fig.3 at its true flat virtual crest and retain physical7 mm observations
inside AWF; hill observations require the certified inverse and correct height origin.

Compute U_b from the complete coupled physical cross-section and T_b=Lx/U_b,
T(z)=Lx/U_0(z), including slowest outer transit Lx/U_0(H). Record conservative
inner diffusion and representative k/epsilon clocks. Proposed Stage B horizon
rule is max(3*T_b,2*T_slow_outer,T_inner_diffusion,T_outer_top_k_over_epsilon),
with terminal monitoring over one slowest outer transit. Its funding, timestep,
application budgets and comparison floors must be frozen in a later native
contract before launch; insufficient duration stays PARTIAL_DURATION/UNRESOLVED.
This finite horizon does not establish stationarity. No tolerance is fitted here.

## Provenance, resources and completion

The original Mac main/appendix hashes are respectively
1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124 and
0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c.
Preserve the quad/\\quad erratum and `(d+a*h)/a^2`. Both local original PDFs
are present: Cao–Tamura SHAadffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df;
Richards–Norris SHA2a8d769dde2d85a6034f2cc64a5da70da31feb9dfdf2c236484cf38520541780.
Richards–Norris printed259/261 (PDF3/5), Eq23/48 and Section3 were visually read.
Do not import that paper's different model constants or optional buffer sources.

New root `/home/jay/build/wj-h13-fullheight-paired-domain-20261113`, healthy
ext4 UUID4c5b928b-0ef9-464e-91cb-3b0b2202e3b0. SSH alias timed out; documented
Tailscale fallback works. New independent600 s allocation:30 fully charged setup,
<=480 launched operation seconds including<=360 science,90 fully reserved/charged
closure/publication/transport. Reserve every operation; retain failed/unknown
charges and monotone halts. One serial worker, six thread controls1 before
imports, host<=72 workers, process-tree RSS<=2 GiB, new output<=256 MiB.

ZERO native governing attempts/CP0s/solver/probe/strict invocations/builds/restart
conversions. Reuse unchanged historical math/native evidence, not executables.
Run only H13-focused profile/flux/geometry/topology/count/initializer/role/guard
fixtures on Linux. Emit artifacts, source hashes, actual commands, tests, numeric
checks, forecasts and closed accounting; verify compact transport and normal-push
coherent commits. Statuses stay separate: preparation/continuum compatibility,
capacity OBSTRUCTED, fresh native initialization/evolution/duration preservation
and experimental comparison NOT_RUN. Artifact success is not preservation PASS.

The immediate next package is mapped-route capacity/backend qualification and
fresh target admission. Domain-duration flat preservation follows; curved physical
open boundaries/initialization and hill evolution follow their own qualification.
Finite-depth/ambient profiles, EP, old measured-rms/Blackadar/20 mm inputs,
retuning, matching-height searches, anisotropic inner stress and old ratio audits
remain outside H13. The all7-mm-outside interval remains empty and no physical
matching height is declared optimal.
