# Neutral ABL URANS implementation proposal

Implement standard k–epsilon and high Reynolds number CLS as active turbulence
models in `IncNavierStokesSolver`, with separate wall policies for a compatible
neutral log law, the original rough-wall analytical wall function, and the
equilibrium correction in the supplied note. Reuse the velocity-correction
solver and scalar field infrastructure, while adding variable diffusion, the
complete Reynolds stress, and pressure boundary conditions consistent with
those stresses.

Use a separate control-volume wall strip coupled to the CG outer domain, as
specified and tested in [Milestone 1](neutral-abl-rans-awf-milestone-1.md).
The note's wall compatibility proof concerns a cell average and one algebraic
dissipation row. Native tests show that its pointwise inner distributions can
have a nonzero higher source moment despite a zero mean. The strip retains the
original integrated row without inventing those additional inner equations.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | plan |
| Origin date | 2026-10-29 |
| Version label | neutral_abl_rans_proposal_v8 |
| Verification status | Milestones 1–7 CLOSED; Milestone 7 closed on 2026-10-30 with a fresh Linux build, eight CTests, 15 public-workflow cases, six standard-library benchmark preparations and rendered documentation; higher-degree 3D and terrain/physics limits remain explicit |
| Mac source baseline | `6522fc9c74edab0e445aade16f7facea3f5987a3` |
| New project branch | `codex/neutral-abl-rans-awf` |
| Linux source at proposal baseline | `812a2c84946d1fdad2b042c265e3e4ed32914138`, clean, on `codex/hit-inflow-generator`; later milestones use their recorded isolated revisions |
| Linux build configuration observed | Release; MPI, regression tests, and unit tests enabled |

The branch includes the existing HIT work. Milestone 1 has parent revision
`d8411938909169c0b30651fb6727a43c47cf0a19`. Its validation uses an isolated Linux
worktree and a separate prototype build against the installed libraries, leaving
the primary checkout, shared build, installation, and ongoing HIT experiments
unchanged. Existing untracked literature, guidance, and figures remain untouched.

[Milestone 2](neutral-abl-rans-awf-milestone-2.md) adds the reusable standard/CLS
closure and original/corrected AWF kernels, with independent tensor/quadrature
checks and a serial/MPI column regression. Its tested parent is `fd82c37`, which
adds the repository workflow instructions to the Milestone 1 branch. Source
hashes and completion evidence are recorded with that milestone.

[Milestone 3](neutral-abl-rans-awf-milestone-3.md) adds the opt-in standard-model
URANS solver, full variable symmetric stress, implicit losses, a coupled weak
pressure iteration, and a compatible log law at a flat standoff boundary.
The independent unsteady manufactured, decay, 3D, serial/MPI, boundary-control
and 30-second periodic ABL tests pass on Linux. This standard boundary does not
add the CLS/AWF conservative strip.

[Milestone 4](neutral-abl-rans-awf-milestone-4.md) adds full cubic CLS transport,
algebraic Yap and the original AWF through the conservative strip, including
its time-dependent inventory, algebraic dissipation and normal velocity.
All 27 CLS/AWF acceptance cases pass, with independent tensor/quadrature
forcing, serial/four-rank field audits, 3D coupling and conservation checks.
Standard transport and kernel regressions also pass. The original baseline
retains the note's source imbalance and evolves away from the CLS reference.
[Milestone 5](neutral-abl-rans-awf-milestone-5.md) couples the corrected
dissipation distribution and seven-sensor switch using the current configured
CLS reference. Its 34-case Linux acceptance checks native assembled residuals,
wall/volume budgets, source/point consistency, exact baseline recovery,
cutoff/regime handling and serial/MPI fields. The coarse 3D reference leaves
the switch plateau because of representation error; degree refinement restores
the plateau with unchanged thresholds. The original-AWF, standard-model and
kernel regressions also pass.

[Milestone 6](neutral-abl-rans-awf-milestone-6.md) completes the inlet/outlet,
restart and longer preservation study. Standard and both corrected CLS policies
meet the 0.1% target over five flows at degree four in 2D and full 3D; the
original AWF leaves that reference. Spatial/temporal, physical-H, perturbation,
numerical-column copy and repartitioned IMEX1/2 continuation results are
recorded with source/binary hashes. An additional degree-five 3D experiment
develops a growing mode and fails the target; it remains an explicit robustness
limit, rather than another preservation pass. Eight CTests and all standard,
original and corrected regression cohorts pass on the final solver.

[Milestone 7](neutral-abl-rans-awf-milestone-7.md) delivers the
[consolidated user handoff](../../solvers/IncNavierStokesSolver/Examples/NeutralABL/README.md),
public case/restart preparation, user/developer documentation and fresh Linux
delivery evidence. Checkpoint, relocated-history, cold-start and repartitioned
continuation workflows pass complete-field checks. Smoke inputs are staged in
the build directory, with no generated solver output in the source checkout.
All 127 recorded numerical C++/header sources remain identical to M6. The
delivery preserves that measured baseline and its higher-degree counterexample.

## Findings at the proposal baseline

At the proposal baseline, no implementation of the requested RANS closures or
wall functions was found in the inspected solver and library sources. The
following existing facilities provide useful integration points.

| Existing code or documentation | Consequence for this project |
| --- | --- |
| [`VelocityCorrectionScheme::v_InitObject`](../../solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp), especially lines 82–109 | All fields before the final `p` field are advanced. The existing `DiffusionCoefficient` function is evaluated once into a scalar per field. |
| `v_EvaluateAdvection_SetPressureBCs`, `v_SetUpViscousForcing`, and `v_SolveViscous` in the same file | There are explicit source and implicit solve hooks, but the current solves divide by constant field diffusivities. Spatially changing eddy viscosity requires a new assembly path. |
| [`VCSWeakPressure`](../../solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionSchemeWeakPressure.cpp) | A useful starting point for weak pressure forcing. Its inherited pressure and outflow machinery still assumes constant molecular viscosity. |
| [`IncBoundaryConditions::Initialize`](../../solvers/IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.cpp) | The custom boundary factory currently receives velocity components and pressure. A wall model that also needs `k`, dissipation, interior samples, and source assembly needs a richer interface. |
| [`IncBaseCondition::AddVisPressureBCs`](../../solvers/IncNavierStokesSolver/BoundaryConditions/IncBaseCondition.cpp) | The viscous pressure contribution is based on constant-viscosity curl–curl. Replacing its viscosity by a local number does not supply the general RANS boundary condition. |
| [`UnsteadyDiffusion::DoImplicitSolve`](../../solvers/ADRSolver/EquationSystems/UnsteadyDiffusion.cpp) and [`ContField`](../../library/MultiRegions/ContField.cpp) | Variable-coefficient Helmholtz solves and weak Neumann/Robin assembly already exist in the framework. |
| [`StdMatrixKey`](../../library/StdRegions/StdMatrixKey.cpp) and `ContField::v_UnsetGlobalLinSys` | Coefficient-dependent matrix keys and matrix cleanup need explicit attention when viscosity changes every step. |
| [`FilterReynoldsStresses`](../../solvers/IncNavierStokesSolver/Filters/FilterReynoldsStresses.cpp) | This measures covariance of resolved velocity fluctuations. The algebraically modelled RANS covariance needs separate output. |
| [Local scalar transport documentation](../user-guide/solvers/incompressible-ns.tex), section “Advecting extra scalar fields” | The established variable order is velocity, additional scalars, pressure. The [published user guide](https://doc.nektar.info/userguide/5.9.0/user-guidese54.html) documents the same convention. |

## Models and reference profiles

### Standard k–epsilon

Use the standard coefficients

\[
C_\mu=0.09,\quad C_{\epsilon1}=1.44,\quad C_{\epsilon2}=1.92,
\quad \sigma_k=1,\quad \sigma_\epsilon=1.3.
\]

Transport total dissipation `eps`. Use the linear Boussinesq covariance and
the usual two transport equations, with
\(\nu_t=C_\mu k^2/\epsilon\). Richards and Norris, section 2.1, give the
consistent high Reynolds number reference

\[
\kappa_{KE}=\sqrt{(C_{\epsilon2}-C_{\epsilon1})\sigma_\epsilon\sqrt{C_\mu}}
\simeq0.43267,
\]
\[
U_0=\frac{u_*}{\kappa_{KE}}\log(z/z_0),\qquad
k_0=\frac{u_*^2}{\sqrt{C_\mu}},\qquad
\epsilon_0=\frac{u_*^3}{\kappa_{KE}z},\qquad
\nu_{t,0}=\kappa_{KE}u_*z.
\]

Compute the implied kappa from the configured model constants. A prescribed
0.40 or 0.41 is a different coefficient choice unless the transport constants
are adjusted consistently. Support the paper's shifted-origin form as an
explicit profile option, with consistent heights in velocity, dissipation,
viscosity, and wall matching.

Provide a compatible rough log-law wall policy as the standard-model control.
For example, its friction scales can use
\(u_k=C_\mu^{1/4}\sqrt{k_P}\) and
\(u_\tau=\kappa_{KE}U_P/\log(z_P/z_0)\), with
\(\tau_w/\rho=u_k u_\tau\). The reference returns \(u_*^2\).
The associated turbulence source assembly must use the same momentum stress
and gradient operators. The finite-volume face formula in Richards and Norris,
section 5, motivates that consistency; it is not directly a spectral-element
formula. Do not import their CFX-specific sand-grain conversion or ground shift
as a universal wall law.

### CLS at high Reynolds number

Implement the complete matrix constitutive relation in the note's section 2,
including all quadratic and cubic terms, the variable \(C(q)\), and the
original algebraic Yap correction using geometric wall distance. In particular:

- Use \(S=\nabla U+(\nabla U)^T\) and
  \(\Omega=\nabla U-(\nabla U)^T\), without factors of one half.
- Use \(C(q)\) in the cubic coefficients as well as in eddy viscosity.
  The equilibrium value of \(C\) is not a replacement for its general function.
- Compute production from the same complete covariance used by momentum:
  \(P_k=-R_{ij}U_{i,j}\). The simple-shear reduction
  \(P_k=\nu_t(U')^2\) is a verification case.
- Keep three-component stresses even in a two-dimensional mean-flow case.
- Retain \(Y_c=(\widetilde\epsilon^2/k)\,0.83(r-1)_+r^2\), with
  \(r=k^{3/2}/(2.5\widetilde\epsilon d_w)\). Replacing \(d_w\) by
  \(d_w+z_0\) would change this model.

The reference constants from the note are

\[
q_*=3.190675081624288,\quad C_{\mu,*}=0.0982278953282113,
\quad\kappa_c=0.440553905186659.
\]

Solve their defining scalar equations from the configured constants, then
check these values for the default configuration. The reference has
\(k_0=q_*u_*^2\), \(\epsilon_0=u_*^3/(\kappa_cz)\), and
\(\nu_{t,0}=\kappa_cu_*z\). Its covariance in physical \((x,y,z)\) order is

\[
R_0/k_0=\begin{pmatrix}
0.92&0&-1/q_*\\0&0.56&0\\-1/q_*&0&0.52
\end{pmatrix}.
\]

Transport `eps_tilde` for CLS and expose total `eps` separately. In the strict
high Reynolds number outer approximation they coincide. If a later mode retains
the molecular split, use
\(\epsilon=\widetilde\epsilon+\nu|\nabla k|^2/(2k)\) consistently in
the k sink and wall prescription. This first project does not require the
complete low Reynolds number 1996 model, its damping, or its second-derivative
near-wall source.

### Molecular terms and the meaning of preservation

Provide two explicitly distinguished resolved-equation modes:

| Mode | Resolved diffusion | Purpose |
| --- | --- | --- |
| High Reynolds number outer equations | Turbulent diffusion; the stated high Reynolds number closure and dissipation convention | Test the analytical references under their derivation assumptions. |
| Molecular plus turbulent diffusion | Molecular viscosity retained with the high Reynolds number closure; any dissipation split must be declared | Practical URANS and assessment of the finite-viscosity departure from the reference. |

The AWF always retains its positive physical molecular viscosity. Merely keeping
it in the AWF does not cancel molecular residuals in the outer equations. A
finite-viscosity column may require a different numerical reference. The new
implicit assembly must also avoid the existing division by `Kinvis`; setting
`Kinvis=0` in the current solver is not an implementation of the first mode.

## Solver architecture

Add an opt-in `RANSVelocityCorrectionScheme`, initially derived from the weak
pressure velocity-correction implementation. Keep the existing solver path
available for regression comparisons. Use fields `u,v,w,k,eps,p` for standard
k–epsilon and `u,v,w,k,eps_tilde,p` for CLS in full 3D.

In a 2D mesh the existing solver uses `u,v` in mesh x,y coordinates. Define a
physical frame explicitly: mesh y can represent the note's vertical z, and
solver v then represents physical W. Embed its gradients in the physical 3 by 3
tensor before evaluating CLS. Full 3D uses the conventional x,y,z frame.

Write
\(R=(2k/3)I-\nu_tS+R^{NL}\), and store
\(\pi=p/\rho+2k/3\), with hydrostatic pressure removed when appropriate.
The velocity equation is then

\[
D_tU=-\nabla\pi+\nabla\cdot[(\nu_b+\nu_t)S-R^{NL}]+f,
\]

where \(\nu_b\) denotes the selected resolved molecular coefficient.
An initial segregated implementation can put
\(\nabla\cdot[(\nu_b+\nu_t)\nabla U_i]\) in the implicit solve and treat
the transpose-gradient contribution and nonlinear stress consistently in the
explicit/Picard part. Variable viscosity makes the transpose contribution
nonzero even for incompressible flow. A scalar effective Laplacian alone is
insufficient.

Use variable Helmholtz coefficients for momentum, k, and dissipation. Add
implicit loss coefficients where useful for the stiff turbulence sinks,
starting with a first-order verification scheme and then a second-order IMEX
scheme. All terms must recover the same stationary balance at a fixed point.
Update turbulence, stress, wall data, and pressure at consistent stages; check
temporal order and use Picard corrections when lagging prevents it. The AWF
auxiliary problem remains the note's quasi-steady problem at each update;
adding a wall-layer acceleration term would be a separately derived variant.

Add a RANS pressure extrapolation path based on the normal momentum balance
with the actual stress split. Include the turbulent and nonlinear stress in
wall, inlet, and outlet pressure treatment. Check reduced versus mechanical
pressure consistently when forming the AWF tangential pressure forcing.

Evaluate nonlinear closure products with adequate quadrature and projection.
Begin equilibrium tests with SVV, GJP, and unrelated forcing disabled. Add
dealiasing convergence checks, bounded matrix-cache lifetime, and iteration and
memory diagnostics. Use configurable positive k and dissipation guards and
report their activation. Any production clipping or realizability limiter must
be an explicit model option; it must remain inactive in the reference tests.

## Wall functions

### Shared original AWF kernel

Implement a field-independent kernel taking physical geometry, molecular
viscosity, current \(K=k_P\), matching tangential velocity, and signed pressure
and convective forcing. Return traction, auxiliary reconstruction, integrated
production, baseline dissipation mean, and baseline point value.

Retain the published roughness correlation, \(\alpha=0.2295\),
\(c_\ell^w=2.55\), the smooth-wall intercept 10.7, and dissipation cutoff
5.1. These are wall-model constants and remain distinct from CLS \(C(q)\).
Handle the smooth limit directly and all viscosity segments, including a
negative virtual viscosity intercept for a fully rough wall.

With the note's integrals, retain

\[
t_A=U_n/J_0-(J_1/J_0)F,\qquad
\overline P_A=H^{-1}\int_0^H\nu_{t,A}|U'_A|^2\,dz.
\]

Use the vector version for 3D tangential velocity and forcing. Preserve signed
reversed traction. Integrate at the viscosity and roughness breakpoints, using
stable logarithms and a nonnegative production evaluation near the molecular
limit. Verify closed forms against independent quadrature.

### Equilibrium correction from the note

Use the current positive K to construct the local counterfactual reference:

\[
u_*(K)=\sqrt{K/q_*},\quad t_0=K/q_*,\quad U_{n,0}=t_0J_0,
\quad A_c=u_*^3/\kappa_c,
\]
\[
P_0=\frac{t_0^2}{H}\int_0^H\frac{\nu_{t,A}}{D_A^2}\,dz,
\quad \Lambda=HP_0/A_c,\quad z_c=H\exp(1-\Lambda).
\]

Require \(z_c\le z_P\), equivalently
\(\Lambda\ge1+\log(H/z_P)\). A requested equilibrium case that fails this
condition must fail preflight. A general run can use an explicitly reported
original-AWF fallback. Clipping the cutoff would destroy the proposed budget.

Compute a single cellwise switch and derive both outputs from the same blend:

\[
\overline\epsilon_{EP}=(1-\chi)\overline\epsilon_A+\chi P_0,
\qquad
\epsilon_{P,EP}=(1-\chi)\epsilon_A(z_P)+\chi A_c/z_P.
\]

Keep original traction, auxiliary velocity, and production for every chi.
Evaluate P0 at zero forcing even in a forced current state. The switch uses
independent resolved information, including the next dissipation sample, and
senses pressure and convection separately. Their sum may cancel. Implement
the note's compact C2 switch with exact plateaus, expose weights and thresholds,
and record them as experimental parameters until calibrated.

Select the inlet velocity offset using the exact matching/interpolation
operator used by the wall assembly. Do not prescribe an independent CLS
aerodynamic z0 as well as h, H, viscosity, and friction velocity. The note's
anchoring generally makes z0 depend on H. Report it and distinguish changing H
from refining the outer discretization at fixed H.

### Spectral element wall coupling

Build a dedicated wall-layer adapter with a control-volume strip on `0 < z < H`
and a CG outer domain beginning at H. Physical P and N sampling locations,
matching, interface fluxes, strip inventories, and algebraic elimination are
fixed in the [Milestone 1 specification](neutral-abl-rans-awf-milestone-1.md).
Use log-coordinate matching for velocity, physical-height matching for k, and
the common auxiliary dissipation distribution for both EP and the outer trace
at H. This last interface prescription is an explicit CV-to-CG extension.

The adapter must satisfy the following contracts before it becomes the default:

1. Apply the AWF tangential momentum flux once, with the outward-normal sign.
   No slip is imposed in the auxiliary problem. Do not also impose strong
   tangential Dirichlet data on the same resolved PDE. Impose impermeability
   consistently with pressure and the selected layer formulation.
2. Account for the separate strip volume once. Supply its k contribution as
   \(\mathcal V(\overline P_A-\overline\epsilon_{EP})\) per unit density,
   the strip's one integrated source. There is no ordinary CG source integration
   over the same layer. Use equal and opposite interface flux contributions in
   the strip and outer weak rows.
3. Express the dissipation prescription as an algebraic evaluation constraint
   at P and replace the corresponding independent transport equation. This
   requires a constrained solve or elimination. The outer trace at H has the
   value from the same auxiliary distribution, replacing that boundary row.
   Do not also impose its flux. Retain strip inertia when eliminating UP/KP.
4. Provide a finite inner representation without evaluating log(z), 1/z, or
   the geometric-distance Yap expression at the wall. Separate physical
   auxiliary values from outer resolved values; the note does not require
   \(U_A(z_P)=U_0(z_P)\).
5. Keep P, N, and H independent of accidental quadrature-point locations as
   polynomial order changes. Validate sample ownership across MPI partitions.

The strip avoids an inner polynomial trial/test space. Its assembly contracts
are verified by the column prototype. Milestone 4 adds multidimensional
tangential transport, strip continuity coupled to pressure and nonlinear
boundary solves; Milestone 5 adds the corrected policy and switch. Milestone 6
verifies full-history restart. A conventional off-wall traction boundary
alone would omit the integrated source and strip inventory.

For the initial deliverable, restrict the wall geometry to flat, aligned
quadrilateral layers in 2D and extruded hexahedral layers in full 3D. General
terrain normals, closest-wall distance, and curved-layer geometry follow after
this assembly is verified. Homogeneous1D also needs separate treatment because
general variable coefficients can couple Fourier modes.

## Proposed source organization and configuration

| Location | Responsibility |
| --- | --- |
| `library/SolverUtils/Turbulence/` | Pure standard and CLS closure kernels, AWF kernels, and neutral-reference construction, separated from the existing HIT code. |
| `solvers/IncNavierStokesSolver/EquationSystems/RANSVelocityCorrectionScheme.*` | Active scalar transport, variable diffusion and source assembly, coupled updates, pressure convention, and residual inspection. |
| Coupled weak pressure implementation in `RANSVelocityCorrectionScheme.*` | Actual normal-velocity trace and full RANS momentum enter the pressure iteration. Flat x+ total-traction outlets are available; general curved/backflow or prescribed-pressure outlets remain outside the validated configuration. |
| `RANSMatrixFree.*`, `RANSHistory.*` | Current weak operator with bounded native preconditioners; complete checked IMEX solution/derivative continuation and MPI repartitioning. |
| `Examples/NeutralABL/` | Public case and continuation preparation with a consolidated user handoff. |
| Solver-local RANS wall-layer classes and boundary adapters | Sampling, geometry, constraints, weak fluxes, layer source replacement, MPI communication, and original versus corrected wall policy. |
| RANS extra field output and solver-local stage/wall CSVs | Modeled covariance, eddy viscosity, production, Yap, wall traction, switch, cutoff, admissibility and balance diagnostics; distinct from existing resolved-fluctuation filters. |
| `library/UnitTests/SolverUtils/` and `solvers/IncNavierStokesSolver/Tests/RANS/` | Algebraic tests, boundary/residual probes, small solver fixtures, and reference comparisons. |
| Existing CMake lists and incompressible solver user/developer guides | Registration, tests, configuration, equations, supported combinations, and limitations. |

The implemented configuration uses `TurbulenceModel=StandardKEpsilon|CLSHighRe`,
`RANSWallModel=None|AWFOriginal|AWFEquilibrium`, and an explicit resolved
molecular-term mode. The standard neutral log wall uses the `RANSLogLaw`
boundary tag with `RANSWallModel=None`. Geometry, physical frame, h or z0,
matching rule, guards, switch parameters and reference parameters must be
validated together. Incompatible closure/wall combinations are rejected.
The [standard k–epsilon AWF extension](standard-kepsilon-equilibrium-awf.md)
derives the standard constants and preservation proof and enables both wall
policies with that closure. The user guide and milestone
specifications give accepted XML rather than these initial design concepts.

Output and full-history restart preserve both transported turbulence fields, the
pressure convention, model constants, wall geometry, and any additional layer
state. Recompute derived fields consistently. Modelled RANS stresses and
covariance of unsteady resolved velocities must have distinct output names.

## Validation plan

The hypothesis is that the note's compatible wall operator, combined with a
consistent discretization of the specified outer equations, can maintain a
neutral constant-stress surface layer. The main response is downstream and
temporal drift of the prescribed profiles. Vary closure, wall policy,
resolution, timestep, matching height, and molecular-term mode separately.

### Algebra and local wall operators

Test zero gradient, simple shear, rigid rotation, pure strain, and rotated
general velocity gradients. Check covariance symmetry, trace for incompressible
states, frame covariance, all cubic terms, production contraction, and the
default equilibrium roots and stresses. Simple shear alone cannot exercise
the cubic terms that cancel in that flow.

For AWF test smooth, transitional, fully rough, and entirely molecular branches;
positive and reversed traction; cutoff equality and failure; switch plateaus;
and exact original-form recovery at chi=0. Check no slip, matching, positivity
of the auxiliary dissipation, the integrated k and mechanical-energy budgets,
and the fixed-input pressure identities. Include cases with cancelling pressure
and convective forcing to test independent sensors.

Reproduce the note's example with
\(u_*=0.4\,\mathrm{m/s}\), \(\nu=1.5\times10^{-5}\,\mathrm{m^2/s}\),
h=0.1 m, H=2 m, and zP=1 m. Its tabulated P0 is about 0.461092609,
whereas the original mean dissipation is about 0.774943967 in m2/s3. Test the
forcing inputs −0.02, 0, and 0.03 m/s2 and the unchanged traction response.
These values are reference targets from the note, not new solver results.
Use nondimensional, scale-aware algebraic tolerances around 1e-11, with tighter
checks where independent closed forms permit them.

### Assembled column residual

Construct a single repeated column using the same expansion, wall layer,
quadrature, interpolation, and boundary assembly as the eventual domain.
Evaluate all unconstrained weak residuals and algebraic constraints on the
projected/sampled analytical state. Inspect momentum, pressure/divergence, k,
and dissipation separately, including every wall-layer moment and the next
resolved layer. Check one complete timestep as well as the spatial residual.

Run this first with an off-wall exact boundary to isolate bulk assembly, then
with the actual wall treatment. A logarithm and 1/z are not generally exactly
represented by a polynomial expansion. Report projection/interpolation error
separately and perform h, p, and quadrature refinement. Do not infer machine
precision preservation from the continuous derivation.

Also solve a stationary discrete column and copy it downstream. Preservation
of this numerical column tests assembly and boundary consistency; agreement
with the analytical column is a separate refinement result. Do not add an
undisclosed reference residual subtraction or repeatedly reset the evolving
fields. An optional method subtracting a discrete equilibrium defect would
need a separate numerical-method specification and comparisons with it off.

### Empty domains and unsteady behavior

| Case | Purpose |
| --- | --- |
| Standard k–epsilon with compatible neutral log wall | Establish the Richards–Norris control and production consistency. |
| CLS with original AWF | Quantify drift from the original wall treatment at otherwise identical inputs. |
| CLS with corrected AWF, fixed chi=1 | Isolate equilibrium correction and admissibility. |
| CLS with corrected AWF and active sensors | Verify reference plateaus and response to controlled departures. |
| Molecular diffusion enabled | Quantify departure from the strict analytical high Reynolds number reference. |
| Deliberately inconsistent kappa or top dissipation flux | Demonstrate that the diagnostic detects boundary/profile drift. |

Use the paper's 300 m by 30 m empty-domain dimensions as one standard-model
benchmark, with U(10 m)=10 m/s and z0=0.02 m. Determine u* using the actual model
kappa. The CLS benchmark uses the note's independently anchored reference;
matching U at one height does not imply identical stress and k across models.

At the top, prescribe reference values or the compatible stress and scalar
fluxes. For the high Reynolds number CLS reference these include

\[
\nu_t U'=u_*^2,\quad \partial_z k=0,\quad
(\nu_t/\sigma_\epsilon)\partial_z\epsilon
=-u_*^4/(\sigma_\epsilon Z).
\]

Use the corresponding standard-model expressions with its kappa. Implement
flux signs and coefficient scaling according to the new weak assembly. A zero
top dissipation gradient is incompatible. Begin with streamwise periodic
repeated columns to isolate the vertical balance, then run inlet/outlet cases
with compatible pressure and scalar conditions. Use periodic spanwise faces
for the full-3D flat case.

Initialize the whole domain from the reference and run for at least five
flow-through times L/Uref, extending if drift has not settled. Sample profiles
at fixed streamwise stations and the same physical heights. Repeat with small
velocity and turbulence perturbations and with controlled forcing so an
algebraic equilibrium is not mistaken for perturbation stability. Verify
second-order temporal convergence on an unsteady manufactured or controlled
case. Test serial/MPI agreement and checkpoint/restart behavior, including
partitions without wall faces and P/N samples owned on different ranks.

Report height-dependent and normalized L2/L-infinity errors in U, k,
dissipation, viscosity, and stress; vertical velocity; divergence; wall
traction; individual production/dissipation/flux budgets; and every guard,
fallback, and switch activation. Normalize near-zero crossflow by u*, rather
than by the reference crossflow. Separate a first wall-layer error from errors
in the next and higher layers.

A proposed engineering target is less than 0.1 percent profile drift over five
flow-through times at the accepted resolution, with observed spatial and
temporal convergence. This is a target to test. A stronger exact-discrete claim
requires normalized residual and timestep defects at the configured linear and
nonlinear solve tolerance; a copied numerical column alone cannot establish
exact retention of analytical samples.

## Implementation milestones

| Milestone | Deliverable and completion condition |
| --- | --- |
| 1 Wall discretization and residual prototype | CLOSED. Conservative strip selected; P/N/H, matching, pressure/no-slip interpretation, integrated source and algebraic elimination fixed. Native column residual, conservation, h/p/quadrature, and MPI checks passed in the Milestone 1 specification. |
| 2 Independent closure and AWF kernels | CLOSED. Default/configured roots, full tensors, original wall branches, corrected budgets, forcing identities, admissibility, switch, and recovery pass independent Linux tests; native serial/MPI column regression passes with the general CLS kernel. See the Milestone 2 specification. |
| 3 Standard model URANS transport | CLOSED. Variable symmetric stress, implicit losses and coupled weak pressure pass unsteady manufactured, decay, 3D, quadrature, gauge and MPI tests. The compatible standard ABL control refines spatially and remains within 0.1% for one nominal periodic flow-through; native kernel and legacy VCS regressions pass. See the Milestone 3 specification. |
| 4 CLS and original AWF integration | CLOSED. Full cubic stress and Yap pass independent unsteady/3D controls and spatial refinement. Conservative original-AWF strip coupling passes point/mean/traction recovery, inventory budgets, nonuniform normal velocity, serial/MPI and invalid-state checks; the original dissipation mismatch remains measurable. See the Milestone 4 specification and its 27-case Linux evidence. |
| 5 Corrected AWF integration | CLOSED. Current configured reference, matching anchor, one integrated source, algebraic EP/EH and the compact switch pass 34 native Linux cases, including assembled residual refinement, independent forcing/sensors, exact original-form recovery, strict preflight/fallback, MPI and 3D sensor refinement. All original/standard/kernel regressions pass. See the Milestone 5 specification and evidence. |
| 6 Preservation and robustness study | CLOSED. Five-flow degree-four standard/fixed/switched profiles meet the 0.1% target in 2D/full 3D. h/p/dt and physical-H studies, perturbation decay, numerical-column copies, MPI and complete IMEX continuation are recorded. Original/inconsistent controls and a degree-five 3D counterexample retain their measured failures. See the Milestone 6 report and closure evidence. |
| 7 Documentation and delivery | CLOSED. Consolidated user/configuration and developer guides, public case and continuation helpers, complete commands and final source/binary evidence are delivered. A fresh Linux build passes eight CTests and 15 workflow cases; six standard-library preparations match M6, and both manuals compile and pass RANS-page visual review. Known numerical and application limits remain explicit. |

Milestones 1 and 2 began with small column and kernel prototypes. Terrain,
thermal AWF, rotating/stratified ABL and complete low-Re CLS extensions require
new implementation and validation beyond this flat neutral baseline. The
AWF paper's dune validation uses a differential length correction; it is not
direct evidence for the note's retained algebraic Yap variant.

## Linux validation and source synchronization

Use the Mac checkout for implementation and Git review. Before any Linux build,
record the branch, HEAD, status, and relevant file hashes in both checkouts.
Transfer only the intended branch commits or a reviewed patch; preserve
uncommitted work and confirm the resulting revision and diff. Each milestone
records its isolated worktree rather than treating the unrelated primary
checkout as already synchronized. Both environments have a GitHub SSH push URL.

Compile and run on Linux using the isolated external-build recipe in the
[user handoff](../../solvers/IncNavierStokesSolver/Examples/NeutralABL/README.md).
It builds the production solver source list against the existing installation
without replacing the shared build or binaries used by other experiments.

The project has a hard limit of 16 simultaneous Linux CPU workers, as directed
on 29 October 2026. The delivered recipe uses four build workers and one
numerical thread per MPI rank. Run focused
new unit tests, new solver regressions, and related legacy scalar, pressure,
wall, and HIT tests. Select explicit CTest concurrency according to each
test's ranks and threads, initially one test at a time. Start solver pilots
with 4 MPI ranks and one OpenMP/BLAS thread each. Count concurrent compiler,
MPI, test-runner, and other threaded workers together: never exceed 16 for this
project or the repository-wide 72-worker total across concurrent server work.

Store each experiment and its TMPDIR under
`/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/<experiment>`.
Keep large meshes, fields, checkpoints, and logs out of the source checkout.
Record commands, revision/diff, constants, mesh and order, timestep, solver
tolerances, MPI/thread allocation, norm definitions, and diagnostic time series
in a run manifest. Monitor residuals, finite values, progress, and resource use;
set an explicit case timeout when each pilot is sized. Small fixtures and
summarized results belong in the repository.

## Source references

- Supplied `cls_neutral_abl_equilibrium_awf_note.tex`, revised 26 October 2026. Sections 2–4 define the closure/reference; sections 5–9 define the wall operator; section 11 states its limits. The tracked [Milestone 1](neutral-abl-rans-awf-milestone-1.md), [Milestone 2](neutral-abl-rans-awf-milestone-2.md) and [Milestone 5](neutral-abl-rans-awf-milestone-5.md) specifications preserve the implemented contracts.
- Craft, Launder, and Suga (1996), [Development and application of a cubic eddy-viscosity model of turbulence](https://doi.org/10.1016/0142-727X(95)00079-6). Constitutive relation, coefficients, transport, and algebraic Yap correction.
- Suga, Craft, and Iacovides (2006), [An analytical wall-function for turbulent flows and heat transfer over rough walls](https://doi.org/10.1016/j.ijheatfluidflow.2006.03.011). Sections 2.1–2.2 and Tables 1–2 define original momentum, production, dissipation and roughness treatment.
- Richards and Norris (2011), [Appropriate boundary conditions for computational wind engineering models revisited](https://doi.org/10.1016/j.jweia.2010.12.008). Standard-model reference and boundary conditions, empty-domain test and production discretization.

The local paper transcriptions/PDFs and supplied note remain optional working
copies under `docs/literature/`; they are not needed to generate or run the
delivered examples and are not included in the Git delivery.
