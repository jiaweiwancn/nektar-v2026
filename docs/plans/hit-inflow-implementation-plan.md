# Homogeneous isotropic turbulence inflow implementation plan

Implement an incompressible, homogeneous isotropic turbulence (HIT) inlet for
Nektar++ using the Fourier synthesis method in Mann, *Wind field simulation*
(1998). Use Batchelor's isotropic spectral tensor with either a von Karman
energy spectrum or a supplied experimental three-dimensional energy spectrum.
Validate the generated inlet first, then simulate its downstream decay and
compare with the grid-turbulence measurements of Comte-Bellot and Corrsin
(1971).

The implementation will use a common three-component synthesis core for full
three-dimensional meshes and `Homogeneous1D` expansions. The first spatial decay
benchmark will use a two-dimensional mesh in x and y, a Fourier expansion in z,
periodic transverse boundaries, and the nonlinear velocity-correction solver.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | plan |
| Origin date | 2026-10-28 |
| Version label | hit_inflow_plan_v3 |
| Verification status | Milestone 2 core build and numerical tests VERIFIED on Linux; inlet integration and downstream validation pending |
| Milestone status | Milestones 1 and 2 CLOSED on 2026-10-28; Milestones 3-8 pending |
| Source baseline | `a3b7efc8fef7e075390f0f59457dfdac30f0b6ed` |
| Project branch | `codex/hit-inflow-generator` |

The binding version 1 interface, solver support, random-number, normalization,
and restart contracts are recorded in the completed
[Milestone 1 specification](hit-inflow-milestone-1.md). Its configuration fixtures
and experimental inputs are checked. The
[Milestone 2 core and its focused numerical checks](hit-inflow-milestone-2.md)
are complete; the runtime inlet and downstream validation remain scheduled work.

## Scope and intended behavior

The first implementation supports a planar inlet normal to the x direction,
a constant mean velocity `Ubar > 0`, all three fluctuating velocity components,
and an isotropic, nonnegative scalar spectrum `E(k)`. Full three-dimensional
meshes and `Homogeneous1D` must use the same spectral normalization and random
realization conventions. A solver with only two physical velocity components
cannot represent the requested three-dimensional HIT.

The inlet supplies a stationary sequence by advecting a frozen synthetic volume
past the inlet plane. The Navier-Stokes solver supplies subsequent nonlinear
evolution, dissipation, and loss of correlation. Prescribing a downstream decay
envelope or continuing turbulence forcing through the measurement region would
invalidate that decay test.

Atmospheric shear, the sheared Mann tensor, spatially varying turbulence
intensity, arbitrary curved inlets, compressible flow, and non-Gaussian
statistics are outside this first implementation. `Homogeneous2D` is a later
extension after its MPI and dealiasing behavior has been assessed. The existing
Synthetic Eddy Method remains a separate forcing option.

## Scientific conventions

### Spectral tensor and Gaussian coefficients

Use angular wavenumbers in radians per metre and velocity fluctuations in
metres per second. Define `k = |k|`, with

\[
\Phi_{ij}(\mathbf{k})
=\frac{E(k)}{4\pi k^2}
\left(\delta_{ij}-\frac{k_i k_j}{k^2}\right).
\]

The energy convention is

\[
K=\frac12\langle u_i'u_i'\rangle
=\int_0^\infty E(k)\,dk,
\qquad
\langle u'^2\rangle=\langle v'^2\rangle=\langle w'^2\rangle=\frac23K
\]

in the isotropic ensemble. Thus `E(k)` has units of cubic metres per square
second. Record this convention in the parser, user documentation, test data,
and diagnostic output.

For a synthetic box with lengths `Lx`, `Ly`, and `Lz`, use the Fourier lattice
`k_i = 2*pi*m_i/L_i` and cell volume `DeltaVk = DeltaKx*DeltaKy*DeltaKz`.
An equivalent form of Mann's Eqs. (46) and (50) is

\[
\widehat{\mathbf u}(\mathbf k)
=\sqrt{\frac{E(k)\Delta V_k}{4\pi k^4}}
\left[\boldsymbol\xi(\mathbf k)\times\mathbf k\right],
\qquad k>0.
\]

For each independent mode, draw
`xi_j = (a_j + i*b_j)/sqrt(2)` with independent standard normal real variables
`a_j` and `b_j`. Generate one member of each conjugate pair, then set
`uhat(-k) = conj(uhat(k))`. Do not independently draw both partners or add an
extra factor of two to their covariance. Set the fluctuating zero mode to zero.
Handle self-conjugate and Nyquist modes explicitly; the initial implementation
may omit modes at a Nyquist limit and use a cutoff below that limit.

This construction gives

\[
\langle\widehat u_i^*\widehat u_j\rangle
=\Phi_{ij}(\mathbf k)\Delta V_k,
\qquad
\mathbf k\cdot\widehat{\mathbf u}=0.
\]

The tensor is positive semidefinite with rank two for nonzero `k`. Its singularity
is physical. Use the explicit transverse factorization rather than a generic
positive-definite Cholesky algorithm that fails on its zero eigenvalue.

### Supported spectra

For the von Karman option, implement Mann's Eq. (9):

\[
E(k)=\alpha\epsilon^{2/3}L^{5/3}
\frac{(Lk)^4}{[1+(Lk)^2]^{17/6}}.
\]

Use an unambiguous amplitude convention. Allow the user to specify either the
amplitude `alpha*epsilon^(2/3)` or a target one-component RMS, with the latter
implemented by a single common scaling of the vector field. Define which
quantity controls normalization and reject conflicting inputs. The spectrum
parameter `L` is distinct from the measured longitudinal integral scale `Lf`
and transverse integral scale `Lg`; document any conversion used during fitting.

The unmodified von Karman model has an inertial tail and no viscous cutoff.
Expose an explicit resolved band and a documented optional smooth radial
high-wavenumber taper. After truncation or tapering, calculate and report the
represented energy and

\[
\epsilon_{\rm represented}=2\nu\int k^2E(k)\,dk
\]

or the corresponding discrete sum. A parameter named `epsilon` in the model
amplitude is not automatically the realized dissipation of the truncated field.

For experimental input, accept a text table of strictly increasing positive
`k` and nonnegative `E(k)`, with an explicit unit convention. Validate finite
values, minimum table length, ordering, and range. Specify interpolation and
both end treatments in the configuration. Prefer interpolation that preserves
nonnegativity and does not overshoot the measured dissipative tail. Treat zero
values separately if logarithmic interpolation is used.

Do not differentiate a noisy longitudinal spectrum inside the solver. The
Comte-Bellot and Corrsin paper already supplies the isotropically reconstructed
three-dimensional spectrum in Table 3. Preserve provenance, interpolation,
low-wavenumber extension, and high-wavenumber extension in a reproducible data
preparation step. Supply only small necessary numeric fixtures in the repository.

### Time dependence and finite representation

Evaluate the inlet from a fixed realization:

\[
\mathbf u_{\rm in}(y,z,t)=U_c\mathbf e_x+
\sum_{\mathbf k}\widehat{\mathbf u}(\mathbf k)
\exp\{i[k_x(x_{\rm in}-U_c(t-t_{\rm ref}))+k_yy+k_zz]\}.
\]

Here `Uc = Ubar` for the benchmark. Frequency follows the streamwise component,
`2*pi*f = Uc*kx`, rather than the magnitude `k`.

Calls at the same absolute time must return the same field, regardless of call
order, rejected steps, or restart history. Never draw new modes per timestep.
The repetition time is `Tbox = Lx/Uc`; require it to exceed the useful recording
interval, or report explicitly when data repeat. Do not splice independent
volumes without a separate continuity and divergence-preservation design.

Match the synthetic transverse periods to the computational periods. Select
spectral cutoffs that the inlet expansion, streamwise mesh, Fourier modes, and
timestep all resolve. Use common radial tapering for all components when an
isotropic cutoff is intended. Record the effect of the finite lattice on
variance, tensor covariance, isotropy, and integral scales.

Mann's Section 6.1 and Figure 13 show low-wavenumber errors in small boxes.
Establish transverse-size convergence rather than treating `8L` as a universal
sufficient bound. Keep the explicit transverse factorization for the first
version. If the tensor convolutions in Eqs. (44)-(47) are added later, verify
their transversality at the actual discrete wave vector before claiming
solenoidality. Averaging tensors at neighbouring wave vectors does not by itself
establish that property; re-projection would also change the target covariance.

## Divergence and mass compatibility

The synthetic three-dimensional Fourier field is analytically divergence-free
because every mode is transverse. A two-dimensional inlet plane is a trace of
that field. It must not be generated by imposing
`d(v')/dy + d(w')/dz = 0`: those terms generally balance `d(u')/dx` in the full
field. Mann's reduced two-dimensional-domain construction is not a substitute
for the required three-component convected-volume construction.

Evaluate the spectral representation at inlet quadrature coordinates. Cache
spatial factors and use FFTs where the sample layout permits them. Establish a
correct reference evaluator before optimizing. Ordinary independent trilinear
interpolation of velocity components must not be assumed divergence-preserving.
Apply any scalar energy normalization to all components together; independently
rescaling components can destroy transversality.

For a complete transverse periodic plane, the streamwise fluctuating mode with
`ky = kz = 0` vanishes for `kx != 0`, and the zero mode is zero. The ideal inlet
therefore has zero plane-integrated fluctuating normal flux. Check this property
after projection using the actual surface quadrature and MPI reductions. Prefer
correcting representation and quadrature errors over introducing an arbitrary
spatial taper or a fluctuating bulk flow.

Distinguish three validation quantities: mode transversality, analytical
three-dimensional divergence, and the solver's discrete divergence. Pressure
projection does not guarantee machine-zero pointwise divergence of sampled
boundary data, and cannot freely alter prescribed Dirichlet traces. Check net
boundary flux and convergence in space and time in addition to the pressure
solver residual.

## Boundary and solver configuration

| Location | Configuration | Purpose |
| --- | --- | --- |
| Inlet | Three velocity Dirichlet conditions supplied by one `HITInflow` instance | Preserve the tensor correlations between components |
| Inlet pressure | Neumann with `USERDEFINEDTYPE="H"` | Include momentum-consistent pressure data and prescribed normal acceleration |
| y boundaries | Paired periodic conditions for `u,v,w,p` | Approximate transverse homogeneity without walls |
| z direction in Homogeneous1D | Fourier periodicity with matching `LZ` | Represent spanwise fluctuations with multiple modes |
| z faces in full 3D | Paired periodic conditions for `u,v,w,p` | Match the same periodic synthetic field |
| Outlet | Existing `HOutflow` family, with an outlet buffer | Let evolved vortical structures leave the domain |

The following inlet syntax is proposed and becomes valid only after the new
boundary type is implemented:

```xml
<REGION REF="0">
    <D VAR="u" USERDEFINEDTYPE="HITInflow" VALUE="Ubar" />
    <D VAR="v" USERDEFINEDTYPE="HITInflow" VALUE="0" />
    <D VAR="w" USERDEFINEDTYPE="HITInflow" VALUE="0" />
    <N VAR="p" USERDEFINEDTYPE="H" VALUE="0" />
</REGION>
```

Define `VALUE` as the mean component and add the generated fluctuations once.
Keep stochastic settings in one region-associated configuration block, with
seed, spectrum, box dimensions, cutoff, time origin, and mean convection speed.
Finalize its exact XML schema against SessionReader conventions in milestone 1.
Do not require the expression parser to evaluate random fields.

`H` requests a computed high-order pressure boundary condition. Its zero XML
value does not request a literal zero pressure gradient. Supply updated velocity
coefficients before the pressure extrapolation reads them.

Begin with the standard semi-implicit velocity-correction scheme and nonlinear
advection, fixed steps, and blank-variant IMEX order 2 with its order-1 startup.
Use the explicit supported-combination matrix in the Milestone 1 specification;
do not treat other integrators or solver variants as supported by default.
Enable appropriate spectral/hp and homogeneous dealiasing for the
resolved decay study. Keep physical viscosity distinct from artificial
stabilization. If spectral vanishing viscosity (SVV) is required, report its
parameters and its contribution to the energy budget, and describe that run as
an underresolved or modeled calculation unless DNS resolution is established.

Compare a simple velocity-Neumann/pressure-Dirichlet outlet with `HOutflow`,
including the convective Robin variant when useful. Move the outlet and extend
the buffer to establish insensitivity at the last measurement station. Keep any
outlet damping outside the measurement region. No-slip and free-slip lateral
walls are inappropriate for the primary periodic HIT benchmark. The experiment
has tunnel walls, so this setup models its central region rather than the full
tunnel and contraction.

## Code architecture and integration points

Paths in this section are relative to the repository root. Names of new classes
and files are proposed and may be adjusted to the existing coding conventions.

| Area | Existing source | Planned change |
| --- | --- | --- |
| Spectrum and synthesis | `library/SolverUtils/` | Add reusable spectrum evaluation, mode generation, deterministic realization, and sampling helpers |
| Inlet class | `solvers/IncNavierStokesSolver/BoundaryConditions/` | Add `HITInflow.h/.cpp` deriving from or fitting the existing boundary abstraction |
| Factory | `BoundaryConditions/IncBoundaryConditions.cpp` | Register the creator and recognized boundary type; associate one object with each inlet region |
| Update scheduling | `EquationSystems/IncNavierStokes.cpp` and boundary manager | Add a prescribed-velocity update phase before pressure-boundary evaluation |
| Pressure coupling | `EquationSystems/VelocityCorrectionScheme.cpp`, `StandardExtrapolate.cpp`, `Extrapolate.cpp` | Reuse existing pressure acceleration treatment and verify current-time inlet coefficients are visible |
| Build | `library/SolverUtils/CMakeLists.txt` and `solvers/IncNavierStokesSolver/CMakeLists.txt` | Register source files, exported headers where needed, and regression cases |
| Tests | `library/UnitTests/` and `solvers/IncNavierStokesSolver/Tests/` | Add focused algebraic, sampling, and solver integration checks |
| User guide | `docs/user-guide/solvers/incompressible-ns.tex` and appropriate XML documentation | Explain configuration, conventions, supported expansions, and an example |

The key scheduling issue is visible at the source baseline:

1. `IncNavierStokes::v_PreIntegrate` calls
   `SetBoundaryConditions(m_time + m_timestep)`.
2. `StandardExtrapolate::v_EvaluatePressureBCs` calls `AddDuDt`.
3. `Extrapolate::IProductNormVelocityBCOnHBC` reads the prescribed velocity
   boundary coefficients.
4. `VelocityCorrectionScheme::v_EvaluateAdvection_SetPressureBCs` currently
   calls the general `IncBoundaryConditions::Update` after
   `EvaluatePressureBCs`.

Adding the turbulence inlet only to the last call would make its acceleration
stale. Add an explicit prescribed-velocity hook at the end of
`SetBoundaryConditions(time)`, preserving pressure-aware post-pressure updates.
For this IMEX variant, the RHS also runs at its completed stage time, so do not
add another inlet advance using `rhsTime+TimeStep` inside pressure evaluation.
Load/broadcast settings on all session ranks before local boundary construction;
ranks with no inlet cannot participate in collectives inside an inlet constructor.
Keep the initial supported integrator restriction until manufactured-inlet tests
establish the numerical behavior. See the completed source audit for anchors.

For full 3D, sample physical coordinates and project onto the face expansions.
For `Homogeneous1D`, follow
`DisContField3DHomogeneous1D::v_EvaluateBoundaryConditions`: set physical boundary
values, call `SetWaveSpace(false)`, then `FwdTransBndConstrained`. The latter
already performs the homogeneous forward transform. Use existing coordinates,
global plane identifiers, and communicator-aware transforms, rather than local
plane indices or a second Fourier transform. Nektar's real Fourier packing and
normalization differ from an unnormalized complex FFT; validate the conversion.

Use a reproducible mapping from the seed and global independent-mode identifier
to Gaussian coefficients. Avoid per-rank random streams whose realization
changes with partition count. Define the RNG and normal-conversion algorithm as
part of the restart contract. Store a configuration fingerprint, seed, time
origin, and any realization data that cannot be reconstructed exactly. Report
the represented spectrum and input-data checksum at startup.

Do not copy the zero-mode-only behavior of `MovingFrameFar` for a fluctuating
inlet. Do not adapt `ForcingIncNSSyntheticEddy` into the new model: its volumetric
source-term algorithm has a different physical and mathematical role.

`Homogeneous2D` is recognized by the solver and has existing regressions, but
`ExpListHomogeneous2D.cpp` rejects homogeneous dealiasing when the column
communicator size exceeds one. Assess available decompositions and diagnostics
before adding that backend. This limitation is not a mathematical reason to
restrict the synthesis core to a particular mesh representation.

## Experimental benchmark definition

### Data and downstream mapping

Use the principal 5.08 cm grid dataset in Comte-Bellot and Corrsin (1971). Its
approach speed is `U0 = 10 m/s`. A contraction of 1.27:1 produces a measurement
region mean speed `Ubar = 12.7 m/s`, explicitly stated in Figure 1 and the
footnote on printed page 289. The local Markdown transcription omits that
footnote, so use the PDF for this parameter.

The station coordinate is

\[
\theta=\frac{U_0t}{M},\qquad
t=\int_0^{x_1}\frac{dx}{\overline U(x)},\qquad M=0.0508\ {\rm m}.
\]

For an inlet representing `theta = 42`, distances in the constant-speed
measurement region follow

\[
x-x_{\rm in}=\overline U\frac{M}{U_0}(\theta-42).
\]

Do not identify `theta` with `x/M` or use `U0` as the actual convective speed.
The table below converts the paper's Table 4 to SI units; distances and elapsed
times after the inlet are calculated from the preceding definitions.

| theta | Time after inlet s | Distance after inlet m | Streamwise RMS m/s | Dissipation m2/s3 | eta m | lambda m | Lg m | Lf m | Re lambda |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 42 | 0 | 0 | 0.222 | 0.474 | 0.00029 | 0.00484 | 0.0127 | 0.0240 | 71.6 |
| 98 | 0.28448 | 3.612896 | 0.128 | 0.0633 | 0.00048 | 0.00764 | 0.0188 | 0.0345 | 65.3 |
| 171 | 0.65532 | 8.322564 | 0.0895 | 0.0174 | 0.00066 | 0.0102 | 0.0228 | 0.0490 | 60.7 |

Table 4 implies approximately `nu = 1.5e-5 m2/s`, both from `u_rms*lambda/Re`
and from `epsilon*lambda^2/(15*u_rms^2)`. Figure 9's caption instead uses
`2*nu = 0.28 cm2/s`. Document this published discrepancy, use the Table 4
consistent value as the initial choice, and assess its effect on dissipation
comparisons rather than asserting the paper specifies one exact viscosity.

The approximate streamwise variance decay law, Eq. (38), is

\[
\frac{U_0^2}{\langle u_1'^2\rangle}
=21(\theta-3.5)^{1.25}.
\]

The transverse coefficient is about 20. The measurements are nearly isotropic,
with residual component and large-scale correlation differences. Compare an
exactly isotropic generator against these data with that distinction preserved.

### Spectrum preparation

Table 3 gives the required three-dimensional `E(k)` at `theta = 42, 98, 171`.
Use the first station to define the experimental-spectrum inlet. Keep the later
stations as downstream prediction targets. Table 2 gives the more directly
measured one-sided longitudinal spectrum `E11^(1)(k1)` for virtual-probe
comparisons.

The unit conversions are

\[
k_{\rm SI}=100k_{\rm cm^{-1}},\qquad
E_{\rm SI}=10^{-6}E_{\rm cm^3/s^2}.
\]

These wavenumbers are angular wavenumbers; no further factor of `2*pi` belongs
in this unit conversion. The one-sided and shell spectra obey

\[
\int_0^\infty E_{11}^{(1)}(k_1)\,dk_1=\langle u_1'^2\rangle,
\qquad
E(k)=\frac12k^3\frac{d}{dk}
\left[\frac1k\frac{d E_{11}^{(1)}}{dk}\right].
\]

The second relation is the paper's Eq. (46) and explains its factor of one half.
Mann uses a two-sided one-dimensional spectrum, which must be doubled for this
positive-axis comparison.

The Table 3 spectrum at `theta = 42` ends at `20 cm^-1`, while its Kolmogorov
wavenumber is about `34 cm^-1`. The unmeasured tail still matters to the
dissipation integral. Prepare and record a smooth nonnegative tail and
low-wavenumber extension, quantify their energy and dissipation contributions,
and assess alternatives. Do not claim that hard truncation plus energy
renormalization reproduces Table 4 dissipation.

Fit the von Karman variant to upstream energy and scale information only.
Report its dissipative-range treatment separately. Avoid adjusting inlet
parameters to fit both downstream stations and then calling those same stations
an independent validation.

### Domain and initialization

Use a rectangular x-y mesh, a periodic y pair, and multiple Fourier modes in z.
Set the main convection speed to `12.7 m/s`. The domain must extend beyond the
8.322564 m station plus an outlet buffer. Choose its transverse dimensions from
the inlet spectrum and the largest downstream integral scales, then establish
size convergence. The experimental tunnel width is not automatically the
required width of a periodic central-region model.

Assess the cost of resolving the inlet `eta` of approximately 0.29 mm across
this long domain before committing to production resolution. Use `kmax*eta`,
represented dissipation, local mesh quality, time resolution, and refinement
results together; a nominal polynomial order alone does not establish DNS.

Use the same synthetic realization at the inlet and, where practical, for a
consistent initial volume. A Gaussian field matches prescribed second-order
statistics but does not start with an established nonlinear cascade. Discard
domain flushing and statistical startup, then assess the persistent inlet
adjustment region. Initial-volume consistency alone does not remove that region.

The direct station-42 inlet experiment provides a controlled first test. If its
adjustment region contaminates the quantitative comparison, add a precursor or
an upstream development region and ensure that the developed reference plane
matches station 42 before predicting later stations. Record the effective
reference location; do not conceal adjustment with a fitted virtual origin.

An optional temporally decaying periodic-box calculation using the same
synthesis core can isolate discretization and initial-spectrum effects before
the long spatial run. It supplements the requested inflow-decay validation.

## Measurements and comparison methods

Use `HistoryPoints` in physical space at several positions on each measurement
plane to obtain virtual hot-wire signals. Record the sample interval and
anti-aliasing treatment. Subtract the local mean before computing velocity
variance or spectra. Use `ReynoldsStresses` or `AverageFields` after startup to
collect mean and second-moment fields, and average over the transverse plane
when testing homogeneity.

For a one-sided fixed-probe PSD `Sf(f)` in frequency, compare with Table 2 using
both the coordinate conversion and its Jacobian:

\[
k_1=\frac{2\pi f}{\overline U},\qquad
E_{11}^{(1)}(k_1)=\frac{\overline U}{2\pi}
S_f\!\left(\frac{\overline U k_1}{2\pi}\right).
\]

Specify the PSD estimator, window, segment length, overlap, normalization, and
frequency resolution. Verify that its integral recovers the signal variance.
Report statistical uncertainty from independent blocks or realizations.

Compute component variances, turbulent kinetic energy, mean shear,
off-diagonal stresses, longitudinal and transverse correlations, integral
scales, Taylor scale, and `Re_lambda`. Compare dissipation from velocity
gradients with energy-decay estimates. In an open spatial domain, account for
transport and boundary terms in the energy balance; do not equate the complete
domain energy change with viscous dissipation alone.

The central two-probe correlation is

\[
R_{11}^{c}(\Delta x)=
\frac{\langle u_1'(x_0,s)u_1'(x_0+\Delta x,s+\Delta x/\overline U)\rangle_s}
{u_{1,\rm rms}'(x_0)u_{1,\rm rms}'(x_0+\Delta x)}.
\]

This follows the mean convection and should lose correlation through physical
evolution. A purely frozen field predicts unity. Use Table 5 and Figures 10-12
for the full-band comparison. Add the matched narrow-band correlations in
Tables 6-7 and Figures 14-16 only after the full-band, energy, and spectrum
checks are satisfactory. Reproduce the filter convention when comparing them.

`ModalEnergy` reports integrated energy per homogeneous mode, and FieldConvert
`powerspectrum` reports homogeneous-direction spectra for `Homogeneous1D`.
Neither is the three-dimensional shell spectrum or a station-local streamwise
spectrum. For mean flow in x, the z spectrum of u is transverse; the z spectrum
of w is longitudinal only by reference to its component direction, and can be
compared with a streamwise longitudinal spectrum under an isotropy assumption.
Sample a suitable local volume and define its window and transform if a shell
spectrum is needed. Report window effects in a streamwise inhomogeneous flow.

## Implementation milestones and completion criteria

### Milestone 1 Finalize interfaces and benchmark inputs

**Status: CLOSED on 2026-10-28.** The
[interface specification](hit-inflow-milestone-1.md), configuration fragments,
source transcription, SI fixture files, and preparation/contract checkers are
saved. The [closure record](hit-inflow-milestone-1-validation.md) and
[Linux receipt](hit-inflow-milestone-1-validation.json) record the tested snapshot.
Both XML examples, the mean-only limit, 26 invalid configurations, all six
spectra, station mapping, units, and data digests passed the one-worker Linux
checks. This closes the interface/input milestone, not the future runtime
parser, stochastic implementation, or solver/numerical milestones.

Confirm branch and checkout states. Finalize the region configuration schema,
RNG contract, Fourier convention, spectrum file format, supported time
integrator, and startup behavior. Prepare the small experimental tables from
the original PDF with units and provenance. Audit the pre-pressure update path
and identify the smallest change that preserves existing boundary classes.

Complete when an example configuration is unambiguous, experimental inputs have
verified units and station mapping, and the exact supported combinations of
solver, expansion, and boundary type are recorded.

### Milestone 2 Implement spectra and transverse mode generation

**Status: CLOSED on 2026-10-28.** The reusable C++ spectrum and Fourier synthesis
core is implemented under `library/SolverUtils/Turbulence/`. All 19 focused
Boost cases passed on Linux at tested revision
`a6790a4a505a78e363759383478599efe2bd837b`, including independent tensor and
derivative checks and two 4096-seed statistical ensembles. The
[closure record](hit-inflow-milestone-2.md) and
[Linux validation receipt](hit-inflow-milestone-2-validation.json) contain the
API, exact commands, numerical diagnostics and worker allocation. The library
and public headers were rebuilt and installed. Solver inlet integration remains
in Milestone 3 onward.

Implement von Karman and tabulated spectra, validation, cutoffs, Gaussian
coefficients, conjugate pairing, zero-mode handling, deterministic mode IDs,
and represented-energy diagnostics. Implement a reference Fourier evaluator
that exposes enough information for analytic derivative checks.

Complete when algebraic tests establish tensor covariance, normalization,
real-valued synthesis, and transversality. Compare Gaussian ensemble statistics
with their sampling confidence intervals, not an exact isotropy requirement on
each finite realization. Invalid inputs must fail with useful messages.

### Milestone 3 Integrate the full three dimensional inlet

Implement the vector boundary class and its pre-pressure update. Register the
factory type, build sources, and documentation entry. Project samples onto a
small periodic-sided full 3D test domain. Verify initialization, mean-only
behavior, normal acceleration, and outlet compatibility.

Complete when a deterministic time-dependent manufactured inlet retains the
solver's expected temporal accuracy and the new type does not disturb existing
boundary regressions. Repeated evaluations at the same time must be identical.

### Milestone 4 Add Homogeneous1D and MPI support

Implement physical boundary sampling with the existing distributed transforms.
Verify real Fourier packing, nonzero-mode content, period matching, all three
components, and absence of a double transform. Exercise both mesh partitioning
and homogeneous partitioning, restart, and FFTW versus matrix transforms where
both are available.

Complete when full 3D and Homogeneous1D samples agree at common coordinates
within projection error, and the inlet realization is independent of MPI
partition count. Require numerical agreement within reduction and solver
tolerances for evolved statistics rather than bitwise equality of MPI solutions.

### Milestone 5 Establish focused regression coverage

Add small reproducible `.xml` and `.tst` cases to the incompressible solver
tests and register them with `ADD_NEKTAR_TEST`. Add core checks to the existing
unit-test structure. Include a mean-only limit, one deterministic transverse
mode with known time dependence, a seeded multimode field, and a tabulated
spectrum case. Include a restart comparison and a small MPI comparison.

For nonzero modes, target normalized transversality below approximately
`1e-12` in double precision. Check analytical derivatives near floating-point
accuracy using a dimensionally normalized divergence measure. Establish solver
divergence and flux thresholds from the manufactured solution, linear-solver
tolerance, and refinement rather than imposing the same algebraic threshold on
the spectral/hp flow solution.

Complete when the new checks pass and relevant existing inflow, outflow,
Homogeneous1D, transform-consistency, and moving-frame boundary regressions pass
on Linux. Select named regressions after inspecting `ctest -N`; avoid running
unrelated heavy suites by default.

### Milestone 6 Run resolution and domain pilots

Begin with inexpensive Linux pilot runs. Refine the inlet representation,
spectral/hp mesh, homogeneous mode count, timestep, transverse periods, and
outlet buffer independently. Measure startup, the adjustment length, recording
duration needed for stable estimates, memory, and runtime. Confirm that the
synthetic repetition time does not affect the record.

A provisional refinement target is less than 2 percent change in component
variances and integrated energy over the comparison region, with overlapping
sampling intervals. Examine spectral shape and dissipation separately so
agreement in an integral cannot hide an unresolved tail. Freeze the production
resolution and uncertainty criteria after these pilots, before comparing or
tuning downstream predictions.

Complete when a feasible production configuration and a documented error budget
exist. If DNS is impractical, report the resolved band and modeled dissipation
and define the result as a limited validation rather than claiming full DNS.

### Milestone 7 Compare downstream decay with experiment

Run the experimental-spectrum case first, then the von Karman alternative using
the same upstream constraints. Sample the reference plane and downstream
stations, compare Table 2 spectra and Table 4 quantities, and evaluate the
full-band convection-delayed correlation. Use more than one independent seed
or demonstrate convergence with independent statistical blocks.

Report discrepancies with numerical and sampling uncertainty. The experiment's
residual anisotropy, finite frequency response, inferred integral scales, and
unmeasured spectrum tails constrain the comparisons. Do not invent a universal
experimental error bar or require exact agreement with all tabulated quantities.
Keep downstream data separate from inlet fitting. If the adjustment region
persists, execute the development-region alternative before a quantitative
physical-validation claim.

Complete when comparisons, uncertainties, represented scales, energy budgets,
and remaining discrepancies are saved reproducibly. An unfavorable comparison
is a result to explain, not a reason to prescribe the desired decay.

### Milestone 8 Deliver documentation and reproducible examples

Provide a small runnable example for each supported expansion, configuration
reference, source data preparation, analysis commands, and a report containing
the tested revision and all solver, MPI, thread, spectrum, and stabilization
settings. Add narrow-band correlation analysis after the primary validation if
the full-band results justify it. Keep large generated data outside the source
checkout and record where it can be found.

Complete when another user can reproduce the small regressions and regenerate
the experimental comparison from the documented Linux environment.

## Linux execution and storage

Make source and Git changes in `/Users/jay/Documents/nektar-v2026`. Use the
Linux checkout `/media/jay/新加卷/nektar-v2026` for compilation, solver execution,
and numerical or performance conclusions, connecting with `ssh linux-usb`.

Before synchronization, inspect `git status`, branch, and `git rev-parse HEAD`
on both machines. Commit intended changes on the project branch, fetch that
branch on Linux, and advance its checkout only when doing so preserves existing
work. If a checkout contains unrelated changes or a conflicting active job, use
an isolated checkout or worktree and a corresponding build; do not reset, clean,
or overwrite that work. Verify revision and local changes again before testing.

Use these existing Linux locations when the checkout and build are compatible:

```text
Source       /media/jay/新加卷/nektar-v2026
Build        /home/jay/build/nektar-v2026-mpi
Installation /home/jay/.local/opt/nektar-v2026-mpi
Environment  /home/jay/.local/opt/nektar-v2026-mpi/env.sh
```

Store each heavier experiment under a separate directory below
`/media/jay/新加卷/nektar-v2026-runs/hit-inflow/`. Put large meshes, field files,
checkpoints, analysis arrays, logs, and temporary files there. Keep small
regression fixtures and analysis scripts in the repository. Store a per-run
manifest containing the revision, tracked diff if any, configuration hashes,
seed, dependency/build information, command, rank decomposition, thread counts,
output locations, and start/end times.

Never exceed 72 total Linux workers, counting simultaneous builds, test
processes, MPI ranks, and threads in every rank. Inspect other active work before
allocating workers. Default to `OMP_NUM_THREADS=1` and
`OPENBLAS_NUM_THREADS=1`; also limit other threaded numerical libraries if used.
Avoid automatic CPU-count parallelism. Run build, tests, and production jobs
sequentially unless their combined allocation has been checked explicitly.

The following are execution templates for later milestones, not commands run by
saving this plan. Confirm test names and case files after implementation.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Use 32 build workers only when at least that allocation is available.
cmake --build /home/jay/build/nektar-v2026-mpi --parallel 32
cmake --install /home/jay/build/nektar-v2026-mpi

# Inspect available tests before selecting the new and related regressions.
ctest --test-dir /home/jay/build/nektar-v2026-mpi -N -R 'HITInflow|KovaFlow_3DH1D'
ctest --test-dir /home/jay/build/nektar-v2026-mpi \
    --output-on-failure --parallel 1 -R 'HITInflow'

# Example pilot allocation: 16 MPI ranks, one thread each.
# Final decomposition and files depend on the prepared case.
mkdir -p '/media/jay/新加卷/nektar-v2026-runs/hit-inflow/pilot-01/tmp'
cd '/media/jay/新加卷/nektar-v2026-runs/hit-inflow/pilot-01'
export TMPDIR="$PWD/tmp"
mpirun -np 16 IncNavierStokesSolver mesh.xml session.xml --npz 4 \
    > solver.log 2>&1
```

Verify that the configured build enables the required solver, MPI, FFTW when
selected, and `NEKTAR_BUILD_TESTS` and `NEKTAR_BUILD_UNIT_TESTS` when their checks
are needed. Do not change a shared build configuration during another job.
Monitor process state, timestep progress, memory, residuals, CFL, mass balance,
and output growth. Set explicit limits appropriate to each pilot and preserve
diagnostics on failure. Numerical conclusions must identify the actual tested
revision and Linux allocation.

## Deliverables and Git workflow

The branch `codex/hit-inflow-generator` contains this plan first. Subsequent
implementation commits should follow the milestones in reviewable units:

1. Spectrum providers and deterministic transverse synthesis with core checks.
2. Full 3D boundary integration and pressure update ordering with a manufactured
   inlet regression.
3. Homogeneous1D transforms, MPI behavior, and restart support with regressions.
4. User documentation, small examples, and experimental data preparation.
5. Benchmark analysis and a report of Linux validation results.

Stage only intended files. Existing untracked literature copies, `AGENTS.md`,
and macOS metadata are separate from the plan commit. Do not add the supplied
PDFs or generated production data merely because they are present locally.
Keep the Linux working tree and unrelated user branches intact during sync.

The Mac supports direct GitHub SSH authentication. Push this branch with
`git@github.com:jiaweiwancn/nektar-v2026.git`; this checkout uses that push URL
while retaining its existing fetch URL. There is no need to relay pushes through
the Linux server.

The final deliverable is the implementation with passing focused checks and a
reproducible downstream comparison, including any documented physical or
resolution limitations. Saving this plan and pushing its commit establishes the
project branch; it does not mark the numerical milestones complete.

## Source references

1. Mann, J. (1998). Wind field simulation. *Probabilistic Engineering Mechanics*,
   13, 269-282. [DTU bibliographic record](https://orbit.dtu.dk/en/publications/wind-field-simulation/).
   Local copies: `docs/literature/Wind field simulation.md` and the corresponding
   PDF. Relevant passages: Section 3, Eqs. (8)-(11); Section 6, Eqs. (39), (44)-
   (50); Section 6.1 and Figures 13-14; Appendix A.
2. Comte-Bellot, G., and Corrsin, S. (1971). Simple Eulerian time correlation of
   full- and narrow-band velocity signals in grid-generated, isotropic
   turbulence. *Journal of Fluid Mechanics*, 48, 273-337.
   [Publisher record](https://doi.org/10.1017/S0022112071001599).
   Local copies: `docs/literature/Simple Eulerian time correlation of full-and narrow-band velocity signals in grid-generated isotropic turbulence.md`
   and the corresponding PDF. Relevant passages: Section 2; printed p.289
   footnote; Eqs. (38)-(53); Tables 2-7; Figures 8-16; Appendices D-E.
3. Repository guidance: `README.md` and the supplied root `AGENTS.md`.
4. Nektar++ documentation: `docs/user-guide/solvers/incompressible-ns.tex`,
   `docs/user-guide/xml/xml-conditions.tex`, `xml-forcing.tex`, `xml-filters.tex`,
   and `docs/user-guide/utilities/fieldconvert.tex`.
5. Solver sources: `solvers/IncNavierStokesSolver/BoundaryConditions/`,
   `EquationSystems/IncNavierStokes.cpp`, `VelocityCorrectionScheme.cpp`,
   `StandardExtrapolate.cpp`, and `Extrapolate.cpp`; forcing implementation
   `Forcing/ForcingIncNSSyntheticEddy.cpp`.
6. Expansion and transform sources:
   `library/MultiRegions/DisContField3DHomogeneous1D.cpp`,
   `ExpListHomogeneous1D.cpp`, `ExpListHomogeneous2D.cpp`, and
   `library/LibUtilities/FFT/NekFFTW.cpp`.
7. Diagnostic sources: `library/SolverUtils/Filters/FilterHistoryPoints.cpp`,
   `FilterModalEnergy.cpp`,
   `solvers/IncNavierStokesSolver/Filters/FilterReynoldsStresses.cpp`, and
   `library/FieldUtils/ProcessModules/ProcessPowerSpectrum.cpp`.
8. Existing solver examples under `solvers/IncNavierStokesSolver/Tests/`:
   `KovaFlow_m8_short_HOBC_3D1H.xml`, `KovaFlow_3DH1D_MVM_FFTW_Consistency.xml`,
   `KovaFlow_3DH2D_xy_FFT.xml`, and `ChanFlow3D_infTurb.xml`.
