# HIT inflow Milestone 1 interface and benchmark specification

This specification closes the interface and input decisions in
[Milestone 1 of the project plan](hit-inflow-implementation-plan.md). It defines
the version 1 contract for implementation in Milestones 2-4. The XML examples
are configuration fixtures, not sessions supported by the current solver.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | plan |
| Origin date | 2026-10-28 |
| Version label | hit_inflow_interface_v1 |
| Verification status | Interface and numeric fixtures checked; solver implementation and numerical behavior remain unverified |
| Audited code revision | `132f4582e438ead0d81f45da671bcd3b813e56b1` |
| Branch | `codex/hit-inflow-generator` |

## Version 1 support contract

| Setting | Initial supported value |
| --- | --- |
| Solver | `SolverType=VelocityCorrectionScheme` |
| Equations | `EQTYPE=UnsteadyNavierStokes`, `EvolutionOperator=Nonlinear` |
| Projection | `Continuous` |
| Variables | Exactly `u,v,w,p`, with pressure last |
| Pressure extrapolation | `Extrapolation=Standard`, including its existing startup policy |
| Integration | `METHOD=IMEX`, empty variant, order 2; order-1 startup supplied by the existing integrator |
| Time step | Fixed positive `TimeStep`; positive integer `NumSteps`; `FinTime=0` |
| Expansions | Full 3D; multiple-mode `Homogeneous1D` with Fourier z |
| Inlet | Stationary planar surface at `x=ORIGIN.X`, fluid on its positive-x side |
| Mean velocity | Constant positive u, zero v and w |
| Inlet velocity BC | Three Dirichlet conditions with the exact `USERDEFINEDTYPE=HITInflow` |
| Inlet pressure BC | Neumann, `USERDEFINEDTYPE=H`, `VALUE=0` |
| Flowrate controller | Disabled, `Flowrate=0` |
| Domain for benchmark | Periodic y and z; open x direction with a downstream buffer |

The configuration parser must reject unsupported combinations before allocating
modes. These include `SingleMode`, `HalfMode`, `Homogeneous2D`, pure 2D flow,
extra transported variables, weak-pressure or implicit VCS, coupled solvers,
substepping, IMEX DIRK, orders 3/4, moving frames, mapping/ALE, adaptive steps,
and other forcing that changes the prescribed mean velocity. A single three-
component inlet region is supported initially; reject multiple `HITInflow`
regions until communicator and sharing tests explicitly add that capability.

Require `CFL=CFLEnd=IO_CheckTime=0` and
`CFLGrowth=TimeIncrementFactor=1` (their defaults when omitted). Use
`IO_CheckSteps` for step-based checkpoints. Reject `FinTime>0` in this first
contract, since finishing at a final time or an output time can shorten a step.
These changes to dt occur in `UnsteadySystem.cpp:289-335` and must not bypass
the fixed-step restriction.

IMEX's name alone is insufficient to select this contract. Check its method,
variant, and order. Orders 3/4 use multistage DIRK startup whereas order 2 starts
with IMEX1. Their boundary timing must be audited separately before support is
extended.

## Configuration grammar

Use exactly one top-level `NEKTAR/HITINFLOW` element alongside `CONDITIONS`.
Tags and attributes below are case-sensitive and uppercase. Reject unknown,
duplicate, or missing elements and attributes. Require `VERSION="1"` and
`UNITS="SI"`. Each custom `REGION REF` is the global boundary-region ID from
`CONDITIONS/BOUNDARYREGIONS`, not the index in a rank's boundary array.

The single region contains one each of `ORIGIN`, `SYNTHESIS`, `SPECTRUM`, and
`CUTOFF`, in that order in the examples. Parsing need not depend on their order.
The examples deliberately use sparse region ID 7, avoiding an accidental
identification of REF with local array index 0.

Ordinary real attributes may contain Nektar parameter expressions, evaluated
once with the session interpreter after `CONDITIONS/PARAMETERS` is loaded.
Preserve the declared spelling of parameter names. Permit constants, parameters,
and pure arithmetic/functions; reject dependence on x, y, z, or t. Check finite
results. For integer grid sizes require exact integral values in the supported
integer range; do not silently round. `SEED`, `VERSION`, and `REF` are strict
unsigned decimal integer text, not floating-point expressions. In particular,
parse `SEED` directly into `uint64_t` without a round trip through `NekDouble`.
Restrict region `REF` and boundary IDs to `0..INT_MAX`, matching Nektar's signed
integer boundary parser rather than the wider seed range.

Mean velocity comes exclusively from the three velocity `VALUE` expressions.
Evaluate them once; require the same constant-parameter restriction and verify
`u>0`, `v=w=0`. Use that u value as the convection speed. There is no independent
speed setting that can silently disagree with the prescribed mean. Each update
overwrites boundary values with mean plus fluctuations; it never adds noise to
values left over from the previous call.

### Origin and synthesis attributes

| Element | Attribute | Meaning and validation |
| --- | --- | --- |
| `ORIGIN` | `X,Y,Z` | Required finite phase origin in metres; X must agree with all inlet points |
| `SYNTHESIS` | `LX,LY,LZ` | Required positive box periods in metres |
| `SYNTHESIS` | `NX,NY,NZ` | Required even grid counts, at least 4, no larger than the implementation's signed indexing limit |
| `SYNTHESIS` | `SEED` | Required decimal integer from 0 through `2^64-1` |
| `SYNTHESIS` | `TIMEORIGIN` | Required finite phase-origin time in seconds, unchanged on restart |
| `SYNTHESIS` | `RECORDDURATION` | Required positive declared duration from TIMEORIGIN, covering the entire intended record including restarts |
| `SYNTHESIS` | `REPEAT` | Required `Error` or `Allow` |
| `CUTOFF` | `KMAX` | Required positive radial cutoff in rad/m |
| `CUTOFF` | `TAPERSTART` | Required number in `[0,KMAX)`; starts a cosine taper of the energy spectrum |

For `REPEAT=Error`, require `RECORDDURATION < LX/Ubar`. `Allow` explicitly permits
periodic reuse and must be reported in the run manifest. Check job start and
end times against the declared record interval with a floating-point time
tolerance. Extending the declared interval on restart is allowed only if it
still meets the period policy; it does not change the realization.

Require `KMAX < min(pi*NX/LX, pi*NY/LY, pi*NZ/LZ)` and
`Ubar*KMAX*TimeStep < pi`. These bounds prevent obvious aliasing; they do not
establish physical resolution or replace refinement. In `Homogeneous1D`, require
the same `LZ` as the homogeneous expansion and additionally
`KMAX < pi*HomModesZ/LZ`, with even `HomModesZ >= 4`. Synthesis `NZ` and solver
`HomModesZ` are independent: the former defines the immutable random lattice;
the latter samples and represents its retained band. Increasing the solver
plane count must not change the realization. This Milestone 4 clarification
replaces the original `NZ=HomModesZ` restriction and allows Fourier refinement
of one fixed inflow. Require the y and z computational periods to match the
synthetic box periods. The full mesh and boundary topology checks are
implemented and tested in Milestones 3-4.

Reject grid products that overflow allocation arithmetic or the RNG counter
range defined below. Log the estimated mode memory before allocation.

### Spectrum attributes

`SPECTRUM TYPE="Tabulated"` requires `FILE`, `INTERPOLATION`, `LOW`, `HIGH`, and
`NORMALIZATION`. `INTERPOLATION` is `LogLinear` or `Linear`. `LOW` is `K4` or
`Zero`; `HIGH` is `Gaussian` or `Zero`. There are no implicit tail choices.
`FILE` is interpreted relative to the process working directory, consistent
with existing file-based boundary inputs, not relative to the XML file.

`SPECTRUM TYPE="VonKarman"` requires positive `LENGTH` in metres and
`NORMALIZATION`. Use the project plan's von Karman formula and name its
amplitude `AMPLITUDE = alpha*epsilon^(2/3)` in SI units. Its value is not a
promise of realized dissipation after cutoff.

`NORMALIZATION` has exactly two values:

- `Preserve`: preserve the supplied spectral amplitude before the explicit
  cutoff. For `VonKarman`, require nonnegative `AMPLITUDE` and forbid `RMS`. For
  `Tabulated`, forbid both `AMPLITUDE` and `RMS`.
- `TargetRMS`: require nonnegative `RMS` in m/s and forbid `AMPLITUDE`. Scale the
  energy spectrum by one common factor so its finite-lattice ensemble energy
  after tapering is `3*RMS^2/2`. For `VonKarman`, start with `AMPLITUDE=1` in SI
  before this normalization. For `Tabulated`, apply the same scaling to the
  supplied spectrum. This option changes its measured amplitude and must be
  explicit in benchmark reports.

Here RMS is the isotropic-equivalent `sqrt(2*expected(K)/3)`. A finite,
anisotropically sampled lattice need not give every component exactly that RMS;
report its component covariance. `RMS=0` or preserved `AMPLITUDE=0` explicitly
selects the mean-only limit. Validate the remaining configuration, then skip
Gaussian draws and mode allocation and return the prescribed mean.

Do not normalize each realization's measured RMS. Such normalization would
alter the Gaussian coefficient distribution. Do not normalize components
independently. If a positive target RMS has zero represented energy, report an
input error instead of dividing by zero.

For interval endpoints `ki,Ei` and `kj,Ej`, define the interpolation fraction
`a=(k-ki)/(kj-ki)` for `Linear`, giving `(1-a)*Ei+a*Ej`. For `LogLinear`, use
`a=log(k/ki)/log(kj/ki)` and `exp((1-a)*log(Ei)+a*log(Ej))`. Require all energies
strictly positive for `LogLinear`; `Linear` permits zero energies.

Below the first table point, `K4` means `E(k)=E1*(k/k1)^4`; `Zero` means zero.
Above the last point, `Gaussian` means

\[
\beta=-\frac{\log(E_n/E_{n-1})}{k_n^2-k_{n-1}^2},\qquad
E(k)=E_n\exp[-\beta(k^2-k_n^2)].
\]

Require the last two energies to be positive and strictly decreasing for that
tail. `Zero` means zero beyond the last point. Return tabulated values exactly
at the endpoints. These extensions are continuous where their formulas imply
it; derivative continuity is not promised. The Gaussian is a declared modeling
choice, not an experimentally measured tail. Later sensitivity studies can
change the choice without rewriting the measurement files.

After spectrum evaluation, multiply its energy by

\[
W(k)=\begin{cases}
1,&k\leq k_t,\\
\tfrac12[1+\cos(\pi(k-k_t)/(k_{\max}-k_t))],&k_t<k<k_{\max},\\
0,&k\geq k_{\max},
\end{cases}
\]

where `kt=TAPERSTART`. Apply normalization next and then the square root when
forming velocity coefficients. Multiplying the velocity amplitude by W instead
would apply W squared to the energy and violate this contract.

### Spectrum file format

The required first nonblank line for an inlet input is

```text
# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2
```

Subsequent `#` comment lines and blank lines are allowed. Every data line has
exactly two whitespace-separated finite decimal numbers: positive increasing
angular k, then nonnegative E. At least two rows are required. Reject duplicate
or decreasing k, NaN/Inf, extra columns, unsupported units, or another KIND.
In particular, `KIND=Longitudinal1D` is a validation spectrum and cannot be
inserted into Batchelor's shell tensor. Files do not contain missing-value
sentinels; individual station files omit absent measurements.

Read and validate the file once on session communicator rank 0 and broadcast
the normalized numeric table and settings. All ranks participate, including
ranks owning no inlet elements. Diagnose read failures collectively before any
dependent collective to avoid hanging other ranks. Record a SHA-256 digest of
the actual file bytes and the resolved numerical settings in checkpoint/run
metadata. Later milestones must supply the runtime digest implementation and
known-answer checks; this milestone's preparation tool already records SHA-256.

### Concrete example files

The two checked configuration fragments are
`solvers/IncNavierStokesSolver/Tests/HITInflow/hit-inflow-tabulated.xml` and
`hit-inflow-vonkarman.xml`. They exercise `Tabulated/Preserve/Homogeneous1D` and
`VonKarman/TargetRMS/full3D`, respectively. Their short record and coarse cutoff
keep their constraints understandable; they are not production benchmark
parameters or claims of resolved turbulence. Neither is registered as a solver
test or runnable until a mesh and complete boundaries are provided.

## Fourier realization and RNG contract

Use the project plan's complex Fourier series, including both members of every
conjugate pair in its mathematical sum. For each coordinate, retain integer
indices `mi` from `-Ni/2+1` through `Ni/2-1`, inclusive, excluding all Nyquist
planes. Remove the zero vector and modes with `k>=KMAX`. The independent member
of a pair is the one whose first nonzero entry in `(mx,my,mz)` is positive.
Generate its vector coefficient, then set the negative partner to its complex
conjugate; do not redraw its Gaussian variables.

Use algorithm identifier `SplitMix64BoxMuller52V1`. The definition below is
normative for future code; no stochastic generator is implemented in this
milestone. Unsigned arithmetic is modulo `2^64`.

1. For the canonical independent member, define
   `p=((mx+NX/2)*NY+(my+NY/2))*NZ+(mz+NZ/2)`.
2. For q from 0 through 5, set `x=SEED+6*p+q`, modulo `2^64`.
3. Define `z=x+0x9e3779b97f4a7c15`, then
   `z=(z xor (z>>30))*0xbf58476d1ce4e5b9`, then
   `z=(z xor (z>>27))*0x94d049bb133111eb`, then
   `h=z xor (z>>31)`.
4. Convert to the open-interval uniform
   `Uq=((h>>12)+0.5)*2^-52`. Taking 52 bits with a half-bin offset avoids both
   zero and a rounded value of one.
5. For j from 0 through 2, use `r=sqrt(-2*log(U[2*j]))` and
   `phi=2*pi*U[2*j+1]`, then
   `xi[j]=(r*cos(phi)+i*r*sin(phi))/sqrt(2)`.
6. Apply the transverse cross-product factor from the project plan. Use
   `DeltaVk=(2*pi)^3/(LX*LY*LZ)`, with no extra factor for half-lattice storage.

Require `6*(NX*NY*NZ-1)+5 <= 2^64-1`, as well as normal allocation bounds.
Canonical mode IDs depend on the grid, not rank number, local point order,
boundary numbering, or traversal order. Changing grid sizes changes the
realization identity. Integer words and uniform inputs are portable; libm
transcendentals may differ slightly between platforms. Require exact integer
known-answer tests and floating-point tolerance tests, not a claim of bitwise
Gaussian equality across different libm implementations.

Use `phase = k dot (x-ORIGIN) - kx*Ubar*(t-TIMEORIGIN)` and evaluate from absolute
time. No rank owns a random stream or advances RNG state on a time callback.
For full 3D and Homogeneous1D, the same coordinates and configuration must
represent the same continuous realization before projection. FFT packing is a
representation adapter; it must not change these coefficients or their
normalization.

For `TargetRMS`, calculate the expected energy over the full retained lattice
as `K0 = sum(DeltaVk*Eeffective(k)/(4*pi*k^2))`. Equivalently, sum independent
members and multiply by two. For a positive target RMS, scale all energies by
`(3*RMS^2/2)/K0`. The velocity coefficients remain Gaussian, while realized energy
is allowed to fluctuate between seeds; only its ensemble expectation is fixed.

## Startup and restart contract

Initialize the boundary manager and all HIT settings before its first prescribed
velocity callback. At the solver's initial or restored time, generate `g(t)`
before initial Dirichlet conditions are imposed. An interior mean-only initial
condition is permitted and is used in the configuration fragments; its
transients must flush before statistics are sampled. There is no automatic
amplitude ramp, random reseeding, or resetting of phase on restart.

Generator restart metadata must contain the schema/algorithm versions, seed,
resolved grid and periods, origin, mean velocity, phase time origin, spectrum
kind/parameters, cutoff/taper, normalization rule and final scalar factor, and
input-file digest plus normalized table when tabulated. Verify these against
the current input and fail on a realization-changing mismatch. Record and check
the numerical environment for exact reconstruction; do not promise bitwise
portability across builds. Reconstruct from the same mode IDs and actual
restored solver time. `RECORDDURATION`, `REPEAT`, output frequencies, and MPI
decomposition do not define the random realization, although their independent
validity checks still apply.

The current solver initializes multistep history on startup. Reproducing the
generator boundary exactly is distinct from restoring all velocity/pressure
history and reproducing an uninterrupted IMEX2 solution bitwise. Milestones 3-5
must test the boundary history, startup convergence, and the accuracy of an
evolved restart separately. This inlet must not silently change existing
pressure-acceleration bootstrap behavior.

## Boundary update audit and minimal change

Add a separate method `UpdatePrescribedVelocity(NekDouble time)` to
`IncBaseCondition` and `IncBoundaryConditions`. The base virtual implementation
does nothing. Only the new inlet overrides it initially; leave its ordinary
post-pressure `v_Update(fields, Adv, params)` empty. Register `HITInflow` in the
factory and recognized-type set. Retain one region object that holds all three
velocity expansions and their original mean expressions.

Load and broadcast the global HIT configuration at the beginning of
`IncBoundaryConditions::Initialize` on all session ranks, before entering its
local boundary factory loop. Store the resulting shared descriptor in the
manager and pass it to local inlet objects. The local class constructor must
not initiate a session-wide collective: ranks without inlet faces never create
that object. Global input validation and local geometry validation are distinct
phases, with collective agreement on errors before stepping.

Invoke the new manager method at the end of
`IncNavierStokes::SetBoundaryConditions(time)`, after generic boundary evaluation
and existing velocity adjustments, guarded for manager availability. Do not move
all existing boundary classes to an earlier phase: `StaticWall` and other
classes also alter pressure. Do not add another inlet update in
`v_EvaluateAdvection_SetPressureBCs` with `rhsTime+TimeStep`.

The relevant call sequence is:

1. Initialization calls `SetBoundaryConditions(t0)` before
   `ImposeDirichletConditions`.
2. The GLM integrator evaluates the initial RHS at t0 using that trace.
3. `v_PreIntegrate` calls `SetBoundaryConditions(tn+dt)` before the implicit
   solve, supplying the required next velocity trace.
4. The one-stage IMEX RHS is evaluated at the completed stage time. Its pressure
   boundary calculation reads the trace already supplied by the common hook.
5. The existing general boundary manager update remains after pressure
   evaluation for existing classes. Its `params["Time"]=rhsTime+dt` is not an
   instruction to advance the new inlet again.

`USERDEFINEDTYPE=HITInflow` does not set `IsTimeDependent()`: the generic boundary
parser recognizes the literal `TimeDependent` for that purpose. The explicit
hook is therefore necessary; changing the flag could trigger generic expression
evaluation and overwrite the new trace.

Establish a mapping between global boundary-region IDs and the local expansion
array order. Do not pass the manager's local `n` as the XML REF without this
mapping. Ranks without inlet faces still need configuration, metadata, and
collective participation. Later tests must cover sparse IDs, no-inlet ranks,
and a partition that distributes the homogeneous planes.

For Homogeneous1D, use `GetCoords`, populate physical samples for all components,
call `SetWaveSpace(false)`, and then `FwdTransBndConstrained`. The latter already
performs the homogeneous forward transform. Preserve its wave-space conventions
for later coefficient reads. Use global plane identifiers and distributed
transforms; do not double-transform or update only the zero plane.

### Verified source anchors

Anchors refer to the audited revision; all paths are repository relative.

| Evidence | Source and line |
| --- | --- |
| Boundary manager constructed during extrapolation setup | `solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp:150` |
| Initial boundary update precedes initial Dirichlet enforcement | `solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp:645` |
| Current generic boundary update path | `solvers/IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp:323` |
| Next-step boundary update | `solvers/IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp:1022` |
| Current general manager update follows pressure evaluation | `solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp:782` |
| Pressure consumes prescribed normal-velocity history | `solvers/IncNavierStokesSolver/EquationSystems/StandardExtrapolate.cpp:69`, `Extrapolate.cpp:72`, `Extrapolate.cpp:545` |
| Initial acceleration bootstrap | `solvers/IncNavierStokesSolver/EquationSystems/StandardExtrapolate.cpp:150` |
| GLM initial explicit derivative and one-stage IMEX timing | `library/LibUtilities/TimeIntegration/TimeIntegrationAlgorithmGLM.cpp:79`, `:626` |
| Order-2 and order-3/4 startup differ | `library/LibUtilities/TimeIntegration/IMEXTimeIntegrationSchemes.h:89` |
| Boundary type allowlist and local index construction | `solvers/IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.cpp:41`, `:101` |
| Global boundary-region definitions | `library/SpatialDomains/Conditions.h:273`, `Conditions.cpp:343` |
| Generic time-dependent flag | `library/SpatialDomains/Conditions.cpp:520` |
| Uppercase custom-element lookup | `library/LibUtilities/BasicUtils/SessionReader.cpp:700` |
| Parameter expression evaluation | `library/LibUtilities/BasicUtils/SessionReader.cpp:1811` |
| Existing working-directory-relative boundary file open | `solvers/IncNavierStokesSolver/EquationSystems/IncNavierStokes.cpp:580` |
| Homogeneous1D physical Dirichlet projection | `library/MultiRegions/DisContField3DHomogeneous1D.cpp:623` |
| Projection includes homogeneous transform | `library/MultiRegions/ExpListHomogeneous1D.cpp:482` |

Milestone 3's manufactured inlet test must confirm numerical temporal accuracy.
The source audit establishes where to integrate the new hook; it does not
substitute for that test.

## Experimental inputs and provenance

The fixture directory
`solvers/IncNavierStokesSolver/Tests/HITInflow/` contains the checked source
transcription, six SI spectrum files, a three-station target file, a deterministic
manifest, a preparation tool, and a contract checker. Its README defines all
units and source pages. The local PDF identity is SHA-256
`e8067700a0ed678c33fa92e5e3f3d737df98d976e5f3b1478612f071536962d3`.

Table 3 at theta 42 is the inlet spectrum, with 19 measured points covering
20-2000 rad/m. The theta 98 and 171 shell spectra remain validation data, with
19 and 18 points. Table 2(a) provides 22 directly measured one-sided longitudinal
points per station. Table 4 provides gross statistical targets. Missing entries
are omitted per station; the source JSON retains them as null. The printed
nonmonotonic Table 3 sequence at theta 98 is preserved.

All files contain only tabulated data converted to SI. Runtime `K4` and
`Gaussian` policies in the example are explicit proposed extensions, not extra
experimental rows. Do not declare their unmeasured low/high tails validated.
The viscosity discrepancy and the finite Table 3 dissipation tail are retained
in the fixture provenance and parent plan.

Station mapping uses `M=0.0508 m`, approach `U0=10 m/s`, and measurement speed
`Ubar=12.7 m/s`. At reference theta 42, theta 98 and 171 correspond to elapsed
times 0.28448 and 0.65532 seconds and downstream distances 3.612896 and 8.322564
metres. The generated target file preserves Table 4's original quantities and
adds only those deterministic relative coordinates.

## Milestone closure evidence

Milestone 1 is complete when this specification, both checked configuration
fragments, original-unit source transcription, SI fixtures, preparation/check
tools, and validation record are committed together on the project branch.
The validation record records exact commands, input digests, checkout states,
and the single-worker Linux check. Runtime parser acceptance, stochastic mode
generation, compiler tests, solver runs, and physical validation belong to the
later milestones and remain open.

The Mac has direct GitHub SSH authentication. Use the SSH push URL
`git@github.com:jiaweiwancn/nektar-v2026.git` for this branch. Linux remains the
authoritative environment for compilation and numerical validation; its source
checkout need not change to validate these small preparation artifacts in an
identified run snapshot.
