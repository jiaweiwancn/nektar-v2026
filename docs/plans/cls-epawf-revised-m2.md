**M2 closure for the revised EP-AWF spectral wall profiles**

M2 is closed on 2 November 2026 on codex/cls-epawf-revised. The tested
code revision is acdacc022d64c72772b46f70a836d51d3f9f874e. The complete
external production-solver build and all 19 Linux CTests pass. Independent
profile/load verification passes; general inner-k positivity fails as
specified by the revised note.

**Material Passport**

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: implementation and native computational verification.
- Origin Date: 2026-11-02.
- Verification Status: VERIFIED for M2 profiles, budgets and dual loads;
  physical inner-k positivity FAILED AS SPECIFIED.
- Version Label: cls-epawf-revised-m2-v1.

The governing
[revised note](/Users/jay/Documents/nektar-v2026/docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex)
and [M0 contract](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-contract-v1.json)
retain their frozen SHA-256 values, respectively
d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f
and 4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a.
The CLS and original-AWF kernels, coefficients and historical solver
selections remain unchanged. The contract's historical availability entries
are preserved; this report records delivery of its six M2 verification IDs.

**Implemented operator**

[EquilibriumAWFProfiles](/Users/jay/Documents/nektar-v2026/library/SolverUtils/Turbulence/EquilibriumAWFProfiles.h)
is an immutable local profile and physical ray-load operator. Construct it
with the current positive sampled K and reconstruct it when those inputs
change. Its physical geometry has independent sampleHeight=r_K and
tailStart=r_m, both allowing the endpoint H. SampleK evaluates a supplied
field at r_K. Neither distance is derived from an FV interpolation point.
The original helper's unused historical query point is internal at H/2.

The operator implements

    p_A0 = nu_t,A * (t0 / D_A)^2       [zero forcing, current sampled K]
    E_c  = A / max(r, r_c)
    E_EP = (1 - chi) * E_A + chi * E_c
    C_k  = chi * (E_c - p_A0)
    S_k  = P_A - (1 - chi) * E_A - chi * p_A0

It assembles S_k in the simplified form once. Evaluate supplies the original
parallel-wall production; AtProduction accepts the caller's actual finite
signed production, including future full-strain patch contractions. Current
forced production remains distinct from p_A0. At the flat reference with
chi=1, both production terms use the same arithmetic path and cancel
pointwise before integration.

The local metric is J=(1+kappa1*r)(1+kappa2*r), with the frozen
B=grad_Gamma(n) convention and positive normal stretches. The cutoff solves
the J-weighted budget using a bracket in log(r_c), retaining 0<r_c<=r_m.
The minimum budget at r_m supplies the admissibility test. Strict requests
reject inadmissibility; the explicit original-AWF fallback records its
reason and zero correction weight. Inactive and low-Re branches do not
solve a cutoff. An unrepresentable cutoff/cap is rejected or recorded as
a fallback; cutoffs and dissipation are not clipped.

Physical Gauss integration splits at roughness, the positive virtual
viscosity intercept, original dissipation cutoff, tail start and corrected
cutoff. Geometric subdivisions resolve the molecular length and reciprocal
tails. Cutoff balance, signed redistribution and loads use that same measure.
RayLoads returns the dual vector B_Q^T W_J s_Q directly. RedistributionFlux
implements chi/J times the integral of J(E_c-p_A0) from zero to r. Its upper
endpoint retains the measured budget residual rather than imposing zero
by assignment.

The scalar C2 weight helper combines the prescribed spatial weight,
seven independently weighted squared indicators, cutoff-margin gate and
high-Re gate. The result is one scalar per ray. Physical sensor extraction,
its derivatives and transition geometry remain M6 responsibilities.

**Verification coverage**

[RANSEquilibriumProfilesProbe.cpp](/Users/jay/Documents/nektar-v2026/library/UnitTests/SolverUtils/RANSEquilibriumProfilesProbe.cpp)
uses actual Nektar StdSegExp bases at degrees 2, 4, 6 and 8. It evaluates
every active basis row and its derivative at the physical split quadrature.
A StdQuadExp degree-4 control checks all 25 tensor-product rows on one
wall ray. Native polynomial interpolation checks the fixed physical K
sampler through p refinement, including r_K=H.

[validate_cls_epawf_m2.py](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m2.py)
uses separate dimensional production/capped-dissipation primitives for
J=1+s1*r+s2*r^2, an independent Brent cutoff root and analytic fluxes.
Its 64-point NumPy Gauss rule uses factor-three physical panels; native
integration uses 16/32-point rules and factor-two panels. SciPy Jacobi
polynomials independently evaluate native Modified_A basis values and
derivatives. Native Ortho_A normalization sqrt((2*j+1)/2) is explicitly
removed when comparing the note's ordinary Legendre moments.

The matrix contains ten CLS fixture families at chi=0, 0.37 and 1, each
with 16/32-point quadrature, plus two standard-model controls: 62 cases.
Families cover the note example, smooth/partially rough profiles, one/two
positive curvatures, mixed curvature signs, smaller/larger current K and
independent r_K/r_m endpoints. The run checks 1,527 native basis rows:
1,488 segment rows, 25 quadrilateral ray rows and 14 orthonormal controls.
The independent driver records 21,202 checks; the native JSON records
71,413 checks and 15 rejected inputs, followed by a final output check.
Separate floating-point
range controls retain a representable subnormal cutoff and reject an
overflowing cap; these augment the 62-case matrix.

| Frozen verification ID | Delivered checks |
| --- | --- |
| M2.scalar_ep_provenance | Unchanged original profiles and forcing response; physical r_K/r_m independent of legacy P/N; fixed-K endpoint controls |
| M2.hp_source_and_moments | Every active segment mode, pointwise reference cancellation, signed supplied production, mean-only nonzero higher moments, source/flux integration by parts, scalar/vector forcing increments for all three chi values |
| M2.physical_weak_load | Positive physical weights/Jacobians, integrated physical volume, native value/derivative samples and independent physical weak integrals |
| M2.metric_cutoff_and_redistribution | Independent analytic metric budgets/root/flux, monotone bracket checks, flat-cutoff rejection on curved data, zero endpoint fluxes, order-16/32 refinement, inactive inadmissible and low-Re branches |
| M2.dual_source_loads | Direct dual-vector values for all active native rows, including 25 quadrilateral ray rows; no second integration or mass conversion |
| M2.sampled_k_positivity_failure | Nonnegative degree-five polynomial with an analytic quartic zero away from r_K, zero ground derivative, native value/derivatives, negative local RHS and local Euler increment |

The general comparison target is 1e-10 after dimensional scaling. The
note example also checks dimensional redistribution identities at 2e-12
and scalar forcing increments at 3e-12. A deliberately corrupted native
source load is rejected by the independent oracle. Pressure and convection
indicators with equal and opposite values cannot cancel in the switch.

**Measured results**

| Quantity | Linux result |
| --- | --- |
| Flat note cutoff | 0.009515269139914396 m |
| Curved note cutoff, kappa1=0.05/m | 0.009458614008737381 m |
| Two-curvature cutoff, kappa1=0.05/m, kappa2=0.08/m | 0.009363561189719717 m |
| Flat note production / minimum budgets | 0.9221852179962751 / 0.24596631258988558 m^3/s^3 |
| Maximum scaled native dual-load versus independent integral error | 3.213e-16 |
| Maximum scaled source/flux integration-by-parts mismatch | 2.258e-17 |
| Maximum scaled metric budget residual | 1.372e-16 |
| Maximum scaled endpoint-flux magnitude | 1.351e-16 |
| Maximum scaled order-16/32 load difference | 1.214e-16 |
| Maximum relative order-16/32 cutoff difference | 0 |
| Maximum note dimensional redistribution / forcing mismatch | 1.909e-17 / 5.656e-17 m^3/s^3 |
| Reference source loads at chi=1 | Exactly zero in every native tested row |
| Representable subnormal cutoff / finite cap control | 1.38929295489488e-309 m / 1.045651790084218e308 m^2/s^3 |
| Overflowing-cap control | Strict rejection; explicit original fallback reason UnrepresentableCutoff |
| Complete Linux regression suite | 19/19 passed |

The mean-only ordinary Legendre control at degrees 0 through 6 gives

    approximately 0,  0.0173551439079370, -0.0220336637626453,
                     0.0220527999939071, -0.0197059614434825,
                     0.0160477076949066, -0.0116765463335830

in m^3/s^3. These match the note's rounded entries within 4.531e-14 in
the stated units. Mean compatibility leaves nonzero higher modes; the
new simplified source cancels all modes at the reference.

For the positivity control, set chi=0 on the smooth-wall ray and use

    k(r) = K * ((r-r0)/(r_K-r0))^4 * (1+4*r/r0)/(1+4*r_K/r0).

This polynomial is nonnegative on [0,H], has k'(0)=0, and has zero value,
first derivative and second derivative at r0=z_v/3. It is represented in
a native degree-6 segment. The measured r0 is 7.4877700787708e-5 m,
r_K=1 m and sampled K=0.5105080130598859 m^2/s^2. At r0, nu_t,A=P_A=0
and the positive original sink is 1335.9909337349875 m^2/s^3. The local
RHS is its negative; dt=1e-6 s gives a local Euler value
-0.0013359909337349878 m^2/s^2.

The analytic zero and native interpolation roundoff are reported separately:
native k(r0)=-3.698e-19 m^2/s^2 and native second derivative=-1.417e-13
in the corresponding units. The negative sink/RHS persists independently
of that roundoff. No source repair, floor or limiter is introduced. This
is a local continuum/source counterexample with native polynomial
evaluation, not a native coupled solver timestep; that map belongs to M4.

**Execution and evidence**

The authoritative run is
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m2-attempt-005
on jay-Super-Server. It used the isolated external Release build compiling
the complete production incompressible-solver source list, native RANS
kernels and source geometry library against the installed MPI libraries.
Both production and external CMake files register the new sources and
RANSEquilibriumProfilesNative / RANSCLSEPAWFM2Oracle tests.

The framework remains 5.9.0 with GCC 11.4.0, CMake 3.22.1, Open MPI 4.1.2,
Python 3.10.12, NumPy 1.21.5 and SciPy 1.8.0. Compilation used 16 workers;
tests ran sequentially pinned to CPU 0. MPI-enabled solver tests used one
MPI process and all numerical thread counts were one. Point-profile probes
do not launch MPI. No MPI partition-equivalence result is claimed.

The probe/kernel rebuild took 4.67 s, the complete target build completed,
the direct oracle took 1.77 s and the full CTest suite took 109.18 s.
The direct oracle, CTest oracle and standalone native test produced
byte-identical profile JSON. Retained M1 native
kernel/constitutive reports are also byte-identical after relinking.
All 37 resolved shared libraries match M1. An explicit geometry-library
LD_LIBRARY_PATH resolves the known external-build runtime-path warning;
ldd found no missing libraries.

The detached Linux source was clean at the tested revision. The unrelated
main checkout remained clean at f4b70baec10016dd098005bc5a5507b87e4a57bf,
and its build/installation were preserved. Source transfer used incremental
Git bundles. Existing Mac untracked files were preserved.

Compact evidence accompanies this report:

- [verification.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m2/verification.json):
  case summaries, errors, controls, physical gate, source/binary and raw-data hashes.
- [execution-receipt.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m2/execution-receipt.json):
  exact tested revision, build/oracle/suite commands, timings, exit codes and log hashes.
- [integration.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m2/integration.json):
  eleven executable hashes, static kernel archive, 37 resolved shared-library hashes,
  compiler flags, complete test inventory and retained M1 identity checks.
- [attempts.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m2/attempts.json):
  all five development attempts, including the initial native size-conversion
  compile failure and the oracle's initial Ortho_A normalization mismatch.

The initial compile issue was fixed with explicit size types. The moment
comparison was corrected to distinguish orthonormal native modes from the
note's Legendre convention. A final range audit replaced the smallest-normal
search bound with the actual representable cutoff/cap bound; the two new
range controls verify the correction. None of these fixes changed the profile law. All raw
reports/logs, earlier attempts, binaries and generated smoke fields remain
on Linux outside the checkout. The following closure commit changes only
documentation/evidence; the code evidence belongs to acdacc0.

**Reproduction and remaining milestones**

Use a clean Linux worktree at acdacc0 and a fresh direct-oracle directory.
The oracle refuses to overwrite a previous run. NumPy/SciPy must be present
at configure time; confirm that both revised-note oracle tests appear.

    source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
    CLS_EP_RUN=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised
    CLS_EP_SOURCE=$CLS_EP_RUN/source
    export TMPDIR=$CLS_EP_RUN/tmp
    cmake -S "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS" \
      -B "$CLS_EP_RUN/build" -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
      -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
    cmake --build "$CLS_EP_RUN/build" --parallel 16
    export LD_LIBRARY_PATH="$CLS_EP_RUN/build/geometry-lib:$LD_LIBRARY_PATH"
    taskset -c 0 python3 "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m2.py" \
      --probe "$CLS_EP_RUN/build/RANSEquilibriumProfilesProbe" \
      --output "$CLS_EP_RUN/m2-replay/verification"
    taskset -c 0 ctest --test-dir "$CLS_EP_RUN/build" --parallel 1 --output-on-failure

This closes the M2 local spectral-profile/load scope. Tangential patch
integration, global multiregion fields, pressure/interface assembly and
coupled full-strain production are M3 work. Actual native time-stage
preservation is M4; curved mesh/tube and flux qualification is M5. Local
metric budgets and tensor-product ray rows do not establish those gates.
The original CLS adverse-diffusion finding from M1 and the inner-k
positivity failure reproduced here remain failed physical properties.
Robust-model and terrain-accuracy qualification remain conditional M8 work.
