**M1 closure for the revised CLS and scalar AWF formulation**

M1 is closed on 2 November 2026. Independent kernel verification passes;
the original CLS constitutive-stability gate fails as specified by the
revised note. The tested code revision is
e5ab5984f00dbbb84c26c0271d353fa98d99ff6a on codex/cls-epawf-revised,
based on 6415398d856323c3116e646eaf7c882c66617602.

**Material Passport**

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: native computational verification.
- Origin Date: 2026-11-02.
- Verification Status: VERIFIED for M1 identities and counterexamples;
  physical constitutive stability FAILED AS SPECIFIED.
- Version Label: cls-epawf-revised-m1-v1.

The governing
[revised note](/Users/jay/Documents/nektar-v2026/docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex)
is committed unchanged, with SHA-256
d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f.
The scope is M1 in the
[approved plan](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-implementation-plan.md).
The CLS/AWF production formulas and original coefficients are unchanged.

**Implemented verification**

[RANSConstitutiveProbe.cpp](/Users/jay/Documents/nektar-v2026/library/UnitTests/SolverUtils/RANSConstitutiveProbe.cpp)
samples the actual CLS kernel, including physical max-invariant evaluation,
the two explicitly selected limiting branches, reference profiles, general
3D gradients and 2D gradients with all three modeled stress components.
It exports stresses and nearby states rather than implementing a second
native constitutive law.

[validate_cls_epawf_m1.py](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m1.py)
compares those samples with a dimensional indexed tensor and differentiates
that independent tensor by complex step. It checks native finite differences
against those actions, with separate truncation/roundoff bounds. Both formal
branches at the invariant tie are distinguished. Smooth physical limiting
sequences approach them without crossing the invariant switch.

The driver also reuses the earlier independent dimensional-index and SciPy
quadrature oracle against RANSKernelProbe: 78 bulk cases, 60 CLS-reference
wall cases and 48 standard-reference wall controls, including configured
coefficients and invalid inputs. The new probe adds 17 constitutive cases.
The native kernel probe completes 19,192 assertions and rejects 30 invalid
inputs; the new driver completes 5,000 checks. Deliberately corrupted stress
and wall-traction outputs are both rejected by the oracle. These counts
describe the executed checks, rather than an all-state proof.

RANSCLSEPAWFM1Oracle and the native probe are registered in both production
and external RANS CMake files. Execution used the external Release build,
which compiles the complete production incompressible-solver source list
against the installed native libraries.

| Revised-note requirement and equation labels | Executed evidence |
| --- | --- |
| Tensor conventions and cubic stress: SO, cubic, coeffs, invariants, C | Dimensional indexed versus native matrix tensor; quadratic/cubic components, all stresses, frame covariance and production; general, rotational, strain, raw-divergence and 2D cases |
| Production and algebraic Yap: production-definition, Yap | Full stress contraction and both scalar sources; reference Yap activity; doubling geometric wall distance changes Yap while leaving stress unchanged |
| Neutral roots/profiles/stress: qroot, kappaeq, reference, Rref | Independent roots and five heights; k, velocity, epsilon, stress, production, epsilon diffusion/source balance and gradient flux |
| Published AWF: DA, brough, J01, Jclosed, J1closed, PA, epsA, epsAmean | Split quadrature for integrals, production, dissipation and reconstruction; smooth, partially/fully rough, h*=70, molecular, weak-turbulent, thin forcing-layer and changed-current-K cases |
| No slip, energy and forcing: UA, energy, pressureresponse, Ppoly | Ground/matching values; reference power balance, signed/reversed traction, finite wall/interface slopes and independent production polynomial |
| Counterexamples: stability-negative-diffusion, stability-transverse-block, stability-inplane-direction | Native transverse block and difference-step convergence, both limiting tangents, nearby smooth 2D adverse direction, covariance and production diagnostics |

**Measured results**

| Quantity | Linux result |
| --- | --- |
| q* | 3.1906750816242875 |
| Cmu* | 0.09822789532821127 |
| kappa_c | 0.4405539051866587 |
| Reference Yap ratio / function | 1.0043449900879187 / 0.003637748895392772 |
| Eigenvalues of reference R/K | 0.3482098773121975, 0.56, 1.0917901226878026 |
| Maximum independent indexed/quadrature scaled error | 3.676129406696578e-14, below the approved 1e-11 target |
| Maximum native production versus independent native quadrature relative error | 1.1940821049274815e-12 |
| Maximum reference wall-energy relative imbalance | 2.1738905623866739e-16 |
| Cases with signed negative wall-traction components | 15 |
| Example wall / matching-surface finite forcing slopes | -0.268331641404983 m / 1.631668358595017 m |
| Clean base suite at 6415398 | 16 of 16 tests passed |
| M1 suite at e5ab598 | 17 of 17 tests passed |

The production-polynomial checks retain independent linear and quadratic
forcing contributions at equal K and matching velocity, including vector
forcing. Current forced production is not substituted into the zero-forcing
reference. The negative virtual viscosity intercept is retained and verified.
Rounded tabulated covariance eigenvalues use a separate 5e-11 comparison
tolerance; finite-difference actions use explicit truncation/roundoff bounds.

The normalized reference transverse flux block from the native physical
stress is

    [-0.9645579065642057  -1.148643029384739]
    [ 1.148643029384739    0.999999999999996]

The independent transverse principal coefficient is -0.9645579065642249.
Native central differences at amplitudes 1e-3, 1e-5 and 1e-7 give
-0.9645559191928263, -0.9645579063654883 and -0.9645579065642057.
The reference covariance remains positive and the first production variation
for this mode is zero.

The nearby smooth strain-dominant 2D state has independent principal form
-0.025967577900340294. Native physical differences at 3e-9, 1e-9 and 3e-10
give -0.025967569347279784, -0.02596769547009381 and
-0.02596739918561919. Every perturbed sample stays on the strain branch.
These small steps are required by its small invariant gap; roundoff
amplification is distinct from the complex-step result.

The formal in-plane strain/rotation reference forms are
-0.025978454472591617 and -0.05298998091898077. Actual smooth limiting
sequences approach the respective values. The molecular diagnostic adds
the expected 0.025 to the transverse coefficient.

These results close M1's requirement to reproduce the failed constitutive
condition. The inner-k positivity counterexample and spectral/hp wall-source
implementation belong to M2. The revised coupled patch and its time map
remain later work; these kernel checks do not establish their implementation
or stable long-time terrain behavior.

**Environment and provenance**

All builds and numerical checks ran on jay-Super-Server under
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised.
The tested repository VERSION is 5.9.0, with GCC 11.4.0, CMake 3.22.1,
Open MPI 4.1.2, Python 3.10.12, NumPy 1.21.5 and SciPy 1.8.0.
The installed Release configuration has NEKTAR_USE_MPI=ON. The note's
version-5.10.0 interface discussion does not establish later wall-patch
API compatibility with this checkout.

Compilation used at most 16 compiler workers. Tests ran sequentially,
pinned to CPU 0, with one MPI process for MPI-enabled solver tests.
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS and BLIS_NUM_THREADS
were one; the independent driver also sets VECLIB_MAXIMUM_THREADS to one.
The point-kernel probes do not launch an MPI job. M1 makes no MPI
partition-equivalence claim.

Git bundles transferred the source into dedicated detached worktrees.
Both tested worktrees were clean. The unrelated main checkout remained
clean at f4b70baec10016dd098005bc5a5507b87e4a57bf; its installed libraries
and build were preserved. Existing local literature and untracked work
were preserved.

Compact evidence is committed with this report:

- [verification.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m1/verification.json):
  inputs, errors, reference values, difference sequences, source/probe
  fingerprints and raw report hashes.
- [baseline-receipt.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m1/baseline-receipt.json):
  exact base build and 16-test commands, times and log hashes.
- [execution-receipt.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m1/execution-receipt.json):
  exact M1 build, direct oracle and 17-test suite, all with exit zero.
  Measured build/oracle/suite wall times were 35.55/1.57/132.93 s.
- [integration.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m1/integration.json):
  framework/configuration/compile flags, three binary hashes and 37 resolved
  shared-library hashes. Direct and CTest runs produced byte-identical
  kernels.json and constitutive.json files.

Raw reports, logs, build outputs and smoke fields remain on Linux. The
primary M1 run is m1-attempt-001; the base run is baseline. The evidence
belongs to e5ab598. The following closure commit updates documentation
and receipts only, preserving the tested code.

**Reproduction**

Use a clean Linux worktree at e5ab598 with the committed revised note.
Choose a fresh direct-oracle output directory; it refuses to overwrite an
existing result. Receipts preserve the executed commands as structured arrays.

    source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
    CLS_EP_RUN=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised
    CLS_EP_SOURCE=$CLS_EP_RUN/source
    export TMPDIR=$CLS_EP_RUN/tmp
    cmake -S "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS" \
      -B "$CLS_EP_RUN/build" -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
      -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
    cmake --build "$CLS_EP_RUN/build" --parallel 16
    taskset -c 0 python3 "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m1.py" \
      --kernel-probe "$CLS_EP_RUN/build/RANSKernelProbe" \
      --constitutive-probe "$CLS_EP_RUN/build/RANSConstitutiveProbe" \
      --output "$CLS_EP_RUN/m1-replay/verification"
    taskset -c 0 ctest --test-dir "$CLS_EP_RUN/build" --parallel 1 --output-on-failure

RANSCLSEPAWFM1Oracle is registered on Linux when NumPy/SciPy are available
at configure time. Confirm it appears in the test inventory; the remaining
smoke tests alone cannot close M1 if those dependencies are absent.
