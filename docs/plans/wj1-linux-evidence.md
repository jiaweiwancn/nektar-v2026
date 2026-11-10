# WJ-1 compiled Linux evidence — 10 November 2026

**Result: 12/12 selected CTest tests passed, zero failures, 18.65 seconds.**
The WJ oracle completed 23,923 comparisons using 544 compiled factory/adapter
queries. The exact reference tangent error is 3.55e-15. The exact rational LDL
certificate proves the incompressible frozen-k/epsilon principal bound
nu_eff/250. Both required adverse witnesses were reproduced without clipping.
These are kernel results, not WJ native coupled equilibrium verification or a
global realizability/stability proof.

## Revisions, source identity and preservation

- Reviewed/fetched parent: `ba0f3068c334bdbe43ce0cadf5ad3ad45340baae`.
  `origin/codex/cls-epawf-pimple-one-sweep` had exactly this head when fetched;
  no advanced head was substituted.
- Mac branch/worktree: `codex/wj-earsm-awf`,
  `/Users/jay/Documents/nektar-v2026/.worktrees/wj-earsm-awf`.
- Production commit: `ec55a770536ffbb2a0a6e61fdda69fef2b690531`.
- Oracle/contract/pressure commit: `70652820bb1f984e5bd403ca523af69796529426`.
- **Final tested source:** `f99755d1b0f4d60ce84ff00088beab3b20273781`;
  Git tree `da96b12e426fd8565d52b1c5031711512c51a9e1`.
  Subsequent evidence-only files do not change the compiled sources.
- Linux run root: `/media/jay/新加卷/nektar-v2026-runs/wj1-20261110`;
  source `source/`, build `build/`, large/temporary data `tmp/` beneath that root.
- USB SSH timed out; the specified `jay@100.105.247.110` fallback was used.
- The original Mac checkout remained at
  `bf381be43b1617763148a412a9dc761d6143e086`; its untracked work and the
  original specification files were preserved. The Linux main checkout was
  observed at `f4b70baec10016dd098005bc5a5507b87e4a57bf` with existing dirty work
  and live CLS campaigns. Neither main working tree was switched or cleaned.
  Git object fetches and isolated worktree registration were the only shared
  repository metadata updates; installed libraries and live executables were
  not replaced.
- NTFS presents 100644 files as executable. `git -c core.filemode=false` was used
  only on individual isolated-worktree commands, without changing repository
  configuration. Hashing every tracked file against its Git blob verified **all
  5,973 file contents**, with zero mismatches and no content changes in the
  tested Linux working tree.
- Full source SHA-256 manifest digest, independently matched on Mac and Linux:
  `f22e265de4a0712e4c5f702ea9706a92f3d8182e2c45cb0a15c2ff7882810c33`.
  The full manifest remains at the Linux run root as `source-content.sha256`.
  The compact changed-file manifest is [changed-source.sha256](wj1-linux/changed-source.sha256).

Actual specification SHA-256 values (Mac originals, Mac copies and Linux copies
all agree with review):

```
1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124  wj_earsm_neutral_abl_ep_awf_note.tex
0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c  wj_earsm_budget_cutoff_awf_appendix.tex
```

The main-note `eq:integrals-linear` `quad`/`\quad` erratum and intended J1
expression are recorded in [the theory-to-code contract](wj1-theory-to-code-contract.md).
The exact specification files retain their reviewed bytes.

## Compiler, dependencies and allocation

Host: `jay-Super-Server`, x86_64 Ubuntu, Linux 6.8.0-138, glibc 2.35.
GCC 11.4.0, CMake 3.22.1, C++17 Release (`-O3 -DNDEBUG`). Kernels build with
`-Wall -Wextra -Wpedantic -Wno-sign-compare -Werror`; the WJ source and probe
also explicitly enable `-Wall -Wextra -Wpedantic -Werror`. NumPy 1.21.5,
SciPy 1.8.0 and SymPy 1.9 were available. The independent polynomial-root,
complex-step, quadrature and exact-arithmetic implementations do not call the
production evaluator to generate expected values.

Read-only installed dependencies:
`/home/jay/.local/opt/nektar-v2026-mpi`, configured Release with
`NEKTAR_USE_MPI=ON`; environment `env.sh` and CMake package
`lib64/nektar++/cmake`. All changed kernel sources were compiled into the isolated
static `RANSKernels` library and source-matched probes. No solver executable or
installation was built/replaced. The shared value/reference types have additions;
future native consumers must rebuild against their matching headers and kernels,
not mix these headers with an older installed kernel ABI.

Every build used `--parallel 1`, every CTest run `-j 1`. No concurrent WJ build
and numerical test jobs, no mpirun multi-rank job, and no long-time/terrain/
optimization/tuning campaign were launched. Serial probes linked to MPI use one
rank. OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS and BLIS_NUM_THREADS
were all explicitly 1. Existing campaigns were inspected before work; the final
prebuild snapshot contains 66 relevant application/observer/runtime threads
(including idle MPI daemon threads), with one additional compilation worker.
Numerical tests add one computational rank at a time. Snapshot files are retained
in the run root; shell command-string matches were excluded from the exact
process-name counts. Allocation stayed within the 72-computational-worker limit.

## Exact commands

Source revisions were transferred by Git bundles and checked out detached in the
isolated Linux worktree after checking content cleanliness. The initial fetch of
the reviewed parent used:

```bash
ssh jay@100.105.247.110 \
  'git -C /media/jay/新加卷/nektar-v2026 fetch --no-tags origin ba0f3068c334bdbe43ce0cadf5ad3ad45340baae'
```

The final run used the following environment and commands on Linux:

```bash
root=/media/jay/新加卷/nektar-v2026-runs/wj1-20261110
export TMPDIR="$root/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh

cmake -S "$root/source/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$root/build" -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_BUILD_SOLVER=OFF

cmake --build "$root/build" --parallel 1 --target \
  WJKernelProbe RANSKernelProbe RANSConstitutiveProbe RANSAcceptedAWFProbe \
  RANSEquilibriumProfilesProbe RANSProfileBranchChecks RANSKBoundQPChecks \
  RANSHistoryRejectionChecks RANSStressTangentAudit

ctest --test-dir "$root/build" --show-only=json-v1 > "$root/test-inventory.json"
ctest --test-dir "$root/build" -j 1 --output-on-failure \
  -R '^(WJKernelContracts|WJIndependentOracle|RANSKernelUnitTests|RANSCLSEPAWFM1Oracle|RANSCLSEPAWFM2Oracle|RANSAcceptedAWFKernel|RANSEquilibriumProfilesNative|RANSProfileBranchChecks|RANSKBoundQPChecks|RANSHistoryRejectionChecks|RANSAcceptedAWFCases|RANSStressTangentOracle)$'

python3 -B "$root/collect_evidence.py" "$root"
```

CTest inventory contains 22 tests. All twelve required selections and all
NumPy/SciPy/SymPy dependencies are present; none of the required tests is missing
or silently skipped. CTest handled selected tests normally with their registered
properties; there were no extra fixture requirements for this selection. The
other ten registered tests were not selected. No new native solver plumbing,
observation or timing component was implemented, and no solver executable was
built; native solver transport/equilibrium smoke campaigns were outside WJ-1.
The unchanged history unit probe was built for its explicitly required regression.

## Final test results

| Test | Result | Seconds |
|---|---|---:|
| RANSKernelUnitTests | PASS | 0.03 |
| RANSCLSEPAWFM1Oracle | PASS | 6.55 |
| RANSCLSEPAWFM2Oracle | PASS | 6.18 |
| RANSAcceptedAWFKernel | PASS | 0.01 |
| RANSEquilibriumProfilesNative | PASS (legacy profile probe) | 1.21 |
| RANSProfileBranchChecks | PASS | 0.01 |
| RANSKBoundQPChecks | PASS | 0.00 |
| RANSHistoryRejectionChecks | PASS | 0.02 |
| RANSAcceptedAWFCases | PASS | 0.18 |
| RANSStressTangentOracle | PASS | 0.09 |
| WJKernelContracts | PASS | 0.00 |
| WJIndependentOracle | PASS | 4.35 |

Raw final log: [ctest-final.log](wj1-linux/ctest-final.log).
Independent results: [wj1-oracle.json](wj1-linux/wj1-oracle.json).
The initial full focused run also passed 12/12 at 7065282 in 15.91 seconds; the
suite was repeated because the final legacy-wall rejection guard changed code.

## Numerical results, tolerances and adverse states

Default WJ parameters are c1=1.8, c2=5/9, Ce1=1.44, Ce2=1.82,
sigma_k=1, sigma_epsilon=1.3 and exact Cd=C*=572/6561. Legacy defaults were
preserved. Scalar diffusion uses Cd, not Ceff; off-reference tests explicitly
exercise unequal scalar and stress viscosities. Nonzero outer molecular requests
and unsupported legacy parameter/variant/bool-branch calls reject.

Base tensor/root comparisons use scaled relative/absolute tolerances about 3e-12;
scaled cubic residual is checked at 1e-14. Directional actions are compared with
independent complex-step derivatives (2e-10 relative, 5e-10 absolute) and centered
differences at steps 1e-3, 3e-4, 1e-4, 3e-5. Worst final normalized centered-
difference error is 2.29e-9, below 1e-6 before roundoff. The complete reference
tangent restricted to all eight trace-free coordinates agrees within 3.55e-15;
the independent exact positive LDL pivots prove nu_eff/250 on every
incompressible rank-one increment. No finite scan is used as that proof.

Finite-nu integrals and source/diffusion weak balances use absolute tolerances
1e-10 or tighter and relative 3e-12 for nonzero scales. Pressure/neutral balances
use about 1e-12. The largest unscaled error across **all** comparisons is 2.33e-10
(on large integral/tensor coefficient scales); it is not an equilibrium residual.
Scalar source modes 0 through 5 use independent adaptive quadrature split at
roughness, viscosity, baseline/corrected cutoffs and the tail start. The first
oracle attempt identified a missing tail-start split in the independent smooth-
wall gate integral; adding the physical breakpoint resolved it without loosening
the tolerance.

Compiled canonical wall values with K=0.8466950351321795, nu=1.5e-5, H=1, h=.1:

| Quantity | Compiled value |
|---|---:|
| c_length_wall | 2.3804029452632984 |
| zv | -0.005393275866854201 |
| a | 0.1909589317733666 |
| J0 | 27.302203755064866 |
| J2 | 4984.456376333159 |
| counterfactual Pbudget | 1.701714806838742 |
| zc | 0.015005743260176297 |
| anchored UH | 6.8255509387662165 |
| z0 | 0.0054420488394214764 |
| budget k_ref(0) | 1.8887084080086478 |
| budget k_ref(.1) | 0.9072068354060853 |
| budget k_ref(.5) | 0.8495674195161087 |
| pressure jump/rho | -0.20906050250177272 |
| molecular power loss | 0.004672927852812337 |

Both compiled Q0 endpoints are zero in this canonical query; the pressure plus
gz+2k_ref/3 is constant to roundoff. EP evaluated at its own returned anchor
cancels pointwise exactly, hence in every scalar mode. Supplying the independently
rounded UH=6.825550938766217 yields at most 4.25e-16 in the six reported loads.
Budget has nonzero nonconstant source modes, balanced by ordinary diffusion,
and no EP redistribution. Full canonical data: [canonical.json](wj1-linux/canonical.json).
The exact stream queries used were:

```
W 0 1 0 0.8466950351321795 1.5e-5 1 .1 .5 0 0 6.825550938766217 0 0 0
W 1 1 0 0.8466950351321795 1.5e-5 1 .1 1 0 0 6.825550938766217 0 0 0
```

They are supplied to `build/WJKernelProbe --stream`; formulation 0 is WJ EP and
1 is WJ budget. The curved test J=(1+.2r)(1-.1r) gives
zc=0.015036191831751875 from the metric root. Frozen-K forced/reversed matching
data leave the counterfactual cutoff unchanged. Smooth molecular-only segments,
chi=0 calibrated baseline, partial chi, failed preservation gates and explicitly
selected direct-baseline fallback are exercised.

Required adverse constitutive witnesses:

- q=3.5: N=4.130462147295278; normalized principal form
  **-0.006333781035394122** in the requested direction. Covariance eigenvalues/k
  are (0.27877891195383653, 0.6666666666666665, 1.0545544213794966), and
  P/epsilon=1.0357609543534572. The independent rational polynomial/sign
  certificate is retained with the exact reference certificate.
- S=diag(2,2,-4), Omega_xy=-Omega_yx=3 sqrt(17)/2: N=3,
  R/k=diag(-2/15,-2/15,34/15), P/epsilon=48/5. The negative variances are
  preserved rather than hidden by PSD clipping or an implicit-3D root.
- An additional smooth budget case (same K,nu,H; h=0, zm=H) passes the cutoff
  gate but gives compiled k_ref(0)=-0.3301275029511697, independently
  -0.330127502951175. This reference is unqualified for positive initialization;
  the API returns it without clipping. The rough canonical reference remains
  positive. Dynamic sampled-K destruction also has no positivity guarantee.

## Binary identities and scope

[identities.json](wj1-linux/identities.json) records all nine probe SHA-256
identities, complete ldd output, 34 resolved linked-library SHA-256 identities,
compiler flags, dependencies, test inventory/properties and thread controls.
In particular:

```
2d08494f7ed0284b8109fa47bbdcffb36216958e747c24a8acd5fd3cc64ddb7e  WJKernelProbe
2485c4bb49bd0f130f218f29118764b5a0ac87109d0dfcfbf306516ec21ffb06  libRANSKernels.a
```

[collect_evidence.py](wj1-linux/collect_evidence.py) reproduces the read-only
identity/content audit. Full build/configuration logs, CTest JSON inventory and
source manifest remain under the Linux run root. No WJ native PDE selection is
available, no coupled native equilibrium preservation is claimed, and neither
reference-local ellipticity nor the wall reference is a global stability,
realizability or evolving-positivity proof.

The bounded WJ-2 checklist is the final section of the
[theory-to-code contract](wj1-theory-to-code-contract.md): complete native
residual/action/diagnostic dispatch; full traction/pressure jump; one outer
epsilon trace; model-aware signatures; zero-step uncancelled reference residual;
two-step versus 1+1 restart; rejection rollback; then spatial/reference and short
temporal tests.
