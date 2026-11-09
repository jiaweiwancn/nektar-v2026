# A93 / A95 / E94E95: first implementation and qualification

Status: implemented and passed through bounded native serial flat-2D integration. Primary numerical revision `1950ded`; guarded legacy replay revision `4ee3dee`.

## Reproducible workspace and specification provenance

The requested branch is `codex/a93-a95-e94e95-awf`, based on reviewed commit
`99112a1cbb9c5448dd3de69f50600bf91465062b`. Fetching the parent on 2026-11-09
returned that same commit, so there are no intervening parent changes.
The original Mac checkout is on `codex/cls-epawf-segregated-picard` at
`bf381be43b1617763148a412a9dc761d6143e086`, with untracked literature and
unrelated plans. None is removed or staged from that checkout.

Mac implementation worktree:
`/Users/jay/Documents/nektar-v2026/.worktrees/a93-a95-e94e95-awf`.
Linux experiment root:
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-awf-20261109`.
Linux source is a separate detached worktree under `source`; `build-baseline`,
`build`, `baseline/runtime`, `runtime-final`, `runtime-qualified`, `baseline`,
`matrix`, and `matrix-qualified` are separate paths beneath this experiment root. The established installation is a read-only
dependency, never an installation destination. Heavy files and TMPDIR stay here.

Direct SSH `linux-usb` timed out; `ssh jay@100.105.247.110` works. The Linux main
checkout is at `f4b70baec10016dd098005bc5a5507b87e4a57bf`, on a different
branch, with many existing changes. It must not be reset or synchronized in place.
Initial conservative campaign process/thread census: 41 threads, including six
serial solvers with three runtime threads each, controllers and observers.
Initial build allocation: eight compiler workers; native tests: one at a time,
MPI ranks 1, OMP/OPENBLAS/MKL/BLIS/NUMEXPR threads 1. Recheck before new jobs.

Both supplied TeX sources were read completely and are preserved byte for byte.
The requested combined filename repeats `E94`; the available filename is
`A93_A95_E94_E95_neutral_ABL_EP_AWF_combined.tex`. No specification was inferred
from the model names. SHA-256:

| Supplied file under docs/literature | SHA-256 |
| --- | --- |
| A93_A95_E94_E95_neutral_ABL_EP_AWF_combined.tex | 056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64 |
| accepted_awf_appendix_A93_A95_E94_E95.tex | 17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40 |

See `a93-a95-e94e95-awf-evidence/specifications.json` for source paths and sizes.
The two `a93_a95_e94e95_note_*.py` scripts in the RANS test directory are
extracted from the combined note's executable listings; only end-of-file
whitespace is normalized. They are mathematical
oracles, not evidence of a native solver implementation.

## Implementation contract

One immutable resolved descriptor supplies closure identity/version, effective
coefficients, total-epsilon semantics, exact equilibrium root, kappa policy,
wall formulation/calibration, outer zero molecular diffusion, inner physical
viscosity, and the named DeviatoricStrainV1 numerical extension. For new models
the native field `eps_tilde` stores total epsilon. Original outer coefficients
are retained except sigma_epsilon for an explicitly prescribed kappa.

Physical tensors are A_ij = partial_j U_i, S=(A+A^T)/2 and Omega=(A-A^T)/2.
The numerical extension evaluates the new constitutive invariants and nonlinear
tensors on S_d=S-tr(A)I/3 and Omega. Production is the complete signed -R:A
with the raw A, including isotropic pressure work. CLS retains its existing
doubled/time-scaled tensors and its existing divergence treatment. Native 2D
uses physical x,z; its second velocity component is z. Keep all nine tensor
entries, the spanwise normal stress, and 2kI/3.

A93: original variable Cmu and full quadratic law, no cubic/Yap/split epsilon.
A95: A0=6.5, original square-root commutator; negative radicand is a domain
failure. E94E95: A0=4, no quadratic stress; epsilon production is
max(.43,eta/(eta+5))*sqrt(2 S_d:S_d)*epsilon, destruction 1.9 epsilon^2/k.
A93/A95 use 1.44 epsilon P/k and 1.92 epsilon^2/k. Preserve signed sources
under the one-sweep numerical source split.

At zero strain and rotation the coefficient limit is unique. Nonzero pure
rotation in A95/E94E95 is rejected because the original scalar viscosity limit
is unspecified. Repeated largest eigenvalues admit forward values; smooth
derivatives are not claimed. A95 zero radicand admits a forward value but not
a general smooth tangent. Unsupported tangent/action routes reject explicitly.

EP-AWF retains SCI constants .09 and 2.55, its declared sampler, constant inner
reference k=K, ray-constant switch, E_pre and D_pre and one redistribution.
The physical source is P_A-D_pre. Evaluate chi=0 directly without requiring the
inactive corrected cutoff. Optional positivity/constitutive modifications are
not primary new-model defaults.

AcceptedBudgetAWF uses endpoint K, compatible c_mu=q^-2 and c_l=kappa*q^1.5,
finite nu in D_A and production, counterfactual p0(K), Lambda>=1 and an unclipped
budget cutoff. It transports the nonconstant reference k initialized from Q0;
source P_A-E_w, zero prescribed redistribution, no chi fade, no local E_w/k
reaction. Its cutoff and sink are frozen during the declared scalar assembly.

Stored pressure is Pi=pi_mechanical+2k/3. Reference pressure offsets and actual
interface traction derive from the complete selected outer covariance, with
one common gauge. Keep common velocity and k traces, conservative k flux, and
only the dependent inner epsilon trace for outer epsilon. No extra inner
epsilon transport or two-sided epsilon flux equation is added. Preserve stage
ordering, histories, and complete rejected-step rollback.

New native support is serial flat 2D PimpleSingleSweep only. Curved geometry,
true 3D, MPI, converged nonlinear/Newton/action routes remain guarded. Original
outer closures and linear inner stress are primary tests. Restart signatures
must include the resolved descriptor and reject model/calibration changes.
Cache admissibility must include velocity-gradient dependence for A95.

## Source/specification mapping and native audit

| Native source | Responsibility / governing specification |
| --- | --- |
| RANSModels.* | Combined note out:gradient through out:epsE; exact roots, compatibility; resolved descriptor |
| AnalyticalWallFunction.* | Original finite-nu SCI kernel and signed forcing/traction; shared infrastructure |
| EquilibriumAWFProfiles.* | Combined note awf:pre-profiles/net-source; appendix mawf:cutoff/k-profile/weak |
| RANSEquilibriumWallEvolution.cpp | Session resolution, native driver, complete histories/restart checks |
| RANSWallEvolutionPatch.* | Fields, reference initialization, raw residual, interfaces and signatures |
| RANSWallEvolutionPimple.cpp / Picard.cpp | One-sweep ordering and frozen source split; guard unsupported alternatives |
| RANSWallEvolutionNonpressure.cpp | Explicit rejection of new-model constitutive actions; no Boolean standard fallback |
| RANSWallEvolutionResolution.cpp | Actual closure epsilonProduction/source diagnostics, including E94E95 |
| RANSWallEvolutionEvaluation.cpp | Exact-step evaluation reuse; gradient-dependent constitutive-domain validation |
| RANSWallEvolutionContact.cpp | Separate specification hashes, resolved wall calibration and model-aware covariance checks |
| RANSWallPhysicsSignature.h | Includes patch signature and session controls; extend without changing compatible CLS identity |

The inherited accepted wall supplied nonconstant ReferenceK, zero redistribution,
frozen prescribed sink and a flat serial one-sweep path. Its former CLS-only
metadata and bound inventory have been extended for the resolved new closures.
The enabled path is RANSEquilibriumWallEvolution and its native patch, not the
older RANSVelocityCorrectionScheme. The standalone CLS resolution replay probe
remains separately guarded for new models; native forward residual/audit
evaluation supports them.

## Qualification matrix and acceptance

Baseline, required before numerical implementation edits: existing relevant
kernel tests and a short native CLS control. Record source-file manifests,
build configuration, executable and every linked-library hash, XML hashes,
environment, commands and worker allocation. A commit label alone is inadequate.

Independent closure checks: symmetry/trace; 3D rotations; A93 independent A-based
and S/Omega-based evaluations; signed production and actual scalar sources;
E source branch; shear roots and epsilon balance; natural/prescribed policies;
zero/singular/domain controls; enabled tangents including coefficient dependence;
expected failures of covariance, strong ellipticity and dynamic positivity.

Independent wall checks: breakpoint-split quadrature; forced/reversed momentum;
finite-nu work and turbulent budgets; Lambda below/equal/above 1; molecular-only
rejection; endpoint dependence; EP plateau and zero-switch baseline; accepted
nonconstant weak tests; smooth negative-reference controls and rough positive
controls; invariant kernel when interior k changes at fixed endpoint.

Native order: E94E95, A95, A93. Each has natural / prescribed .41 crossed with
EP / accepted walls (12 rows). For each: admissible reference and documented
polynomial representation; unanchored initial spatial residual; consistent
histories; 1,2,8,16 accepted steps; drift, stage constraints, interface jumps,
admissibility, signed/absolute correction budgets. Add mild perturbations,
spatial/quadrature and dt refinement. Never cancel a reference residual or reset
an evolving field to the reference. Distinguish continuum, projection, spatial
and temporal errors. Report separate momentum, continuity, pressure/gauge,
inner/outer k, outer epsilon, dependent scalar and interface residuals.

Regression gates: baseline CLS and standard kernels; compatible restart;
changed model/calibration rejection; rollback; representative cache on/off;
unsupported routes. Large data remain in Linux run directories. Compact JSON
results and minimal failure reproducers belong in this evidence directory.

## Results and numerical interpretation

The reviewed baseline and the isolated implementation build pass on Linux.
At `597b59b`, all twelve primary cases reach 16 accepted one-sweep steps at
dt=1e-4, T=0.0016, without bound projection, constitutive clipping, anchoring
or reference resetting. The final qualification revision is `1950ded`; its
complete source manifests, executable/library identities and final receipts
are recorded separately below. No long-time, physical, complete temporal-order
or curved/3D/MPI qualification is claimed.

The two unchanged combined-note executable listings pass on Linux: all interval/exact-sign certificates, and 873 wall assertions over 22 cases. The isolated baseline external build succeeds with explicit Nektar++_DIR=`/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake`, Release, `-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS`, eight workers. Existing kernel oracle, accepted-AWF probe and profile-branch checks pass before numerical edits. The new Linux worktree has NTFS executable-bit differences only: `git -c core.filemode=false status --porcelain` is empty. No content is reset to hide these differences.

Baseline native CLS control passes two serial steps at dt=0.000625, UniformEPFlat, diagnostic R2 guard, with its existing SaturatedR2/regularized controls unchanged. The private ELF/library snapshot and baseline.json are under `baseline/`. Numerical implementation edits began only after this result was observed.

The independent production-kernel suite at `ff915d3` passed 1326 checks; the
expanded final suite at `1950ded` passes 1346 checks on Linux.
Its first run identified a test classification error: a pure strain A95 state
with tau=1.2/.7 lies outside the printed square-root domain. The corrected
harness requires that rejection, instead of changing the closure. Both original
and corrected reports remain in the Linux experiment. All stated adverse
covariance/ellipticity/production and smooth-reference controls retain their
prescribed mathematical signs.

Native new-model audit schema reuses the established `accepted_awf_*` keys for
model-aware scalar, budget and covariance checks, including EP wall cases. The
`resolved_physics` and `wall_formulation` metadata identify the actual model;
these keys do not select the accepted wall equations. New diagnostics partition
the actual common-CG k equations into inner-only basis support, outer-only
support and shared interface rows. A shared weak conservation row is counted
once; reconstructed pointwise traction/flux jumps remain separate diagnostics.
All three restricted dual norms use the full k volume and K0*u*/H scaling.
The ordinary `coupled_norms` order in 2D is u, vertical velocity, k, outer total
epsilon, inner continuity, outer continuity, pressure gauge, constraint maximum.
CLS-specific nonlinear/source-forensic and counterfactual-regularization routes
are explicitly rejected for the new closures; primary native audits use actual
closure source outputs and complete 3D covariances.


The initial reference is generated by the native descriptor and projected into
continuous p=3 velocity/k/epsilon and p=1 pressure spaces. The flat tensor mesh
has x in [-1,1], H=2, upper height 4, and radial panel endpoints
0,.0025,.005,.01,.02,.04,.1,.25,.5,1,2,2.25,2.5,3,3.5,4.
Native GLL reconstruction uses eight points per coordinate; wall integration
uses normal order 32 and tangential order 8. The accepted cutoff is about .015,
resolved by the graded panels. Physical nu=1.5e-5, roughness=.1 and u*=.4.
The primary EP reference uses its declared point sampler at H and chi=1;
the accepted wall uses mandatory endpoint sampling, its model-compatible wall
constants and its nonconstant initialized k. IMEX1 native history starts from
the complete represented initial state. Every accepted step and checkpoint is
retained on Linux.

These polynomial references are not exact discrete equilibria. For the E94E95
natural-policy examples at `597b59b`, the unanchored initial momentum norm is
about .4904; k is .01920 for EP and .37129 for accepted; outer epsilon is .001700.
Continuity, pressure gauge and imposed constraints are near floating-point
roundoff. Halving radial panels and increasing wall normal quadrature to 48
reduces momentum residual by a factor .418, accepted k by .261 and epsilon by
.241. The accepted 16-step velocity drift decreases from 1.068e-4 to 2.257e-5,
k drift from 2.182e-5 to 4.579e-6 and epsilon drift from 6.479e-7 to 1.021e-7.
The EP refinement shows the same momentum/epsilon behavior and reduces its k
residual by .459. This is evidence of representation/spatial error, not a
roundoff preservation claim.

At fixed mesh, halving dt and doubling steps to the same T changes the physical
velocity L2 by about 2.42e-7, accepted k by 9.34e-8 and epsilon by 1.05e-9.
One dt refinement distinguishes a large implementation error from splitting
error; it does not establish an order of accuracy. Single-sweep refreshed
nonlinear source/trace defects are reported rather than repaired after the k
solve. A weak shared interface conservation row and a reconstructed pointwise
traction jump measure different quantities.

The first continuation harness at `597b59b` correctly exercised all three
incompatible restart rejections, but its report summarizer indexed an empty
initialization-state list. Its rollback comparison also hashed JSON integer
`-0` against a binary negative zero: all numerical values were equal, but the
reader had lost the sign. Raw failures remain in `qualification/`. Revision
`1950ded` handles empty startup receipts and writes negative zero as `-0.0`;
it does not modify rollback mathematics. The final rerun checks binary64
state and native history including signed zero.

Mathematical and architectural discrepancies found and resolved:

* The original uncapped EP scalar assembly omitted its prescribed wall source;
  it now assembles P_A-D_pre exactly once, explicitly, with no manufactured
  local sink rate. The legacy capped route is unchanged.
* The inherited Boolean non-CLS source diagnostic would have reconstructed
  standard epsilon production for the new models. Diagnostics now consume the
  closure's actual epsilonProduction, including the distinct E94E95 term.
* Scalar-only admissibility reuse is insufficient for A95. New-model checks
  evaluate the current one-sided native velocity gradient; CLS reuse remains.
* The new full tensor evaluator is distinct from doubled-time-scaled CLS
  tensors. DeviatoricStrainV1 retains the raw-divergence production term and
  reports raw divergence; it is an explicit numerical extension of the notes.
* The published adverse controls show no global covariance, ellipticity or
  dynamic positivity guarantee. They are retained as negative controls, not
  used to tune or regularize the original closures.

The two wall formulations remain separate. EP uses original SCI wall constants
and its redistribution; accepted uses model-compatible constants, an unclipped
budget cutoff, frozen sink and zero redistribution/reaction. EP also admits
its declared flat Upstream switch. A small-dt structural switch control is
separate from the twelve equilibrium-reference runs: the original chi=0 sink
has a stiff finite-nu molecular cutoff and is not an equilibrium initial state.


## Twelve-configuration bounded integration receipt

Production solver and primary kernel probe: `1950dedf4e04a903a3e919b64bc3ce7a25c25358`.
Later commits add evidence collection and the separately guarded legacy replay
probe; production closure/wall/native solver sources are unchanged.
All rows run serial, IMEX1, dt=1e-4, T=.0016, eta=0, bound projection Off.
The drift columns are the native scaled velocity, k and epsilon L2 departures
from the represented initial reference after 16 steps.

| Closure | Kappa policy | Wall | Accepted steps | Velocity drift | k drift | Epsilon drift |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| E94E95 | Natural | EP | 16 | 0.0001064 | 4.457e-06 | 5.997e-07 |
| E94E95 | Natural | Accepted | 16 | 0.0001068 | 2.182e-05 | 6.479e-07 |
| E94E95 | Prescribed | EP | 16 | 0.0001064 | 4.457e-06 | 6.01e-07 |
| E94E95 | Prescribed | Accepted | 16 | 0.0001069 | 2.199e-05 | 6.503e-07 |
| A95 | Natural | EP | 16 | 0.0001046 | 3.621e-06 | 5.256e-07 |
| A95 | Natural | Accepted | 16 | 0.0001071 | 1.962e-05 | 5.789e-07 |
| A95 | Prescribed | EP | 16 | 0.0001046 | 3.621e-06 | 5.215e-07 |
| A95 | Prescribed | Accepted | 16 | 0.0001068 | 1.92e-05 | 5.723e-07 |
| A93 | Natural | EP | 16 | 0.0001062 | 4.308e-06 | 6.301e-07 |
| A93 | Natural | Accepted | 16 | 0.0001059 | 2.015e-05 | 6.647e-07 |
| A93 | Prescribed | EP | 16 | 0.0001062 | 4.308e-06 | 6.38e-07 |
| A93 | Prescribed | Accepted | 16 | 0.0001069 | 2.151e-05 | 6.8e-07 |

All twelve initial pressure offsets match the independent full-stress oracle
to 4.2e-17 absolute. Final full native tensor witnesses match the independent
physical tensor oracle to 1.82e-16 relative; both include all three normal
stresses. Minimum sampled covariance eigenvalue across the primary trajectories
is .17727, minimum transported k .53186 and minimum outer total epsilon .03735.
These are finite sampled checks, not global realizability proofs.
Minimum wall budget ratio is 5.315 and maximum relative cutoff budget defect
is 1.89e-16. Bound correction calls and signed/absolute correction integrals
are zero for every primary row.

The weak continuity/gauge equations are near roundoff, while the reconstructed
raw divergence reaches about 2.3e-4 of the local reference shear. The outer
constitutive law therefore uses DeviatoricStrainV1, and production still
contracts the complete R with this raw gradient. The compact receipt includes
dimensional raw divergence at the inner and outer covariance witnesses.

The imposed frozen-stage epsilon trace is checked to 1e-12, but the refreshed
dependent trace after the subsequent k solve has a measured splitting defect,
up to 1.36e-7 in epsilon units. Common velocity and k traces remain at roundoff.
The defect is reported in the final dependent-scalar residual and is not
repaired after the scalar solve. Single-sweep source-storage identities are
checked to 1e-9 and continuity/gauge to 1e-10 in the native scaled norms.
The full JSON also reports actual inner/outer k partitions, outer epsilon,
pressure/gauge, interfaces, milestone drift, and pointwise one-sided traction
and flux jumps. Those pointwise jumps are not substituted for shared weak
interface conservation.


## Expected negative controls

These are specified outcomes of the original models, and pass only when the
failure or adverse sign is retained. No coefficient clipping, realizability
repair or constitutive regularization is added to make them positive.

| Control | Expected outcome, independently verified |
| --- | --- |
| A93 pure rotation, tau omega=8 | R_yy=-.144: covariance positivity fails |
| A93 simple shear q=8 | Normalized acoustic value -51301/296100, although covariance is positive |
| A95 simple shear q=1 | Normalized transverse acoustic value -.044533122876222, with positive radicand and covariance |
| E94E95 strain diag(10,9,-19) | Normalized acoustic value -19/68, with positive covariance |
| A93 mixed strain/rotation example | Signed production -.3114625443 with positive covariance; no clipping |
| A95 strain diag(2,1.8,-3.8), k=epsilon=1 | Radicand -.248768: constitutive-domain rejection |
| A95/E94E95 nonzero pure rotation | Rejection: no uniquely specified original scalar-viscosity limit |
| Accepted smooth wall with budget Lambda=4 | Admissible cutoff but negative reference ground k; native positivity is not implied by cutoff admissibility |
| Accepted rough-wall primary controls | Positive nonconstant reference and correct weak diffusion/source cancellation |
| Accepted zero actual velocity at positive fixed K | Negative prescribed net source: no general dynamic positivity guarantee |
| Lambda<1 or molecular-only accepted kernel | Domain rejection, no cutoff clipping |
| EP chi=0 with inadmissible inactive correction | Direct original baseline succeeds |
| New-model eddy viscosity underflow | Explicit range rejection, no coefficient floor |

The trace-free mathematical certificates also verify positive neutral acoustic
margins. They do not imply global ellipticity, full solver stability or smooth
constitutive tangents at repeated largest eigenvalues or A95's square-root
boundary. Every new constitutive tangent/action route remains explicitly
unsupported; enabled native one-sweep solves use forward values and frozen
coefficients rather than an unverified new-model derivative.


For E94E95 natural kappa at fixed T, the refreshed epsilon trace defect decreases
from 1.02634e-7 to 5.12501e-8 with dt halving for accepted AWF, and from
7.16726e-8 to 3.57872e-8 for EP. Spatial/quadrature refinement decreases those
defects to 7.44802e-9 and 5.17928e-9 respectively. The inner/outer/shared k
stage residuals also decrease by approximately one half under dt halving.
This identifies a measured splitting contribution without inserting a trace
repair or claiming a complete temporal-order campaign.


## Regression, continuation and provenance results

All 28 checks in `a93-a95-e94e95-awf-evidence/matrix-results.json` pass:
compatible restart, complete rejected-step rollback, cache on/off, mild inner-k
perturbation, spatial/quadrature and dt refinement for both walls; changed
model, kappa/calibration and wall restarts are rejected. A95 has an additional
bitwise cache comparison using its current velocity gradient. The accepted
interior-k perturbation changes the transported field while preserving the
endpoint kernel and cutoff exactly. A flat EP Upstream test spans chi=1,
transition and chi=0, and accepts two steps at dt=1e-6. It is a structural
switch/source control, not a continuum equilibrium or production campaign.

A genuine checkpoint generated by the reviewed `99112a1` binary continues with
the new binary to step four and reproduces its complete binary64 coupled state,
native histories and history times bitwise. Independent legacy resolution
replay gives exactly the same native block norms for the reviewed whole run
and the restarted new-binary run (maximum absolute difference zero). New EP
and accepted packets are explicitly rejected by that legacy replay tool,
rather than evaluated under a default CLS descriptor. Its tested revision is
`4ee3dee7749224670d4f459454db7731b8dc4604`; production solver sources and the
frozen solver executable are unchanged from `1950ded`.

The existing 5000-check standard/CLS kernel verification passes, the accepted
wall kernel passes 57 checks, and profile-branch checks pass. Its inherited CLS
physical-stability gate remains `failed_as_specified`; retaining that negative
control is a regression pass and does not qualify global CLS stability.
The baseline two-step CLS native control is repeated successfully. Seven
serial CTest regressions pass: history rejection, standard transport, CLS
transport, original AWF, equilibrium AWF, standard equilibrium AWF and the
independent transverse stress-tangent oracle. No unexpected numerical failure
remains in the exercised scope.

The prototype orchestration completed its numerical checks but initially failed
its final provenance lookup (`KeyError: frozen`). The original log remains in
`qualified-driver.log`. The committed fresh-run driver fixes that lookup;
`finish-linux.py` completed the recorded run's replay checks and final identity
receipt without rerunning or changing any trajectory or checkpoint. This is a
receipt repair, not a dissipation-trace or numerical-state repair.

Compact evidence (relative to this document):

* `a93-a95-e94e95-awf-evidence/matrix-results.json`: all 12 rows, 28 controls,
  unanchored initial/stage residual groups, actual interface diagnostics,
  independent full stress/pressure checks, drift and correction budgets.
* `a93-a95-e94e95-awf-evidence/kernel-results.json`: 1346 independent checks and
  all equilibrium roots, effective sigma_epsilon, cutoff/budget reference rows.
* `a93-a95-e94e95-awf-evidence/baseline.json`: pre-edit source, runtime and control
  identities, including the original source-file manifest.
* `a93-a95-e94e95-awf-evidence/source-build.json`: exact source-file hashes,
  compiler, CMake cache, frozen solver/probe and linked-library identities,
  actual commands, final immutable-runtime check and receipt-helper hashes.
* `a93-a95-e94e95-awf-evidence/actual-loaded-libraries.json`: 70 actually mapped
  native libraries/plugins, executable/cwd and controlled environment from
  the first primary run. This supplements ldd and captures dynamic MPI plugins.
* `a93-a95-e94e95-awf-evidence/regression-summary.json`: standard/CLS regression
  commands, launch-time runtime identities and genuine reviewed-checkpoint
  continuation provenance.
* `a93-a95-e94e95-awf-evidence/replay-build-provenance.json`,
  `resolution-runtime-identity.json`, `environment.json`: the guarded replay
  probe, unchanged production-source proof, compiler/Python/NumPy/SciPy,
  MPI and thread environment identities.

Primary native source-file manifest SHA-256:
`fa16b66c7ac0d9610aa03c6782410bdc8113d84ca8aa6722b6fe435ab692574a`.
Frozen solver ELF SHA-256:
`22e85bab007d2a871a286be1c3c490ae5ba48ba471ecea68a6ec85eb5bf12c08`.
CMake cache SHA-256:
`3f150b0dce850e1c7c1f6e6b0f5b2211aa536b6dca45282685b9a20202f5c098`.

The runtime is MPI-enabled Nektar++ 5.9, used with one rank and one numerical
library thread per run; the unsupported-MPI guard alone is exercised with two
ranks. GCC 11.4.0, CMake 3.22.1, Open MPI 4.1.2, Python 3.10.12, NumPy 1.21.5,
and SciPy 1.8.0 are recorded. Build jobs use eight workers. Conservative process
snapshots plus each allocation remain below 72; runtime MPI progress threads
are counted in the census. Large reports, XMLs, checkpoints, source bundles,
build logs and all trial failures remain under the Linux experiment root.
No established installation or active CLS runtime is replaced.

## Reproduction commands

For the recorded experiment, use the individual commands/paths in
`source-build.json`; every XML and full native report is hashed in the compact
matrix. The executable is
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-awf-20261109/runtime-qualified/IncNavierStokesSolver`.
The complete primary matrix can be regenerated in a new output directory:

```bash
export A93_RUN=/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-awf-20261109
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$A93_RUN/tmp"
export LD_LIBRARY_PATH="$A93_RUN/runtime-qualified/lib"
python3 -B "$A93_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_a93_a95_e94e95.py" \
  --probe "$A93_RUN/runtime-probe-qualified/RANSNewModelsProbe" \
  --output "$A93_RUN/kernel-reproduction.json"
python3 -B "$A93_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/run_a93_a95_e94e95_matrix.py" \
  --solver "$A93_RUN/runtime-qualified/IncNavierStokesSolver" \
  --descriptor-probe "$A93_RUN/runtime-probe-qualified/RANSNewModelsProbe" \
  --output "$A93_RUN/matrix-reproduction" --steps 16
python3 -B "$A93_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/qualify_a93_a95_e94e95_native.py" \
  --solver "$A93_RUN/runtime-qualified/IncNavierStokesSolver" \
  --descriptor-probe "$A93_RUN/runtime-probe-qualified/RANSNewModelsProbe" \
  --matrix "$A93_RUN/matrix-reproduction/matrix.json" \
  --output "$A93_RUN/continuation-reproduction"
```

These commands refuse existing destinations. The source revision reported by
new case generation is the checkout used for that generation; the frozen
executable's build identity remains `1950ded` and is recorded independently.

For a fresh isolated build, use a fresh experiment root with separate
`source` (this branch) and `baseline-source` (reviewed `99112a1`) worktrees.
Set A93_RUN to that fresh root and export the thread/TMPDIR settings above.
The capture helper from this branch identifies the reviewed source and binary
explicitly. Check available worker budget before each build:

```bash
mkdir -p "$A93_RUN/tmp"
cmake -S "$A93_RUN/baseline-source/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$A93_RUN/build-baseline" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS='-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$A93_RUN/build-baseline" --parallel 8 --target \
  IncNavierStokesSolver RANSKernelProbe RANSConstitutiveProbe \
  RANSAcceptedAWFProbe RANSProfileBranchChecks
python3 -B "$A93_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/capture_a93_baseline.py" \
  --source "$A93_RUN/baseline-source" --build "$A93_RUN/build-baseline" \
  --output "$A93_RUN/baseline"
cmake -S "$A93_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$A93_RUN/build" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS='-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$A93_RUN/build" --parallel 8 --target \
  IncNavierStokesSolver RANSNewModelsProbe RANSKernelProbe RANSConstitutiveProbe \
  RANSAcceptedAWFProbe RANSProfileBranchChecks RANSStressTangentAudit \
  RANSHistoryRejectionChecks RANSWallResolutionProbe
python3 -B "$A93_RUN/source/docs/plans/a93-a95-e94e95-awf-evidence/qualify-linux.py" \
  --root "$A93_RUN"
```

Without a prebuilt SpatialDomains override this builds a private geometry
library too. The recorded incremental build instead uses the private baseline
copy via `-DRANS_PREBUILT_SPATIAL_LIBRARY=$A93_RUN/dependency/libSpatialDomains.so.5.9.0`;
its hash is included in every runtime identity. It never installs into the
shared prefix. The complete fresh-run driver includes the additional A95
cache, reviewed CLS checkpoint, replay and unsupported-replay checks as well
as the standard seven CTest controls. The separate `finish-linux.py` is only
for completing the preserved prototype receipts, not needed for a fresh run.

## Remaining validation gate and scope boundary

Implemented and passed: the three original outer closures, both kappa policies,
both distinct finite-nu walls, model-aware native serial flat-2D one-sweep
integration, full 3D stresses/covariances, independent kernel checks, the twelve
bounded trajectories, representative perturbation/refinement, restart,
rollback, cache and regression checks above.

Implemented but failed a required unexpected test: none remains. Original-model
covariance/ellipticity/domain/dynamic-positivity negative controls retain their
specified adverse outcomes. Coarse polynomial references have measurable
spatial residual and splitting drift; they are not called exact discrete
equilibria.

Next unresolved gate: an asymptotic preservation study before production or
physical-validation claims. The smallest next experiment is another
cutoff-aligned radial/p refinement of the E94E95 natural/accepted case, at the
same bounded T, crossed with one further dt halving. Compare the unanchored
initial momentum/k/epsilon residuals, the refreshed dependent-epsilon trace,
and all three k stage partitions against the existing coarse/fine cases.
This would separate remaining profile representation/quadrature effects from
the measured one-sweep splitting defect. It does not require a long campaign.

Long-time stability, a complete temporal-order campaign, physical/terrain
benchmarks, general curved or true-3D integration, MPI qualification and a
complete nonlinear/accepted-wall Newton method are outside this first task
and have not been attempted. Unsupported execution/action routes reject
explicitly; no evidence here establishes global realizability or positivity.
