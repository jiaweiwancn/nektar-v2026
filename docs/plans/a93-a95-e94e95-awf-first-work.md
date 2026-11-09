# A93 / A95 / E94E95: first implementation and qualification

Status: preparation; no new-model qualification claimed yet.

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
`build`, `runtime-baseline`, `runtime`, `baseline`, and `matrix` are separate
paths beneath this experiment root. The established installation is a read-only
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
The two `a93_a95_e94e95_note_*.py` scripts in the RANS test directory are exact
extractions of the combined note's executable listings. They are mathematical
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
| RANSWallEvolutionNonpressure.cpp | Separate Boolean CLS/non-CLS tangent dispatch; must not default to standard for new enum |
| RANSWallEvolutionResolution.cpp | Source diagnostics currently reconstruct standard epsilon production separately; use actual outputs |
| RANSWallEvolutionEvaluation.cpp | Exact-step evaluation reuse; gradient-dependent constitutive-domain validation |
| RANSWallEvolutionContact.cpp | Wall metadata currently fingerprints the CLS appendix; model-aware identity and covariance checks |
| RANSWallPhysicsSignature.h | Includes patch signature and session controls; extend without changing compatible CLS identity |

The inherited accepted wall already supplies nonconstant ReferenceK, zero
redistribution, frozen prescribed sink, and a flat serial one-sweep path.
It currently requires CLS, CLS-specific appendix metadata and bound inventory.
Extending only RANSVelocityCorrectionScheme is insufficient.

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

## Results and next gate

Workspace/specification audit completed. Linux baseline build and short control
are the current gate; results and exact reproduction commands will be added
before numerical implementation edits. No new-model build, numerical pass,
long-time stability, terrain/physical, 3D/MPI or temporal-order claim is made.
