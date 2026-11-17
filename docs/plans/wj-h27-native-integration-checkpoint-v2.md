# H27 approved-v2 native integration checkpoint

This is a substantive implementation checkpoint, not H27-B/C/D closure.
`TASK3_READY=NO`. The user approved the isotropic lower-wall amendment in
`wj-h27-model-amendment-proposal.md`. No further physical-model amendment is
included here.

The normal `IncNavierStokesSolver` factory now registers the opt-in
`WJNativeCoupledV2` system. It advances current mixed velocity/modified pressure,
full-domain k and outer-only epsilon, with dependent interface dissipation and
cutoff, conservative scalar transport, complete sampled AWF forcing,
component-dependent physical boundaries, actual pressure-gauge determination,
and BDF2/EX2 history with BE startup. Current constitutive, scalar and interface
data are refreshed to convergence. Historical default systems are unchanged.
Full-state checkpoint/restart and fresh-field transfer have distinct identities;
steady candidates are not accepted as qualified production references.

## Tested identities and execution

The Linux source-v9 capsule was assembled from parent
`ae134d6787553bf7036af86965313422ea11f595` plus the 90-file overlay recorded in
`wj-h27-linux/native-integration-v2b/source-v9.json`. Those file hashes were
checked against this Mac worktree before committing. The commit containing this
document therefore binds the tested source contents; the Linux build did not
originate from an already committed revision. Runtime binary, loader, shared
libraries, compiler commands, native inputs and individual receipts are recorded
in the same evidence directory. Large raw fields remain in the isolated Linux
run namespace `/home/jay/build/wj-h27-coupling-r1-20261117`.

All builds and numerical conclusions below are Linux results. Build parallelism
was four. Solver runs were serial, one MPI rank, with OMP, OpenBLAS, MKL, BLIS,
VECLIB and NumExpr threads set to one. The shared-host watchdog enforced at most
72 actual workers. Source-v9 build took 208.716 s.

## Independent results and remaining gates

* Mixed closed/slip/outlet action tests pass, including the actual gauge and
  transported-Q3-k action. Independent small-action errors are below 1e-10.
  Scalar operators, conservative advection, changing diffusivity/reaction and
  compatible essential traces pass their independent checks.
* The reference table, full inlet and cutoff are independently qualified over
  `[0.9,1.1] K_in`; the table SHA256 is
  `b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d`.
  Out-of-range states reject rather than clamp.
* Continuous coupled MMS passes with current full WJ, evolving scalars,
  pressure/AWF/known-k actions and separately activated body-force terms. Both
  Historical and PolynomialScalarsV1 variants have independent strong checks.
* A genuinely perturbed polynomial steady solve takes four Newton iterations,
  37 Krylov iterations and 47 Jacobian-vector products. Its final normalized
  momentum/k/epsilon residuals are 1.0490e-11 / 3.5373e-11 / 2.5463e-11 in the
  independent current-equation evaluator. The independently constructed,
  exactly representable polynomial reference passes the 1e-8 field gate and
  split-quadrature equation checks. The earlier optional undamped block-Picard
  reference algorithm failed its range guard; that failed outcome is retained,
  and is not represented as successful verification.
* The preliminary temporal study measures U orders 2.113/2.053, k
  2.015/2.011, epsilon 2.019/2.013, and pressure 1.991/1.987. Its U mesh
  separation is 13.27% of the finest temporal error, exceeding the prospective
  10% requirement. The study therefore remains NOT_MET and needs a finer base
  mesh. Quadrature and tolerance separation pass. Physical-gradient integrated
  budgets improve on refinement; arithmetic reaction accounts are separate.
* On the physical 700-cell case the first BE step, dt=0.005, passes current
  native and independent equations. Independent momentum/k/epsilon residuals
  are 9.0367e-11 / 6.48e-13 / 2.05e-12; continuity is 6.65e-16. Full principal
  margin divided by current mu is 0.004046 and covariance/k margin is 0.2814.
  AA depth four needs 44 nonlinear iterations and 46.485 s; the unaccelerated
  control needs 79 iterations and 45.979 s. This is evidence of admissibility
  and complete coupling, not an affordable production forecast.
* The bounded physical steady diagnostic reduces the normalized residual from
  0.23487 to 0.0030378 in three accepted Newton updates, but reaches its 240 s
  native deadline. Its candidate is explicitly unqualified and accepted history
  is preserved. No developed reference or physical-development envelope is yet
  qualified. The fixed envelope and numerical thresholds are unchanged.
* Full time-history/EX cache actions agree independently at arithmetic scale.
  Max-iteration and incompatible-viscosity negative outcomes are retained.
  Actual cold restart, three recovery rotations, terminal scientific retention,
  longer baseline/startup/perturbation pilots, refined physical comparisons and
  complete cost gates remain open.

## Resource and continuation status

The published OPEN account snapshot totals 3998.5351 aggregate operation
seconds, including 703.6786 s carried from the closed original H27 feasibility
account and 120 s of bootstrap/control. Failed operations are included. This is
one continuing 43200 s allowance with 3600 s protected for evidence, not a new
allowance. The first publication attempt failed due to a missing summary
directory; its receipt is preserved alongside the successful create-only
publication. The evidence archive was 1480634 bytes, SHA256
`9b122fc69fb9fc384214cf6b5a34122eb876ef5eeb565ac1a205bbd880990c7f`.

Continue with equivalent positive directional stress splitting, direct current
flux residual actions, and reducing redundant principal audits in uncommitted
Newton finite-difference probes. Preserve full current-state validation on
accepted nonlinear states and final steps, qualify against the assembled route,
then complete the finer MMS study and physical steady reference before admitting
the longer physical pilot matrix. Full Task3 execution remains outside H27.
