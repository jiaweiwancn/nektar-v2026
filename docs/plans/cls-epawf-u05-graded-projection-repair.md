# H=.5 graded projection repair

2026-11-08, branch `codex/cls-epawf-pimple-one-sweep`. The inspected local and
remote starting commit was `a0ccc1c7febc3ed18efb5f4e5fc6d1116778f22f`.
Implementation and native validation use
`faef649474071996a621028650a4c88fbb8b6900`.
Compact [native evidence](cls-epawf-u05-graded-projection-repair-evidence/native.json)
and the [artifact index](cls-epawf-u05-graded-projection-repair-evidence/artifact-index.json)
bind the immutable Linux inputs, results, runtime closures and dated status.

**The saved attempt-65 candidate is repaired and hard-audited with the existing
900-second allowance. Native restart reproduced the exact original raw system,
accepted step 65, and continued past the former failed time.** The original
terminal experiment remains unchanged. A separately identified repaired
campaign is supervised toward T=.05, .1 and later integer-step targets.

The repair resolves the accumulated finite QP before invoking moving physical
separation again. It changes the numerical optimization path, not the physical
floor, reserve, metric, equations or V3 coverage requirements. It does not make
the unforced Upstream case an equilibrium or establish temporal accuracy.

## Failure and focused change

The exact saved raw candidate, pre-k state, actual callback history, external
data and accepted checkpoint were used. The replay constructs no new h0, calls
no Lift and performs zero governing solves. It reconstructs **PicardKMass M_ff**
on the saved free rows. The archived matrix is transport **A_k**, not M_ff.

| Immutable input | SHA256 |
|---|---|
| Accepted step-64 checkpoint | `4a953092ce4e85c6a6fec9c172ce5eee83bf820c8b0600762637afd106a9abe0` |
| Attempt-65 packet | `54740a789b17757d5bc74d42d47018b1a7267056d93fda0293ec9c65215b5222` |
| Raw coupled candidate | `17a7a306373898d129f06e7c7a0815ce99608d708fb02d3f07f2d78f983ae659` |
| Actual callback history | `b7dab6b2f21a94089bae47f9ebda24b91ec90813f8e0a32ccea7384c0f0d055b` |
| Frozen transport A_k | `2f262970f272bfc477ea72a229ed48955db7a06cf8aca237b5d0a222c0396e9f` |

The frozen original engine reproduced the wall-time failure in 901.72 seconds:
759.62 seconds rebuilding inventories, versus 18.57 in projection algebra.
It exhausted neither 2000 work units per Add nor 20000 total work units.
Its last working set had tiny stationarity, but the accumulated finite pool
still had normalized primal violation 2.79e-13 and physical floor defect
1.05e-11. Thus a certificate for its current working set was insufficient.
The earlier archived terminal has different stop counts because wall time,
rather than a deterministic iteration count, triggers termination.

`RANSWallPimpleKBoundSolver=WorkingSet|RetainedPool` selects the implementation;
**WorkingSet remains the default**. `RetainedPool` reports version
`retained-pool-v1` and adds three operations to the existing mass-whitened,
normalized, twice-reorthogonalized working-set solver:

1. Reuse only exactly identical physical rows **and** right-hand sides.
   Distinct nearly parallel inequalities remain in the pool.
2. After a working-set solve, scan every accumulated inequality and re-enter
   the most violated inactive cut. Finish this finite QP before returning to
   state-dependent inventory reconstruction. Per-Add and total work limits
   continue to count the internal re-entry work; they are not reset.
3. Preserve the latest completed primal with its corresponding solved dual
   if a later inner re-entry exhausts its budget. Unfinished working data and
   incomplete full-pool feasibility remain explicitly labelled.

Every subproblem minimizes half delta_k^T M_ff delta_k **relative to the original
raw candidate**, with prescribed coefficients fixed. All current and retained
physical cuts remain subject to their original floor/reserve requirements.
The algebraic pool tolerance is 64 machine epsilons times its whitened scale;
it does not replace or relax the final physical floor, complete inventory,
continuous-ground or constitutive checks.

## Exact-candidate Linux results

| Quantity | Frozen original replay | RetainedPool replay |
|---|---:|---:|
| outcome | wall-time rejection | repaired and hard-audited |
| correction wall seconds | 901.72 | 196.85 |
| algebra / inventory seconds | 18.57 / 759.62 | 7.16 / 160.19 |
| physical separation passes | 1297 | 492 |
| inventory rebuilds | 569 | 182 |
| completed active-set work | 2551 | 1444 |
| KKT solves | 2551 | 1423 |
| finite-pool scans / internal re-entries | not performed | 743 / 252 |
| full-pool normalized primal violation | 2.79e-13 | 2.17e-19 |
| physical floor / continuous-ground floor defect | 1.05e-11 / 0 | 0 / 0 |
| independently checked stationarity | 3.52e-19, working set only | 4.37e-19 |
| complementarity | see evidence | 1.95e-22 |
| dual feasibility defect | 0 | 0 |
| physical mass objective | 9.5663264683e-7 | 9.5657766933e-7 |
| primal-dual gap, repaired full pool | — | roundoff |

The independent repaired stationarity divided by the whitened correction norm
is 3.16e-16; the independently formed physical primal-dual gap is -8.47e-22
(roundoff at objective 9.57e-7). One mass factorization and 1234 factor
applications are separately recorded. Neither the working-set optimality
numbers nor the finite-pool tolerance waive the full physical checks.

The relative floor remains 1e-8, physical floor 5.1050801306e-9, arithmetic
reserve 4.6430433309e-13. All four limits remain **2000 per Add, 20000 total,
2000 physical separations, 900 seconds**. The complete hard audit includes the
existing scalar/sampler inventory, covariance and normalization conditions;
it is not merely “no exception occurred.” Narrow R2 qualification is reported
separately and remains unsatisfied.

Non-k and prescribed k coefficients are bitwise unchanged by the corrector;
shared CG values remain shared. The repaired state passes static, pre-k,
current/prospective V3 and retained-cut evaluations, continuous ground and
the final constitutive audit. This is the existing finite-inventory guarantee,
not positivity throughout every volume polynomial.

The raw linear residual is 1.6653e-16. After correction, the frozen k residual
is **.14669247** in coefficient maximum norm. Independent accounting verifies

    r_k + A_k (Delta_raw + delta_bound)
      = A_k delta_bound + (r_k + A_k Delta_raw).

The added frozen residual has physical mass-dual norm 10.00132, whereas the
distinct storage reaction M_ff delta_bound/dt has norm 2.21307. The signed
physical correction integral is .00020583115. The corrected field is not
described as a roundoff solution of the original frozen transport system.

### Cut-path dependence is retained as a limitation

Resolving each complete finite pool changes subsequent moving-cut discovery.
The legacy failed replay accumulated 1297 cuts; the successful replay uses a
different, smaller pool. Their returned corrections differ by 7.9004e-6 in
physical mass norm. The new vector is **not** feasible for every transient cut
generated by the separate failed legacy search: the cross-pool normalized
violation is 6.4233e-7. Conversely the old unfinished vector does not satisfy
the successful replay's pool.

No cut is discarded within the repaired call. Each generated pool is retained
and checked, including inactive inequalities. The failed call's internal
optimizer state is not native checkpoint history and is not imported into the
new stateless call. This is an explicitly versioned algorithmic repair, not a
claim of bitwise equivalence to the old incomplete correction or of an
algorithm-independent global volume-positivity projection. Both executions
use the same raw center, mass metric and declared physical constraint policy.

## Native verification and continuation

The strict native reader loads an unmodified copy of checkpoint 64 with the
same physical signature, mesh, H=.5, z_m=.5, ColumnAverage K, IMEX1 dt=.000625,
unit relaxation, SaturatedR2/capped destruction, zero eta and forcing. The
stateless finite-QP implementation is a separately recorded numerical control;
the floor, metric, V3 inventory, sampler, spaces and physical parameters retain
their strict checkpoint identities. Campaign resume additionally checks the
recorded backend. No checkpoint byte/signature or native history is edited.

The native startup bound is inactive and preserves the coupled state and all
native histories. Actual native attempt 65 exactly reproduces the archived
committed, post-flow, pre-k, raw, history, free rows, transport residual,
increment and A_k bytes. Its corrected state is bitwise identical to the
successful saved-candidate replay. The accepted step-65 checkpoint SHA256 is
`d282f3e5812f363e22085d134b61bccf5a55bdfb7e376717e591256294de2830`.
The native verification costs 253.75 wall seconds, including startup,
governing solves, audits and output.

Each accepted step still performs exactly **one vector momentum predictor,
one pressure/gauge solve, one outer-epsilon solve, one k transport solve and
one logical final bound call**. There is one implicit callback and zero extra
governing solves. Pool scans, re-entry, mass-factor applications and KKT solves
are separately counted. No extra outer sweep, field damping, source change or
epsilon-trace restoration is introduced.

Focused checks:

- Linux CTest passed RANSKBoundQPChecks and RANSHistoryRejectionChecks (2/2).
  Existing no-op, consistent-mass, dependent, between-node ground and inactive
  constraint cases remain; new tests cover full-pool re-entry, exact reuse,
  preservation of near-parallel inequalities and budget exhaustion.
- A forced replay-only stop at total work 75, during re-entry in Add 57,
  retains the latest coherent completed primal/dual. Its returned-state mass
  objective agrees exactly with the saved whitened primal. Full-pool status is
  incomplete, no step is accepted, and zero governing solves occur.
- Earlier active step 9 is bitwise unchanged against the pinned legacy solver
  using WorkingSet. RetainedPool also reproduces that earlier accepted state
  exactly. This does not imply identical late cut-generation paths.
- Native step-65 identity and startup/history preservation passed as above.
- Active corrected checkpoint/restart is bitwise identical to an uninterrupted
  two-step run through step 66. A forced step-66 rejection after an activated
  correction preserves the accepted state, time and every native history
  exactly. Startup at corrected checkpoint 65 is also an exact no-op.

The production continuation is **separate**:

```
<run-root>/repaired-campaign/trajectory
```

It starts from accepted step 64, preserves the old terminal experiment, and
has accepted step 69 at **T=.043125**, beyond the failed time .040625.
The checkpoint SHA256 is
`13f60231facc84786015ae485e0fd9a32f08941971f5dcd0d3f375e5d2c13bd6`.
The dated compact evidence gives the subsequently observed frontier and
checkpoint. Targets are 65,80,160,800,1600,3200,8000 steps
(T=.040625,.05,.1,.5,1,2,5), with four-step blocks and explicit 7200-second
native block timeouts. The prefix correction budget .005024866793258327 is
carried once, not re-added at every restart. Numerical rejection remains
terminal; the supervisor only recovers infrastructure interruptions.

At the 2026-11-08 02:08 China-time capture, the new driver, supervisor and
native-child identities and inherited campaign lock were valid. Step 65 and
the next four-step block completed successfully. The next automatic milestone
is step 80 (T=.05), followed by step 160 (T=.1); those targets are not yet
reported complete. Supervisor PID 175843 and driver PID 175844 are recorded
with Linux process start tokens; native PIDs change at block boundaries.

| Accepted step | T | Bound seconds | Correction/raw-step mass-norm ratio | Signed correction integral |
|---:|---:|---:|---:|---:|
| 65 | .040625 | 207.46 | .5070 | .000205831 |
| 66 | .041250 | 224.70 | .5107 | .000208312 |
| 67 | .041875 | 169.80 | .5141 | .000210794 |
| 68 | .042500 | 587.32 | .5174 | .000212602 |
| 69 | .043125 | 71.22 | .5205 | .000214388 |

All five corrections pass the retained/current finite inventory and continuous
ground checks. Correction remains persistent and the cost is uneven; this is
a demonstrated repair of attempt 65, not a guarantee that later attempts will
fit the same bounded work allowance. The completed four-step block took
1220.79 seconds (305 seconds per accepted step including all work). At that
recent average, reaching step 80 would take about one more hour, subject to
changing projection cost. The supervisor remains durable and checkpointed.

The cumulative signed correction integral is .00607679439210, including the
prefix once; .00105192759884 was supplied in steps 65--69. These are physical
integrals, whereas the .51--.52 ratios above compare mass norms of increments,
not fractions of the total turbulent energy. Positive and negative
redistribution and fixed z<=.01/.1 bands are retained in the compact evidence.

At T=.043125, mean wall stress is .0173379. Original refreshed alpha=0
dimensional mass-dual residual norms (u,w,k,outer-epsilon) are
(1.94349,.211547,10.69028,.198592), and the refreshed matching defect is
.00112618. The solved epsilon trace was retained during this readback.
These are materially nonzero residuals. The short repaired prefix remains
unresolved transient evolution with substantial bounding, and supplies no
steady-convergence or temporal-accuracy claim. Narrow reference qualification
is false while the existing hard numerical checks pass.

The independent Upstream one-element and Uniform G/H campaigns retain their
existing runtime closures and physical policies. Their dated process,
checkpoint and executable/library checks are included with the new campaign's
status. No running executable or library was rebuilt or replaced in place.

| Independent experiment at the same capture | Accepted step | T | Status |
|---|---:|---:|---|
| Original U05-graded terminal | 64 | .04 | numerically rejected; preserved |
| Existing Upstream U05-one | 446 | .278750 | running |
| Existing Uniform G | 1271 | .794375 | running |
| Existing Uniform H | 1623 | 1.014375 | running |

Uniform G's cumulative bound integral is zero at this snapshot. Uniform H has
activated bounding; its current cumulative budget is in the dated evidence.
They retain their original WorkingSet runtime and are
independent controls. Their progress is not evidence for the repaired graded
case's later-time behavior.

This repair establishes numerical completion of the captured projection and
native passage beyond its failure. It does not establish constitutive
qualification, accuracy of persistent correction, original-equation steady
convergence or a fully converged M8 stage. The physical Upstream depletion
mechanism remains unchanged.

## Reproduction and artifacts

Heavy artifacts are under:

```
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-u05-graded-projection-repair-20261108
```

All source edits/Git work used the existing Mac worktree. Linux validation
used the documented Tailscale fallback after `linux-usb` timed out. The shared
Linux checkout, dirty at revision `f4b70bae`, was left untouched. A new source
clone/build and private frozen executable/library closures were used. Initial
build parallelism was four, incremental builds two; numerical runs use one
rank and all computational thread controls equal to one. Concurrent workloads
were inspected and stayed within the 72-worker allocation.

Frozen solver SHA256:
`c016c85792ee2d41a72ade53ed88afabff93623c434d01bf701f1db168dd6525`.
Frozen probe SHA256:
`75d2b8979b2a4f25232f68c367f5a7eed1ed4b1a2f61821119efaa3fd5e90125`.
The artifact index includes their complete linked-library identities and
verification of compiled source contents against commit faef649.

Let R be this run root, O be the preceding
`cls-epawf-upstream-branches-20261107` root, and S the tests directory below
R/source. Use distinct output names when reproducing:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR="$R/tmp"
S="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
cmake -S "$S" -B "$R/build" -DRANS_BUILD_SOLVER=ON \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS='-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$R/build" --parallel 4 --target IncNavierStokesSolver \
  RANSWallResolutionProbe RANSKBoundQPChecks RANSHistoryRejectionChecks
ctest --test-dir "$R/build" -j1 --output-on-failure \
  -R 'RANS(KBoundQP|HistoryRejection)Checks'

python3 -B "$S/replay_cls_epawf_graded_timeout.py" \
  --archive "$O/terminal-U05-graded" \
  --probe "$R/pool-v1/probe/RANSWallResolutionProbe" \
  --output "$R/reproduce-replay" --finite-solver RetainedPool \
  --per-add 2000 --total 20000 --separations 2000 --seconds 900 --full-trace
python3 -B "$S/test_cls_epawf_graded_repair.py" \
  --solver "$R/pool-v1/solver/IncNavierStokesSolver" --old-root "$O" \
  --replay "$R/replay-pool-v1" --output "$R/reproduce-first" --phase first
# Repeat with --phase transactions or earlier in separate output directories.
python3 -B "$S/analyse_cls_epawf_graded_repair.py" \
  --root "$R" --old-root "$O" --output "$R/reproduce-analysis.json"
```

The durable supervisor's recorded command uses `--upstream-case U05-graded
--wall-policy Upstream --inventory ConstitutiveV3 --bound-solver RetainedPool`,
the verified `restart64` copy, original prefix report/ledger, and all four
explicit work controls. The active campaign owns an inherited lock; its exact
PID/start/executable/cwd and command are in the dated status evidence.

`replay-original` and `replay-pool-v1` retain full traces, raw/corrected vectors,
actual M_ff and normalized pools. `replay-budget-stop` retains the coherent
intentional failure. `native-{first,earlier,transactions}` contain independent
verification trajectories; these do not add solves to the production path.
`repaired-campaign` retains durable native checkpoints and the idempotent
accepted-step ledger. The compact evidence/index records paths and hashes;
large fields, matrices and logs remain on Linux.
