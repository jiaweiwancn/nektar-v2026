# Segregated CLS/EP-AWF Picard R3: reduced search implementation

**Status: implemented, native numerical validation blocked by the unmounted
Linux data volume. The original M8 convergence problem is unresolved.**

This phase continues the reviewed `1f2b67db0260ccf0e09271cc7ab7d75424c4e6c7`
on `codex/cls-epawf-segregated-picard`. The new implementation is commit
`819c135` (the full producing commit is in the evidence manifest). R1/R2
implementation, reports, evidence, legacy selectors, and the unrelated Mac
files were preserved. The shared Linux checkout and installation were not
modified. Serial, flat, tensor-cubic ground support remains explicit.

No new native M8 solve was run in this phase. The last numerical result is
the **historical R2 original-start result**, not an R3 result:

| Quantity | R2 result | Original requirement |
| --- | ---: | ---: |
| Complete original residual | 2.9088338887423277 | <= 2e-12 |
| Continuous ground minimum k | about 2.14e-16 | arithmetic admissibility |
| Maximum normalized outer-gradient departure | 0.00382640524 | <= 0.001 |
| Minimum outer-epsilon ratio | 0.9744548514 | >= 0.98 |
| Accepted physical steps / final time | 0 / 0 | qualified first stage, then .0075 |

There is no new evidence of reduced stationarity or nonexistence. A failed
subspace model, failed non-k solve, or future exhausted line search must not
be interpreted as such a certificate.

## Direction construction

`RANSWallPicardReduced=On` selects the new path. It is opt-in; the conditional
k method remains a reproducible control. The old optional settling code is
retained, but it does not define a reduced residual for the new method.

The reusable `EquilibrateNonK` operation holds **all independent k
coefficients** fixed. It uses the extracted native velocity predictor and
pressure Schur correction followed by the original epsilon transport and
nonlinear matching-trace solve. Its merit excludes k and includes momentum,
both continuity equations, common gauge, outer epsilon and trace. Accuracy
and update budget are arguments. Its return value includes converged/failed
status, attained block norms, residual, work counts and failure reason.

The input state and physical history are const. Failure returns the incoming
state, never the partially equilibrated iterate. Expected factorization,
domain and line-search failures permit retries. Dimension errors and the
fixed-k invariant throw programming errors separately. The numerical trial
state is stored only for rejected-state diagnostics; it never advances the
native physical history.

Every response and trial requires explicit H convergence. A column is a
one-sided difference of **converged** non-k solutions from the same base
warm start. Each candidate response is checked at h and h/2, with a tighter
non-k tolerance for the latter. Columns exceeding 3% relative residual-action
error or 5% maximum normalized observation-action error are discarded.
Insufficient usable columns trigger tighter accuracy and an expanded model,
not silent use of partial solves. `RANSWallReducedVerify=On` additionally
checks repeated k evaluation, a disturbed warm start, and deliberate zero-
budget failure/rollback. These native checks are implemented but **unrun**.

The response basis includes conditional Newton and rate-Picard seeds,
ground-inward modes, sampled-K and normalization modes on both sides of the
support, independent epsilon-bound witness seeds, spatially distinct signed
qualification witnesses, and the previous displacement. Feasible one-sided
probes are constructed using the physical ground polynomial. Two-pass
rank-revealing orthogonalization rejects dependent seed columns. Failed
models can expand the requested basis by three directions, up to 32.

The reduced action B enters the small direction optimization **before**
choosing the correction. This differs from choosing a conditional Jkk
direction and then optionally settling the other fields. No full all-field
Jacobian is assembled or solved. The verified free-k analytic matrix and
exact alpha mass remain the conditional seed, not a substitute for B.

## Constraints, restoration and actual trials

The small quadratic model uses the original free-k mass-dual residual norm.
The displacement scaling is the Euclidean free-k coefficient norm divided
by reference K. Regularization is 1e-3 times its squared norm plus a 1e-10
coefficient diagonal for nearly dependent response coordinates.

Both signs of every `OuterBoundsR5` observation enter direction selection:
individual gradient components, local outer k, sampled K, actual outer
epsilon, and projected matching trace. Each observation receives the
equilibrated full-state or projected-trace response. It is not a derivative
of the single worst violation. Ground constraints start with all continuous
ground extrema and add cuts at actual candidate extrema. Positivity cuts
also cover the complete current scalar inventory. At moving-ray locations
these constrain the current points; the actual candidate rebuilds and scans
its own moving points separately. This is not a full-element positivity
theorem.

Outside qualification, a minimax slack restoration objective reduces the
largest normalized signed violation. It may increase residual merit. A
filter tracks residual and qualification separately. Inside the envelope,
the model targets reduced residual decrease subject to the signed bounds.
Continuous ground admissibility and native evaluation-domain validity
remain enforced throughout; no clipping, floor, forcing or envelope change
is introduced.

After cut handling, the **combined actual correction** is bounded by the
radius, and all predictions are recomputed. Scaling this constrained model
solution is not claimed to preserve its affine KKT optimality. Logged
`unscaled_affine_kkt` values describe only the small unscaled model, not
nonlinear reduced stationarity.

Every assessed fraction constructs k, converges H to the same required
accuracy, rebuilds dependent quantities, evaluates the original full and
reduced residuals, scans continuous ground extrema and the complete signed
qualification inventory, and compares actual versus predicted progress.
Radius and accuracy adapt to agreement and failures. The log records
response rank, basis identities, active cuts, work, H errors, actual
displacement, model errors, signed observation changes, full-inventory
witnesses, ground minima and original qualification margins. Large model
discrepancy triggers native inner-destruction decomposition. The R8/R9
reports informed this diagnostic; their coupled optimizer was not imported.

Final return requires the original complete residual tolerance and
`CompleteRestrictedAudit`, including covariance and normalization-lambda
checks. Eta continuation was not attempted: no R3 native numerical evidence
exists yet to justify choosing a continuation experiment.

## Checks completed and their scope

Compact Linux compiler commands, source hashes, executable identity,
runtime-library identities, worker allocation and output are in
[`linux-checks.json`](cls-epawf-segregated-picard-r3-evidence/linux-checks.json).
Only small temporary artifacts were placed under
`/tmp/nektar-picard-r3-syntax` and `/tmp/nektar-picard-r3-model` while the
heavy-test volume was unavailable. One compiler worker and one OpenMP/BLAS
thread were used, below the 72-worker limit.

The standalone `RANSReducedModelChecks` test covers unconstrained and binding
solutions, duplicate observations, discovery of newly active constraints,
opposing-witness restoration, a separate lower-bound condition, a known
feedback-dependent correction, the radius on the final combined step,
infeasible starts, exhausted budgets, factorization failures and programming
shape errors. It uses address/undefined sanitizers and bounds assertions on
Linux. It is explicitly labelled **model checks, not M8 convergence**.

The new native sources received syntax-only checks on Linux and Mac. The
Python harness/checker parses successfully; fixture preparation preserves
the original physical XML semantically after removing only the new solver
controls. `git diff --check` passes.

The native post-run checker `test_cls_epawf_reduced_response.py` verifies
repeatability, fixed-k/rollback checks, usable two-scale actions, used-basis
checks, final radius, H accuracy, and native accepted-state/history retention.
It has **not** been executed because it needs a native run. No small-model
check is evidence that sequential H solves converge on this fixture, that
reduced predictions agree there, or that qualification is restored.

## External prerequisite and resumption

The server reboot left `/dev/sda2` unmounted. It is the 3.6 TB NTFS volume
labelled `新加卷`, UUID `5AB8A669B8A642FD`. The required source, original
saved endpoints, and designated heavy-run root are therefore absent.
Tailscale SSH works; `linux-usb` timed out. UDisks mounting requires desktop
authorization, and passwordless sudo is unavailable. The user elected to
mount it manually with their password and request resumption afterward.
The supplied [`mount-volume.sh`](cls-epawf-segregated-picard-r3-evidence/mount-volume.sh)
uses the verified UUID and installed NTFS-3G driver; it formats or repairs
nothing.

After that prerequisite, the intended isolated experiment root is
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard-r3`.
Check the shared and isolated revisions/working trees before synchronizing.
Build only the isolated source, using the installed MPI environment,
`-D_GLIBCXX_ASSERTIONS -fno-fast-math -ffp-contract=off`, and at most eight
compiler workers. Record actual executable/library identities. Do not reuse
a binary merely because this syntax check passed.

The immediate native sequence remains outstanding:

1. Revalidate the shared predictor extraction with
   `test_cls_epawf_picard_verification.py` and relevant existing regressions.
2. Run the new method with `--reduced --reduced-verify` at the original h0 and
   the R2 endpoint packet
   `.../cls-epawf-segregated-picard-r2/reduced-radius/target/readback-input.json`.
   Check H convergence, branch repeatability and two-scale signed responses
   before increasing outer budgets. Exercise both gradient witnesses and
   the epsilon lower bound, and inspect complete-inventory witness migration.
3. Correct numerical issues exposed by those tests and run the unchanged
   original-start target (2165 unknowns, alpha400, t.0025, tolerance2e-12,
   zero prescribed forcing). Retain independent saved-state readback.
4. Only after a qualified first stage, run three accepted native steps to
   .0075, admissible guess perturbations, rejected-state/history checks and
   independent final residual reconstruction.

Thus the controlled response and constrained direction have been implemented,
but their native accuracy and usefulness are not established yet. Original-
start R3 residual/qualification, physical evolution and any mathematical
obstruction assessment remain open; the current blocker is external access
to the data volume, not a numerical nonexistence result.
