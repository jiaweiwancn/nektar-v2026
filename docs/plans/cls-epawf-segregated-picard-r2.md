# Segregated CLS/EP-AWF Picard R2

The predictor verification defect was repaired and revalidated. The
restricted nonlinear k correction reduced the full residual of the
unchanged M8 problem from the previous 3.1777 plateau to 2.90883389.
A qualifying first stage remains unresolved. All failed runs retained
the initialized state, four native history slots, and time labels.

Work continued on `codex/cls-epawf-segregated-picard` from the reviewed
`866f24e69301f99b4c7e2e760ff2039220fd77fc`. The separate R8 worktree on
`codex/cls-epawf-revised` was preserved. The producing revisions are
`3cd7e2c44c6694fd558947d85e25cf5d86b33abd` for the verification repair,
`43ed2d9f5849be80d748667d605dabf655501f2b` for the conditional k solve,
`d5a013436550b6cd128696e2bf44cf75d8fe02d4` for fixed-k feedback
measurements, and `3f942fdbe864cac673bbf37e3f3f6a78134c1453` for the
reduced trial. The original Picard and legacy selectors retain their
defaults.

## Verification repair

The predictor contains 980 velocity increments. The old verification
applied it to 1,065 flow rows, including 85 pressure and gauge rows.
The predictor diagnostic now applies its increment only to the velocity
rows. The generic block update requires an exact shape match and uses
bounds-checked state access.

`test_cls_epawf_picard_verification.py` used a Linux build with
`-D_GLIBCXX_ASSERTIONS`. Two accepted sweep hashes and full norm vectors
were bitwise identical with verification on and off. Both runs reached
their deliberate two-sweep stop with the initialized physical state
intact. The corrected predictor's inner and outer continuity norms were
(8.54316\times 10^{-4}) and (2.53250\times 10^{-4}). Its gauge norm
was (1.02830\times 10^{-16}). The corrected-flow values were
((6.87347\times 10^{-16}, 8.68001\times 10^{-16}, 6.21722\times
10^{-17})). The Schur and direct saddle solutions differed by at most
(3.61\times 10^{-16}) in coefficient units. The nonzero predictor
gauge previously reported is superseded by this checked observation.
The unverified historical ground-projection endpoint provides a separate
record of the k plateau.

## Restricted k and the saved endpoint

The original M8 input has 2,165 unknowns, physical alpha 400, time .0025,
zero prescribed forcing, and tolerance (2\times 10^{-12}). All runs
used the original mesh, reference, boundary, closure, quadrature, and
initialized history. The saved higher-relaxation endpoint is supplied
through `RANSWallPicardGuessFile`. The native callback checks its shape,
history values, alpha, time and external data. The packet replaces only
the nonlinear guess.

The analytic assembler scatters only the 504 free k rows and columns.
It adds the physical `alpha Mff` exactly. It retains the derivatives of
sampled K, moving cutoffs, quadrature weights and test evaluations, stress
and production, limiter and complete-ray normalization, and transition
sensors including the mechanical-pressure contribution from k. The
conditional trial calls the original native residual at the changed k
and retains all independent velocity, pressure, epsilon and matching
trace coefficients. The epsilon block then resolves its own trace
constraint. No full all-field matrix is assembled by this path.

Six independent directions were checked against actual conditional
residual changes at both h0 and the saved endpoint. The check uses a
range from (10^{-4}) to (10^{-8}), central differences where both
signs are admissible, and one-sided differences near contact. It
includes a ground-inward direction, a sampled-K direction, a near-ground
normalization-ray direction and three near-tangent directions. It
screens each trial against the actual continuous polynomial and scalar
domain before evaluating the constitutive model. The maximum of the best
relative action errors was (9.54\times 10^{-11}) at h0 and
(9.79\times 10^{-7}) at the endpoint. The endpoint's frozen and true
near-tangent actions differ by as much as 98% of the true action norm.
The sampled-K direction changes moving panels and the normalization-ray
direction changes the complete-ray budget at the endpoint. The base ray
cap remains inactive.

The old projected endpoint has an original k norm of 3.1777167530. The
raw frozen affine solution has a continuous ground minimum of
`-0.0142171`. The model is evaluated at the old state and the affine
candidate is inspected as a polynomial. At k relaxation .75 the largest
frozen coefficient increment shrinks from `0.190477` before projection
to `3.85e-6` after it. Seven cut additions produce a stationarity
defect of `1.12e-13`, zero measured active primal defect and
`1.34e-16` complementarity. The actual k norm drops by
`1.53e-9` on that trial. The frozen projected update is nearly stationary
with a nonzero native residual.

Using the restricted Jacobian in the same mass-dual projection gives a
predicted k norm of 2.94504 and an actual conditional norm of 3.02829.
At one eighth of the step these values are 3.13482 and 3.13437. The true
conditional model therefore has feasible descent at the saved state. The
old frozen projection does not establish stationarity for the nonlinear
conditional problem. The heavy diagnostics retain the free mass, both
matrices, both increments, active rows and multipliers, and all continuous
ground extrema.

## Sequential feedback

The optional safeguarded k correction bounds the raw Newton proposal,
projects against continuous ground constraints, and uses the exact
conditional k norm for its line search. A sweep can fall back to the
original implicit-rate direction if the conditional trial cannot be
accepted. The original full stage residual and qualification remain the
only physical acceptance gates.

The non-k diagnostic holds k fixed and equilibrates momentum, continuity,
gauge and epsilon including the matching trace using separate native
flow and epsilon solves. Its merit excludes k. It always stops before
physical history advancement. With the old projected endpoint and
perturbations of (10^{-4}) and (10^{-5}), the non-k tolerance was
(2\times10^{-9}). The base equilibration changed k's norm by
approximately (2.2\times10^{-9}). The equilibrated action differed
from the conditional action by about 2.2% for the inward ground direction,
1.8% for the sampled-K direction and 12% for the tested near-tangent
direction. The two perturbation sizes give the same sign in the signed
residual projections. These are targeted reduced responses, with the
other governing equations equilibrated separately.

The gradient-restoring sampled-K direction was also tested at amplitudes
.005 and .015 with H re-equilibrated. At .005 the maximum gradient
departure only falls to .00278565 as the witness moves to the opposite
side of the support at ((-2.604649609,0,2)). The sampled K ratio reaches
1.00974 and the k residual rises to 4.41291. At .015 the new gradient
witness is still outside the bound at .00278622, and the sampled K and
epsilon ratios also exceed their upper bounds. The k residual rises to
9.57532. This direction cannot restore qualification in the tested
range, and the witness migration limits extrapolation of the local
derivative measured at x=1.60465.

The reduced trial optionally re-equilibrates flow and epsilon at fixed
trial k when a full sweep would otherwise be rejected. It accepts inner
updates using the non-k merit and resubmits the resulting state to the
complete residual check. At the old warm start, this can remove the
momentum defect while increasing k. For example, on one rejected sweep,
the first momentum norm fell from 0.31729 to (6.66\times10^{-5}),
while the k norm rose from 3.07296 to 3.52100. The feedback can therefore
undo a useful conditional correction. Sequential equilibration of that
trial is insufficient on its own. This test does not imply that every
coupled direction near the ground constraint is blocked.

## Linux reproduction

Run storage is
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard-r2`.
The shared checkout remained at `f4b70baec10016dd098005bc5a5507b87e4a57bf`.
Only the isolated `source` and `checked-build` were changed. The
existing MPI-enabled installation supplied the dependency libraries.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
run='/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard-r2'
old='/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard'
src="$run/source"
build="$run/checked-build"
prefix=/home/jay/.local/opt/nektar-v2026-mpi
export TMPDIR="$run/tmp"
harness="$src/solvers/IncNavierStokesSolver/Tests/RANS"
cmake -S "$harness" -B "$build" -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$prefix" "-DNektar++_DIR=$prefix/lib64/nektar++/cmake" \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$build" --parallel 8 --target IncNavierStokesSolver
python3 "$harness/test_cls_epawf_picard_verification.py" \
  "$build/IncNavierStokesSolver" "$run/verification-replay"
python3 "$harness/run_cls_epawf_picard.py" --solver "$build/IncNavierStokesSolver" \
  --resolution-probe "$old/build/RANSWallResolutionProbe" --ground-projection \
  --k-newton --equilibrate-after-k --k .75 --epsilon 1 --k-radius .40 \
  --sweeps 90 --timeout 870 --output "$run/reduced-radius-replay"
```

The conditional runs used `--k-diagnostic --ground-projection --k-newton
--k .75 --epsilon 1`, with 90 sweeps for h0 and 8 for the saved packet.
The feedback script `measure_cls_epawf_picard_feedback.py` uses the saved
packet and its k diagnostic, with `--workers 3`. The reduced comparisons
use `--equilibrate-after-k`. The default k proposal radius is .15 times
reference k, with .40 tested separately. Every solve uses one rank and
one thread per runtime. Builds use eight jobs. The focused regression
tests use `ctest --parallel 1`. Per-run manifests and `ldd` files record
the images and resolved dependencies.

All five native regressions passed. The Linux constitutive probe passed
28,434 checks. A separate legacy R3 action regression reproduced all
bytes of the full 2,165 by 2,165 stage matrix at the identical h0 input,
using the pre-change R7 probe as the reference. Its SHA256 is
`0ee7bd325c28d36fe5bbf5636a0edd6f82d72295e01bda996714decd754c0f1e`.

## Stage results

| Run | Initial guess | Completed sweeps | Original max norm | Outer gradient departure | Solver seconds | Stop |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Conditional k | h0 | 90 | 2.90925243 | .003826398 | 334.77 | Sweep limit |
| Conditional k | Saved endpoint | 2 | 2.99944187 | .003048721 | 83.67 | No admissible decrease in the allocated trials |
| Reduced trial | h0 | 110 | 2.90883705 | .003826405 | 385.52 | Sweep limit |
| Reduced trial | Saved endpoint | 2 | 2.99944187 | .003048721 | 74.03 | No admissible decrease in the allocated trials |
| Reduced trial, radius .40 | h0 | 39 | 2.90883389 | .003826405 | 167.78 | No admissible decrease in the allocated trials |

Every fresh readback reproduced the logged terminal full norm vector.
The best original-start vector in native order was
((2.08144\times10^{-4},1.33589\times10^{-6},2.90883389,
4.74807\times10^{-8},6.59346\times10^{-16},7.81750\times10^{-16},
5.59042\times10^{-18},4.68252\times10^{-12})).
Continuous ground k had a minimum of (2.14\times10^{-16}).
The independent qualification audit found a gradient witness at
((1.666666667,0,2)) on the outer side of element 34 and a minimum
outer epsilon ratio of .974455, below .98. All five endpoints failed
qualification.

At the best original-start state the true k projection used ten cut
additions and had stationarity defect (2.02\times10^{-12}).
Its largest projected coefficient increment was (5.75\times10^{-5}).
The predicted and actual conditional reductions were both
approximately (4.11\times10^{-8}). At the earlier 90-sweep
conditional endpoint the corresponding actual reduction was still
(3.94\times10^{-4}). The subsequent runs therefore reached a much
flatter conditional plateau than the earlier endpoint. The affine
subproblem and actual trial agree closely there. This is local
evidence about the tested ground-constrained update.

The first accepted conditional sweep from h0 exceeded the original
outer-gradient bound. The departure grew from `0.000500195` to
`0.00237553`, with the witness at ((1.604649609,0,2)) on the outer
side of element 34. The bound is .001. The k-only action leaves the
velocity gradient fixed. In the sequential response test, the
sampled-K direction changes the re-equilibrated witness by
`-1.31e-5` at amplitude (10^{-4}), while the base stays at
`0.00288904`. A local directional response does not establish
compatibility or incompatibility of the entire qualification envelope
with a root.

Zero physical steps were accepted and the final native time was 0 in
every run. The native three-step trajectory and perturbed-guess runs
remain conditional on a qualified first stage. The next numerical test
is a ground- and qualification-aware reduced search with independent
directions at both gradient witnesses. The measured response to a
single sampled-K direction and the settled trial results motivate
checking both constraints together before attempting a longer continuation.

The [compact validation evidence](cls-epawf-segregated-picard-r2-evidence/linux-results.json)
contains sweep histories, sampled qualification, projection and action
checks, and the producing executable hashes. Its
[collector](cls-epawf-segregated-picard-r2-evidence/collect-linux.py)
validates the paired state histories and terminal readbacks. The Linux
run root retains the raw matrices, full active rows, native solver logs,
readbacks, terminal-audit inputs, regression log and dependency lists.
