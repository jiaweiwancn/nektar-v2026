# R3 durable exact-target entry

The reviewed baseline is `e4c9d44758f08109ba5ea8c2b909790a02b3acc8`.
The R3 additions are opt-in. `Legacy` and `ScalarFeasibleR2` retain their
existing branches and output. R3 uses `RANSWallNonlinearPolicy=ScalarFeasibleR3`
paired with `RANSWallActionPolicy=VerifiedR3Action`; the selected stress remains
`SaturatedR2`, with every physical coefficient, field, source, quadrature and
finite inventory fixed by the R2 contract.

The new EquationSystem selectors are `JacobianBenchmarkR3` and `ExactTargetR3`.
Both retain original h0, time zero and all native history slots. They evaluate
only alpha `400` at time `.0025`, with eta zero and no prescribed body/scalar
forcing. The benchmark assembles one complete current Jacobian and performs
one factorization, then returns without a Newton update or target attempt.
The target performs one coupled search. Neither entry advances native history.
There is no R3 `3200,1600,800` solve sequence.

`RANSWallR3PredictorFile` names a versioned manifest. The preparation helper
binds the exact saved failed SaturatedR2 auxiliary packet and its binary
payloads, extraction ledger and historical trace. It compares protected XML
from the historical Sat session, current initialization and target, permitting
only the declared numerical selectors and maximum-Newton-update choice.
It preserves the original h0 hash
`f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca`.
The native entry rechecks the manifest's state/h0 binaries, current executable,
exact target XML and current native physics signature. It evaluates the
current scalar inventory, constraints, exact target residual, covariance and
fixed-region applicability. Predictor coefficients are not lifted, clipped or
substituted for history. The historical auxiliary remains nonconverged and is
not a physical time sample.

The R3 complete Jacobian combines the new nonpressure spatial action, the
already qualified complete pressure/gauge action, and exact alpha mass on
unconstrained differential rows. Current dependent epsilon-H constraint
derivatives remain in the matrix. There are no complete-residual coordinate
differences in this branch. The action implementation and its selected
one-sided/semismooth branches require independent verification before use;
this protocol does not supply that gate by naming the policy.

The initial implementation retains the declared R2 globalization: scalar
domain screening before residual evaluation, maximum fraction from the same
actual finite inventory and dependent trace action, `.995` only when the bound
is below one, Armijo coefficient `1e-4`, at most 32 backtracks, collapse below
`2^-40`, and stagnation after three accepted relative improvements below
`1e-4`. The coupled tolerance stays `2e-12`. The stage entry allows an explicitly
chosen maximum of 1–16 Newton updates, recorded in the protocol and budget.
No field or residual tolerance is relaxed.

## Durable publication

`RANSWallR3ProgressFile` defaults to `case.r3-progress.jsonl`. R3 never creates
the old large per-solve R2 JSONL record. Each completed internal iterate has
separate immutable state, exact residual and constraint-pair binary files,
with the unchanged history payload shared within the solve. Payloads use
little-endian IEEE754 binary64 and SHA256 identities. A complete immutable
snapshot and an atomic `.last.json` snapshot are published before its completed
event is appended and before a later Jacobian starts. Data files, event writes
and directory publications are fsynced. A timeout can leave an unpublished
event tail; it does not invalidate an already complete atomic snapshot.

Every Jacobian gets an `incomplete` status marker and flushed start event
before assembly. A matrix is published with a SHA256 identity only after the
whole row-major file is complete. A later marker/event identifies completion;
factorization has a separate completion event. A kill during assembly or
factorization therefore leaves the last completed state/residual and the
appropriate incomplete phase. Trial states and evaluated residuals are small
separate binary payloads, with compact classification, domain, fraction and
branch metadata. Accepted updates are saved immediately, including one that
triggers a subsequent stagnation stop.

Profiling records outer solver domain/residual/constraint/norm calls, complete
action assembly, factorization, linear solve and line search. Complete action
metadata separately records wrapper domain/constraint/pressure/mass work and
backend point/ray/kernel/moving-rule work. Assembly time includes its internal
components; line-search time includes its domain/residual calls. These nested
timings must not be added as disjoint totals.

`collect-progress.py` validates all published binary/snapshot/matrix identities,
the original history, event ordering and complete prefix. It preserves a
fully published last snapshot even if its event was interrupted, and exposes
an incomplete Jacobian or missing terminal result explicitly. Native residual
convergence is copied separately from independent accuracy, arithmetic/sign,
constitutive and physical qualification. Internal states are never labelled
physical trajectory samples.

## Bounded execution and focused checks

All compilation and execution use the authoritative Linux host over Tailscale,
external run directories, one MPI rank and the six numerical thread variables
set to one. Compilation uses an explicitly allocated maximum of 16 workers;
the global limit remains 72. `build-epoch.py` creates a fresh isolated Release
build and records exact commands, flags and source/image pins. It does not
replace the shared installation. Runtime fingerprints and strict floating
point/build checks remain separate required gates.

`run-stage.py` requires a current-source/executable-bound independent action
gate. Target mode also requires the completed current-J benchmark receipt,
an explicit seconds budget and a frozen budget note. There is no default
four-hour target budget or automatic auxiliary/retry campaign. The root
supervisor decides whether the measured cost supports a complete bounded
target; an impractical benchmark requires a focused numerical decision first.
The runner uses disconnect-safe progress output, always waits for a stopped
child before taking its after-runtime fingerprint, and collects the durable
state even after a timeout. Fresh outputs prevent silent retries.

The tiny standalone `durable-probe.cpp` tests SHA256 known vectors and actual
payload publication without a constitutive/PDE calculation. On an allocated
Linux worker:

```sh
c++ -std=c++17 -Wall -Wextra -Werror -I "$SOURCE/solvers" \
  "$TOOLS/durable-probe.cpp" -o "$RUN/durable-probe"
python3 -B "$TOOLS/test-durable.py" \
  --probe "$RUN/durable-probe" --output "$RUN/durable-tests-001"
```

The focused test closes the metadata supervisor's stdout before its first
write, times out a child held immediately after its incomplete-J marker,
waits for that child, and verifies that collection still returns a complete
last state/residual. It also checks unpublished event tails and changed binary
payload rejection. No Nektar solver or scientific operator runs in this test.

An accepted target still requires independent reconstruction at the unchanged
coupled tolerance, scalar/covariance/fixed-region and actual trace/constraint
gates, and the specified source, energy, continuity and traction accounting.
Only after those pass may the separately authorized actual native startup
step verify integration bookkeeping. The stage/benchmark entries themselves
provide no such history advance. All broader M8 and physical terrain gates
remain open.
