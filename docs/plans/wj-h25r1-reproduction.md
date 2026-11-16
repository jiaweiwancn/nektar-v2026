# H25R1 executed recipes and evidence

The selected configuration is an **unperturbed flat control**, not an admitted
production strategy. It uses original full3D WJ/DevStrainRawWorkV1, finite-nu
budget AWF, eta0, unchanged boundary/interface laws, IMEX1, and one
flow→k/bound→epsilon sweep. The 700-cell Q3/Q1 space and governing quadrature
12/16 are retained. The .005 time step is qualified only on the recorded short
unperturbed C/M/F comparisons. The prescribed finite seed fails initial principal
quality and is never advanced. Neither the passing 630-cell nor 8/12 initial
screen licenses evolved h/q substitution.

Use the branch-owning Mac worktree for Git. Linux alone supplied compilation,
solver, numerical and timing results. SSH `linux-usb` timed out; the executed
fallback was `ssh jay@100.105.247.110`. The unrelated general Linux checkout and
other jobs were preserved. The H25 capsules use the already authorized ext4
exception under `/home/jay/build`, since the old external mount alias is absent.

```sh
R=/home/jay/build/wj-h25-production-20261116
B=$R/resumption-r1
S=$B/source-v3/solvers/IncNavierStokesSolver/Tests/RANS
T=$B/tools-v4/solvers/IncNavierStokesSolver/Tests/RANS
L=$B/tools-v6/solvers/IncNavierStokesSolver/Tests/RANS
F=$B/tools-v7/solvers/IncNavierStokesSolver/Tests/RANS
```

Each operation's `request.json`/`receipt.json` is the exact executed argument
vector, input/script hash, deadline, worker admission, thread environment, RSS
and charge. These are historical commands: labels and output directories are
write-once and the final accounts are closed. A repeat needs a new bounded
experiment/journal and output paths; never reopen these accounts or overwrite
an accepted checkpoint. Reading/verifying the retained artifacts needs no rerun.

The selected policies are:

```text
RANSWallNativeLookup=TensorExactV1
RANSWallScalarEvaluation=Off
RANSWallAssemblyScatter=ElementExactV1
RANSWallEvaluationReuse=ExactStepV1
RANSWallImmutableRowReuse=Off
RANSWallRetention=BoundedV1
RANSWallMonitoringPolicy=CheckpointAuditV1
RANSWallPimpleSnapshots=Off
RANSWallPimpleFlowVectors=Off
RANSWallPimpleArithmetic=Off
RANSWallPimpleVerify=Off
```

The governing runtime is selected by `science/implementation-selected-v1.json`,
which pins `runtime-build-scatter-v1.json` and its wrapper/ELF. It comes from
source-v2, commit `39ff5e9110dc0202dd0af9072a180e979a46f2b4`, over source-v1 and
the immutable historical source-v8. The memo source-v3/build-memo-v1 was tested
and rejected for performance. The seed producer used separately built source-v1
at `955600e1569feea0317814f6b15fea4d5d3ce9aa`. Source and tool manifests bind every
capsule file. Documentation HEAD is never substituted for executable provenance.
All wrappers retain the pinned H22 private glibc loader/geometry libraries and
the installed MPI dependencies; operation native-exec maps verify actual binding.

Builds were Release/MPI-enabled, with eight jobs for the seed and four for the
scatter/memo builds. Science is one rank, with OMP, OpenBLAS, MKL, BLIS, NumExpr
and vecLib thread controls all one. Total counted host scientific/helper workers
plus reservations stay at most 72; the initial 131-worker obstruction admitted
no native work. Science RSS 8 GiB/build 16 GiB, new code 8 GiB/new science 8 GiB and
free space 15 GiB limits are explicit. The later native reservation includes ten
workers for driver, observer, rank/helper and bounded stdout pump.

The actual operational command family is:

```sh
python3 -B "$S/run_wj_h25r1.py" --root "$R" \
  --label C-wide-r1-v1 --seconds 900 --phase native --attempts 8 -- \
  python3 -B "$T/advance_wj_h25r1.py" --root "$R" --mode wide
python3 -B "$S/run_wj_h25r1.py" --root "$R" \
  --label C-temporal-r1-v1 --seconds 300 --phase passive -- \
  python3 -B "$T/temporal_wj_h25r1.py" --root "$R"
```

The wide control imports genuine H24 C4/.02, advances to C12/.06 and checks
exact parity with retained C5--8 before comparing C10/M20 at .05 and C12/M24
at .06. Both compared callbacks are completed; no interpolation is used.
Sparse/dense modes import newly qualified C8/.04 and reach C12/.06 independently,
using scientific strides4/1 and recovery stride1. Their terminal states/history
are compared with uninterrupted wide control, with a fresh startup and each
scheduled scientific assessment. `strict --label C-sparse-restart-v1` performs
the ordinary QualifiedV1 terminal read with genuine producing callback/history,
zero projection/history initialization and zero governing solves. `rotate` then
deletes only owned, superseded, off-cadence recovery payloads while keeping the
latest two, all scientific endpoints, completion metadata and the source anchor.
The actual receipts, not this command description, establish which modes passed.
The final post-rotation read uses `tools-v13/.../terminal_wj_h25r1.py` under
operation `C-terminal-read-r1-v1`; it reopens the retained terminal after deletion,
checks genuine alpha200 and scalar admission, and retains compact weak/mass-dual
residuals. Actual final costs use `tools-v7/cost_wj_h25r1.py`; final identity/account
closure uses `tools-v14/close_wj_h25r1.py`. These later tool identities are separate
from the unchanged selected native runtime.

The pilot log writer retains three 256 KiB native stdout chunks and a 16 KiB first
prefix, records total bytes/rotations and a streaming full digest, and limits
observer trace to 1 MiB. Full scientific state/failure payloads are independent
of stdout retention. Two focused Linux tests cover byte order, digest, bounded
retention and preservation of an existing log. The earlier wide logs remain.

`fresh_wj_h25r1.py` executes genuine NumSteps0 native projection/lift/reference/
defect/FlatInletTrace and own-dt history, compares all initial payload bits with
original C CP0, then independently assesses current physical state. It measures
fresh construction rather than relabelling a no-projection seed as initialization.
The complete-cost recipe counts this initial archive stage, two advancing native
processes (one planned restart), every native.step, all scheduled audit/wait/IO,
terminal exports and a separate strict read. Fixed process costs exclude all
native.step intervals. The largest observed block mean is used separately for
ordinary, scientific and fixed costs. The fixed horizon/window are unchanged;
ceil(T/.005) reaches 505.025 and is not relabelled 505.02466521278404.

Compact R0/R1/R2/C evidence is under `wj-h25-linux/resumption-r1/`. Each transport
manifest lists original and stored hashes/sizes; decode gzip before verifying
the original digest. Heavy arrays, checkpoints, reports, logs, binaries and
libraries remain at their exact indexed Linux paths. Final validation rehashes
retained provenance and explicitly distinguishes the intentionally rotated owned
recoveries. Closed account artifacts retain original ceilings, every failed
charge, all fixed reserves and separate incremental/cumulative costs.
