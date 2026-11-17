# Reproduce H27's constitutive feasibility result

This is a diagnostic reproduction, not a Task3 run. The account and outputs
under `/home/jay/build/wj-h27-coupled-20261117` are closed/create-only. Allocate
a fresh isolated namespace and execution allowance, respecting the original
aggregate H27 ceiling if continuing this goal. Do not reuse labels or reopen
that account. Mac owns source/Git; execute these numerical commands on Linux.

The exact executed source capsule was `source-v3`, whose compiled inputs are
byte-bound to `529eb9337bde9eb1be3606160f0ee1617ca4b54a`. The subsequent
`86b8f0d2b3493a147b93ebeee475e7e313eca63e` adds only the proposal and its
separately labelled tensor screen. The final publication checkpoint adds the
evidence and reporting tool. Consult the preserved requests for exact argv.

Use `ssh linux-usb`, or the verified fallback `ssh jay@100.105.247.110`.
Activate `/home/jay/.local/opt/nektar-v2026-mpi/env.sh`. Reserve at most four
compiler workers and one analysis/science process initially, within the live
72-worker total. Set all numerical-library threads explicitly:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
```

In the commands below, set `SRC` to a fresh Linux checkout/capsule with the
pinned source, `NEW` to a fresh writable result namespace on healthy ext4,
and `E` to the preserved publication directory (or its verified Git copy).
Supply an appropriate admitted runner around every command. The new runner
must carry prior aggregate charges for a continuation, and protect the
3600-second evidence reserve. The original exact commands/limits/identities
are in `E/operations/*/request.json`.

```bash
T="$SRC/solvers/IncNavierStokesSolver/Tests/RANS/h27"
python3 -B "$T/reference_oracle.py" \
  --k-factors 1 --scan-points 41 --principal-directions 128 \
  --output "$NEW/reference-kin.json"
python3 -B "$T/reference_oracle.py" \
  --principal-audit-input "$NEW/reference-kin.json" \
  --output "$NEW/planar-bound.json"
python3 -B "$T/inlet_principal.py" \
  --reference-json "$NEW/reference-kin.json" \
  --output "$NEW/fixed-inlet.json"
python3 -B "$T/outer_reference.py" --diagnostic-seconds 140 \
  --output "$NEW/outer-reference.json"
cmake -S "$T" -B "$NEW/build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS='-ffp-contract=off -fno-fast-math' \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$NEW/build" --target WJH27PrincipalProbe --parallel 4
```

The host's installed libraries require the pinned private loader arrangement
recorded in `raw/runtime-v3.json`. To invoke the diagnostic directly with the
preserved independent witness inputs:

```bash
L=/home/jay/build/wj-h22-early-phase-replay-20261115/runtime
"$L/glibc/ld-linux-x86-64.so.2" \
  --library-path "$L/glibc:$L/geometry:/home/jay/.local/opt/nektar-v2026-mpi/lib64" \
  "$NEW/build/WJH27PrincipalProbe" "$E/raw/native-states-v2.tsv" \
  > "$NEW/native-principal.json"
```

Use an unused output name; shell redirection is shown only inside the fresh
namespace. Preserve the new request, elapsed cost and all failure outcomes.
The native program reports the exact directional action, 2048-point planar
angular search/refinement and four actual-adapter FD levels for each input.
Compare it with `E/raw/native-principal-raw-v2.json` and
`E/raw/native-comparison-v2.json`. The independent root budget can pass while
the selected reference is scientifically rejected; execution success must
not be relabelled as model qualification.

The proposed alternate tensor has a separate, explicitly proposal-only
command. Its positive sampled result neither adopts the model nor qualifies
the coupled solver:

```bash
python3 -B "$T/proposal_principal.py" \
  --reference-json "$NEW/reference-kin.json" \
  --output "$NEW/proposal-only.json"
```

There is **no Task3 launch command**: a registered, admissible, verified,
accuracy-qualified coupled U/Pi/k/epsilon recipe has not been produced.
`task3-launch-hold.json` records null executable/state/recipe/forecast/commands
and abort reasons. It preserves T=505.02466521278404 and terminal window
[447.8483832973273,505.02466521278404]. Do not substitute the diagnostic
executable, an H26 prescribed-coefficient driver or old H25 histories.
