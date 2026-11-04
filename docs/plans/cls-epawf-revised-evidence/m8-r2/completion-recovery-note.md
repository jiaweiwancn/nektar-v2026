# Additive recovery of the frozen-scalar completion

The original `epoch-002/frozen-scalar-001/receipt.json` remains failed and
unchanged (SHA256
`3656d28e44cf878a0177653e7e7bb21dbd3ca0aeec16b738acdb627c477ccfc6`).
Its supervisor encountered `BrokenPipeError` while publishing progress. The
native child continued and later published a complete report and JSONL record.
The supervisor's `finally` fingerprint ran before native completion, and a
second broken progress write prevented its final provenance fields. That early
fingerprint cannot serve as a post-completion fingerprint.

The lost native process exit status and elapsed time remain **null**. Solver
timing text stays in the native log and does not replace process timing. The
supervisor wall cap was not verified after the stdout failure. An additive
recovery does not change the failed supervisor status or claim its normal
completion.

`completion-check-001/receipt.json` records that the complete report existed
before a fresh fingerprint, no producing solver process was observed, source
revision `9f20b3053c303e8805f5183f668a8ec468d9656b` was clean, and the original
receipt, XML, log, report and trace stayed unchanged. Its fresh
`runtime-postcompletion.json` has SHA256
`9aefcf2e4b1295b3f65c29b067d8f9033f463f80f033e9d97a5cb4eef0bd7a43`.

`recover-stage-completion.py` consumes these retained files without executing a
native solver, fingerprint tool, scientific oracle or Git command. It rehashes
the artifacts, complete JSONL record and every matrix, verifies matrix byte
lengths, and matches the final report's state/history/norms to the complete
record. It requires actual unchanged native history and the exact protected
preflight XML. It also rechecks the original matched prerequisite and its
historical authorization assessment, including the absence of positive coupled
target convergence evidence. The matched receipt and assessment SHA256 values
are respectively
`35d554bf555fa582dd2cf57717d933e81c513baf8e3d5d58bd15c1d9f09eafeb` and
`5d07cf84294b3b0204095767de99632a54b96f6497317d6938fab9d45ef32d38`.

The collector produces a **separate** recovery receipt. The assessor's optional
`--recovery` path recomputes those facts and requires an exact canonical match;
its normal receipt path remains strict. The recovered family status is
`recovered_native_artifacts`, with the original supervisor status retained as
`failed`. Missing process status/duration and unverified wall-cap enforcement
stay explicit. Complete native nonconvergence is retained as an observation;
metadata recovery supplies no derivative, path, physical-trajectory or root
qualification.

`curate-stages.py --recovery` independently recomputes the same recovery and
assessor gates before producing the compact native summary. Its `runtime_after`
and source-end observation come from the verified fresh post-completion pin.
It preserves `original_supervisor_status: failed`, null native exit/duration,
and unverified wall-cap enforcement. This is an explicitly recovered-artifact
status, without relabelling the original supervisor receipt.

The complete-record iteration count comes from the JSONL `iterations` length.
The final report's `iterations` field is not substituted for that count. This
metadata recovery does not classify trial paths or infer a failure mechanism.
The separate terminal accounting audit retains the finite measurements near
the scalar boundary; measured positive values below their recorded arithmetic
roundoff scale do not establish a certified positive margin.

The root supervisor executes the following commands on the authoritative Linux
host, with one explicitly allocated CPU and the six numerical thread variables
set to `1`. `TOOLS` is a clean immutable analysis checkout containing the same
collector and assessor versions; `R2` is the external run root. Each output path
must be fresh.

```sh
taskset -c 50 python3 -B "$TOOLS/recover-stage-completion.py" \
  --receipt "$R2/epoch-002/frozen-scalar-001/receipt.json" \
  --completion-check "$R2/epoch-002/frozen-scalar-001/completion-check-001/receipt.json" \
  --expected-source 9f20b3053c303e8805f5183f668a8ec468d9656b \
  --output "$R2/epoch-002/frozen-scalar-001/completion-recovery-001.json"
taskset -c 50 python3 -B "$TOOLS/assess-stage-evidence.py" \
  --receipt "$R2/epoch-002/frozen-scalar-001/receipt.json" \
  --recovery "$R2/epoch-002/frozen-scalar-001/completion-recovery-001.json" \
  --expected-source 9f20b3053c303e8805f5183f668a8ec468d9656b \
  --output "$R2/epoch-002/frozen-scalar-001/recovered-assessment-001.json"
taskset -c 50 python3 -B "$TOOLS/curate-stages.py" \
  --run "$R2/epoch-002/frozen-scalar-001" \
  --recovery "$R2/epoch-002/frozen-scalar-001/completion-recovery-001.json" \
  --output "$R2/epoch-002/frozen-scalar-001/recovered-native-summary-001.json"
```

The focused stdlib test recipes are:

```sh
python3 -B test-recover-stage-completion.py
python3 -B test-assess-stage-evidence.py
python3 -B test-stage-progress.py
```

They exercise complete recovery, unchanged normal-path rejection of the failed
raw receipt, missing/stale fresh fingerprints, altered raw receipts or native
artifacts, incomplete or mismatched complete-stage evidence, and invented
integer native exit status in the completion check, raw receipt or additive
receipt. All fixtures are synthetic metadata; no native or scientific code is
called.

The future `run-stages.py` progress emitter catches only broken-pipe errors,
disables later progress writes, and redirects broken actual stdout to
`/dev/null` (or detaches a closed stream) so interpreter shutdown cannot exit
with a final-flush EPIPE. Child polling, process groups, timeout enforcement,
native commands and receipt management are unchanged. The focused tests also
close a real metadata-only Python child's stdout before its first write and
verify continued simulated polling/final receipt work and exit code `0`.
