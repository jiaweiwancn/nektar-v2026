# H22 compact evidence and reproduction

H22 is CLOSED/OBSTRUCTED_OBSERVER_AUDIT_ROLE. The requested seven-attempt
diagnostic is incomplete: C1 and M1/M2 were captured; F1–F4 were not launched.
M2 remains completed-but-unqualified. The passive correction is not an executed
native monitor qualification. Read the evidence report before using a PASS
inside a partial-result table.

Authoritative Linux root:
`/home/jay/build/wj-h22-early-phase-replay-20261115`.
`authoritative-index-v1.json` binds `publication-v1` and separately CLOSED
`closure-v1`. Publication accounting deliberately remains its earlier OPEN
snapshot. `local-storage-aliases.json` maps only hash-verified copied files;
other historical/native/science paths remain remote identities, not local files.
Large arrays, checkpoints, vector binaries and sparse matrices remain on Linux.
Compressed JSON copies are lossless and retain original/copy identities in
`publication-v1/copy-bindings.json`.

Useful starting records:

- `publication-v1/component-gates.json`, `timeline.json`, `online-counts.json`;
- `publication-v1/attempt-bindings.json.gz`, `roles.json.gz`, phase/window tables;
- `publication-v1/passive-correction/result.json`, exact retained-prefix proofs;
- `publication-v1/M/first-replay-capture-failure.json`, failure ACK and unchanged cursor;
- `publication-v1/C/step-1-physical.json.gz` and complete signed moments;
- `publication-v1/controls/runtime-build-v1.json.gz`, compiled-source/object pins,
  runtime dependency and actual-exec-map records;
- `closure-v1/execution-counts.json`, `verification.json`, `accounting-final.json`;
- `closure-v1/operation-receipts/` for post-publication operations and all receipt copies.

Native source capsule is057e6b7465660cac5ec8617abf5c876fa8bff71f.
Governor ELF SHA-256 is
`1c1df27b32ca5b6925ee467431f5291e25996e010c6c0026bebd12130e4d2e1d`;
reader ELF SHA-256 is
`8ca139e140637bfe14ae875b1b5815b50fbb131229ac7020b3d9168791daa739`.
The runtime checker reuses `study_wj_h17.runtime_qualify`; its inherited scope
text still says H17. The actual operation paths, selected new ELFs, PID samples
and library pins in that record are H22. This label is not a claim that H22 used
an old executable. No historical manifest was edited.

`publication-v1/source-versions/executed` preserves the actual failing Python;
`passive-code-v3` is the separately pinned retained-data correction;
`closure-v1/final-source` includes the later common sealed-index selector/test.
The current Git observer must not be mistaken for the observer executed in the
failed M block. The latter is bound by execution-config-v2 and its per-file pins.

Commands below document what ran. The root and ledger are closed: rerunning a
governor or resetting reservations is prohibited. Exact argv/thread/allocation,
time, source and log pins are in each operation receipt. Set all six
OMP/OpenBLAS/MKL/BLIS/NumExpr/VECLIB thread controls to1 before numerical imports.

```bash
# The frozen runner and explicit interpreter admitted the actual C/M commands.
R=/home/jay/build/wj-h22-early-phase-replay-20261115
P="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
/usr/bin/python3 -B "$P/run_wj_h22.py" --root "$R" \
  --label target-C-actual-v1 --seconds 600 --phase target \
  --kind target --workers 2 --attempts 1 -- \
  /usr/bin/python3 -B "$P/study_wj_h22.py" --root "$R" --mode target --level C
# M was the corresponding target-M-v1 reservation, --attempts 2 --level M.
# F had no reservation or launch. Final strict commands selected C1 and M0.

# Read-only correction/partial analysis; preserved under its own charged receipt.
/usr/bin/python3 -B "$R/passive-code-v3/repair_wj_h22_roles.py" --root "$R"
# Its existing output is immutable; this is the executed command, not permission
# to overwrite/repeat the closed operation.
```

The build command/flags and object-reuse proof are in the affected-build receipt
and controls. The native serializer fixture was actually executed through the
private loader, not merely registered in CTest. Initial Python tests are listed
in `publication-v1/operations/focused-tests-v1/command.log`; the additional
session/cursor/audit-role/index fixtures have their own receipts. There were24
distinct Python methods and one synthetic native component, not a new small CFD
trajectory. No kernel or saved governing-matrix operation ran.

The two figures were visually checked on Mac after exact archive/file hash
verification. They preserve signed contributions and distinguish the exact
upper/lower normal owners at10.2430071555mm. Their C/M-only coverage is explicit;
they supply no F or later-window stages. No duration forecast is derived from
this diagnostic capture mode.
