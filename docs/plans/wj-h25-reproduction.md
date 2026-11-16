# H25 reproduction and retained evidence

Source/Git is the branch-owning Mac worktree. Numerical/build commands run on
Linux, with all six numerical thread controls1 and the qualified private loader.
The finite ledger is `/home/jay/build/wj-h25-production-20261116/accounting.json`.
Closed/failed operation labels are immutable: these commands document executed
work and must not be relaunched into an existing reservation.

```sh
R=/home/jay/build/wj-h25-production-20261116
P=$R/source-v4/solvers/IncNavierStokesSolver/Tests/RANS
python3 -B "$P/run_wj_h25.py" --root "$R" \
  --label profiling-build-v2 --seconds 900 --phase build -- \
  python3 -B "$P/build_wj_h25.py" --root "$R"
python3 -B "$P/run_wj_h25.py" --root "$R" \
  --label W0-profile-four-v1 --seconds 600 --phase native --attempts 4 -- \
  python3 -B "$P/profile_wj_h25.py" --root "$R"
```

The original eight-step W0 uses immutable `source-v2/study_wj_h25.py --mode
baseline`. Successful stock S0 uses `tools-v4/stock_opt_wj_h25.py` and the explicit
NoCollection `.opt`, not the failed initial tuning commands. Each operation's
receipt pins its exact executable command, source/script, threads, measured
duration and RSS. Newer governor maps are sampled after native ELF execution;
the original W0 receipt has no fresh map capture and uses the separately verified
historical closure. Do not confuse pre-exec Python host maps with native mixing.

`source-v5`/`build-lookup-v1`/`runtime-lookup-v1` contain the first exact-lookup
candidate. Baseline command `profile_wj_h25.py --variant lookup-off` retains Off;
`--variant lookup-on` is the separately admitted candidate comparison. The
independent numerical evaluator/cache remains unchanged with explicit new-path
compatibility records. The cache's historical pins are never edited.

Compact copied Linux evidence is in `wj-h25-linux/checkpoint-a-v1/`; its transport
manifest distinguishes original and compressed hashes. Coefficient arrays,
checkpoints and large reports remain at their exact Linux paths. The open
accounting snapshot is not a CLOSED final receipt. Completion will add separately
sealed candidate decisions, pilot/restart/storage and final accounting evidence.
