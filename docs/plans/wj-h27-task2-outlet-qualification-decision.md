# Flat WJ–AWF outlet qualification decision

Canonical spatial qualification remains open for `ZeroTransverseV1`. Bn2 and Bf pass their independent current and physical audits. Their direct comparison exceeds the stress, shear and gradient limits in 12 of 690 spatial rows. Every row in x∈[−2.5,8] meets its spatial limit. The outlet Rxr and shear maxima reach .003147 and .003313 tau0 at x=15, r=.06. The independent terminal x witness meets the spatial difference limits and exceeds the independent momentum limit. The pressure-action audit now uses exact cell integrals at its original tolerance. [Signed spatial evidence and independent gates](wj-h27-task2-outlet-qualification-cp2.md)

`FullTractionV1` remains the default. It historically missed 16 rows in the paired corner mesh screen and 180 rows in the isolated radial collar comparison. `ZeroTransverseV1` remains the explicit flat candidate. It prescribes zero transverse velocity and the retained normal traction at x=15, with free axial velocity away from the shared ground vertex and a measured transverse reaction. The tests use Q3/Q2, 6×8 native quadrature, 1e−12 linear tolerance and 1e−10 coupled tolerance. They retain the exact finite-viscosity reference, ScalarCdV2, sampled AWF, current-K cutoff, full production and the implicit positive stress split with its complete remainder. All runs stay at zero simulated time. [Pinned model, meshes and limits](wj-h27-task2-outlet-qualification-contract.md) · [Historical A and B comparison](wj-h27-task2-outlet-decision.md)

## Outlet location

The move to x=17.5 is deferred pending an independently qualified canonical mesh and a complete matched-recipe forecast. The historical 118×69 test on B2 carried a 2,464 s and 9.92 GiB forecast. A selected radial grid and numerical recipe determine the new cost. The comparison preserves the upstream x edges through 14.375 and measures at x=15 as an interior station. Native parity on the retained x=15 domain precedes a matched extension. [Frozen extent criteria and disposition](wj-h27-linux/task2-outlet-qual/cp3/outlet-qual-outlet-location-disposition.json)

The retained window is x∈[−2.5,8]. The outlet-location RMS/MAX limits for the extension-minus-canonical differences are .01/.02 for U and W, .05/.10 for k/Kin, .10/.20 for local outer epsilon and .05/.10 for ground shear. Stress, Pi and traction carry signed effects and within-B numerical uncertainty. Their outlet-location physical acceptance thresholds were unspecified. Physical departures from the inlet reference are case-specific. Flat outlet-location sensitivity remains unmeasured. Hill wake sensitivity remains a separate decision. [Physical and numerical observation contract](wj-h27-linux/task2-outlet-qual/cp3/outlet-qual-outlet-location-disposition.json)

## Resource and launch disposition

The smallest next targeted experiment splits the two Bf children of [.06,.06875], giving 111×94 cells. The full forecast is 2,783 charged seconds and 12.21 GiB process-tree RSS. The exact lower-band axes, limiting signed shear budget and time and memory anchors are sealed. The upper blend and flanks, streamwise station coverage, one-sided H convergence and the separate quadrature and tighter-algebraic controls remain open. The latest state and mesh hashes are recorded as audited observations. [Next test and cost](wj-h27-linux/task2-outlet-qual/cp2/outlet-qual-next-experiment.json) · [Launch hold and state identities](wj-h27-linux/task2-outlet-qual/cp3/outlet-qual-task3-launch-hold.json)

**TASK2_COMPLETE=NO. TASK3_READY=NO.** Physical timestep accuracy, cold restart, matched stock timing and complete advancing cost are future gates after the flat recipe qualifies. The one-day production criterion remains `1.25*C_total ≤ 86,400 s` per trajectory.

## Reproduction

The Release OpenMPI runtime ran on `jay@100.88.70.60` on the ext4 experiment volume. The solver source identity is `54b2f4268a5903fcf894f032b968bbfcb42a3e999a9c0409d032e999114894a8`, the application SHA-256 is `3e8db421575122676808adc6c9783f0c897136e9d9d283b9d0cb200c5ec2df45`, and the exact-inner table SHA-256 is `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`. Each named solve uses one MPI rank and six single-thread library settings. The execution lock and watchdog carry the 80-worker host limit and case-specific RSS checks. The Bf retry receipts record the actual permitted CPU masks. [Per-case commands, resource receipts and Linux archive](wj-h27-linux/task2-outlet-qual/cp2/operations.json)

```bash
W=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
R="$W/runs/wj-h27-task2r-20261118"
source "$W/env.sh"
S="$R/tools/outlet-control-v12/task2r_run.py"
D="$R/tools/outlet-python-v27/h27"
export PYTHONPATH="$R/tools/outlet-control-v12:$D:$R/tools"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python3 -B "$S" "$R" balance
python3 -B - "$R" "$S" <<'PY'
import json
from pathlib import Path
import shlex
import sys

root, runner = sys.argv[1:]
for label in (
    "outlet-qual-bf-recovery-forecast", "outlet-qual-bf-admission-1",
    "outlet-qual-bf-continue-1-inputs", "outlet-qual-bf-continue-1-steady",
    "outlet-qual-bf-continue-2-inputs", "outlet-qual-bf-continue-2-steady",
    "outlet-qual-bf-audit", "outlet-qual-bf-compare",
):
    path = Path(root) / "operations" / label / "receipt.json"
    request = json.loads(path.read_text())["request"]
    print(shlex.join(["python3", "-B", runner, root, "run", label, request["phase"],
                      str(request["cap_seconds"]), "--", *request["command"]]))
PY
```

The printed commands identify the sealed run. A fresh run requires new admitted output paths and labels.

## Settlement and checkpoints

The final Linux settlement records 59,341.329 cumulative charged seconds against 64,800 s. It leaves 1,858.671 nonprotected seconds, 3,600 protected seconds and zero reserved seconds. The forecast shortfall for the 111×94 repeat is 924.329 s. The case-specific 12 GiB cap also falls below its 12.21 GiB estimated peak. The SHA-256 of the sealed CP3 archive is `d5849a9dd346ce4d920e7df120c8d979700f707e1e4752d67220b8b08462464f`. [Settled ledger](wj-h27-linux/task2-outlet-qual/cp3/accounting.json) · [Final publication receipt](wj-h27-linux/task2-outlet-qual/cp3/operations.json)

| Verified checkpoint | Git SHA |
| --- | --- |
| CP0 protocol and resource amendment | `bf34514ef0a348304403f9bb122c6bc3576a476d` |
| CP1 replay and B3 challenge | `5999a4a8e4322212f91946cf0778a83a72c5f23e` |
| CP2 spatial qualification disposition | `59c0fff2b9249a71a037cc60565dec29553e6b41` |
