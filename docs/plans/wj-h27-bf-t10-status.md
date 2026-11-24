# Bf advancement through T=10

Both Newton variants and the two implicit fixed-budget Picard runs completed step 500 at \(\Delta t=0.02\) on the retained-inlet Bf mesh. All used outlet `ZeroTransverseV1`, one MPI rank, and one compute thread on `jay@100.88.70.60`. The [qualified steady Bf solution](wj-h27-task2-outlet-campaign-cp2.md) supplies the common field reference.

## Newton endpoints

The pressure-corrected and Sparse-LU runs both reached \(t=0.32\) with BDF2/EX2. The next step hit the outlet stress-covariance guard during Newton backtracking. The diagnostic recovery used current implicit advection, a forcing multiplier of four, and a Krylov target capped at 0.3. Each run resumed from its own verified history and kept its established full-step tolerance.

| Variant | Full-step tolerance | Native maximum at T=10 | Independent momentum | Independent continuity | Independent k | Independent outer epsilon |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Sparse-LU | 1e-10 | 7.902e-12 | 2.619e-10 | 1.872e-14 | 1.451e-10 | 1.119e-10 |
| One pressure corrector | 3e-6 | 4.013e-7 | 4.013e-7 | 5.078e-10 | 7.073e-8 | 4.958e-8 |

Both native stops passed. The independent Gauss-20/24 audit passed all four equation limits for Sparse-LU. For pressure correction, continuity passed and momentum, k, and outer epsilon exceeded their independent limits. The independent limits are 1e-9 for momentum and continuity and 1e-8 for each scalar. Both audits passed the quadrature, history cache, boundary history, and accounting checks. [Sparse-LU endpoint evidence](wj-h27-linux/bf-t10/sparse-t10.json) · [Pressure-corrected endpoint evidence](wj-h27-linux/bf-t10/pressure-t10.json)

The field RMS values use area weighting on the shared signed Bf maps, with outer epsilon integrated over its transported domain. The T=0.08 column is the previous four-step retained-inlet Sparse-LU endpoint.

| Field | T=0.08 RMS from steady | Sparse-LU T=10 RMS from steady | Pressure-corrected T=10 RMS from steady | Direct T=10 RMS difference |
| --- | ---: | ---: | ---: | ---: |
| Axial velocity | 2.609e-4 | 1.601e-4 | 1.600e-4 | 6.780e-7 |
| Transverse velocity | 3.582e-5 | 2.547e-5 | 2.549e-5 | 4.066e-7 |
| Modified pressure | 2.389e-5 | 1.734e-5 | 1.736e-5 | 1.296e-7 |
| k | 1.554e-5 | 1.515e-5 | 1.515e-5 | 4.806e-8 |
| Outer epsilon | 1.149e-6 | 7.966e-7 | 7.964e-7 | 2.124e-9 |

The Newton endpoint distances from steady agree within 0.11% across fields. Axial velocity still changed by 2.31e-7 RMS in the last Sparse-LU step and 2.37e-7 in the pressure-corrected step.

The tested Newton diagnostic used base revision `2f28cce00cb1284e9a5736b8a92d57f599a301ae`, source identity `4e893ab954e7132967a9ddb260b89ab9a869cd2e46fe674c4072abeb3db5b31e`, and application SHA-256 `c6a88ef1b19a2a872007df3a7870bfbe5a9afce525d1360dbea2e8e301446c15`. The [Newton diagnostic patch](wj-h27-linux/bf-t10/newton-diagnostic.patch) applies to that base with `git apply --unidiff-zero`. The [audit reader patch](wj-h27-linux/bf-t10/newton-audit-scheme.patch) applies to this snapshot's transient oracle and selects implicit advection from the diagnostic report marker.

## Fixed-budget Picard endpoints

The checked-in Picard changes add optional implicit advection and the velocity-tangent stress linearization. The Linux build using the matching source files completed both the three-sweep and four-sweep cases with one pressure corrector per sweep. Their native step-500 residuals were as follows.

| Fixed budget | Momentum | Continuity | k | Outer epsilon |
| --- | ---: | ---: | ---: | ---: |
| Three sweeps | 2.873e-11 | 6.253e-11 | 3.326e-10 | 2.185e-11 |
| Four sweeps | 4.930e-11 | 6.312e-11 | 1.612e-11 | 2.185e-11 |

Both final Picard checkpoint histories verified. The Linux analytic controls passed for explicit and implicit advection at both startup and BDF2. The Picard rows describe native endpoint checks. [Picard runtime, checkpoints, and analytic controls](wj-h27-linux/bf-t10/picard-t10.json)

## Timing comparison

Steps 201–500 cover \(4<t\le 10\) at \(\Delta t=0.02\). This interval is inside the final uninterrupted segment of all four runs. Complete step wall time includes the advance, journal, checkpoint, and scientific output for each step. The receipts record one MPI rank and one compute thread per case, pinned to the CPUs below, and no supervisor pauses in these segments.

| Solver | Pinned CPU | Total wall (h) | Mean (s/step) | Median (s/step) | P10–P90 (s/step) | Time relative to Sparse-LU |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Coupled Newton, Sparse-LU | 105 | 11.37 | 136.5 | 136.4 | 135.9–137.1 | 1.00× |
| Picard, 4 sweeps × 1 corrector | 93 | 12.48 | 149.8 | 149.1 | 146.6–151.5 | 1.10× |
| Coupled Newton, 1 pressure corrector | 14 | 12.74 | 152.9 | 141.9 | 92.2–244.8 | 1.12× |
| Picard, 3 sweeps × 1 corrector | 3 | 14.27 | 171.3 | 171.1 | 159.2–182.0 | 1.26× |

The two Newton variants used full-step residual tolerances of `1e-10` and `3e-6`, respectively. Both Picard variants used implicit advection, the velocity-tangent stress linearization, a fixed sweep count, and a final continuity tolerance of `1e-10`. The timings compare the archived configurations with these acceptance criteria and their own restart histories.

### Where the time went

- Sparse-LU spent 82.5 s per step on rebuilds, about 61% of native step time. Pressure correction lowered that cost to 42.2 s. It averaged 36.5 Jacobian actions per step versus 13.0 for Sparse-LU, with 27.3 s more assembly and 11.3 s more constitutive work per step. Its step cost varied widely.
- Four-sweep Picard took 12.6% less wall time than three-sweep Picard over the matched interval. The pressure solves averaged 86.7 versus 114.3 s per step, with about 185 aggregate pressure Krylov iterations per step in each run. Rank CPU and complete step wall means agreed within 0.3 s for all four cases. The runs used different pinned cores on a shared host, which limits attribution of the timing difference to sweep count.

### Time cost across the full march

The accepted native step totals sum the `seconds` field for steps 1–500. The journals form continuous state and advection histories across the restarts. Launcher totals sum receipts on the retained histories and include startup, output, pauses, and incomplete work in those launches. Time between launches and separate test and recovery probes are outside those totals.

| Solver | Accepted native step time (h) | Launcher wall time (h) |
| --- | ---: | ---: |
| Coupled Newton, Sparse-LU | 19.96 | 21.24 |
| Picard, 4 sweeps × 1 corrector | 21.47 | 21.69 |
| Coupled Newton, 1 pressure corrector | 25.75 | 26.43 |
| Picard, 3 sweeps × 1 corrector | 23.54 | 23.76 |

The early Newton histories include explicit advection through step 16 and recovery to implicit advection at step 17. Each physical time band below contains 100 accepted steps. The values are median native seconds per step, including any restart effects on the accepted steps in that band.

| Solver | \(0<t\le 2\) | \(2<t\le 4\) | \(4<t\le 6\) | \(6<t\le 8\) | \(8<t\le 10\) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Coupled Newton, Sparse-LU | 137.2 | 136.2 | 136.3 | 136.3 | 136.0 |
| Picard, 4 sweeps × 1 corrector | 166.0 | 154.5 | 150.3 | 148.4 | 146.9 |
| Coupled Newton, 1 pressure corrector | 226.6 | 151.5 | 142.9 | 132.7 | 147.1 |
| Picard, 3 sweeps × 1 corrector | 159.8 | 164.8 | 160.8 | 172.8 | 174.7 |

Four-sweep Picard continued to get cheaper, with a 2.3% drop in median cost from \((4,6]\) to \((8,10]\). Three-sweep Picard rose 8.7% across those bands. Sparse-LU stayed near 136 s per step, and pressure-corrected Newton slowed again in the final band.

The [timing evidence](wj-h27-linux/bf-t10/timing.json) records the segment paths, tested source identities, and timing breakdowns. To regenerate it from the retained Linux archives, run the following command from the repository root. [Collector](wj-h27-linux/bf-t10/summarize-timing.py)

```bash
ssh jay@100.88.70.60 \
  'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 - /media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs' \
  < docs/plans/wj-h27-linux/bf-t10/summarize-timing.py
```
