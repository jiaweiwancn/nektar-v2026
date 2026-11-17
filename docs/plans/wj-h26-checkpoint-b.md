# WJ-H26-B — Working verified native core

The compiled candidate passes the required numerical gates on Linux, including the corrected nonlinear flow with explicitly zero continuity forcing. The final selected-binary timing matrix remains pending, so integration readiness is NOT_READY at B. Full coupled model and production remain NOT_QUALIFIED.

Tested source: `8aeb00a2f37437cd17e043a9ba6f86a6f208e1d6`; verification tools are included in the same tested revision. Native source/build/runtime are isolated under `/home/jay/build/wj-h26-native-flow-20261117` (`source-v8`, `build-v8`, `runtime-v8`).

| Executed gate | Result |
|---|---:|
| Independent operator / total-traction discrepancy | 5.41011e-15 |
| Tiny independent direct-saddle agreement | 3.06655e-12 |
| Independent momentum / continuity residuals | 4.41083e-14 / 1.34076e-14 |
| BDF2 velocity orders | 2.0373, 2.0163 |
| Pressure temporal orders | 2.0713, 2.0330 |
| Current nonlinear flow residual / iterations | 6.69454e-12 / 6 |
| True-history restart | All four state/history vectors identical at t=.4 |
| Changed-H restart | Rejected before advancement |

| Geometry | U L2 order | U H1 order | Pi L2 order | Natural weak-traction order | WJ collar weak-traction order |
|---|---:|---:|---:|---:|---:|
| Flat | 4.2924 | 3.0922 | 5.3812 | 3.5730 | 2.9448 |
| Curved | 4.3029 | 3.0924 | 5.3687 | 3.5716 | 2.9444 |

The WJ collar test retained the initial n=4->8 rate near2.74, which missed the preregistered2.8 target; n=16 resolves the asymptotic rate near2.94. No C1 constraint or realizability clipping was used. The selected R spectra stay positive. The complete blended-tensor weak-load identity passes, and the finite molecular term is retained separately from the pure high-Re WJ helper.

The material20% viscosity test holds the body load fixed, changes the solution, and passes the current-system independent direct comparison with the same preconditioner. Mutation tests detect missing transpose strain, wrong AWF pressure sign and stale viscosity. The independent SymPy continuous forcing checks both polynomial and trigonometric problems. Actual stock/native constant-nu and constant-mu limits show convergent velocity and pressure errors; all reported orders and errors are in the JSON.

The v3 slow prototype averaged0.703085529s/step. Exact affine operator/source caches and a gauge-projected sparse preconditioner reduced the v4 complete-step pilot to0.0704072744s/step. The full predecessor v6 matrix is retained as preliminary profile evidence; its largest required mean is0.1300996484s/step. The final v8 matrix must be rerun after the restart-identity repair. These are prescribed/frozen flow measurements, not transported-turbulence or production timings.

Reproduce the selected numerical suite and constant limits with:

```bash
R=/home/jay/build/wj-h26-native-flow-20261117
T=$R/source-v8/solvers/IncNavierStokesSolver/Tests/RANS/h26
python3 -B "$T/campaign.py" --root "$R" --version 8 --mode verify
python3 -B "$T/campaign.py" --root "$R" --version 8 --mode constant_limit
```

Existing passing operations are verified and reused; a failed operation requires a new explicit label. To reproduce from scratch use a new isolated namespace/capsule and the same pinned inputs. All solver stages use one MPI rank and OMP/OpenBLAS/MKL/BLIS/NUMEXPR/VECLIB=1; the runner checks total live host workers<=72, RSS/storage/free-space limits and the protected campaign reserve.

Evidence: [numerical summary](wj-h26-linux/checkpoint-b/summary.json), [source/runtime identities](wj-h26-linux/checkpoint-b/identities.json), [Linux raw/input/receipt index](wj-h26-linux/checkpoint-b/linux-index.json), [preliminary complete-step matrix](wj-h26-linux/checkpoint-b/preliminary-v6.json). Large operator/field/checkpoint artifacts remain at the exact indexed Linux paths.

The earlier curved nonlinear fixture manufactured a nonzero continuity datum from an interpolated target. That result is retained on Linux as a predecessor, but does not close the zero-continuity flow gate. The selected fixture first projects a native weakly incompressible velocity and explicitly sets the nonlinear continuity forcing to zero; its final true continuity residual is 2.10010e-16.
