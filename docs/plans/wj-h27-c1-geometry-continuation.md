# C1 physical hill geometry continuation

The continuation starts from a developed, stationary discrete solution at height fraction .002 (0.08 mm), using the retained finite-viscosity inlet and the exact-inner reference table. The full experimental hill has height 40 mm, half-length 100 mm, and height fraction 1. The effective wall, independent outer mesh, fixed roughness displacement, fixed $h$, sampled forces, boundary laws, and cutoff family are defined in the [C1 governing and verification report](wj-h27-c1-normal-hill.md).

## Code and acceptance

| Operation | Implementation |
| --- | --- |
| Typed physical seed import | `WJNativeCoupled.cpp`, `WJNativeInitialFields.cpp` and `Coupling::InitialiseFromFields` in `WJNativeCoupling.cpp`. `WJGeometrySeedFromHeight` binds the declared source geometry, source revision, reference table, outlet policy, signed basis and coefficient hash. A positive-time packet needs a residual-controlled implicit history. A same-height, same-grid continuation uses its ordinary cold checkpoint. A strict nested grid can transfer the Q3/Q2 polynomials and start fresh. |
| Velocity transfer | `Coupling::RemapVelocity` carries contravariant volume flux between wall-layer charts. A tapered outer correction removes the jump in the Piola tangential limit at $H$ while preserving the mapped normal trace. `ProjectContinuity` restores weak divergence on the target chart after the velocity and scalar lifts are rebuilt. `WJGeometrySeedTaper` controls the radial extent of the seed correction. |
| Time advancement | `FreshHistory` constructs BE and BDF2 histories for each changed height or nested grid. `WJNativeTimeNewton.cpp` solves the complete new-time residual. The physical pilot clock counts only accepted time on its own geometry. |
| Independent checks | `c1_oracle.py` reconstructs physical derivatives, cutoff panels, outer quadrature and budgets from exported coefficients. `c1_margins.py` uses independent WJ companion roots, two centered directional differences, and planar angular stationary values at the finer volume rule and both sides of matching. `c1_continue.py` requires both audits before a state seeds its successor. |

The advancing target is $10^{-10}$ for the native refreshed coupled norm. The independent scaled limits for momentum, continuity, k and outer epsilon are $10^{-9},10^{-9},10^{-8},10^{-8}$. The dependent epsilon projection limit is $10^{-10}$. The pointwise polynomial representation error is recorded separately. The independent principal scan requires a margin greater than $10^{-8}$ with the derivative-difference allowance, and covariance over k greater than $10^{-10}$. This is a finite spatial scan. The whole-interface K extrema are taken from the signed cubic trace, including interior extrema.

## Build and run identities

The authoritative workspace is `/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace` on the Linux host `jay` at Tailscale address `100.88.70.60`. When the peer key expired during the last run, the pinned SSH host key and ext4 mount both verified on `jay.local`. The workspace and running jobs stayed on that host. The inspected live WJ branch tip was `57845c53a29ba32cca49c31a89e0f4ec4da4f300`. The starting local governing implementation and report were at `c0d5962f38e67e827c0467c501db3ea6d1dcab05`. The final standard-main ELF SHA-256 is `2ad3f3e5ac9a9b988393078288d0358a0c42e2ae349372c1974bfdf61f86f058`. Its embedded numerical source identity is `4dea57192d813744bffd938ee67907e37dc5ec9b7cf470d4587842e2000123c5`. The audit that additionally integrates the separate Cartesian forces is at `7810225fd6c1d30f1c1d3ac64af3e400b1d73658`. Numerical revisions and alternate executable hashes are listed in the [Linux evidence manifest](wj-h27-c1-geometry-continuation-evidence.json).

The Release builds use GCC 11.4 and the workspace MPI installation. The standard-main solver runs with one MPI rank, no MPI core binding, and one thread in each numerical library. Eight build workers were used at most for this work. The shared 72-worker limit includes concurrent workloads. Large fields, checkpoints, rejected states, XML inputs, logs and scripts are under `runs/wj-c1-geometry.uZz3Uh7t` in the workspace. The manifest binds the primary session, seed, audit and restart files by size and SHA-256. The retained physical calibration uses $L_0=.04\,\mathrm m$, $H=.15$, $d_0=.0875$, $h=.11541061001735275$, $K_{\rm in}=.008809015145515192$ and $c(K_{\rm in})=.0243960321782528$. The reference table SHA-256 is `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`.

```bash
ssh jay@100.88.70.60
# On the local network, use the already pinned key if the Tailscale key has expired.
# ssh -4 -o HostKeyAlias=100.88.70.60 -o StrictHostKeyChecking=yes jay@jay.local
ws=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
source "$ws/env.sh"
archive="$NEKTAR_RUNS/wj-c1-geometry.uZz3Uh7t"
RUN=$(mktemp -d "$NEKTAR_RUNS/wj-c1-geometry-replay.XXXXXXXX")
mkdir "$RUN/tmp"
tests="$NEKTAR_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/h27"
table="$NEKTAR_SOURCE/docs/plans/wj-h27-linux/task2r-r1/tables-exact-inner-v1.json"
export TMPDIR="$RUN/tmp" PYTHONPATH="$tests" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cmake -S "$tests" -B "$RUN/build" -DCMAKE_BUILD_TYPE=Release \
  "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
cmake --build "$RUN/build" --parallel 8 --target IncNavierStokesSolver WJH27HillGeometryProbe
solver="$RUN/build/IncNavierStokesSolver"
run_case() {
  (cd "$1" && /usr/bin/time -v mpirun --stdin none --bind-to none -np 1 "$solver" \
    --use-opt-file run.xml.opt run.xml >run.log 2>&1)
}
```

For the stationary anchor, reuse the certified geometry-bound packet under `nested-developed-002` with the final numerical source `7c3d6737378a5eaaeed16580a502ccbb8f54ea26`. To rebuild the stationary state, use the bound `taper-developed-002/state.result.json` packet with its explicitly declared source fingerprint as a steady-preparation guess. The prior coarse stages and their inputs are in `flat-coarse-with-feet`, `flat-coarse-continue-fd4`, `coarse-shallow-00025`, `coarse-shallow-00025-q8x12`, and `continuation-track`.

```bash
prior="$archive/taper-developed-002/state.result.json"
source_id=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["state"]["source_identity"])' "$prior")
python3 "$archive/make-coarse-mapped.py" hill "$RUN/developed" --height .002 \
  --quadrature 10 16 --initial-fields "$prior" --initial-source "$source_id" \
  --steady-seconds 120 --step-seconds 180
run_case "$RUN/developed"
python3 "$archive/make-local-refined.py" hill "$RUN/streamwise" --height .002 \
  --from-height .002 --initial-fields "$RUN/developed/state.result.json" \
  --quadrature 12 20 --seed-taper .30 --dt .02 --steps 2 --science 1 --step-seconds 180
run_case "$RUN/streamwise"
for test in c1_oracle c1_margins; do
  python3 "$tests/$test.py" --table "$table" --native "$RUN/streamwise/state.result.json" \
    --order 24 --seconds 300 --out "$RUN/streamwise/$test.json"
done
python3 "$tests/c1_continue.py" --solver "$solver" --fixture "$archive/make-local-refined.py" \
  --table "$table" --start "$RUN/streamwise/state.result.json" \
  --start-audit "$RUN/streamwise/c1_oracle.json" --start-margins "$RUN/streamwise/c1_margins.json" \
  --root "$RUN/streamwise-track" --quadrature 12 20 --increment .00315 --minimum .00025 \
  --seed-taper .30 --attempts 24 --target 1
```

The coarse x edges are `[-10, -2.5, 0, 2.5, 15]`. The nested streamwise edges are `[-10, -5, -2.5, -1.25, 0, 1.25, 2.5, 6, 15]`. Both retain the 92 radial cells and the explicit hill feet. The radial witness bisects radial cells intersecting $[.80,1.36]$ in chart radius using `make-local-foot-radial.py`. For that witness, use the streamwise result as the seed with `--from-height .002 --height .002`, run two steps and audit again before substituting the radial fixture into the ladder.

## Observations

### Developed starting state and quadrature

The direct coarse Newton preparation stalled in its Krylov solve at $1.18\times10^{-9}$. Increasing the finite-difference perturbation to $4\times10^{-6}$ and restarting from the saved candidate converged at $2.31\times10^{-11}$. The .00025 hill solved at 6×8, but the independent outer and inner reintegration found a scaled k residual of $2.07\times10^{-7}$. Re-solving at 8×12 reduced it to $1.89\times10^{-11}$. At height fraction .002 the 8×12 independent momentum residual was $1.22\times10^{-9}$. Re-solving that state at 10×16 gave the stationary anchor shown below. The settled outer logarithmic profile and known physical residual remain in the retained inlet construction. Stationarity here refers to the developed discrete solution.

A full-depth Piola transfer had initially rejected a point well above the hill at $r=5.52$. Tapering the correction localized the perturbation. Removing the two one-sided Piola maps' tangential jump at matching let a .00515 seed clear the initial guard. Independent flux and trace checks at 129 stations found agreement within $1.6\times10^{-16}$ for the tested map pairs. A steady solve still approached the principal-response floor in the upstream collar. Two physical BE/BDF2 steps passed after the change in height, supporting continuation with explicitly reset histories.

| State | Native quadrature | Independent quadrature | Momentum | Continuity | k | Outer epsilon |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Developed .002 | 10×16 | 20×20 | 2.61e-11 | 3.59e-16 | 9.09e-11 | 1.22e-11 |
| Nested streamwise .002 at .04 | 12×20 | 24×24 | 8.77e-14 | 5.94e-16 | 5.71e-12 | 3.80e-13 |
| Nested radial .002 at .04 | 12×20 | 24×24 | 7.92e-14 | 5.72e-16 | 5.71e-12 | 3.37e-13 |
| Independent limits | | | 1e-9 | 1e-9 | 1e-8 | 1e-8 |

At .002 the developed state has $K/K_{\rm in}\in[.9598766966,1.0283625561]$, margins .0598766966 and .0716374439 to the certified table edges, and a minimum native principal margin of $2.19\times10^{-3}$. The 20×20 independent outer principal margin is $1.83\times10^{-3}$. The change to a strict nested streamwise grid held height fixed and reduced its weak continuity residual from $7.11\times10^{-4}$ to $1.50\times10^{-12}$. The projected velocity, dependent epsilon trace and rebuilt physical history then passed the independent 24×24 evaluations. Wrong source height, generic positive-time import, same-grid same-height reseeding, and coarsening were exercised as rejection controls.

### Recovery and location of the outer loss

A 10×16 coarse state at height fraction .009686 returned a positive native principal margin of $3.23\times10^{-5}$. Evaluating the same polynomials on 20×20 gave $-2.28\times10^{-4}$ in the outer layer. At the nearby .009749 trial, the final guard rejects a nonlinear solve with maximum refreshed residual $3.51\times10^{-11}$. The independent directional differences recover the sign and magnitude of the reported $-1.56\times10^{-5}$ response, with positive covariance. Finer native quadrature moved the failed height downward, and the independent scan excludes the .009686 result as a seed.

The coarse .00767 state passed a same-signature cold continuation. At physical time .16 the independent residuals were `[2.58e-11, 2.61e-16, 9.34e-11, 9.16e-13]`. The independent outer margin declined from $3.36\times10^{-4}$ at .04 to $1.81\times10^{-4}$ at .16, $7.44\times10^{-5}$ at .24 and $-2.88\times10^{-5}$ at .32. Height and numerical signature were fixed. The loss therefore occurs in the represented evolution even without a new height transfer.

The .30-code-unit outer seed taper raised the highest coarse independently admissible height to .01019 (0.4076 mm). A .55 taper also reached .01019, with a smaller margin. The 2.50 taper reached .00767. The .30 taper's 24×24 outer scan at .01019 follows the same pattern at two timesteps on the same native rule. Both use a new BE startup from the same remapped preceding geometry. The .02 history after its first two steps resumes by ordinary cold restart.

| Time at .01019 | $\Delta t=.02$ outer principal | $\Delta t=.01$ outer principal |
| ---: | ---: | ---: |
| .04 | 1.309e-4 | 1.316e-4 |
| .12 | 1.929e-5 | 1.909e-5 |
| .16 | -3.860e-5 | -3.888e-5 |

Near loss, the witness is at chart $a\simeq2.49$, $r\simeq1.044$ in the transported outer region, just upstream of the downstream hill foot. Its physical elevation is roughly 45 mm above the original floor. The whole-face reference interval and scalar positivity retain margin. The decrease survives these changes to the timestep and integration rule. The independent sampled response tests the current represented gradients and tensor. The comparison leaves full physical spatial and timestep convergence for C2.

The current numerical source reproduces byte-identical uninterrupted and cold flat and manufactured-curved checkpoints within each test, at $t=.06$. Independent curved residuals are approximately `[8.7e-16, 2.6e-17, 4.3e-14, 3.1e-14]`. The new full-height geometry probe gives minimum normal-tube stretch .81248, effective $J$ .75997 and an outer determinant lower bound .76960. The independently reconstructed face measure agrees within $8.89\times10^{-16}$. Both curvature limits and the rigid-rotation checks remain covered by the [governing verification](wj-h27-c1-normal-hill.md).

### Refined continuation and full-height witness

| Grid and protocol | Height fraction | Principal margin at independent volume rule |
| --- | ---: | ---: |
| Four columns, .30 taper | .01019 | 1.309e-4 |
| Eight columns, .30 taper | .01019 | 1.369e-3 |
| Eight columns with local radial splits, .30 taper | .01019 | 1.567e-3 |
| Eight columns, .30 taper | .01523 | 3.940e-5 |
| Eight columns with local radial splits, .30 taper | .0158348 | 8.350e-5 |

On the eight-column grid, shortening the taper from .30 to .05 at the same .018254 preflight improved the sampled principal response from $-4.01\times10^{-4}$ to $-3.48\times10^{-4}$. Both preflights rejected before the first step. For the eight-column path, the .015482 candidate failed its second physical step at the downstream foot with maximum refreshed residual $1.28\times10^{-12}$. The native principal response is $-3.06153\times10^{-5}$. Independent fourth-order directional differences at two spacings give $-3.06153\times10^{-5}$ in both checks. The first BE state and its recovery file were preserved. The fixed-height .014222 cold run later reached the same outer location and rejected at .34, with a sampled margin of $-7.21\times10^{-6}$ after advancing from the .04 history. The outer principal scan was $1.04\times10^{-4}$ at .16, $-2.01\times10^{-5}$ at .24 and $-9.48\times10^{-5}$ at .32.

The locally radially refined .016238 trial reached .02 and rejected the next step at a refreshed residual of $1.24\times10^{-12}$. Its principal response was $-1.07257\times10^{-5}$, confirmed by both independent directional differences. The failed witness had moved from chart radius 1.044 to 1.117 near the same downstream foot. Both refined rejections preserved the accepted BE checkpoint byte for byte. The nested 14-column transfer at fixed .01019 reached .02 and rejected the next step in the downstream-foot outer region. Its maximum refreshed residual was $1.15\times10^{-12}$ and its principal response was $-3.58\times10^{-5}$, confirmed independently. The combined 14-column and local radial ladder passed through .0083, with independent outer principal margin $3.56\times10^{-5}$. Its .008615 trial failed, bringing the next increment below the declared minimum after eight attempts. The combined 14-column and local radial mesh prepared from .002 also rejected its .01019 BDF2 candidate at the downstream foot. Its refreshed residual was $1.80\times10^{-11}$ and its principal response was $-5.78\times10^{-5}$. Together these controls distinguish an evolving local stress response from the nonlinear stopping rule. The changed location and margin with refinement limit the interpretation to the computed states.

The actual 40 mm preflight used a separate large jump from the audited .014222 radially refined state. It reduced weak continuity from $5.64\times10^{-2}$ to $7.05\times10^{-13}$. The inner-collar guard then failed at chart $(a,r)=(.0115246,.0906210)$ near the crown. The native principal response was $-5.99109225\times10^{-3}$, and the independent fixed-state directional differences agree to $3\times10^{-12}$. Covariance over k was .36118. No 40 mm physical time was accepted by this continuation. The full-height travel target remains 1.4741618623 code units at the retained approach velocity .6783515607, or $t=1.48$ on the .02 timestep grid. The prior direct full-height .02 BE witness is documented separately in the governing report. To reproduce the continuation preflight, use the radial fixture with `--height 1 --from-height .014222 --steps 0 --quadrature 12 20 --seed-taper .30`. The source packet is `nested-radial-track/04-height-0.0142220-q12x20/state.result.json`.

The archived radially refined .0158348 state gives the following short rejection replay with the 7c3d6737 executable. The expected solver status is 1. The held .02 state and the failed BDF2 residual are written in the case directory. In a new directory the replay reproduced the failure point, all six refreshed residuals and the held checkpoint byte for byte.

```bash
accepted="$archive/nested-radial-track/07-height-0.01583480-q12x20/state.result.json"
python3 "$archive/make-local-foot-radial.py" hill "$RUN/reject" --height .016238 \
  --from-height .0158348 --initial-fields "$accepted" --quadrature 12 20 \
  --seed-taper .30 --dt .02 --steps 2 --science 2 --step-seconds 180
status=0
run_case "$RUN/reject" || status=$?
test "$status" -eq 1
python3 "$tests/c1_assess.py" --table "$table" --guard "$RUN/reject/run.log" \
  --out "$RUN/reject/independent-guard.json"
```

### Refined physical fluxes, budgets, and cost

Across the 17 ladder states accepted by both independent checks, the whole-interface $K/K_{\rm in}$ interval is $[0.95957035,1.02840883]$. The minimum lower and upper table-edge margins are 0.05957035 and 0.07159117, and the observed cutoff range is $[0.02423070,0.02464546]$. This range stays within the calibrated table.

The .0158348 case has 760 cells, two accepted physical steps to .04, and maximum native refreshed residual $1.15\times10^{-12}$. Its 24×24 scaled equation residuals are `[7.23e-14, 5.61e-16, 5.60e-12, 3.18e-13]`. The dependent trace projection defect is $1.97\times10^{-16}$. The pointwise epsilon trace representation error is $1.81\times10^{-7}$. The whole-interface $K/K_{\rm in}$ bounds are $[.9607621,1.0221991]$ and the margins to the certified edges are .0607621 and .0778009. The cutoff lies in $[.02426617,.02463785]$. Minimum k and outer epsilon are .0069250 and $1.74427\times10^{-5}$. The native covariance and principal margins are .271326 and $1.57048\times10^{-4}$. The independent outer principal margin is $8.35005\times10^{-5}$. The independent outlet minimum velocity is outward at $1.70\times10^{-6}$. A same-signature ordinary cold restart advanced this BDF2 history to .06 and reproduced the predecessor coefficient hash. Its 24×24 equation residuals are `[7.10e-14, 5.53e-16, 5.70e-12, 3.13e-13]`, with an outer principal margin of $3.69\times10^{-5}$.

Both matching owners give signed mass flux $-6.17783\times10^{-4}$ toward the outer side. Their raw k fluxes are $-6.22921\times10^{-5}$ and $-6.22956\times10^{-5}$. The largest weak free matching row is $2.41\times10^{-16}$, and the inner-layer weak k budget defect is $4.82\times10^{-16}$. The exact discrete mean-power accounting defect is $7.25\times10^{-16}$. It includes the independently measured weak divergence work $3.12\times10^{-9}$ and the BDF kinetic-work correction $-5.79\times10^{-10}$. Pressure, k and convection sampled work are $-1.57708\times10^{-5}$, $2.32394\times10^{-6}$ and $-6.17757\times10^{-6}$. The integrated Cartesian pressure, k and convection forces are respectively $(-4.72821\times10^{-5},-7.92638\times10^{-6})$, $(8.23434\times10^{-6},-2.37416\times10^{-7})$ and $(-3.91870\times10^{-5},1.24710\times10^{-5})$. The manifest contains both raw Cartesian momentum fluxes, the traction jump, the one-sided outer epsilon reaction, molecular dissipation, production and physical boundary work.

The .0158348 run took 47.66 wall seconds including construction, initialization, both solves, recovery and scientific output. Peak RSS was 766,792 KiB. Startup and BDF2 complete step times were 16.82 s and 20.98 s. The measured setup and initialization were 3.34 s and 5.54 s. Recovery took .117 s, scientific frames .134 s and final export .079 s. These are included in the complete run time.

| Disjoint work in final BDF2 step | Wall seconds |
| --- | ---: |
| Constitutive evaluation | 4.127 |
| Wall evaluation | .974 |
| Quadrature and assembly | 8.900 |
| Factor rebuilds | 1.324 |
| Coupled momentum and pressure application | .066 |
| Scalar application | .014 |
| Constitutive guards | 4.540 |
| Refreshed residual | .796 |

A same-signature cold extension of the .014222 case, through its failed attempt, took 311.41 s and peaked at 788,120 KiB. These timings characterize the sampled shallow grids and include rejected work where indicated. The production allowance applies to complete measured trajectories at the established horizon.

## Gate decision and next package

| Flag | Status |
| --- | --- |
| `CLEAN_IMPLICIT_BASELINE` | PASS |
| `CURRENT_CURVED_GOVERNING_VERIFIED` | PASS |
| `REAL_HILL_PILOT_PASSED` | BLOCKED |
| `C2_READY` | NO |
| `LONG_FLAT_QUALIFIED` | OPEN |
| `PRODUCTION_AFFORDABLE` | OPEN |
| `EXPERIMENT_VALIDATION_COMPLETE` | OPEN |

The obstruction is a signed local principal-response loss on the computed paths. Its location, timing and height depend on representation and the starting fields. Existence of another admissible full-height trajectory remains open. The next C1 work should test a divergence-compatible predictor that controls the resolved velocity-gradient change near the downstream foot as height changes. Use the replayed .016238 BDF2 rejection and the full-height crown-collar preflight as discriminating controls. Keep the same constitutive equation and admissibility floor. A constitutive amendment would require evidence that the loss persists for a spatially resolved state on an independently prepared path.

After C1 closes, C2 should measure physical spatial and temporal sensitivity at common times, test top and outlet truncation, and qualify the developed whole-interface table range. The developed, fresh-start and signed-perturbation flat classes retain the qualification horizon 26.13350407707765. The production horizon is 505.02466521278404, with the terminal window beginning at 447.8483832973273. Complete trajectory cost must satisfy $1.25 C_{\rm total}\le86{,}400$ wall seconds with rejected work, recovery, output and assessment included. H1 applies the qualified physical hill to the experimental comparison.
