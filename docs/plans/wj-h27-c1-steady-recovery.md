# WJ–EARSM/AWF C1 steady recovery

This campaign tests whether an admissible developed state of the original WJ–AWF equations can bypass the failed physical startup. The branch started at `69b6a8ca2c0e1da5abc0ac886e1c67e7cb428417`, which matched the fetched remote at inspection. Evidence includes complete rejected fields, a bounded search of the original steady equations, and a separate diffusion homotopy. Height fraction `g=1` means the 40 mm cosine-squared hill. The inherited flat and curved governing verification retains its scope. [C1 governing report](wj-h27-c1-normal-hill.md) · [Machine-readable campaign evidence](wj-h27-c1-steady-recovery-evidence.json)

## Model and code

| Implementation | Treatment in this campaign |
|---|---|
| `WJNativeCoupling`, `WJNativeTimeNewton` | Retain the finite-viscosity WJ covariance, lower wall and collar, `ScalarCdV2`, outer-only epsilon, current momentum and scalar advection, and all three sampled AWF actions. Keep the physical $G_h-C_\Pi$ pressure block, common matching measure, and the $A$K$/H$ trace lift. |
| `WJNativeCoupled`, `WJNativeInitialFields` | Compile an explicit `DiagnosticSteady` route from coefficients to free unknowns without constructing a history. Typed packets carry role, attempted time or nonlinear iteration, numerical and geometry identities, coefficients, and the input coefficient or time-history identities. Diagnostic import requires the packet digest. Physical restart accepts the existing signed history format. A packet is published atomically and its destination cannot be overwritten. |
| `WJNativeSteadyKrylov` | Shares feasible Newton and flexible GMRES with the coupled time solver. In diagnostic mode all nonlinear merit and convergence checks use the freshly evaluated steady equations. The principal-response threshold is observed for the returned diagnostic field. Root/domain, covariance, scalar positivity, whole matching $K$ range, geometry, and outlet guards apply to the diagnostic search. |
| `WJNativeScalar` | Optional dyadic Bernstein subdivision gives a tighter lower bound for represented diagnostic scalars. The ordinary physical guard still uses the one-cell bounds. Both bounds are reported for final diagnostic candidates. |
| `c1_oracle.py`, `c1_margins.py`, `c1_rejected_points.py` | Reintegrate signed Q3/Q2 fields with independent chart metrics and companion WJ roots. Candidate inputs need explicit diagnostic mode. The point replay uses four-point directional differences at three spacings and samples both sides of selected faces. |
| `c1_steady_case.py`, `c1_steady_continue.py` | Pin input and session identities. Predict with the existing normal-chart velocity remap and weak continuity projection, rebuild the dependent trace, and correct at fixed geometry. Reintegrate equations and response before choosing the next diagnostic seed. Geometry preparation advances zero physical time. |

The correction equation uses $J_F + M/\Delta\tau + \beta L_+$. Here $M$ is the constrained consistent Galerkin mass on the differential variables, including the variation of the dependent epsilon lift. Pressure, continuity and pressure gauge rows have zero additive action. $L_+$ is the unit-coefficient symmetric velocity-gradient form with homogeneous essential corrections and a single coefficient across matching. With $R=\lVert F(X_m)\rVert$ and $R_0$ its value at the first nonlinear base, the frozen reciprocal pseudo-time receives a factor $\mathrm{clamp}(R/R_0,10^{-8},1)$. The correction diffusion receives the square root of the same factor. Arnoldi, the restarted residual and the freshly evaluated linear checks use the additive action. `WJSteadyPreconditionerMass` remains an inverse approximation control. [Kelley and Keyes](https://epubs.siam.org/doi/10.1137/S0036142996304796), [PETSc pseudo-transient form](https://petsc.org/release/manualpages/TS/TSPSEUDO/)

The numerical initialization problem $F_\eta$ adds $\eta\nu_{\mathrm{reg}}(G+G^T)$ only to the conservative momentum stress. The profile is prescribed throughout each solve. With $x$ normalized between the exterior planes and $z$ normalized from the flat effective datum to the physical top, its shape is $[4x(1-x)]^2(1-z^8)^2$. A separately declared normal-radius ramp can suppress the coefficient in the highly sheared lower layer. This ramp reaches one below matching. The Cartesian part of the profile agrees on both sides, and the exterior inlet, outlet and top receive zero added traction. In the ramped case the coefficient also vanishes at the essential ground. The matching traction contains the extra stress. An original-model audit on the same coefficients supplies the original traction defect. The scalar equations retain $P_0=-R_0:G$. Mean power includes the extra loss $\eta\int \nu_{\mathrm{reg}}|G+G^T|^2/2$. Independent $F_0$ and $F_\eta$ audits use the same coefficients and label the tested residual. A rank-one incompressible planar perturbation gains exactly $\eta\nu_{\mathrm{reg}}$ in raw directional response at fixed profile. Normalization always uses the original $\mu$.

The reference remains `tables-exact-inner-v1.json` with SHA-256 `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`. The same physical top, `ZeroTransverseV1` outlet at $a=15$, finite-viscosity inlet, $H=.15$, fixed $h=.11541061001735275$, and certified whole $K/K_{\rm in}\in[.9,1.1]$ apply. The complete governing-note identity and physical coordinate map are recorded in the C1 governing report. [C1 governing report](wj-h27-c1-normal-hill.md)

## Rejected fields

The archived executable reproduced the $g=.016238$, $t=.04$ checkpoint byte for byte. Its SHA-256 is `3adfa6e43e4d9f2a012446a5b6f5af4b9046f116f21e1aded3f800d7598a884b`. The instrumented build produced the same coefficient hash `e662cb07872410869a2aa12c136483e63bed78b08830dd5464a1e99171e5f387` under its own numerical signature. Its cold BDF2 attempt at $t=.06$ stored candidate `c9de1957a40475ba2368aaca490ec17d8831abd0d8d94a954cee53f7db1a8654` and preserved its $t=.04$ checkpoint byte for byte. The recorded refreshed rejected residual is 1.03039e-12. The order-32 independent momentum, continuity, $k$ and outer epsilon residuals are 7.07e-14, 6.83e-16, 5.75e-12, and 3.16e-13. The accepted order-32 outer principal margin is +2.29833e-5.

| Control and sample | Chart $a,r$ | Principal divided by original $\mu$ | Raw response |
|---|---|---:|---:|
| Rejected BDF2, first native sample | (2.4969920,1.1078989) | -3.46787e-6 | -7.52685e-8 |
| Rejected BDF2, hill-side foot limit | (2.5,1.1078989) | -3.11163e-5 | -6.75415e-7 |
| Rejected BDF2, flat-side foot limit | (2.5,1.1078989) | +2.62955e-3 | +5.62438e-5 |
| Rejected BDF2, 49-point-per-axis local scan | (2.5,1.0968308) | -3.40997e-5 | -7.32816e-7 |
| Full-height crown preflight, first native sample | (0.0115246,0.0906210) | -5.99109e-3 | -6.62477e-6 |
| Full-height crown preflight, 49-point-per-axis local scan | (0,0.0917969) | -3.34975e-2 | -3.57772e-5 |

The first foot sample has $\chi=1$, zero physical auxiliary viscosity and a covariance-to-$k$ ratio of 0.277125. At that point, keeping the accepted turbulence scalars and inserting the rejected gradient gives -6.73228e-6. Keeping the accepted gradient and inserting the rejected scalars gives +4.06161e-5. The turbulence time scale $k/\varepsilon$ changes from 28.013972 to 28.013671. The gradient change is sufficient for the sign change at this point. The three-spacing spread for the rejected point is 1.12e-13 in normalized response. The one-sided and radial series, physical coordinates, covariance, WJ invariants, and minimizing directions are in the retained point reports.

The crown preflight has $\chi=.653245$ and positive covariance at its first sample. The complete order-20 scan finds (-.388262) in the inner region and (-.361006) just outside matching. These are properties of the rejected initial polynomial. No time step was attempted. A targeted refinement of the foot startup changes the first failing point to (1.87350, .606067) before time advancement, with principal ratio (-.00161988). This refinement exposes sensitivity of the projected seed. The sign of the refined $t=.06$ BDF2 field remains unresolved.

The ordinary `WJRestart` control rejected the diagnostic packet at checkpoint parsing. The explicit diagnostic coefficient importer used the same BDF2 candidate in a steady search. The side-by-side JSON uses distinct schemas for the attempted coefficients and the accepted history. Empty residuals in a preflight packet mean no equation solve was attempted.

## Numerical attempts and gates

The independent residual limits are $10^{-9}$ for momentum and continuity, $10^{-8}$ for $k$ and outer epsilon, and $10^{-10}$ for the projected dependent trace. The fresh native steady maximum is limited to $10^{-10}$. Native physical acceptance retains strictly positive whole-cell scalar bounds and a covariance-to-$k$ ratio greater than $10^{-10}$. The original planar response divided by $\mu$ must exceed $10^{-8}$ after allowance for the independent directional spacing uncertainty.

The flat and curved zero-shift diagnostic fixtures reproduced all coefficient hashes from the ordinary steady reference. Shifting the correction mass and diffusion gave original-model maxima 1.67e-11 and 5.56e-13. The curved action check gave a dependent trace derivative relative discrepancy 3.04e-11. Its four independent equation residuals are 5.68e-14, 2.41e-17, 3.49e-13, and 5.56e-13. The analytic scalar control gives an accumulated correction of $1/3$ with a zero shifted linear residual and original nonlinear residual $2/3$. Its independent direction action is 3. The Bernstein control certifies the positive quadratic after subdivision and rejects its sign-changing counterpart.

The final clean build retained the curved action result, and the three-step manufactured trajectory and its cold restart from step two give identical final histories with SHA-256 `08afd9d2e75c1eb7a5c23d4c731d3254fa6f1823385e9d843a95f71e827701c1`. The physical advancement in the foot replay uses `CoupledJFNKV1` with the Sparse-LU mixed flow inverse and new-time advection from BE startup. The T=10 flat Sparse-LU run used the same coupled Newton family. Its first 16 steps used extrapolated advection, followed by an implicit suffix implemented with the archived patch. The steady search here solves the zero-time-derivative equations. [T=10 solver record](wj-h27-bf-t10-status.md)

### Original-residual path

| Geometry and changed control | Fresh original maximum | Independent original minimum principal ratio | Observed result |
|---|---:|---:|---|
| $g=.002$, developed anchor re-evaluated on 10×16 | 2.16e-11 | +1.83e-3 at order 20 | Stationary original anchor. The four independent equations pass. |
| $g=.004$, unshifted from anchor | 5.72e-5 | Native −1.31e-2 | The line search rejected 125 scalar bounds. The one-cell lower bound on $k$ was 9.94e-10 while a 41×41 per-cell scan found minimum $k=0.00761$. |
| $g=.004$, added correction mass with initial pseudo-time 5 | 7.41e-5 | Native +1.56e-3 | The second correction exhausted feasible Armijo trials. |
| $g=.0024$, smaller height change at 10×16 | 1.15e-12 | −4.61e-3 at order 20 | Original steady root. The upstream collar at $a\approx-8.17$ sets the adverse sign. |
| $g=.0024$, same mesh solved at 16×24 | 3.60e-11 | −4.27e-3 at order 24 | The reintegrated equations pass and the upstream sign persists. |
| $g=.004$, same seed as the unshifted attempt, with Bernstein subdivision | 1.75e-12 | −2.88e-2 at order 20 | Original steady root with a positive subdivided scalar certificate. Its one-cell $k$ bound is −.00315. |
| $g=.0046$, original ladder at 10×16 | 1.10e-12 | −3.70e-2 at order 20 | Independent momentum is 1.19e-9 against 1e-9. |
| $g=.0046$, re-solved at 16×24 | 3.46e-11 | −3.65e-2 at order 24 | Independent equations pass. The negative inner response persists. |
| $g=.0052$, original ladder resumed on 16×24 | 7.43e-13 | −4.42e-2 at order 24 | Independent equations pass. Inner and outer samples include negative response. |
| $g=.00571$, same-height diagnostic resume | 3.32e-11 | −5.04e-2 at order 24 | Independent equations pass. |
| $g=.00622$, same-height diagnostic resume | 1.07e-12 | −5.63e-2 at order 24 | Independent equations pass. The upper $K/K_{\rm in}$ margin is .00033815. |
| $g=.016238$, rejected BDF2 candidate as an explicit diagnostic guess | 2.12e-3 | Native −2.74e-2 | Correction and finite-difference probes encounter the table range. This is a steady residual on a formerly time-dependent field. |

The $g=.0024$ targeted mesh adds abscissae at −8.125, −6.25 and −4.375 and bisects radial cells through the adverse collar. The original steady residual fell from $1.65\times10^{-3}$ to $2.61\times10^{-4}$, then to $1.20\times10^{-4}$ after one bounded restart. The constrained correction reduced it to $3.28\times10^{-6}$. Restarting with the last effective shift failed to improve the base norm in the declared GMRES cycles. Restoring a stronger correction mass gave $7.28\times10^{-7}$ after two accepted corrections and exhausted the wall budget during the third linear solve.

At that refined candidate, order-24 momentum, continuity, $k$ and outer epsilon norms are $8.10\times10^{-8}$, $4.50\times10^{-13}$, $7.28\times10^{-7}$ and $4.73\times10^{-7}$. The projected trace defect is $2.24\times10^{-16}$. The original inner and outer principal minima are +.005226 and +.001724, with directional-spacing allowance below $5.90\times10^{-11}$. Whole-face $K/K_{\rm in}$ lies in $[.96571458,1.09056757]$, with cutoff $[.02389441,.02460642]$. The native whole-cell lower bounds for $k$ and outer epsilon are .001413 and $1.74425\times10^{-5}$. These finite measurements make it the best positive-response refined candidate. Its equation residuals leave steady convergence unresolved.

The $g=.00622$ original root has $K/K_{\rm in}\in[.91716802,1.09966185]$, with the maximum at chart $a=.557466$, physical $(x,z)/L_0=(.558064,.242987)$. A $g=.00634$ probe reduces the residual to $4.11\times10^{-6}$ and halts under the wall cap. Its whole-face upper margin is $8.21\times10^{-7}$. The line search records 36 table-range rejections and decreases the last accepted step below .0005. The independent audit corroborates the residual. This bounds the tested coarse original branch search near the table edge. The branch already fails the principal guard before that edge.

### Temporary momentum diffusion

| Geometry and dose | Fresh $F_\eta$ maximum | Fresh $F_0$ maximum at the same field | Minimum independent original / modified principal ratios |
|---|---:|---:|---:|
| $g=.002$, $\eta=1$ | 9.98e-4 | 9.98e-4 | −.04763 / −.01511 |
| $g=.002$, $\eta=.05$ | 7.44e-13 | 7.07e-5 | +.00184 / +.00199 |
| $g=.002$, released to $\eta=0$ | 9.53e-13 | 9.53e-13 | +.00183 / same |
| $g=.0024$, $\eta=.10$ | 3.21e-11 | 1.37e-4 | +9.23e-4 / +1.24e-3 |
| $g=.0024$, $\eta=.075$ | 9.48e-13 | 1.04e-4 | −2.22e-4 / +2.29e-5 |
| $g=.0024$, $\eta=.025$ | 1.32e-12 | 3.59e-5 | −2.87e-3 / −2.82e-3 |
| $g=.0024$, literal $\eta=0$ | 2.18e-11 | 2.18e-11 | −4.27e-3 / same |
| $g=.0024$, $\eta=.10$ with lower-layer ramp | 4.25e-11 | 9.79e-6 | −5.27e-3 / −5.04e-3 |
| $g=.00308$, $\eta=.25$ | 4.58e-11 | 3.14e-4 | −.03275 / −.02549 |

The initial $\eta=1$ path stalled in covariance, scalar and table trial guards. The smaller $\eta=.05$ dose converged at the developed anchor. Fixed-dose continuation attained $g=.004236$ with passing modified equation audits and negative response. At $g=.0056232$ the native solve converged, with independent quadrature still unresolved. The dose ramp at $g=.0024$ suppressed forcing at the lower wall. Its modified root converged while retaining the negative sign. The global-profile release at $g=.0024$ recovered the adverse upstream root. On that tested path, the independent modified margin crosses zero between $\eta=.025$ and $.075$, and the original margin crosses between $\eta=.075$ and $.10$. At $\eta=.10$ the native extra mean-energy loss is $2.4483495186916756\times10^{-4}$. Independent integration agrees within $3.3\times10^{-17}$ and its power accounting discrepancy is $2.24\times10^{-16}$.

At full height the original steady solve exhausted its first linear-solve budget at maximum residual .05380. A separate full-height attempt added the radial ramp over $[.06,.089453125]$. The independent preflight prescription uses a 30% reserve on measured raw deficits. Dose 8 covers the crown and matching samples. Dose 50 covers the next interior deficit of approximately −.00199 at physical $(x,z)=(-1.363,2.069)$. The exterior-zero profile still leaves adverse outer samples near the inlet. The diagnostic dose scan showed that scaling this profile enough to correct the near-inlet witness would require $\eta>1.65\times10^7$. The subsequent scan still found a negative sample on the same coefficients. The bounded solve used dose 50. Four accepted corrections had steps of order $10^{-6}$. It rejected 64 trial states that left the certified matching $K$ range and had an independent modified maximum of .05380.

The stronger correction mass at dose 50 reduced the modified residual to $.0504903$ and exhausted feasible Armijo trials. The independently evaluated modified momentum, continuity, $k$ and outer epsilon norms are $.0186588$, $6.26\times10^{-9}$, $.0504903$ and $.0397574$. Its whole-face $K/K_{\rm in}$ range is $[.93821014,1.099999995]$. The limiting trial face is geometry 540 on chart $0\le a\le1.25$. The retained maximum is at $a=.782017$, physical $(x,z)/L_0=(.892046,.987860)$. The candidate is $5.0\times10^{-9}$ from the certified upper edge. Sampled scalars and covariance are positive, the whole-cell $k$ certificate is negative, and the independent modified principal margin is −.04939. The binding table edge is evidence about the tested modified problem. The original $F_0$ branch may have different trace support.

At the accepted shallow $t=.04$ state, the order-32 projected epsilon trace defect is $2.75\times10^{-16}$ and its pointwise approximation error is $1.80\times10^{-7}$. The inner weak $k$ budget defect is $4.53\times10^{-16}$ and the power accounting discrepancy is $2.65\times10^{-16}$. Raw one-sided traction and $k$ diffusion jumps are $2.85\times10^{-8}$ and $6.53\times10^{-10}$ in RMS. The attempt manifest gives the corresponding values on each audited candidate, including the modified and original stress tractions on the same field. A small power accounting discrepancy checks the identity even for a field with a large equation residual.

| Complete Linux process | Wall seconds | Peak RSS, KiB |
|---|---:|---:|
| Instrumented shallow BE and BDF2 through .04 | 108.02 | 1,231,416 |
| Instrumented cold BDF2 rejection including separate packet | 77.16 | 1,171,216 |
| Original developed anchor | 8.08 | 237,932 |
| Original $.004$ subdivision attempt | 175.35 | 408,992 |
| Original refined mesh unshifted resume | 550.35 | 762,536 |
| Original refined mesh bounded correction shift | 430.71 | 759,288 |
| Original refined mesh stronger mass at the last viable state | 365.54 | 739,212 |
| Full-height $F_0$ steady diagnostic | 167.62 | 834,508 |
| Full-height $F_{50}$ steady diagnostic | 433.43 | 964,444 |
| Full-height $F_{50}$ stronger correction mass | 117.95 | 964,716 |

The replay spent 7.75 s in construction, 11.76 s in initialization, 39.40 s on BE and 48.03 s on the accepted BDF2 step. The BDF2 step attributes 9.90 s to constitutive evaluation, 2.09 s to wall work, 21.55 s to assembly, 1.29 s to inverse rebuilds, .063 s to mixed solves, .014 s to scalar applications, 10.89 s to guards and 1.98 s to residual work. Recovery and scientific frames over both steps cost .42 s within the complete step timings. Final export took .095 s. The steady ledger records residual and validity evaluations, preconditioner preparations and applications, line trials and Arnoldi time. The initial continuation budget was 3,600 elapsed seconds from its launch after the archived replay. Measured process times sum to 6,409.27 s across 55 solver invocations, 2,914.81 s for the recorded independent audits and 548.06 s for the isolated builds. These totals include overlapping processes and do not estimate production-trajectory time. The original physical cost requirement remains $1.25C_{\mathrm{total}}\le 86400$ wall seconds per qualified production trajectory. The full-height target remains $t=1.48$ under the original equations. Preparation and steady iterations advance zero physical time.

## Decision

| Gate | Result |
|---|---|
| `REJECTED_CANDIDATE_REPLAY_VERIFIED` | PASS |
| `ORIGINAL_RESIDUAL_STEADY_SEARCH_VERIFIED` | PASS |
| `FULL_HEIGHT_ORIGINAL_STEADY_STATE_FOUND` | BLOCKED |
| `ZERO_REGULARIZATION_ENDPOINT_VERIFIED` | BLOCKED at full height |
| `FULL_HEIGHT_ORIGINAL_PHYSICAL_PILOT_PASSED` | BLOCKED |
| `C1_CLOSED` | NOT MET |
| `C2_READY` | NOT MET |

The best independently qualified developed original state is still $g=.002$ (0.08 mm). The best qualified dynamic state remains $g=.016238$ (0.64952 mm), $t=.04$. Its later cold rejection is reproduced. The original steady path establishes that the loss of response can occur at an independently reintegrated steady root, and therefore cannot be attributed solely to startup or timestep error. The upstream location, the difference from the transient foot witness and the positive-response refined candidate leave spatial representation and alternative steady branches open. The full-height modified candidate remains far from a root. There is no demonstrated admissible full-height state from which to release the required physical pilot, signed disturbances or timestep comparison.

The next C1 work is to improve the coupled correction inverse on the saved positive-response refined field and qualify that original root on both its native and enriched integration rules. The last 355 s of diagonal-mass correction yielded two accepted steps and a final reduction factor of .85. Repeating the same recipe at that rate gives poor leverage for the remaining four orders of magnitude. A verified coupled scalar and momentum correction including the sampled loads and dependent trace would test whether the remaining limitation is linear conditioning. Any proposed extension to the cutoff family should follow an admissible original branch that actually needs the coverage. The inherited CLEAN_IMPLICIT_BASELINE and CURRENT_CURVED_GOVERNING_VERIFIED gates retain their verified scope. LONG_FLAT_QUALIFIED, PRODUCTION_AFFORDABLE and EXPERIMENT_VALIDATION_COMPLETE remain open.

## Reproduction and custody

The campaign is under the ext4 workspace `/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-steady-recovery.LOp8mfob` on the authorized Linux server. The manifest records paths, digests and wall time with peak resident memory. Relative payload paths are relative to that run root. Large fields, checkpoints, compile logs and individual independent audit reports stay there. Every solve used one rank, one numerical thread and an explicit timeout. The isolated builds used eight workers. The final executable was compiled with GCC 11.4.0 and CMake 3.22.1. Its Git revision is `e9b1020183c3c631866313e04c2d382e019e20a6`, source fingerprint `6a8ff82511eb8efaac4998027aabcbe7e2110771846bfe400720aca1aeeca091`, executable SHA-256 `c5e3a8f65e2ae08d7598b17b200d920cba3a393e560dbba22899841fe07ca5a6`.

```bash
source /media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh
evidence="$NEKTAR_RUNS/wj-c1-steady-recovery.LOp8mfob"
work=$(mktemp -d "$NEKTAR_RUNS/wj-c1-steady-rerun.XXXXXXXX")
mkdir "$work/tmp"
export TMPDIR="$work/tmp"
table="$NEKTAR_SOURCE/docs/plans/wj-h27-linux/task2r-r1/tables-exact-inner-v1.json"
archive="$NEKTAR_RUNS/wj-c1-seed-gradient.BoAhhRpl"
old="$NEKTAR_WORKSPACE/build/wj-c1-seed-gradient-reach.BoAhhRpl/IncNavierStokesSolver"
bash "$evidence/replay.sh" "$work" "$archive" "$old"

git -C "$evidence/source-v4" worktree add --detach "$work/source" e9b1020183c3c631866313e04c2d382e019e20a6
tests="$work/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
cmake -S "$tests" -B "$work/build" -DCMAKE_BUILD_TYPE=Release \
  "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
cmake --build "$work/build" --parallel 8 --target IncNavierStokesSolver WJH27SteadyKrylovProbe WJH27ScalarProbe
solver="$work/build/IncNavierStokesSolver"

anchor="$evidence/anchor-original"
python3 "$tests/c1_steady_case.py" "$anchor/run.xml" "$work/coarse-original-004" \
  --seed "$anchor/state.diagnostic-steady.json" --height .004 --remap \
  --scalar-subdivision 3 --seconds 210
bash "$evidence/run-capped.sh" "$work/coarse-original-004" "$solver" 245
python3 "$tests/c1_oracle.py" --table "$table" \
  --native "$work/coarse-original-004/state.diagnostic-steady.json" --diagnostic \
  --out "$work/coarse-original-004/equations-q20.json" --order 20
python3 "$tests/c1_margins.py" --table "$table" \
  --native "$work/coarse-original-004/state.diagnostic-steady.json" --diagnostic \
  --out "$work/coarse-original-004/margins-q20.json" --order 20
```

Use the matching archived executable for byte replay. A rebuilt executable carries its current source identity. A fresh instrumented replay makes the $t=.04$ state with that executable before passing its binary history to the same executable for $t=.06$. The manifest identifies the exact original and instrumented sessions and coefficient packets.

The following two inputs exercise the original steady release and the limiting full-height table trial. Their source candidates are the immutable files from the preceding attempts.

```bash
last="$evidence/eta-release-0024/eta-.025"
python3 "$tests/c1_steady_case.py" "$last/run.xml" "$work/zero-at-0024" \
  --seed "$last/state.diagnostic-steady.json" --scalar-subdivision 3 \
  --seconds 190 --eta 0 --regularization-viscosity 6.423173803526447e-5
bash "$evidence/run-capped.sh" "$work/zero-at-0024" "$solver" 225

last="$evidence/full-height-eta50"
python3 "$tests/c1_steady_case.py" "$last/run.xml" "$work/full-height-table-control" \
  --seed "$last/state.diagnostic-steady.json" --scalar-subdivision 3 --seconds 420 \
  --eta 50 --regularization-viscosity 6.423173803526447e-5 \
  --regularization-radius .06 .089453125 --pseudo-time .5
bash "$evidence/run-capped.sh" "$work/full-height-table-control" "$solver" 465
```

The continuation journals record the exact commands for geometry rungs and their independent audits. The saved lower-layer and exterior stress profile study is reproducible with `choose-homotopy-dose.py` under the run root. Its manifest binds the original order-20 margin scan, the complete preflight packet, and each attempted dose.
