# WJ-C1R3 original height progress

## Decision

The original steady equations have an independently admissible root at $g=.003$, a 0.120 mm hill on the 864-cell nested mesh. Both coupled correction policies reached the root from the exact saved coefficients. A fresh original trajectory and two signed downstream velocity disturbances reached $t=.06$. The quiet BDF2 cold restart reproduced the final history byte for byte.

The first larger step, $g=.0042$ or 0.168 mm, converged in the original equations and passed independent integration of the equations. Its one-sided principal response at the forcing-support boundary was negative. The full 40 mm physical pilot remains open. Every solve here used literal $\eta=0$.

## Frozen base and correction

The archived packet, session, and last accepted coefficients have SHA-256 prefixes `6ac39b56633a`, `ac4ec53ed2a9`, and `a6e302be8db4`. The full digests and numerical signatures are in the [evidence manifest](wj-h27-c1-original-height-progress-evidence.json). Both the archived binary and the new binary reproduced the coefficient digest, native principal ratio $+.0039032243$, and these refreshed 10×16 residuals.

| Equation | Frozen residual |
|---|---:|
| Momentum | $3.101160530149\times10^{-6}$ |
| Continuity | $1.79866\times10^{-13}$ |
| $k$ | $7.238399815114\times10^{-6}$ |
| Outer $\varepsilon$ | $1.508549529425\times10^{-5}$ |

The canonical signed inventory $[u,w,\Pi,k,\varepsilon]$ is $[8092,8092,3667,8092,2212]$. The respective free counts are $[7776,7462,3667,7803,2106]$, totaling 28,814. Total outlet traction fixes pressure. The replay of the qualified $g=.0024$ nested control gave $9.02204\times10^{-12}$ on the 16×24 rule. The frozen $g=.003$ packet's whole matching $K/K_{\rm in}$ range was $[.9697797942,1.0367564993]$.

`WJNativeSteadyKrylov` keeps the current algebraic base, scales, safeguards and true residual checks. With `ColoredSparseV1`, only the right inverse may reuse its LU for up to four subsequently accepted bases. A reused solve has a 12-iteration work cap and one cold retry from the same base if its current true linear check fails. Each invocation starts cold. `WJNativeCoupling` prepares the current dependent-trace derivative and correction-only shift before every correction. The mass action is on constrained differential fields, including the lifted outer epsilon direction. Pressure and continuity stay algebraic. The current $J_{F_0}$ and $F_0$ use the original sampled loads, boundary conditions, scalar transport and finite-viscosity WJ/AWF stress. The direction model uses the current Jacobian action.

The curved constrained control used correction-only $\Delta\tau=1000$ and $\beta=2\times10^{-5}$. Rebuild and reuse reached $2.79\times10^{-11}$ and $5.60\times10^{-11}$ respectively. A one-iteration reuse cap detected a $1.06\times10^{-9}$ current linear residual and refreshed once. The cold retry gave $6.11\times10^{-15}$ and reproduced the rebuild coefficients. The C++ probe verified a new cold invocation, known additive action and the original direction model. The Python driver tests separately exercise budget, audit, nonlinear failure and pending custody decisions.

The driver retains the last admitted root, the actual last accepted iterate and any better uncommitted trial as separate records. An improving native deadline resumes the accepted packet at the same geometry. An unfinished audit holds the stage. Reported linear residuals are null with an explanatory reason when the operation did not reach a linear measurement. Physical release continues to require an original root and three digest-bound audit receipts.

The native curved deadline control tested 16 s first and finished both corrections. At 11 s it finished one correction, retained $9.88264\times10^{-8}$ and stopped during the next LU preparation. Its second linear result is null. The driver classified the exit as `BUDGET_EXHAUSTED_WITH_PROGRESS`. Fixed-height import used the exact accepted coefficient hash `464f9d9c026832471a79292c330dbe9240ce3aab9b74ad28efac5e85b26f255c` and converged to $1.16420\times10^{-11}$. A separate analytic deadline probe retained a better evaluated trial distinct from the accepted vector and verified restart from the latter.

### Matched cost at $g=.003$

Each current arm starts from `a6e302be8db4`, with the same 10×16 original residual, fixed norms, finite difference setting, allowable work, merit, direction measurements and pressure constraints. The archived 1,800 s continuation from that exact packet took 518.41 s complete and converged after three full corrections. The current rebuild arm reproduced its coefficient hash.

| Work or criterion | Rebuild each base | Bounded reuse |
|---|---:|---:|
| First accepted $F_0\le10^{-6}$, solver s | 352.72 | 189.62 |
| First accepted $F_0\le10^{-10}$, solver s | 523.10 | 205.05 |
| Complete native process, wall s | 535.60 | 217.56 |
| 16×24 recheck and three independent audit processes, sum of wall s | 257.85 | 266.83 |
| Native solve plus those qualification processes, sum of wall s | 793.45 | 484.39 |
| Maximum native solve RSS, KiB | 1,083,408 | 992,440 |
| Colored preparations / assembly residual evaluations | 3 / 1,290 | 1 / 430 |
| Current Jacobian actions / accepted corrections | 6 / 3 | 11 / 4 |
| Refreshed native maximum at 10×16 | $8.07223\times10^{-13}$ | $8.41597\times10^{-13}$ |

The complete native solve was 2.46 times faster with reuse. The summed cost through qualification was 1.64 times faster, including the richer native and independent checks. Pattern work, factor application, guard callbacks and independent operator times are recorded separately in the manifest. Both routes qualified without a richer-rule Newton update. Absolute canonical differences were $1.65\times10^{-13}$ for the mixed flow vector, $1.93\times10^{-14}$ for $k$, and $1.59\times10^{-15}$ for $\varepsilon$. These timings apply to diagnostic steady preparation on this case.

### Measured nonlinear direction

The prediction uses the freshly checked accumulated direction $d$. It records $F_0+\alpha J_{F_0}d$, the actual refreshed $F_0(X+\alpha d)$, and their difference in the same scaled equation blocks. When a correction shift is active its exact additive action is subtracted from the current linear action. The shifted curved control exercises this subtraction.

At the frozen $g=.003$ base, the first full step predicted a maximum of $1.334\times10^{-8}$ and produced $1.14301\times10^{-5}$. Its actual-to-predicted decrease ratio was 0.243. The dominant outer-$\varepsilon$ free row 78 was shared by cells 632, 633, 641 and 642 near $(a,r)=(-6.25,.15232)$. The measured term remainders there were $3.997\times10^{-8}$ advection, $1.14214\times10^{-5}$ diffusion, $6.739\times10^{-8}$ production and $1.060\times10^{-8}$ destruction. Their signed reconstruction differed from the full remainder by $2.30\times10^{-13}$. Halving the derivative spacing changed the comparison by $1.78\times10^{-10}$. Holding the dependent lift fixed changed the remainder comparison by $2.75\times10^{-5}$. That comparison changes overlapping operator terms and is recorded separately. The lift coefficient remainder was $4.93\times10^{-8}$. Later step ratios were between .997 and 1.000 to the reported precision. Ordinary halving accepted all measured full steps on this base and was retained.

## Scientific admission at $g=.003$

| Measurement | Rebuilt root | Reused root | Limit |
|---|---:|---:|---:|
| Fresh 16×24 native original maximum | $3.96725\times10^{-11}$ | $3.96840\times10^{-11}$ | $10^{-10}$ |
| Independent 28×28 momentum | $5.98245\times10^{-13}$ | $5.94477\times10^{-13}$ | $10^{-9}$ |
| Independent continuity | $6.8230\times10^{-16}$ | $6.8164\times10^{-16}$ | $10^{-9}$ |
| Independent $k$ | $4.00110\times10^{-11}$ | $3.99738\times10^{-11}$ | $10^{-8}$ |
| Independent outer $\varepsilon$ | $9.15456\times10^{-14}$ | $9.19894\times10^{-14}$ | $10^{-8}$ |
| Projected trace defect | $5.40406\times10^{-16}$ | $6.42923\times10^{-16}$ | $10^{-10}$ |
| Whole-cell $k,\varepsilon$ lower bounds | $4.54546\times10^{-3},1.74421\times10^{-5}$ | same at shown precision | positive |
| Native covariance divided by $k$ | .263488 | .263488 | $>10^{-10}$ |
| Minimum native principal response over $\mu$ | $+.00306512$ | $+.00306512$ | $>10^{-8}$ |
| Minimum local signed response after derivative uncertainty | $+.0026786065$ | $+.0026786065$ | $>10^{-8}$ |
| Whole matching $K/K_{\rm in}$ | [.96712285, 1.04047016] | same at shown precision | [.9, 1.1] |

The one-sided, volume and support scans include both cell owners, both matching owners and the moving cutoff. Three-spacing directional checks use independent constitutive differencing. The limiting sample is at the inner forcing-support line, $a=1.875,\ r=h$. The independent outer face minimum is $+.00385685$. Sampling gives finite coverage. On the rebuilt root the minimum Jacobian is .99943743, the minimum outlet normal velocity is $+1.21\times10^{-6}$, and the whole cutoff range is [.02416281,.02459753]. The trace pointwise representation error is $8.83\times10^{-7}$. Discrete mean-power and mass identity defects are $7.79\times10^{-16}$ and $2.64\times10^{-15}$; the maximum weak free interface residual is $4.85\times10^{-16}$. Raw one-sided traction and $k$ diffusion jumps have RMS $3.60\times10^{-8}$ and $9.94\times10^{-10}$ respectively. The manifest includes the force, flux, production and budget terms.

### Scoped original physical release

The native release rebuilt ordinary advection and startup history at $t=0$ from the rebuilt qualified packet. The timestep solver was `CoupledJFNKV1`, with BE startup, implicit BDF2 thereafter and $dt=.02$. The velocity disturbance uses the native Q3 curl of a $C^1$ streamfunction in the flat downstream chart. Its support is $a\in[2.5,15]$ and $r\in[1.04495,1.35626]$. It shares the interface and essential traces. The two signed $10^{-4}$ proposals failed the initial physical principal guard in outer cell 764 without advancing time. Both original rejected diagnostic packets are retained. Reducing both peaks to $2.5\times10^{-5}$ left a sampled speed change $2.45338\times10^{-5}$, gradient change .00181639 and scaled weak continuity defect $3.44\times10^{-22}$.

| Trajectory at $g=.003$ | Time reached | Independent maximum at $t=.06$ | Smallest independent sampled principal | Native wall s |
|---|---:|---:|---:|---:|
| Quiet | .06 | $4.00110\times10^{-11}$ | $+.0026786065$ | 53.05 |
| Cold from $t=.04$ | .06 | same bytes as quiet | same bytes as quiet | 34.66 |
| Positive velocity pulse | .06 | $6.14154\times10^{-12}$ | $+.0026786053$ | 116.85 |
| Negative velocity pulse | .06 | $6.30732\times10^{-12}$ | $+.0026786076$ | 118.74 |

The final quiet and recovered checkpoint hash is `d0112523446a8e6bc5638a9725e5be0ffac5250ab8182341f00e16119eaf0a27`. The disturbed branches remain measurably distinct from the quiet branch at this common time, with maximum mixed coefficient differences $6.50721\times10^{-5}$ and $6.51273\times10^{-5}$. This release tests the nested representation locally over three timesteps.

## Continuation and localized response

The constrained recent-root remap at $g=.0042$ had original maximum $5.10950\times10^{-5}$ and native response $+.00367689$. An explicitly typed secant guess between $g=.0024$ and .003 was freshly evaluated at the same target. Its maximum was $7.04806\times10^{-4}$ and its native response was $-.0107779$. Selection retained the recent-root remap.

The 864-cell 10×16 steady solve reached $4.78923\times10^{-11}$ in 15 accepted corrections and three factors. The independent 28×28 residuals were $[4.644\times10^{-11},3.442\times10^{-15},2.076\times10^{-9},4.794\times10^{-11}]$. The native minimum was $-.00896091$. At $(1.875,h)$ the two streamwise owners returned $-.00974956$ and $-.00913749$, with derivative spread below $4\times10^{-12}$ at the first owner. Its $G_{xz}$ increased from 1.30170 at the qualified .003 root to 1.37966, and $k/\varepsilon$ decreased from 2.87995 to 2.84728. The independent outer face minimum remained $+.00351028$. The adverse result is in the inner support region, with positive scalars and covariance and whole $K/K_{\rm in}\in[.92854339,1.07242687]$.

| $g$, height | Cells and native rule | Original maximum | Angular minimum after uncertainty | Custody |
|---|---|---:|---:|---|
| .0024, 0.096 mm | 970, 10×16 | $1.01\times10^{-12}$ | $-.007206$ | Adverse diagnostic root after same-height resume |
| .0030, 0.120 mm | 864, 16×24 | $3.97\times10^{-11}$ | $+.002679$ | Qualified and released |
| .0030, 0.120 mm | 970, 10×16 | $1.44\times10^{-12}$ | $-.021181$ | Adverse diagnostic root |
| .0032, 0.128 mm | 864, 10×16 | $7.36\times10^{-13}$ | $+.00030267$ | Independently qualified |
| .0032, 0.128 mm | 864, 16×24 | $6.96\times10^{-11}$ | $+.00030267$ | Independently qualified with identical coefficients |
| .0032, 0.128 mm | 970, 10×16 | $2.28\times10^{-12}$ | $-.025810$ | Adverse diagnostic root |
| .0036, 0.144 mm | 864, 10×16 | $4.52\times10^{-11}$ | $-.004074$ | Adverse diagnostic root |
| .0042, 0.168 mm | 864, 10×16 | $4.79\times10^{-11}$ | $-.009750$ | Adverse diagnostic root |
| .0042, 0.168 mm | 864, 16×24 | $8.03\times10^{-13}$ | $-.009750$ | Adverse diagnostic root |
| .0042, 0.168 mm | 970, 10×16 | $2.17\times10^{-11}$ | $-.048640$ | Adverse diagnostic root |

The 970-cell alternative started from the qualified $g=.003$ field, with a native margin of $+.003781$ after geometry transfer and weak continuity projection. Its converged original maximum at $g=.0042$ was $1.0205\times10^{-12}$, with the same adverse region and a native minimum of $-.04256723$. Maximum canonical differences between the two 970-cell solutions were $2.23\times10^{-10}$ in flow, $7.32\times10^{-11}$ in $k$ and $3.66\times10^{-12}$ in $\varepsilon$.

The richer 864-cell solution preserves the sign and location. The 970-cell mesh splits the streamwise interval at $a=1.40625$ and the radial interval below $h$ at $r=.1146829575130125$. Its converged root independently integrates below $2.26\times10^{-9}$ in every equation block, within the individual original limits. The negative minimum moves to $(1.40625,h)$ on a refined owner and $G_{xz}$ rises to 1.52180 there. The raw traction and $k$ diffusion RMS jumps fall from $1.17\times10^{-7}$ and $2.12\times10^{-9}$ to $7.52\times10^{-8}$ and $1.43\times10^{-9}$. Integrated tangential wall force changes from $-.0657630$ to $-.0658388$. The localized magnitude and location are sensitive to spatial representation. The negative sign is reproducible across these solved representations.

The coefficient digest of the admitted $g=.0032$ root is `ce61d6e9a921a661c518e8154b4b63b7a9e2dae425b5d4c4508fa6a6cb15d810`. The whole matching interval is $[.96227319,1.04498690]$. The independent limiting derivative uncertainty is below $1.0\times10^{-11}$. The next higher 864-cell root, $g=.0036$, had a negative response at the same physical support location after a converged original solve.

At $g=.003$ and $g=.0032$ the independently equation-converged 970-cell roots lose the positive sign seen on 864 cells. Both refined checks limit at $(1.40625,h)$, with $G_{xz}=1.37752$ and 1.40032, respectively. On the 864-cell roots the limiting values were 1.30170 and 1.31556. Both 970-cell roots retain strictly positive scalar polynomials, covariance and complete table coverage. The advance to 0.128 mm is admitted on 864 cells and its richer integration rule. Further use as a height predictor depends on resolving the spatial sensitivity.

The lower $g=.0024$ control started from the independently qualified 864-cell developed anchor. The 970-cell solve fell from $1.34973\times10^{-3}$ to $1.21314\times10^{-9}$, then ran out of time while preparing a new factor. The one bounded extension reused the exact accepted coefficient hash `3e67638b8b3528d376a1a05923dcd1f088e077c448091194ac83b8943fdaac4b` with identical geometry and 10×16 rule. A fresh current correction reached $1.01413\times10^{-12}$. Its native principal ratio was $-.00435269$.

At $g=.0024$ or 0.096 mm the independent order-28 maximum was $5.39\times10^{-11}$ and each equation met its limit. The angular minimum was $-.00720601$ after uncertainty, again at $(1.40625,h)$. At that witness $G_{xz}=1.31199$ and differentiation spread was $5.48\times10^{-12}$. This solved lower control shows that the representation-dependent sign precedes the new height increment on the 970-cell grid. The independently audited results establish an obstruction on the tested continuation roots. The continuum sign and existence of other admissible roots remain open.

The Linux campaign used 7,475 elapsed seconds including builds, setup, solver runs and independent audits. The lower-grid correction and audit used 275 seconds beyond the initial envelope, within the single 540-second extension. The recorded native process wall-time sum is 6,002.84 s and the independent audit process sum is 4,084.35 s. These sums include overlapping runs and rejected work. The peak solver RSS was 1,304,132 KiB. Builds used eight explicit jobs. Solver and audit processes used one rank and numerical-library teams of one, with at most five concurrent run or audit lanes.

## Reproduction and custody

Use the authoritative Linux run volume. The recorded scripts, exact sessions, cost files and all coefficient fields sit beneath `$run` below. Their digests are in the manifest. Use a new output directory for reproduction. The archived source fingerprint is `1e446e90b3486cd0128286030a6333e4210622d0c0fa6887948dad033a3e888f`. The tested current numerical Git SHA is `485d6393ffd91f4bfeb1e4b44edd148b8e64e33d`, its source fingerprint is `2729c99891f4c7526927d73dbd0c8f7abd9a6248c2be5fdbbf4850b017f89037`, and its executable SHA-256 is `94b96c153f104f6e0dd428f25b8cc876635767bcb039410213c1d5b6e2c06251`. The helper tree for the height path is `de31973d56ff54eac7b0793893d568fc21644f1f`. The fixed-$h$ reference SHA-256 is `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`.

```bash
ws=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
source "$ws/env.sh"
run="$NEKTAR_RUNS/wj-c1-original-height.dIppuznG"
old="$NEKTAR_RUNS/wj-c1-refined-steady.kS72Dxve"
frozen="$old/nested-small-height/001-height-0.003"
runner="$NEKTAR_RUNS/wj-c1-steady-recovery.LOp8mfob/run-capped.sh"
suite=solvers/IncNavierStokesSolver/Tests/RANS/h27
helper="$run/route/$suite"
table="$run/route/docs/plans/wj-h27-linux/task2r-r1/tables-exact-inner-v1.json"
new=$(mktemp -d "$NEKTAR_RUNS/wj-c1r3-repeat.XXXXXXXX")
mkdir "$new/tmp"
export TMPDIR="$new/tmp" PYTHONPYCACHEPREFIX="$new/tmp/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cmake -S "$run/source/$suite" -B "$new/build" -DCMAKE_BUILD_TYPE=Release \
  "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
cmake --build "$new/build" --parallel 8 --target IncNavierStokesSolver WJH27SteadyKrylovProbe
exe="$new/build/IncNavierStokesSolver"
sha256sum "$frozen/run.xml" "$frozen/state.diagnostic-steady.json" "$table" "$exe"
python3 "$helper/c1_steady_case.py" "$frozen/run.xml" "$new/replay" \
  --seed "$frozen/state.diagnostic-steady.json" --evaluate-only --seconds 90
bash "$runner" "$new/replay" "$exe" 140
for arm in rebuild reuse; do
  opt=()
  if [ "$arm" = reuse ]; then opt=(--maximum-inverse-age 4 --reuse-linear-iterations 12); fi
  q="$new/$arm"
  python3 "$helper/c1_steady_case.py" "$frozen/run.xml" "$q" \
    --seed "$frozen/state.diagnostic-steady.json" --seconds 1500 \
    --eta 0 --inverse ColoredSparseV1 --direction-model "${opt[@]}"
  bash "$runner" "$q" "$exe" 1560
  python3 "$helper/c1_steady_case.py" "$q/run.xml" "$q-rich" \
    --seed "$q/state.diagnostic-steady.json" --quadrature 16 24 --seconds 900 \
    --eta 0 --inverse ColoredSparseV1 "${opt[@]}"
  bash "$runner" "$q-rich" "$exe" 960
  bash "$run/audit-qualified.sh" "$helper" "$table" "$q-rich"
done
q="$new/rebuild-rich"
release=(--packet "$q/state.diagnostic-steady.json" --equations "$q/equations-q28.json"
         --response "$q/margins-q28.json" --points "$q/points-q28.json")
python3 "$helper/c1_release_original.py" "$q/run.xml" "$new/quiet" "${release[@]}" --steps 3
bash "$runner" "$new/quiet" "$exe" 240
python3 "$helper/c1_release_original.py" "$new/quiet/run.xml" "$new/cold" \
  --restart "$new/quiet/state.recovery-2.bin" --steps 1
bash "$runner" "$new/cold" "$exe" 180
cmp "$new/quiet/state.final.bin" "$new/cold/state.final.bin"
for peak in .000025 -.000025; do
  p="$new/pulse-$peak"
  python3 "$helper/c1_release_original.py" "$q/run.xml" "$p" "${release[@]}" \
    --steps 3 --velocity-radius 1.1 --velocity-peak="$peak"
  bash "$runner" "$p" "$exe" 300
done
```

The actual standalone physical-path audits are in `$run/physical-release.sh`. `$run/start-height.sh`, `next-height.sh` and the stage journals record the exact handoff and driver invocations, including caps and path digests. The archived frozen packet and the two source trees remain immutable.

For the localized adverse case, the stored `g0042-quadrature` session pins the 864-cell original root's coefficients and 16×24 integration. The following command sequence recreates the spatial mesh, transfers the original coarse root, and solves and audits it independently. Substitute the admitted `release-root` packet for the alternative preparation.

```bash
p="$run/height-track/001-height-0.0042"
python3 "$helper/c1_refined_witness.py" nested "$p/run.xml" "$new/support-mesh" \
  --streamwise 1.40625 --radial .1146829575130125
q="$new/support-original"
python3 "$helper/c1_steady_case.py" "$new/support-mesh/run.xml" "$q" \
  --seed "$p/state.diagnostic-steady.json" --height .0042 --remap \
  --quadrature 10 16 --seconds 1100 --eta 0 --inverse ColoredSparseV1 \
  --maximum-inverse-age 4 --reuse-linear-iterations 12
bash "$runner" "$q" "$exe" 1135
bash "$run/audit-qualified.sh" "$helper" "$table" "$q"
```

For a rule-only reevaluation, inherit `run.xml` and omit the height override. The attempted decimal override from `0.0031999999999999997` to `0.0032` correctly failed the fixed-geometry import. Keeping the existing geometry bits resolved the setup failure and preserved the coefficient digest.

The bounded lower-control continuation used the following fixed-height command from its actual accepted packet.

```bash
first="$run/g0024-spatial"
q="$new/lower-control-resume"
python3 "$helper/c1_steady_case.py" "$first/run.xml" "$q" \
  --seed "$first/state.diagnostic-steady.json" --quadrature 10 16 --seconds 325 \
  --eta 0 --inverse ColoredSparseV1 --maximum-inverse-age 4 --reuse-linear-iterations 12
bash "$runner" "$q" "$exe" 360
bash "$run/audit-qualified.sh" "$helper" "$table" "$q"
```

The last qualified packet is `$run/height-bracket-rich-bound/state.diagnostic-steady.json` on 864 cells. The two preparations of the adverse 0.168 mm result and the resumed adverse 0.096 mm control are retained under their diagnostic roles. All launched nonlinear solves needed for the decision have a converged endpoint.

The next discriminating test is a solved original $g=.0024$ lower anchor on an independently placed support refinement. It would keep $h$, the chart and the exterior laws fixed and move the added streamwise and radial cuts away from the 970-cell alignment. Compare both owners at the same physical support points, the new searched minimum, $G$, $k/\varepsilon$ and wall force. Resume qualified height continuation only from a common positive-response result across the targeted representations.

## Gates

| Gate | Status | Scope |
|---|---|---|
| `FROZEN_G003_CASE_REPRODUCED` | PASS | Exact coefficient SHA and original residuals on both executables. |
| `DEADLINE_AWARE_CONTINUATION_VERIFIED` | PASS | Native fixed-height resume and separate analytic uncommitted-trial check. |
| `CURRENT_DIRECTION_MODEL_VERIFIED` | PASS | Current true linear checks and measured unshifted residual predictions, including the shifted fixture. |
| `COUPLED_INVERSE_REUSE_VERIFIED` | PASS | Same frozen base, curved constraint, same-base refresh and cold entry. |
| `QUALIFIED_SOLVE_COST_COMPARISON_COMPLETE` | PASS | Both arms met the same richer and independent gates. |
| `NESTED_G003_ORIGINAL_ROOT_QUALIFIED` | PASS | Exact requested 864-cell representation. |
| `NESTED_ORIGINAL_PHYSICAL_RELEASE_VERIFIED` | PASS | Quiet, signed downstream and exact cold recovery through $t=.06$ at $g=.003$. |
| `ORIGINAL_HEIGHT_BRANCH_ADVANCED` | PASS | Independently admitted through $g=.0032$ on 864 cells. |
| `FULL_HEIGHT_ORIGINAL_PHYSICAL_PILOT_PASSED` | BLOCKED | Adverse sign on converged continuation roots. Full-height physical time is zero. |
| `C1_CLOSED` | OPEN | Full-height original root and physical pilot remain outstanding. |
| `C2_READY` | OPEN | Depends on C1. |
| `LONG_FLAT_PRODUCTION_COST_CLOSED` | OPEN | The later limit remains $1.25C_{\rm total}\le86400$ s per qualified trajectory. |
| `REMOTE_PUSH_VERIFIED` | PASS | `git ls-remote origin refs/heads/codex/wj-earsm-awf` matched the evidence head `768a65ba164a775109c051cbaf62f95745b28f8b`. |
