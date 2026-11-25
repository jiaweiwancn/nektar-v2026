# C1 refined steady recovery at fixed shallow hill geometry

The campaign solves the original finite-viscosity equations at `eta=0`. On the 658-cell representation at hill fraction $g=.0024$, a full coupled Newton correction converged the saved positive-response field. The inlet, reference table, normal and outer charts, matching trace, sampled AWF loads, and `ZeroTransverseV1` outlet are those of C1R. Constitutive homotopy was NOT_NEEDED for this campaign. [Governing implementation and reference](wj-h27-c1-normal-hill.md) · [C1R custody and original equations](wj-h27-c1-steady-recovery.md)

## Fixed comparison

The archive at `wj-c1-steady-recovery.LOp8mfob` is read-only for this comparison. Native replay on the clean Linux correction build reproduced the archived coefficient digests and refreshed residuals for all three controls. Independent quadrature of order 24 reproduced their residuals and signed response.

| Archived preparation | $g$ | Native order | Refreshed $F_0$ maximum | Independent momentum / continuity / $k$ / outer epsilon | Independent inner / outer principal divided by $\mu$ |
|---|---:|---:|---:|---|---|
| Developed original anchor | .002 | 10×16 | 2.15568e-11 | 2.60637e-11 / 3.59e-16 / 9.09247e-11 / 1.21937e-11 | +.00235848 / +.00178326 |
| Negative coarse root | .0024 | 16×24 | 3.59855e-11 | 7.74290e-12 / 2.57e-16 / 3.63957e-11 / 8.13712e-12 | −.00427141 / +.00096505 |
| Positive refined field A | .0024 | 10×16 | 7.27866e-7 | 8.09506e-8 / 4.50e-13 / 7.27855e-7 / 4.72693e-7 | +.00522588 / +.00172401 |

The fixed refined grid has seven chart columns and 94 radial rows with 175 outer cells. The signed Modified_A coefficient spaces and essential constraints are as follows. Outlet traction fixes pressure, so there is no gauge row on this physical mesh.

| Block | Canonical coefficients | Free unknowns | Essential coefficients |
|---|---:|---:|---:|
| Cartesian $u$ | 6,226 | 5,922 | 304 |
| Cartesian $w$ | 6,226 | 5,620 | 606 |
| Modified pressure $\Pi$ | 2,835 | 2,835 | 0 |
| $k$ | 6,226 | 5,943 | 283 |
| Outer epsilon | 1,672 | 1,575 | 97 |
| Total | 23,185 | 21,895 | 1,290 |

Preparation B is the coarse root projected to A's signed inventory. The continuity projection reduced the scaled mass defect from 3.10272e-4 to 2.594e-13. The whole trace $K/K_{\rm in}$ remained [0.95525011, 1.03496777], with the dependent epsilon lift rebuilt from it. The initial refined $F_0$ maximum became 1.64931e-3. The native sampled principal minimum became −.00497789. The changed polynomial tests can detect defects that were outside the coarse test space.

## Coupled correction

| Code | Responsibility |
|---|---|
| `WJNativeCoupling::SearchSteady` and `MakeSteadySparseInverse` | Pack the original coupled equations and pressure constraints, supply signed dependency maps and dependent lift derivatives, freeze correction terms, and refresh acceptance metrics. |
| `WJNativeSteadyCorrection` and `WJNativeSteadyKrylov` | Assemble and factor the scaled colored inverse, apply original or shifted Newton actions, record checked linear residuals by equation, and globalize with the original norm. |
| `WJNativeInitialFields::ValidateOriginalSteadyRelease` and `Coupling::ReleaseOriginalSteady` | Validate the signed state and audit receipts, check the original steady state, and build a fresh physical history. |
| `c1_steady_case.py`, `c1_refined_witness.py`, `c1_release_original.py` and `c1_steady_continue.py` | Pin exact preparations, rebuild a nested chart, independently inspect both sides, and drive original fixed-height solves and tests. |

`ColoredSparseV1` prepares a right inverse of the full scaled Newton correction. The row-disjoint coloring uses actual signed velocity and scalar cell maps. Inner wall data and moving cutoff panels share matching $K$ columns. The outer matching owner supplies the epsilon lift chain. Momentum in the forced layer also sees the matching velocity and sampled $k$ columns. Exact columns carry $G_h-C_{\Pi}$, divergence and the gauge if present. Centered differences of the original coupled residual refresh production, current transport, scalar reaction and diffusivity, dependent traces, and WJ stress together. The shifted preparation, when selected, adds the frozen correction mass and positive stiffness only to the linear system. The original unshifted equations supply globalization and the final acceptance test.

Both inverses used the same 21,895 free variables, original equations, row and variable scales, finite-difference setting $4\times10^{-6}$, and initial forcing factor 0.3. Its physical dual target was 2.18360e-7. The comparison used the same base coefficient digest `932de8ae565a8b7cebd09bd33ebc220061c23a2b1b8f7e7f6bdc3d854cc6582b`.

| Common-base correction | Krylov iterations | Best true momentum / continuity / $k$ / epsilon | Accepted nonlinear $F_0$ | Complete process / peak RSS |
|---|---:|---|---:|---:|
| `BlockV1` | 160 | 7.81954e-8 / 4.32019e-13 / 7.24829e-7 / 4.58454e-7 | 7.27866e-7, no accepted correction | 154.95 s / 647,676 KiB |
| `ColoredSparseV1` | 1 | 1.88610e-10 / 3.32382e-14 / 3.39097e-10 / 1.41811e-12 | 7.25606e-9, full step | 165.72 s / 735,988 KiB |

The sparse inverse used 212 colors, 424 residual calls, 2,077,368 structural entries, and 28,505,757 factor entries. Assembly took 131.54 s and sparse LU took 2.57 s. The separate direction audit took 14.20 s and 48 calls. Krylov residual callbacks took 1.37 s, and applying the new inverse took .027 s. The old inverse spent 132.39 s in residual callbacks and 6.75 s in factor preparation and application. Both base corrections had zero infeasible derivative trials. The old solve also failed before globalization. At this base the linear coupling was the limiting operation.

The frozen audit used four displacement lengths 1.038e-4, 2.594e-5, 6.485e-6 and 1.621e-6 in each scaled direction. The velocity and mixed Taylor remainder fell by factors near 16 when the displacement decreased by four. The scaled centered error reached 1.69e-5 for velocity and 2.88e-6 for the mixed direction before approaching the assembly difference floor. Bulk $k$, matching $k$ and outer epsilon reached 1.06e-11, 3.34e-11 and 1.81e-11. Pressure was linear with a Taylor remainder at 1e-15. The largest scaled action outside the declared support was 1.02e-11. Actual true linear checks continued to use uncolored full-residual differences.

The fresh shifted probe and the probe after changing the base returned the independently predicted action 3. The curved additive-action fixture used $\Delta\tau=1000$ and $\beta=2\times10^{-5}$. The maximum discrepancy in the scaled shift action was 2.29e-10 with dependent trace derivative discrepancy 3.04e-11. Its freshly evaluated original residual was 2.79e-11. The out-of-table rejection control forced a matching vertex to $1.12K_{\rm in}$, retained every coefficient byte under SHA-256 `871845170dc15504ba5b9f27a67b6fa3c62e83f831f942c6dcbce61c87c24385`, exited 1, and recorded a null cutoff range and dependent lift. The packet remains a rejected initialization with no physical history.

## Original roots and refinement

Preparation A accepted two full corrections in total. B accepted six corrections before exhausting the initial 900-second native budget during preparation of the seventh inverse. Its last accepted original maximum fell from 1.64931e-3 to 7.63974e-7, with positive margin and room inside the table. A digest-bound extension with 650 native seconds reached the gate in one correction. The two endpoints agree within [2.70e-10, 7.86e-11, 1.20e-10, 7.71e-11, 5.62e-12] in the maximum absolute canonical coefficients of [u, w, Pi, k, epsilon]. They use the same signed Modified_A inventory and numerical signature. This supports the same root from the two preparations at the measured precision. The negative coarse result is representation-sensitive at this height.

| Root at $g=.0024$ | Native order | Fresh native maximum, limit 1e-10 | Independent order-28 momentum / continuity / $k$ / epsilon | Projected trace defect, limit 1e-10 | Smallest sampled principal divided by $\mu$ |
|---|---:|---:|---|---:|---:|
| A converged | 10×16 | 1.16393e-12 | 7.65557e-10 / 4.70e-16 / 1.64790e-10 / 6.57938e-12 | 4.23e-16 | +.00163373 |
| B converged | 10×16 | 4.93348e-11 | 7.70403e-10 / 4.79e-16 / 1.66978e-10 / 4.95852e-11 | 2.80e-16 | +.00163373 |
| A re-solved | 16×24 | 8.05524e-13 | 8.08388e-13 / 4.42e-16 / 5.82576e-12 / 3.20383e-14 | 3.21e-16 | +.00163372 |
| A on targeted nested mesh | 10×16 | 8.39232e-13 | 4.00488e-13 / 6.81e-16 / 1.03838e-11 / 1.15252e-13 | 4.96e-16 | +.00385157 |
| Independent limits | | | 1e-9 / 1e-9 / 1e-8 / 1e-8 | 1e-10 | >1e-8 after uncertainty |

The rich rule changed the native residual of the old coefficients to 7.65531e-10. One further correction was required. The richer root's whole $K/K_{\rm in}$ interval was [0.96573794, 1.09056075] with margins 0.06573794 and 0.00943925. Independent order 28 scanned 450,268 inner and 137,200 outer volume points, plus both owners of cell and matching faces and cutoff samples. The maximum estimated differentiated form uncertainty across the regions was 5.97e-11. The 25 and 49 point local grids and the independently minimized angular checks gave a least normalized response of +.00163372347 at the downstream hill-side face, $(a,r)=(2.5,.79431570)$.

The rich root's common trace and matching-measure discrepancies were each 1.11e-16, and its largest weak matching row was 6.28e-16. The inner weak k balance was 4.00e-16 and the discrete power error relative to free-test work was −1.31e-15. The raw one-sided traction and k-diffusion differences were 1.70483e-7 and 1.82087e-9 RMS. The pointwise approximation to the nonlinear epsilon trace was 1.85618e-6 with projection defect 3.21e-16. Tangential wall force on the fluid was −.06614380667 in code units. The manifest binds the corresponding velocity, gradient, divergence, turbulence time scale and root diagnostics for upstream, crown, matching and both downstream foot limits.

The nested mesh adds streamwise cuts at $a=-8.25,1.875$ and radius cuts at $r=.1102,1.10$. It uses 864 cells with the same geometry, inlet, traction and outer epsilon ownership. Refinement initially exposed a weak continuity defect of 8.48e-4. Projection reduced it to 6.24e-13, and the original residual of that projected seed was 4.92105e-3 with principal ratio −.00093617. Five coupled corrections restored an admissible original steady solution. The nested whole $K/K_{\rm in}$ range was [0.97387561, 1.02869095], the cutoff range was [0.02422909, 0.02455521], and the independent local minimum was +.00385157275 at $(2.5,9.91641)$. The nested inner and outer volume minima were +.01061814 and +.00386257.

| Fixed observation | 658-cell richer root | 864-cell nested root |
|---|---:|---:|
| Upstream $u$ at $(-8.29,.1111)$ | .40651004 | .39670437 |
| Upstream $G_{xz}$, divergence and $k/\varepsilon$ there | .902153 / +.00936135 / 2.90014 | 1.18491 / −.00408783 / 2.82452 |
| Downstream-foot $u$ at $r=1.107899$ | .70359269 | .70446771 |
| Downstream-foot $k/\varepsilon$ there | 27.98816 | 28.01314 |
| Downstream-foot hill and flat one-sided principal ratios there | +.00225252 / +.00296368 | +.00729590 / +.00701882 |
| Tangential wall force on fluid | −.06614381 | −.06578459 |
| Raw matching traction and k-diffusion defects, RMS | 1.70483e-7 / 1.82087e-9 | 2.00315e-8 / 6.40769e-10 |
| Inner weak k balance and discrete power defect | +4.00e-16 / −1.31e-15 | −1.68e-17 / −2.32e-16 |

These comparisons resolve the positive sign on the tested shallow steady branch with finite spatial sampling and a targeted refinement. The upstream collar still has a measurable representation dependence in longitudinal velocity, shear and divergence. The previously rejected positive-time foot and full-height crown histories remain distinct controls from C1R.

On the nested mesh the 16×24 native rule evaluated the root at 9.02204e-12, already inside the original tolerance, with native principal ratio +.00386940. The independently integrated 28×28 residual and face checks above use the same coefficients. The established flat shifted control gave the same declared additive action within 2.63e-10 and converged its original residual to 1.94e-11 in four full corrections. Both preparation fixtures kept pressure constraints algebraic.

## Release of the shallow original root

The fresh release option verifies the zero-regularization packet and numerical signature, and the SHA-256 digests of independent equation, response and local angular receipts. It evaluates the original steady residual again, reconstructs the dependent trace and runs the usual admissibility checks before building the BE startup history. The regular coefficient import and cold restart remain separate entry points. The release build re-evaluated the richer root with identical coefficients and new source-bound receipts. The wrong-receipt-digest control failed before constructing a history.

Physical advancement uses ordinary `CoupledJFNKV1` with new-time advection from BE startup, followed by BDF2. The original steady preparations use `ColoredSparseV1` in diagnostic mode. The archived T=10 Sparse-LU Newton trajectory used extrapolated advection through step 16 and then an implicit suffix. [Archived T=10 comparison](wj-h27-bf-t10-status.md)

The undisturbed physical run advanced to $t=.06$. Restarting from $t=.04$ with the same numerical signature reproduced the final complete history byte for byte. Its SHA-256 is `1b8edf71c057d581f13d8840a2696a786f2e46037aaad59cef96dfda19ddf1b9`. The signed perturbations use a single Q3 interior $k$ bubble near the downstream foot, with peaks of $\pm.002K_{\rm in}$ or about $\pm1.76\times10^{-5}$. The velocity, discrete divergence and matching and essential traces at release are shared by the two starts. All three histories began at physical time zero with the same original equations.

Both signed branches reached $t=.06$. They took three Newton corrections per timestep and the maximum accepted original timestep residuals were 1.73e-12 and 1.76e-12. The final native principal margins were +.00183309 and +.00183298. Both independent order-28 equation audits passed, with face margins no lower than +.00163351. The ordinary restart control rejected the diagnostic packet at the coupled checkpoint checksum, while the digest-bound release rebuilt its physical history. This is a short disturbed original trajectory at .096 mm hill height. The full-height horizon and common-time timestep comparison remain open.

## Original height continuation

The fixed comparison justified a bounded original geometry continuation. A rebuild of the chart and constrained velocity transfer preceded every fixed-height solve. A nonzero correction pseudo-time was tested at the same nested $g=.008$ base after unshifted corrections made little progress. All these attempts retain literal $\eta=0$ and use $F_0$ for line search and acceptance.

On the 658-cell mesh, the attempt at $g=.0036$ reduced $F_0$ from 4.82983e-5 to 3.92287e-5. Larger trials exceeded the whole matching-trace table near $a\in[-10,-8.125]$. The last accepted $K/K_{\rm in}$ maximum was 1.09904360 and the native principal margin was +.00109634. The halved-height attempt at $g=.003$ reached 1.81142e-5, with maximum 1.09501687. Both stopped at their bounded preparation timers.

On the 864-cell mesh, the unshifted $g=.008$ attempt reduced $F_0$ from 2.38466e-4 to 2.15394e-4. The last accepted $K/K_{\rm in}$ range [0.95841790, 1.04884434] had room at both edges. The independent signed response was −.00161334 at the one-sided chart point $(1.875,h)$ in cell 582. The same point on the converged nested anchor gave +.01050773. Its shear $G_{xz}$ rose from 1.25885 to 1.32735, and $k/\varepsilon$ fell from 2.89459 to 2.86932. Here $\chi=1$. Substituting just the candidate gradient into the anchor state gives −.00390233, and substituting just the candidate scalars gives +.01287571. The three-spacing variation at the adverse point is below 7.4e-12. This locates a shear-driven response loss in an unconverged iterate near the sampled-load support boundary. The earlier coarse upstream collar and positive-time outer foot failures occupy other locations.

At the exact same $g=.008$ coefficient base and equation identity, initial pseudo-time 160 adds the differential mass with initial factor .00625. Its correction reduced $F_0$ to 2.25979e-4. The independent face minimum worsened to −.00444389. The unshifted and shifted trials both encountered merit and scalar or table feasibility rejections, despite reserve in their last accepted polynomial traces. No constitutive regularization entered either residual. Above the qualified root, these outcomes show where globalization and response loss occur for the attempted paths. They leave the sign and existence of a converged taller root unresolved.

The smaller nested step to $g=.003$ used a ninefold smaller height increment than the first nested attempt. It produced the closest higher accepted candidate at $F_0=1.50855{\times}10^{-5}$, governed by the outer epsilon residual. The signed native margin was +.00390322 and $K/K_{\rm in}$ lay in [0.96977979, 1.03675650]. Its three true linear residual maxima were 2.80e-8, 2.56e-8 and 2.65e-8. Merit rejected five larger steps, with accepted fractions .25, .25 and .5. The fourth correction had no measured linear residual at the deadline. This same-base evidence gives a concrete next problem for original-residual globalization. The `nested-small-height/001-height-0.003` session and digest-bound packet retain its exact input and last accepted coefficients.

The independent order-24 evaluation of that candidate gives momentum 3.10116e-6, continuity 1.799e-13, $k$ 7.23840e-6 and outer epsilon 1.50855e-5. Its independent smallest face margin is +.00386230. The next discriminating problem is a fixed-$g=.003$ globalization test from this packet. Measure the nonlinear outer-epsilon remainder and the actual-to-predicted decrease for the coupled correction, then adjust the trust in that direction against the original $F_0$. Its three measured linear solves already meet their targets. The cutoff table has reserve on the nested candidate and the represented scalars are positive. A converged original result there would warrant the next geometry step.

The intermediate nested step at $g=.0052$ reduced the refreshed original maximum from 1.19203e-4 to 9.91481e-5 under the unshifted correction. Its final whole-trace interval was [0.96353328, 1.04413192], with principal ratios +.00147455 native and +.00052591 on the independent scan. The outer epsilon residual remained dominant and scalar polynomial and merit rejections limited the accepted updates. It is a diagnostic candidate at zero physical time.

Elapsed Linux campaign time was 5,529.56 s, including a 129.56 s tail beyond the initial 5,400 s budget for the final bounded audit. Recorded process wall-time sums were 7,867.70 s for native runs and 2,406.30 s for independent audits, with overlapping cases. Peak per-process solver RSS was 1,123,348 KiB. The A remainder solve, B initial plus extended solve, richer A solve and nested solve took 164.23 s, 1,078.42 s, 381.93 s and 943.21 s respectively, excluding setup and audits already charged to their own measurements. The signed shallow physical runs took 107.06 s and 95.21 s end to end. These campaign measurements have separate scopes from the later production cost estimate.

## Reproduction and custody

All large files are under the two Linux run roots below. The archive packet, coefficient, source, executable, table, geometry and session identities are pinned in the [machine-readable manifest](wj-h27-c1-refined-steady-closure-evidence.json). The Linux builds used GCC 11.4.0 and CMake 3.22.1 with one rank per solver and one thread per numerical library. The Linux limit is 72 workers and the campaign used at most eight compiler jobs and three concurrent solver or audit processes. Existing suspended jobs were left in place. All equivalent reruns should use a new output directory.

| Purpose | Clean numerical git revision | Native source fingerprint | IncNavierStokesSolver SHA-256 |
|---|---|---|---|
| Archived preparation A and first correction replay | `e9b1020183c3c631866313e04c2d382e019e20a6` | `6a8ff82511eb8efaac4998027aabcbe7e2110771846bfe400720aca1aeeca091` | `c5e3a8f65e2ae08d7598b17b200d920cba3a393e560dbba22899841fe07ca5a6` |
| Fixed-geometry correction and roots | `e50ec928586852214d921804d9580e2b4ea2117f` | `32c84d55c046ae835eff79e2da90b7ff96dea24e01ad2f3ddcbc6bd530fa0b30` | `714204ac8bb183d3cc47218fa3f2506878b09859e505ae0b089c015daed31ac8` |
| Receipt, shallow physical runs and height checks | `29c5a6a8d7ca9b7cf55bec25ce9c5821a2df2893` | `1e446e90b3486cd0128286030a6333e4210622d0c0fa6887948dad033a3e888f` | `4b7afeb7d324ab71a48d8a5d36b23a66f1a539bdc3ac938e1c65473c091a70f1` |

The independent helper revisions are `83bf1f2ee5832bfcbb98c6589fb7d21f4a2262ed` for the nested mesh and one-sided scans and `169362d798b27a192aac864bb8b4d003bb57ee96` for the guarded original height driver. The final report revision adds evidence to these tested sources. The table SHA-256 is `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`.

```bash
ws=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
source "$ws/env.sh"
old="$NEKTAR_RUNS/wj-c1-steady-recovery.LOp8mfob"
campaign="$NEKTAR_RUNS/wj-c1-refined-steady.kS72Dxve"
out="$(mktemp -d "$NEKTAR_RUNS/wj-refined-reproduce.XXXXXXXX")"
mkdir "$out/tmp"
export TMPDIR="$out/tmp" PYTHONPYCACHEPREFIX="$out/tmp/pycache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
path=solvers/IncNavierStokesSolver/Tests/RANS/h27
for tag in correction release; do
  if [ "$tag" = correction ]; then
    revision=e50ec928586852214d921804d9580e2b4ea2117f
  else
    revision=29c5a6a8d7ca9b7cf55bec25ce9c5821a2df2893
  fi
  git -C "$NEKTAR_SOURCE" worktree add --detach "$out/$tag-src" "$revision"
  cmake -S "$out/$tag-src/$path" -B "$out/$tag-build" -DCMAKE_BUILD_TYPE=Release \
    "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
  /usr/bin/time -v -o "$out/$tag-build.time" timeout --kill-after=15s 420s \
    cmake --build "$out/$tag-build" --parallel 8 --target IncNavierStokesSolver WJH27SteadyKrylovProbe
  "$out/$tag-build/WJH27SteadyKrylovProbe" "$out/$tag-probe.json"
done
tests="$out/release-src/$path"
solver="$out/correction-build/IncNavierStokesSolver"
table="$out/release-src/docs/plans/wj-h27-linux/task2r-r1/tables-exact-inner-v1.json"
a="$old/original-0024-refined-mass-reset"
b="$old/original-0024-richer"
prepare() { python3 "$tests/c1_steady_case.py" "$@"; }
run() { bash "$old/run-capped.sh" "$1" "$solver" "$2"; }
audit() {
  local case="$1" order="$2" packet="$1/state.diagnostic-steady.json"
  python3 "$tests/c1_oracle.py" --table "$table" --native "$packet" --diagnostic \
    --order "$order" --seconds 210 --out "$case/equations-q$order.json"
  python3 "$tests/c1_margins.py" --table "$table" --native "$packet" --diagnostic \
    --order "$order" --seconds 210 --one-sided-faces --out "$case/margins-q$order.json"
  python3 "$tests/c1_refined_witness.py" inspect --table "$table" --native "$packet" \
    --order "$order" --seconds 140 --out "$case/points-q$order.json" \
    --point -8.29 .1111 --point 2.5 1.107899 --point 0 .0917969 --point 0 .15 \
    --minima inner "$case/margins-q$order.json" --minima outer "$case/margins-q$order.json" \
    --minima faces_inner "$case/margins-q$order.json" --minima faces_outer "$case/margins-q$order.json"
}
for inv in BlockV1 ColoredSparseV1; do
  args=()
  if [ "$inv" = ColoredSparseV1 ]; then args=(--correction-audit); fi
  prepare "$a/run.xml" "$out/$inv" --seed "$a/state.diagnostic-steady.json" \
    --seconds 300 --newton 1 --inverse "$inv" "${args[@]}"
  run "$out/$inv" 360 || test -s "$out/$inv/state.steady.json"
done
prepare "$a/run.xml" "$out/A" --seed "$out/ColoredSparseV1/state.diagnostic-steady.json" \
  --seconds 900 --inverse ColoredSparseV1
run "$out/A" 980
prepare "$a/run.xml" "$out/B" --seed "$b/state.diagnostic-steady.json" --remap \
  --seconds 900 --inverse ColoredSparseV1
run "$out/B" 980 || test -s "$out/B/state.steady.json"
# Extend from B's accepted digest if the improving case exhausted its timer.
prepare "$out/B/run.xml" "$out/B-finish" --seed "$out/B/state.diagnostic-steady.json" \
  --seconds 650 --inverse ColoredSparseV1
run "$out/B-finish" 720
audit "$out/A" 28
audit "$out/B-finish" 28
prepare "$out/A/run.xml" "$out/richer" --seed "$out/A/state.diagnostic-steady.json" \
  --quadrature 16 24 --seconds 1060 --inverse ColoredSparseV1
run "$out/richer" 1120
audit "$out/richer" 28
python3 "$tests/c1_refined_witness.py" nested "$out/richer/run.xml" "$out/nested-template" \
  --streamwise -8.25 1.875 --radial .1102 1.10
prepare "$out/nested-template/run.xml" "$out/nested" --seed "$out/richer/state.diagnostic-steady.json" \
  --remap --quadrature 10 16 --seconds 1200 --inverse ColoredSparseV1
run "$out/nested" 1260
audit "$out/nested" 28
```

The launch scripts in the manifest retain the exact command lines for the replay controls, receipt-bound release, signed trajectories, and geometry stages, along with the environment location and caps. The corresponding session XML includes the actual cutoff table path and boundary laws. A finished process's `identities.sha256` is taken immediately before execution. The digest-rejection case changes the receipt hash after generating its initial input manifest, and the executed session digest records that change.

```bash
solver="$out/release-build/IncNavierStokesSolver"
prepare "$out/richer/run.xml" "$out/qualify" --seed "$out/richer/state.diagnostic-steady.json" \
  --seconds 120 --evaluate-only
run "$out/qualify" 160
audit "$out/qualify" 28
q="$out/qualify"
receipt=(--packet "$q/state.diagnostic-steady.json" --equations "$q/equations-q28.json" \
  --response "$q/margins-q28.json" --points "$q/points-q28.json")
python3 "$tests/c1_release_original.py" "$q/run.xml" "$out/quiet" "${receipt[@]}" --steps 3
run "$out/quiet" 180
python3 "$tests/c1_release_original.py" "$out/quiet/run.xml" "$out/cold" \
  --restart "$out/quiet/state.recovery-2.bin" --steps 1
run "$out/cold" 160
cmp "$out/quiet/state.final.bin" "$out/cold/state.final.bin"
for sign in + -; do
  python3 "$tests/c1_release_original.py" "$q/run.xml" "$out/pulse$sign" "${receipt[@]}" \
    --steps 3 --pulse-cell 586 --k-peak "${sign}0.002"
  run "$out/pulse$sign" 220
done
```

## Gate decision

| Gate | Result |
|---|---|
| `FROZEN_CASE_REPRODUCED` | PASS |
| `COUPLED_CORRECTION_VERIFIED` | PASS |
| `REFINED_ORIGINAL_ROOT_CONVERGED` | PASS |
| `SAME_MESH_PREPARATION_COMPARISON_COMPLETE` | PASS |
| `TARGETED_REFINEMENT_CONCLUSION_SUPPORTED` | PASS |
| `REFINED_ORIGINAL_ROOT_ADMISSIBLE` | PASS at $g=.0024$ |
| `DIAGNOSTIC_TO_PHYSICAL_RELEASE_VERIFIED` | PASS at $g=.0024$ |
| `FULL_HEIGHT_ORIGINAL_PHYSICAL_PILOT_PASSED` | OPEN |
| `C1_CLOSED` | OPEN |
| `C2_READY` | OPEN |
| `REMOTE_PUSH_VERIFIED` | Verify the published commit against `git ls-remote origin refs/heads/codex/wj-earsm-awf` |

The stronger coupled solve and resolution witnesses qualify a shallow original branch, and the higher continuation attempts delimit its present evidence. Full-height original admissibility, the disturbed $t\ge1.48$ pilot and the common-time timestep comparison remain C1 work. Production, long-flat qualification and experimental validation retain their own gates. The later production limit is $1.25C_{\rm total}\le86400$ wall seconds per qualified production trajectory.
