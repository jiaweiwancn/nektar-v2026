# Divergence-compatible hill seed test

The foot fit reached height fraction .016238 (0.64952 mm) through its first BDF2 step with both independent audits passing at $t=.04$. This run uses the verified BE predecessor on the same .0158348 source geometry. The original .04 source packet also benefited from the fit, but its finer constitutive scan remained below the acceptance floor. The 40 mm preflights failed. The earlier continuation, current equations and definitions of the C1 gates are documented in the [geometry continuation report](wj-h27-c1-geometry-continuation.md).

## Discrete seed

`Coupling::ProjectContinuity` solves a constrained least-squares problem over the native Q3 velocity coefficients. With $u_P$ the tapered Piola velocity and $G_s$ the source Cartesian gradient at the same chart label, the objective, up to a constant and a common factor, is

$$
\int_{\Omega_t} \left(|u-u_P|^2 + (.4H)^2 |\nabla_t (u-u_P)|^2
+ \ell_f^2 w_f |\nabla_t u-G_s|^2
+ \ell_c^2 w_c |\nabla_t u-QG_sQ^T|^2\right)\,d\Omega_t.
$$

The native target continuity matrix constrains the Q2 weak divergence to zero with the prescribed velocity traces. The source derivatives use the source inverse metric. $Q$ rotates the source wall frame to the target frame in the crown fit. The projected velocity receives a fresh BE/BDF2 history. The pressure multiplier of the projection is discarded. The time stages use the existing momentum, scalar and constitutive equations.

The foot weight is full within .25 chart units of the downstream foot and between $5.5H$ and $9H$ in radius. It fades to zero at .80 from the foot, $4H$, and $11H$. The crown weight is full within .2 hill half-lengths of the crown and between $.4H$ and $.8H$. It fades to zero at .6 half-lengths, $.15H$, and the specified radial reach. The reach defaults to $H$. The reflected quintic evaluates each decreasing transition directly. This corrected a small negative roundoff weight in the first collar assembly. Both fit lengths default to zero.

`WJGeometrySeedPacketSHA256` and `WJInitialSourceIdentity` bind the archived positive-time coefficient snapshot for the matched comparisons. The native loader verifies the declared model, geometry, signed layout, coefficient digest and accepted physical state metadata. The incorrect-digest and unpinned cross-revision controls rejected. Ordinary cold recovery still checks the complete numerical signature.

`c1_seed_compare.py` reconstructs signed velocity and pressure basis actions independently using the source and target physical charts. The first BE packet carries the exact prepared coefficients in its input history. The comparison verifies their hash, integrates the patch gradient errors, and computes a separate pressure-mass norm of the Q2 weak divergence. `c1_oracle.py` checks the discrete equations and budgets. `c1_margins.py` checks the principal response and covariance on the finer rule. The principal floor is $10^{-8}$ with the directional-difference allowance.

## Matched controls and strength sweep

The original Release executable reproduced the archived foot BDF2 rejection. Its first accepted checkpoint matches the preserved one byte for byte. At the BDF2 rejection, its residual was $1.24\times10^{-12}$ and its principal response was $-1.07257073\times10^{-5}$ at $(a,r)=(2.4884754,1.1168174)$. The archived full-height preflight reproduced $-5.99109225\times10^{-3}$ at $(.0115246,.0906210)$. The new executable with zero fit reproduced both unaltered gradient witnesses. The foot comparison retains the 760-cell grid, .30 Piola taper, and .02 timestep.

| Foot fit length with the archived .04 parent | Native time reached | Independent 24×24 outer principal at that time |
| ---: | ---: | ---: |
| 0 | .02 | -5.402e-5 |
| .3 | .04 | -7.101e-5 |
| 1 | .04 | -2.748e-5 |
| 3 | .04 | -9.278e-6 |
| 10 | .04 | -9.819e-6 |
| 30 | .04 | -1.233e-5 |

On the .04 parent, the independently measured foot gradient RMS at initialization is $4.53118\times10^{-5}$ with zero fit and $7.44471\times10^{-6}$ with length 3. The weak continuity norms are $3.33\times10^{-13}$ and $3.91\times10^{-13}$. The imposed velocity coefficients agree exactly. Further strengthening the fit at lengths 10 and 30 reduced the prepared gradient RMS to $3.83\times10^{-6}$ and $1.49\times10^{-6}$. The post-BDF2 margins from these cases remained negative. Halving the timestep at length 3 left the common-time outer margin negative at $-1.04228\times10^{-5}$. A finer 24×24 native replay kept the original parent's accepted BE state and rejected its BDF2 candidate with principal response $-6.11\times10^{-6}$.

## Earlier parent and recovery

The parent replay advances .014222 to .0158348 for one BE step with the archived standard-main solver. The coefficient SHA-256 matches the .04 parent's prior history exactly. The replay has independent outer margin $1.32300\times10^{-4}$ at .02. Its independent equation residuals are `[4.27e-12, 2.40e-15, 5.89e-11, 2.02e-11]` in momentum, continuity, k and outer epsilon order.

| Seed from the verified .02 parent at height .016238 | Native rule | Independent rule | Time | Outer principal |
| --- | --- | --- | ---: | ---: |
| Zero fit | 12×20 | 24×24 | .04 | -5.546e-5 |
| Foot length 3 | 12×20 | 24×24 | .02 | 7.742e-5 |
| Foot length 3 | 12×20 | 24×24 | .04 | 3.441e-5 |
| Foot length 3, final numerical revision | 24×24 | 32×32 | .04 | 2.298e-5 |

From the .02 parent, the foot gradient RMS drops from $4.53106\times10^{-5}$ to $7.44447\times10^{-6}$, and the fitted seed has independently integrated weak continuity $4.01\times10^{-13}$. The positive .04 state has an independently sampled upslope face-limit principal margin $7.66\times10^{-6}$ at the downstream foot. The downstream-side limit is $2.65\times10^{-3}$. The one-sided scan includes radial element ends and 64 interior Gauss points per sampled edge element.

On the final build with native 24×24 and independent 32×32, the scaled equation residuals at .04 are `[7.18e-14, 6.20e-16, 5.83e-12, 3.22e-13]` against `[1e-9, 1e-9, 1e-8, 1e-8]`. Whole-interface $K/K_{\rm in}$ is $[.96101178,1.02185123]$. The one-rank case took 102.36 s including startup, both steps and recovery. Peak RSS was 1,227,808 KiB. The original 12×20 cold continuation from .04 has independently negative outer response at .06 and .08, at $-6.62\times10^{-6}$ and $-4.54\times10^{-5}$. The same-signature 24×24 cold run rejects the .06 candidate at residual $1.03\times10^{-12}$. Its principal response is $-3.46787\times10^{-6}$ at $(2.4969920,1.1078989)$, and independent reconstruction gives $-3.46787\times10^{-6}$. The held .04 checkpoint matches the source checkpoint byte for byte.

## Full-height preflight

| Foot length | Crown length | Crown radial reach | Principal at first rejected point | Chart point $(a,r)$ |
| ---: | ---: | ---: | ---: | --- |
| 0 | 0 | .15 | -5.991e-3 | (.01152, .09062) |
| 1 | 0 | .15 | -5.873e-3 | (.01152, .09062) |
| 1 | .12 | .15 | -1.928e-3 | (.05993, .09414) |
| 1 | .4 | .15 | -1.028e-3 | (1.26152, .11540) |
| 1 | .8 | .15 | -1.984e-3 | (-.01152, .13962) |
| 1 | 1.2 | .15 | -2.081e-1 | (-.01152, .14480) |
| 1 | .4 | .30 | -1.348e-2 | (1.23848, .11248) |
| 1 | .8 | .30 | -1.057e-3 | (1.26152, .11248) |
| 1 | .4 | .45 | -1.560e-2 | (1.23848, .11248) |
| 1 | .8 | .45 | -1.308e-3 | (1.23848, .11248) |

The stronger fits with reach $H$ moved the rejection into the last part of the radial fade. Extending that fade through matching moved the first rejection to the lee side of the crown patch. The fixed-state independent directional checks recover the sign in every row. Each value is the first rejected sample on that case's native rule. Full-height physical advancement remains blocked before BE.

## Linux reproduction and artifacts

The source and runs use `/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace` on the requested Linux host. The verified LAN SSH route `jay.local` uses the pinned key for `100.88.70.60`. The inspected GitHub WJ branch was `57845c53a29ba32cca49c31a89e0f4ec4da4f300` and the local starting branch was `f24ee6f9e1fc612206b50c9fdb57641936c4e986`. The intermediate reflected build is from `f4e179b7d5b57568941d4c473f66c985295a2bb2`. The final numerical revision is `8036cabae4a7ee267ca9fbc09ea9a52bb6621ea0` and its executable SHA-256 is `6746f9f205296174515fe9984f14b6fd4a884721bf1245b08e008fbdc4beb49e`. The run directory retains the fixtures, sessions, field packets, rejection logs, recovery histories and audits. An intermediate finer-rule solver run completed before its shell wrapper encountered a concurrent script edit. Its two independent audits were recovered from the completed packet in a separate invocation. [Evidence manifest](wj-h27-c1-seed-gradient-evidence.json)

```bash
ws=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
source "$ws/env.sh"
a="$NEKTAR_RUNS/wj-c1-geometry.uZz3Uh7t"
study="$NEKTAR_RUNS/wj-c1-seed-gradient.BoAhhRpl"
RUN=$(mktemp -d "$NEKTAR_RUNS/wj-c1-seed-gradient-replay.XXXXXXXX")
mkdir "$RUN/tmp"
t="$NEKTAR_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/h27"
table="$NEKTAR_SOURCE/docs/plans/wj-h27-linux/task2r-r1/tables-exact-inner-v1.json"
export TMPDIR="$RUN/tmp" PYTHONPATH="$t" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cmake -S "$t" -B "$RUN/build" -DCMAKE_BUILD_TYPE=Release \
  "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
cmake --build "$RUN/build" --parallel 8 --target IncNavierStokesSolver
solver="$RUN/build/IncNavierStokesSolver"
run_case() {
  (cd "$1" && /usr/bin/time -v mpirun --stdin none --bind-to none -np 1 \
    "$solver" --use-opt-file run.xml.opt run.xml >run.log 2>&1)
}
seed="$study/replay-source-be/state.result.json"
identity=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["state"]["source_identity"])' "$seed")
packet=$(sha256sum "$seed" | cut -d' ' -f1)
python3 "$a/make-local-foot-radial.py" hill "$RUN/foot" --height .016238 \
  --from-height .0158348 --initial-fields "$seed" --initial-source "$identity" \
  --seed-packet-sha256 "$packet" --seed-taper .30 --seed-foot-length 3 \
  --quadrature 24 24 --dt .02 --steps 2 --science 1 --step-seconds 180
run_case "$RUN/foot"
python3 "$t/c1_oracle.py" --table "$table" --native "$RUN/foot/state.result.json" \
  --order 32 --seconds 300 --out "$RUN/foot/independent.json"
python3 "$t/c1_margins.py" --table "$table" --native "$RUN/foot/state.result.json" \
  --order 32 --seconds 300 --out "$RUN/foot/independent-margins.json"
python3 "$t/c1_seed_compare.py" --source "$seed" --target "$RUN/foot/state.science-1.json" \
  --prepared --order 24 --out "$RUN/foot/prepared-gradients.json"
```

For the unrefined-rule pair, use `--quadrature 12 20` and audit at order 24, with `--seed-foot-length 0` or `3`. For the full-height preflight, use `nested-radial-track/04-height-0.0142220-q12x20/state.result.json` from `$a` as the source, recompute both bindings, and use `--from-height .014222 --height 1 --steps 0 --quadrature 12 20 --seed-foot-length 1 --seed-collar-length .8 --seed-collar-reach .30`. The expected process status is 1. `c1_assess.py --guard` reproduces the failed point from the run log.

The independently qualified height bracket reaches .016238 through .04. The later loss during ordinary cold recovery bounds this result in time. The full-height physical pilot is blocked and C2 remains closed.
