# B3 direction accuracy and original nonlinear curvature

**B3 remains pending.** The original maximum fell from \(6.44120254493\times10^{-6}\) to \(6.43842147203\times10^{-6}\), a 0.04318% decrease over seven updates. The target is \(10^{-10}\). The fixed-rule tangent improved accuracy on the tested vectors. The coupled trial remedy missed its budgets, and the fine hybrid action did not restore useful nonlinear contraction.

The case retains the signed 5,238-cell mesh at \(g=.0024\) (0.096 mm), \(H=.15\), \(d_0=.0875\) and reference \(h=.11541061001735275\). The original WJ law, molecular viscosity, AWF, scalar transport and production, dependent \(\epsilon_H(K)\), sampled loads, \(G_h-C_\Pi\), ZeroTransverseV1, Q3/Q2 spaces and 12×12 GLL expansion remain in \(F_0\). Correction mass, correction diffusion and residual regularization were zero. Inverse-only mass was \(.009539982703780335\). All iterations advanced zero physical time.

## Exact starting field and retained diagnostics

Linux verified the four archived input and equation-receipt hashes from R8. The primary packet `B3-seeded-trials/state.accepted-000001.json`, packet SHA-256 `30034e13516180631745c573427108001ad9f7ec8196550cde763fb03238e41e`, is under `wj-c1-b3-unshifted.DNPSiHDh`. The coefficient digest `710989059963482483be5f097f1c1c73641c0b26e89a6b1423471242c41b337a`, signed free inventory, essential lifts, geometric and equation settings and all six original residual blocks matched exactly. Each source-bound equation identity verified against its own signature. The native momentum, continuity, \(k\), and outer \(\epsilon\) norms were \(3.28902979639\times10^{-7}\), \(2.37785886136\times10^{-12}\), \(6.44120254493\times10^{-6}\), and \(2.35457631093\times10^{-6}\).

The audit saves one initial usable direction, up to two spacing rejects after successful linear checks, and one feasible Armijo rejection. Names bind base, proposal and sequence. The receipts preserve raw primary and half-spacing actions, original models, signed free directions, actual stencil displacements, both norms and inverse provenance. Their nonaccepted schema cannot supply an accepted checkpoint. The focused control reproduced a nonzero \(-.8\) correction with passing linear checks and failed sensitivity, and a separate feasible rejection with remainder \(5\). Neither created an accepted observation. The failed-linear Arnoldi trigger remains separate.

## Action accuracy and baseline progress

| Fixed direction | Physical model error | Scaled model error | Physical allowance | Scaled allowance |
| --- | ---: | ---: | ---: | ---: |
| First baseline at 512e−6 | 7.195e−7 | 2.113e−5 | 5.301e−7 | 7.302e−6 |
| Late baseline at 512e−6 | 1.689e−6 | 2.491e−5 | 3.099e−6 | 4.069e−5 |
| Same late vector at 2048e−6 | 1.781e−6 | 2.549e−5 | 7.632e−7 | 1.267e−5 |

The late case's forcing targets were \(1.932\times10^{-6}\) physical and \(4.283\times10^{-5}\) scaled, with allowance limits at half those values. At 512e−6 the direction passed both true linear tests and failed both sensitivity tests. At 1024e−6 its physical allowance was \(1.476\times10^{-6}\). The single wider check at 2048e−6 passed both tests on the same base, vector and scales. The sampled stencils were central with zero feasibility reductions. The repeated actions were bitwise equal. The late parameter step was \(1.043\times10^{-6}\) at 512e−6 and \(4.171\times10^{-6}\) at 2048e−6. Taylor remainders decreased by about four on halving. The first current action came from a cold factor at the latest field and differs from R8's reused-factor direction.

The baseline accepted five fractions from \(3.860\times10^{-5}\) to \(3.088\times10^{-4}\) and reduced the original merit by **0.0431%**, to \(6.43842382363\times10^{-6}\). The maximum scaled direction norm fell from 210.8 to 48.0 across the successful bases, then rose to 1154.6 at the captured reject. Accepted displacements ranged from .00814 to .02493 in the fixed scaled infinity norm.

The first feasible rejected fraction .00014500887 used the same direction and factor as the accepted fraction \(3.85993\times10^{-5}\). The \(k\) remainders were \(1.018\times10^{-7}\) and \(7.211\times10^{-9}\), and the outer-\(\epsilon\) remainders were \(5.321\times10^{-7}\) and \(3.770\times10^{-8}\). Both ratios were about 14.1, matching the squared fraction ratio. The active original norm stayed \(k\). This repeated R8's curvature-limited behavior with a checked current action.

## Matched coupled trial experiment

At a feasible rejected plain trial the opt-in method forms the full signed \(q=F_0(X+\alpha d)-R-\alpha J_0d\). It checks the blocking block against twice the primary action allowance at that fraction. The auxiliary RHS is \(-q\), and one-sided FD uses the actual \(F_0(X)\) as its base. The solve retains the same inverse and the unshifted current action. It can try one auxiliary correction at each of at most two feasible rejected fractions in a primary search. The shared product limit counts every auxiliary action, including the fresh and half-spacing checks.

The block error plus its action allowance must fit half the available slack from the primary prediction to the original Armijo threshold after the primary uncertainty. The scaled budget is the smaller of a quarter of the scaled defect and half the verified scaled model slack at the fraction. The auxiliary size limit is half the actual scaled plain displacement. Any corrected state rebuilds its constraints, passes the full guards, and uses the original Armijo test at the primary \(\alpha\). The controller keeps the primary model as the denominator of its reduction ratio. The product, size and error controls enter the numerical-policy fingerprint.

The two-component control changed a plain norm of 1.1000 to 0.9205 with four products and a true scaled error of \(2.79\times10^{-12}\) using a second-order one-sided stencil. The checkpoint observed the expected \((-0.3,-5.25)\) to numerical precision and kept the primary .1 fraction and .795 ratio. Its resumed and continuous states and controller ledgers agreed. The .15 comparison gave the predicted defect ratio 2.25. A smaller size cap preserved accepted custody. Native restart also reproduced the exact accepted coefficients, all residual blocks and the controller before any factor preparation. A one-product change in the auxiliary cap under the same restart entry failed the fingerprint check.

| B3 fraction at the same base, direction and factor | Plain norm | \(q_k\) / \(q_\epsilon\) | Auxiliary scaled error / budget | Auxiliary size / limit |
| --- | ---: | ---: | ---: | ---: |
| 6.17589e−4 | 2.92824e−4 | 5.586e−5 / 2.907e−4 | 1.421e−6 / 3.729e−8 | 417.81 / .35708 |
| 6.17589e−5 | 6.47873e−6 | 5.558e−7 / 2.905e−6 | 1.111e−8 / 3.729e−9 | 4.0223 / .03571 |

Both coupled defects changed by about 100 at the tenfold fraction. Both auxiliary solves exceeded their error budgets before a corrected PDE state could be evaluated. Their measured candidate sizes were also 1,170 and 113 times the respective limits. The fallback accepted \(6.91430\times10^{-7}\), with scaled displacement .000800 and original decrease \(2.207\times10^{-12}\). The next direction had scaled norm 4316 and failed the sensitivity allowance at 2048e−6. This run saved the actual rejected direction. The tested curvature method did not improve original reduction per work.

Both approximate auxiliary corrections had a cosine above .9999999993 with the primary correction in the fixed scaled coordinates. Their largest components were pressure and \(k\). The measured scaled action per unit \(\ell_2\) primary direction was about \(6.94\times10^{-9}\). This explains why the attempted correction required large displacements for such small residual changes.

## Fixed-rule diffusion discriminator

The late baseline action gap was dominated by \(k\) diffusion at cells 47 and 101 in \(x=[4.0625,5.625]\) and \(r=[0,6.25\times10^{-5}]\). Both use fixed integration. The physical \(k\) differences for diffusion, source and advection were \(1.54925\times10^{-6}\), \(9.63\times10^{-10}\), and \(2.30\times10^{-8}\), respectively. Their signed action total closed the measured \(k\) difference to \(4.71\times10^{-13}\). At identical spacing the factored product difference reduced the \(k\) difference to \(1.386\times10^{-6}\), leaving the sensitivity above target.

The exact tangent uses \(D_k\nabla\delta k+\delta D_k\nabla k\) on 5,130 fixed-rule cells. The derivative includes the transported fields and the projected dependent trace in the outer region, and the reference \(K\)-derivative of inner \(\epsilon\), including its cutoff-dependent plateau. All 108 moving-rule cells and the other original terms retain full current differences. The physical residual and factor assembly at \(4\times10^{-6}\) are unchanged. On the frozen late baseline vector the diagnostic tangent substitution reduced the maximum physical difference from \(1.549\times10^{-6}\) to \(4.636\times10^{-8}\) and the scaled difference from \(2.035\times10^{-5}\) to \(7.009\times10^{-7}\). This localizes the measured numerical effect to the fixed-rule contribution.

The integrated tangent at 512e−6 matched the earlier diagnostic replacement on the same vector. The fresh physical/scaled model errors were \(1.73846\times10^{-6}\) and \(2.48102\times10^{-5}\). The empirical allowances were \(9.27124\times10^{-8}\) and \(1.40180\times10^{-6}\). The state coefficients and all six \(F_0\) blocks reproduced exactly. The focused analytic-split control retained the remaining terms in a nontrivial one-sided stencil.

At the newer rejected vector from the curvature run, the 2048e−6 hybrid action still failed. Its physical model error was \(2.99232\times10^{-6}\) and physical allowance \(3.37717\times10^{-6}\). At 512e−6 with a 256e−6 witness the allowances fell to \(2.750995\times10^{-7}\) physical and \(4.317827\times10^{-6}\) scaled, within the original limits. The old vector still missed both fresh model error targets. The retained analytic-substitution and integrated-action differences agree. The dominant remaining scaled \(k\) row is shared by upstream cells 1136 and 1190 over \(x=[-8.125,-6.25]\). Cell 1136 intersects the certified cutoff range and retains moving quadrature. Cell 1190 uses fixed quadrature. The signed residual after the tangent substitution was dominated by the remaining \(k\) diffusion difference at that row. The finer hybrid attempt therefore began a new Krylov construction at 512e−6 with the factor held at 4e−6 and coupled trial correction disabled.

The first rebuilt fine hybrid direction passed both fresh residual checks and both spacing allowances. Its physical model error was \(9.116\times10^{-7}\) and physical allowance \(2.697\times10^{-7}\). Its scaled infinity norm was 4,342. A feasible trial at \(1.45009\times10^{-4}\) returned \(2.27792\times10^{-4}\). The model predicted \(6.43754\times10^{-6}\). Reducing the fraction tenfold reduced its \(k\) and outer-\(\epsilon\) defects by about 100, while the original norm still increased. After four feasible rejections the solver accepted \(4.72853\times10^{-8}\), with scaled displacement .00020533 and original decrease \(1.443\times10^{-13}\).

The next direction grew to 17,053 in the scaled infinity norm. Its physical and scaled true errors were \(6.466\times10^{-7}\) and \(1.272\times10^{-5}\). The measured physical allowance \(1.1064\times10^{-6}\) exceeded \(9.6576\times10^{-7}\) and the scaled allowance passed. The late derivative parameter displacement was \(7.0614\times10^{-8}\). The 1024e−6 / 512e−6 pair was limited by \(k\). The 512e−6 / 256e−6 pair was limited by momentum. Further halving increased the discrepancy. This run stopped before trying a nonlinear step on the unresolved direction.

The terminal no-factor split reproduced both captured actions bitwise. The physical momentum change was \(5.5320\times10^{-7}\). Its constitutive term had norm \(5.5308\times10^{-7}\). Removing ordinary pressure, sampled pressure and sampled \(k\) left \(5.5319\times10^{-7}\). The signed term split closed to \(2.21\times10^{-11}\). The largest scaled row was in the streamwise velocity at upper-domain cells 5211 and 5212 near \(x=1.484375\), \(r=[16.1625,19.9125]\). This is the specific remaining action limitation on the captured vector.

The baseline, curvature fallback and fine hybrid run together accepted seven updates and reduced the original norm by **0.04318%** from the R8 checkpoint. The case meets the operational criterion for changing the remedy. The latest two attempts delivered no trend toward larger accepted fractions or useful contraction. The next discriminating numerical experiment should measure the achievable original linear error versus correction size on a frozen accepted base, focusing on the amplified pressure/\(k\) direction. That experiment tests whether both original forcing targets can be met inside a useful scaled radius.

## Qualification, spatial decision and cost

**B3 original qualification remains open.** The best accepted coefficients are `ced15ecef5037f0e638664999a6ec82f46e02585a95895f907aa14ed3704c600`. The native 10×16 blocks are \(3.28660\times10^{-7}\) momentum, \(2.37866\times10^{-12}\) continuity, \(6.43842147203\times10^{-6}\) \(k\), and \(2.85784440493\times10^{-6}\) outer \(\epsilon\). Whole-cell scalar bounds and native covariance passed, and the sampled principal response was \(+.00400912\). The whole-interface \(K/K_{in}\) range was [.95679514, 1.03736834].

The 16×24 replay retained the exact coefficients and passed the direct/assembled native action comparisons. Independent order-28 blocks were \(3.28660\times10^{-7}\), \(2.37914\times10^{-12}\), \(6.43841856811\times10^{-6}\), and \(2.85784437207\times10^{-6}\). The projected trace defect was \(6.29\times10^{-16}\), and the pointwise trace approximation was \(1.1014\times10^{-6}\). The mass identity and algebraic power defects were below \(1.7\times10^{-15}\). The native equation failure keeps the independent margin and point-owner qualification open. The retained field remains a pending diagnostic checkpoint. Spatial outcome **D** remains unresolved and the third pair is pending.

The retained 2,175-location, four-approach core has SHA-256 `73c264cf3d237311049fec71f07c32fa44a5b38c2cd2cd308c85103c9f62ab98`. The existing A2/A3 difference is .02551944209, about 33.6 times the .0007594305446050912 ceiling. The 5,292-cell mesh remains an unsolved generator control. This diagnostic work leaves the full-height physical and `1.25*C_total <= 86400` second production qualification open.

The numerical campaign began at 20:08:28 Asia/Shanghai on 2026-11-26. Linux evidence collection reached 22:16:20 after 7,672 elapsed seconds of the 10,800-second allocation. The 3,128-second balance included the 1,800-second qualification reserve. The 1,250-second estimate for another cold factor plus output reserve would have used almost all of the nonreserved balance. The observed direction sizes and contraction provided no basis for repeating that allocation. The extension remained unused.

The largest observed aggregate at launch was 66 counted workers, including other jobs. The 47-rank external job showed one compute thread per rank. Its blocked MPI service threads were recorded separately. Builds used eight workers and numerical-library threading was one. The native solve lane was serial. Independent audits and diagnostics overlapped within the allocation.

| Native attempt | Process wall (s) | Assembly / factor (s) | Remaining solver work (s) | Original decrease per remaining solver second | Peak RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 512e−6 baseline | 1,447.59 | 913.16 / 74.76 | 407.28 | 6.823e-12 | 7.886 |
| Coupled curvature and fallback | 1,586.74 | 1,035.14 / 95.99 | 399.38 | 5.527e-15 | 7.888 |
| Fine hybrid tangent | 1,528.78 | 1,097.45 / 101.99 | 265.60 | 5.432e-16 | 7.886 |

The rate denominator is native solver time less assembly and factorization. It includes warm corrections, residual checks, trials and associated diagnostics. Main processes sum to 4,563.11 s. Separate native work sums to 1,599.41 s and independent equation processes to 1,495.80 s. Build processes sum to 1,459.38 s, including the first attempt that stopped on a duplicate test identifier. Configure and focused probe timers add 3.88 s. These overlapping process sums have separate scopes from elapsed campaign time.

| Tested scope | Source revision | Linux solver SHA-256 |
| --- | --- | --- |
| Capture controls and baseline | `8bde00dfd086131ebf351d937ca778a20ec317ef` | `639fd531f615864b11ab09440458576bd11076f92379b2726490b498b5e947a0` |
| Matched curvature and scalar term audit | `fcbf08698a74a9766605a5c09d20027f3f6e45ff` | `a6744203e838dfbc844343642bc0a11a9f6bfeb2e287d605da2b806e0af8a3de` |
| Integrated tangent and final accepted field | `5be879349f5d8581bdaf2c0e15524163e4903742` | `f8d8cb6724651a803196ef0951b304c26822e4cc604c40ba3e74b4977bb0308e` |
| Terminal momentum replay | `76b74ec07567811364c244650dd07175e07c0ab4` | `c94625eff672026beb59813dae00b612493e961138e5f03715762073707dd165` |

The remaining diagnostic builds and source fingerprints are pinned in the manifest. The terminal momentum build added an opt-in observation and replayed the final field and both rejected actions bitwise. The nonlinear trajectory and its unchanged-policy restart use the v6 identity.

## Gates

| Gate | Status | Scope |
| --- | --- | --- |
| `BEST_B3_STATE_REPRODUCED` | PASS | Exact signed coefficients, equation bindings and six original blocks in the current-source replay. |
| `MERIT_REJECTED_DIRECTION_CAPTURE_VERIFIED` | PASS | Focused custody controls and three actual saved ACTION_UNRESOLVED directions. Representative feasible rejected-trial identities verified. |
| `CURRENT_DIRECTION_ACCURACY_RESOLVED` | OPEN | Scoped frozen repair and first rebuilt hybrid direction passed. Terminal hybrid physical allowance 1.1064e-6 exceeds 9.6576e-7. Constitutive momentum sensitivity is localized. |
| `ORIGINAL_PROGRESS_RATE_ASSESSED` | PASS | Seven actual original updates give 0.04318% overall improvement. Last two decreases were 2.207e-12 and 1.443e-13. |
| `COUPLED_CURVATURE_REMEDY_ASSESSED` | PASS | One-sided, size and restart controls pass. Both matched B3 auxiliary solves were ineffective under the declared error and size budgets. |
| `B3_ORIGINAL_ROOT_QUALIFIED` | OPEN | Native maximum 6.43842e-6 versus 1e-10. Independent k exceeds its limit. |
| `PAIRED_AXIAL_LEVEL_THREE_COMPLETE` | OPEN | A3 retained qualified, B3 pending. |
| `NEXT_SPATIAL_DISCRIMINATOR_COMPLETED` | OPEN | Qualified third pair needed to select the next spatial solve. |
| `POSITIVE_ORIGINAL_ANCHOR_RESOLUTION_SUPPORTED` | OPEN | Paired and parent-child spatial evidence pending. |
| `ADVERSE_ORIGINAL_BRANCH_RESOLUTION_SUPPORTED` | OPEN | Retained sampled response positive. Converged independent branch evidence pending. |
| `HEIGHT_CONTINUATION_RELEASED` | OPEN | Fixed-height spatial admission pending. |
| `FULL_HEIGHT_ORIGINAL_PHYSICAL_PILOT_PASSED` | OPEN | Physical time advanced zero. |
| `C1_CLOSED` | OPEN | Full-height independent qualification and pilot pending. |
| `C2_READY` | OPEN | C1 closure pending. |
| `LONG_FLAT_PRODUCTION_COST_CLOSED` | OPEN | Separate production criterion remains unqualified. |
| `REMOTE_PUSH_VERIFIED` | OPEN | Publication check pending. |

## Reproduction and continuation

The Linux evidence is under `/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM`. The manifest pins archive inputs, signed coefficients, session/equation identities, executable/source identities, exact launch scripts, accepted lineage, term and direction payloads, and timing receipts. [Machine-readable evidence](wj-h27-c1-b3-direction-accuracy-curvature-evidence.json)

Use a fresh output directory with the existing tested binaries. These commands produce the baseline, coupled and fine hybrid comparisons from their pinned inputs. The native wrapper records exit status `1` for a bounded pending solve. The original acceptance criteria remain active.

```bash
base=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
record="$base/runs/wj-c1-b3-dir-curvature.SfYckSXM"
old="$base/runs/wj-c1-b3-unshifted.DNPSiHDh"
r=$(mktemp -d "$base/runs/wj-c1-b3-direction-reproduce.XXXXXXXX")
source "$base/env.sh"
mkdir "$r/tmp"
source "$record/threads.sh"
launcher="$NEKTAR_RUNS/wj-c1-steady-recovery.LOp8mfob/run-capped.sh"
prefix=solvers/IncNavierStokesSolver/Tests/RANS/h27
common=(--eta 0 --quadrature 10 16 --factor-fd 4e-6 --newton 6
  --inverse ColoredSparseV1 --inverse-mass .009539982703780335
  --maximum-inverse-age 6 --reuse-linear-iterations 30 --action-gap-audit
  --direction-model --globalization ModelBacktrackV1 --steady-forcing .3 .001 .3
  --merit-mass 0 --merit-reanchor 0 --merit-entry Reanchor --merit-work 96)

python3 "$record/source/$prefix/c1_steady_case.py" "$old/B3-seeded-trials/run.xml" "$r/baseline" \
  --seed "$old/B3-seeded-trials/state.accepted-000001.json" "${common[@]}" \
  --fd .000512 --seconds 1890 --merit-trial-fraction 1.4500886948443608e-4 \
  --accepted-checkpoint 2100 55 1250
bash "$launcher" "$r/baseline" "$record/build/IncNavierStokesSolver" 2100

python3 "$record/v3/source/$prefix/c1_steady_case.py" "$record/baseline/run.xml" "$r/curvature" \
  --seed "$record/baseline/state.accepted-000005.json" "${common[@]}" \
  --fd .002048 --seconds 2190 --merit-trial-fraction .0006175886032317118 \
  --coupled-curvature 18 .5 .5 --accepted-checkpoint 2400 55 1250
bash "$launcher" "$r/curvature" "$record/v3/build/IncNavierStokesSolver" 2400

python3 "$record/v6/source/$prefix/c1_steady_case.py" "$record/curvature/run.xml" "$r/hybrid-fine" \
  --seed "$record/curvature/state.accepted-000001.json" "${common[@]}" \
  --fd .000512 --seconds 2150 --merit-trial-fraction 1.4500886948443608e-4 \
  --fixed-k-tangent --accepted-checkpoint 2360 55 1250
bash "$launcher" "$r/hybrid-fine" "$record/v6/build/IncNavierStokesSolver" 2360
```

The exact probes and guarded launches are retained as `run.sh`, `late-panel.sh`, `curvature.sh`, `stable.sh`, `tangent.sh`, `hybrid.sh`, `hybrid-fine.sh`, and `terminal.sh`. `pending-counted.sh` runs the exact-coefficient 16×24 replay and independent order-28 audit. `resume-check.sh` contains the native policy control. Use the worker snapshot before launching numerical work. The final v6 Resume command was verified for exact coefficients, residual and controller before factor preparation.

The best accepted packet is `hybrid-fine/state.accepted-000001.json`, SHA-256 `84dae45ef585f24321137adb4b131ffc3157c15644ed1cd16ac7221b34caabb6`. An unchanged-policy restart uses the tested v6 executable and the saved controller. Use a fresh allocation. The following reproduces continuation from those coefficients with a cold inverse. The measured contraction motivates the frozen-base direction-size experiment before allocating another ordinary factor lifetime.

```bash
base=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
record="$base/runs/wj-c1-b3-dir-curvature.SfYckSXM"
r=$(mktemp -d "$base/runs/wj-c1-b3-direction-resume.XXXXXXXX")
source "$base/env.sh"
mkdir "$r/tmp"
source "$record/threads.sh"
bash "$record/workers.sh" "$r/workers.ps" 68
seed="$record/hybrid-fine/state.accepted-000001.json"
printf '%s  %s\n' 84dae45ef585f24321137adb4b131ffc3157c15644ed1cd16ac7221b34caabb6 "$seed" | sha256sum -c
suite="$record/v6/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
python3 "$suite/c1_steady_case.py" "$record/hybrid-fine/run.xml" "$r/case" \
  --seed "$seed" --eta 0 --quadrature 10 16 --fd .000512 --factor-fd 4e-6 \
  --seconds 1890 --newton 6 --inverse ColoredSparseV1 --inverse-mass .009539982703780335 \
  --maximum-inverse-age 6 --reuse-linear-iterations 30 --action-gap-audit --fixed-k-tangent \
  --direction-model --globalization ModelBacktrackV1 --steady-forcing .3 .001 .3 \
  --merit-mass 0 --merit-reanchor 0 --merit-entry Resume --merit-work 96 \
  --merit-trial-fraction 1.4500886948443608e-4 --accepted-checkpoint 2100 55 1250
bash "$NEKTAR_RUNS/wj-c1-steady-recovery.LOp8mfob/run-capped.sh" \
  "$r/case" "$record/v6/build/IncNavierStokesSolver" 2100
```
