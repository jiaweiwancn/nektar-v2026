# WJ-H26-C — Qualified native flow kernel

H26_NUMERICS = **PASS**  
H26_COST = **PASS**  
H26_INTEGRATION_READINESS = **READY**  
FULL_COUPLED_MODEL_AND_PRODUCTION = **NOT_QUALIFIED**

The opt-in serial native Q3/global-continuous-Q2 flow kernel meets the H26 numerical, engineering-screen and evidence gates. The selected physical parameters, equations, solution tolerances and common 700-cell comparison mesh were retained. H27 scalar/wall-budget coupling and H28 production qualification remain separate tasks.

Tested native source: `8aeb00a2f37437cd17e043a9ba6f86a6f208e1d6`. Final publication tools: `6165f8bd94c4e845d1d5a2f7ca7484643cd3ed0c`. Subsequent publication commits change reports/Python evidence tools only; all compiled native/stock source/header/CMake inputs were checked identical to the tested revision.

| Checkpoint | Commit |
|---|---|
| H26-A: prospective contract/fixtures | `98a9e6182646a30da90c59d0896cad141b81e5f5` |
| H26-B: working verified core/profile | `80aced5edeaf73437de6df7a4e8a4370cbf63bff` |
| H26-C: qualified candidate/evidence | Commit containing this report; exact remote SHA is reported with delivery |

A and B were pushed and exact remote refs verified before continuing. C is committed/pushed and its exact remote ref verified as the final publication step. No force push is used.

## Current equations and executed numerical gates

The current operator is alpha*M + K_sym(mu), with full symmetric strain and physical molecular nu throughout the new domain. The pressure action is G-C_AWF, including the physical tangential pressure trace sampled at H and its roughness-layer extension. Pure high-Re WJ remains called with outerNu=0; the new PDE supplies molecular terms. The complete blended tensor enters the weak stress load. The nonlinear fixture uses current WJ coefficient/remainder, current pressure coupling and explicitly zero continuity forcing.

Geometry/topology, native quadrature/physical derivatives, signed continuous maps, sparse mass/Riesz operators and useful sparse preconditioners are retained. The physical mean-pressure operator remains current; only the preconditioner factors a pressure-pinned equivalent and projects back to the mean gauge. There is no legacy custom dense pressure/mass machinery in the ordinary 700-cell candidate. Tiny independent dense direct solves are correctness oracles only.

| Gate | Executed result | Target |
|---|---:|---:|
| Independent operator / natural total-traction discrepancy | 5.4101077e-15 | <=1e-10 |
| Independent tiny direct saddle agreement | 3.0665451e-12 | <=1e-8 |
| Independent momentum residual | 4.4108342e-14 | <=1e-9 |
| Independent continuity residual | 1.3407581e-14 | <=1e-9 |
| Complete blended-stress weak identity | 3.7554156e-17 | <=1e-10 |
| Zero-forcing nonlinear momentum / continuity | 6.6945404e-12 / 2.1000963e-16 | <=1e-9 each |
| Nonlinear target agreement / refreshes | 6.8474345e-12 / 6 | <=1e-8 |
| Adapter dimensionless identity error | 3.9385185e-16 | <=1e-10 |
| True-history restart at t=.4 | All four serialized state/history vectors identical | same completed time |
| Changed-H restart | Rejected before advancement | scientific identity must match |

The four independent cases cover flat/curved geometry with closed essential boundaries and a natural total-traction top. Viscosity is nonconstant and positive; both velocity components, cross strain, nonzero pressure and nonzero Dirichlet lifts are active. Closed cases have exactly one unaugmented pressure null mode; natural-traction cases have none. Constant pressure has zero C_AWF action. Work/dissipation and the known -2*C_AWF*k/3 action pass independent checks.

| Geometry / gauge | Operator/traction max | Direct error | Momentum | Continuity | 20% mu change: iterations / rebuilds |
|---|---:|---:|---:|---:|---:|
| curve=0.0, gauge=True | 1.58731e-15 | 3.06655e-12 | 3.26475e-14 | 1.34076e-14 | 5 / 0 |
| curve=0.0, gauge=False | 1.58731e-15 | 9.90992e-13 | 2.01273e-14 | 9.54995e-15 | 6 / 0 |
| curve=0.1, gauge=True | 5.04209e-15 | 1.73527e-12 | 4.28190e-14 | 1.30874e-14 | 5 / 0 |
| curve=0.1, gauge=False | 5.41011e-15 | 2.96315e-13 | 4.41083e-14 | 4.08777e-15 | 6 / 0 |

The20% coefficient-change tests hold the load fixed and genuinely change the solution (about1.08–1.27%). They retain the same preconditioner and solve the changed physical operator. Missing transpose strain, flipped AWF sign and stale viscosity mutations are detected. Independent Gauss quadrature and analytic graph geometry/Jacobi basis derivatives are used by the oracle. Norms are the preregistered physical reference and mass/dual norms, with no dimensionful max(1,rhs) floor.

The adapter passes the independent planar physical-root/C(q)/d(q) reduction, trace2k, general rotational covariance, raw stress work, full-WJ collar recovery, finite-nu physical total traction and selected realizability checks. Its minimum unit-fixture R eigenvalue/K is 0.5488038. No clipping is used. A full-WJ collar is essential to the declared extension; without it, nonzero-shear endpoint matching additionally requires nu_A(H)=C(q_H)K²/epsilon_H.

## Continuous manufactured convergence and history

Forcing is independently differentiated from the physical continuous PDE. Compiled polynomial and trigonometric formulas agree with independent SymPy expressions. It is verification-only forcing, not discrete A_h*u_h balance forcing or log-law preservation forcing. The advancing scheme uses genuine BE/EX1 startup, then BDF2/EX2, new-time coefficient/stress/forcing/BC data and true accepted histories.

| dt to t=.4 | Velocity L2 error | Pressure L2 error | Velocity order | Pressure order |
|---:|---:|---:|---:|---:|
| 0.04 | 2.07295966e-05 | 5.46029736e-06 | — | — |
| 0.02 | 5.05009180e-06 | 1.29923655e-06 | 2.037311 | 2.071315 |
| 0.01 | 1.24837605e-06 | 3.17457850e-07 | 2.016257 | 2.033027 |

| Geometry / n | Velocity L2 | Velocity H1 seminorm | Pressure L2 | Natural traction L2 | Weak traction moments | Strong divergence |
|---|---:|---:|---:|---:|---:|---:|
| 0.0 / 2 | 1.3096063e-01 | 1.8089778e+00 | 4.2339724e-01 | 1.0925117e-01 | 9.3306187e-02 | 2.5040105e-01 |
| 0.0 / 4 | 4.6127098e-03 | 1.8959146e-01 | 9.2724518e-03 | 1.1502582e-02 | 7.8455132e-03 | 8.7695789e-02 |
| 0.0 / 8 | 2.3540437e-04 | 2.2231103e-02 | 2.2247710e-04 | 1.3849780e-03 | 6.5922964e-04 | 1.0597193e-02 |
| 0.1 / 2 | 1.3131878e-01 | 1.8111553e+00 | 4.2260102e-01 | 1.1075015e-01 | 9.2620462e-02 | 2.5797660e-01 |
| 0.1 / 4 | 4.6766346e-03 | 1.9025714e-01 | 9.3477840e-03 | 1.1443857e-02 | 7.6486675e-03 | 8.7840264e-02 |
| 0.1 / 8 | 2.3693971e-04 | 2.2307283e-02 | 2.2624097e-04 | 1.3868216e-03 | 6.4330726e-04 | 1.0644619e-02 |

| Geometry | Final U L2 order | U H1 order | Pi L2 order | Natural weak traction order | WJ-collar weak traction order |
|---|---:|---:|---:|---:|---:|
| 0.0 | 4.292402 | 3.092243 | 5.381222 | 3.573015 | 2.944789 |
| 0.1 | 4.302879 | 3.092364 | 5.368692 | 3.571629 | 2.944406 |

Spatial tests use dt=.0001 to t=.01. Weak continuity passes the1e-9 target; strong divergence is reported separately and decreases under refinement. A weakly constrained C0 velocity is not expected to have pointwise zero divergence.

The dedicated WJ-collar trace fixture retains distinct native element gradients at H, common positive K/epsilon and full physical stress. Pressure traces agree exactly in these fixtures. Initial n=4->8 traction rates near2.74 missed the preregistered2.8 finite-grid target; n=16 resolves rates near2.944, with positive selected R spectra. The initial miss is retained in the raw evidence and the threshold was not lowered.

| Curve / n | WJ traction L2 jump/K | Weak moments/K | Minimum R eigenvalue/K |
|---|---:|---:|---:|
| 0.0 / 2 | 1.34419193e-02 | 4.41904282e-03 | 0.3045484 |
| 0.0 / 4 | 6.03203804e-03 | 2.06376919e-03 | 0.3313225 |
| 0.0 / 8 | 9.03633649e-04 | 3.08416761e-04 | 0.3340104 |
| 0.0 / 16 | 1.17374235e-04 | 4.00560655e-05 | 0.3341714 |
| 0.1 / 2 | 1.35193783e-02 | 4.90159087e-03 | 0.3057015 |
| 0.1 / 4 | 6.11982315e-03 | 2.53551324e-03 | 0.3316830 |
| 0.1 / 8 | 9.17216578e-04 | 3.80709550e-04 | 0.3348620 |
| 0.1 / 16 | 1.19142679e-04 | 4.94583135e-05 | 0.3351521 |

Tightening linear tolerance from1e-11 to1e-13 changes velocity/pressure L2 by 2.297304e-13 / 5.146823e-14, respectively 1.840e-07 / 1.621e-07 of their temporal discretization errors. Algebraic error is subordinate.

| Constant limit / solver | Velocity errors at dt=.02,.01,.005 | Last velocity order | Last pressure order |
|---|---|---:|---:|
| nu / S | 1.2313134e-05, 3.0803222e-06, 7.7036597e-07 | 1.999465 | 2.022915 |
| nu / N | 1.2326581e-05, 3.0819693e-06, 7.7055861e-07 | 1.999876 | 2.021794 |
| constant / S | 1.4750666e-06, 2.7480716e-07, 5.9544668e-08 | 2.206375 | 2.255082 |
| constant / N | 1.0571132e-06, 2.3437223e-07, 5.5704413e-08 | 2.072938 | 2.016464 |

S denotes actual stock VCSWeakPressure; N denotes the new mixed full-strain core. Both converge to the same continuous constant-coefficient problem. Pressure splitting, weak grad(div U), and boundary treatment prevent a general bitwise trajectory claim. The stock physical target-difference and weak-continuity measurements are published as such, not relabelled as candidate algebraic residuals.

The flow-only nonlinear fixture uses prescribed positive k/epsilon, current adapter data and zero continuity forcing. Six refreshes take 0.0648032s in the retained small case; constitutive work takes 0.0544947s. This includes its native incompressible target projection/setup and does not bound700-cell or transported-scalar nonlinear cost.

## Matched complete-step timing matrix

Linux host jay-Super-Server, AMD EPYC7763; one MPI rank, all six numerical-library thread controls1. Each block advances420 accepted steps at dt=.005 over[0,2.1], retaining400 warmed steps over[.1,2.1] after20 warm-up advances. Three repetitions are interleaved Snu,Nnu,Smu,Nmu,Nspace,Ntime. Every warmed block lasts more than3s. Native quad coordinates/topology are the unchanged materialized700 physical cells; all cases use velocityQ3, global continuous pressureQ2 and common GLL8x8 quadrature, startup and boundary traces.

Actual stock uses unchanged VCSWeakPressure/native IMEX(2), direct pressure Poisson and two component Helmholtz solves. Its prescribed forcing callback is analytically de-extrapolated so BE/EX2 produces the same new-time continuous forcing as the candidate. The scalar weight identity is verified in the pinned matrix configuration. Native uses converged full-strain nonsymmetric G-C coupling. Natural total traction is qualified separately. These are declared differences, not identical-discretization or fixed-ratio claims.

Physical nu=6.4231738035264469e-05; representative mu_ref=0.19399018791383052. The prescribed project-based spatial mu range is[0.000114909711292,0.38791614409]; the10% turbulent modulation has envelope[0.000109841913967,0.426701335325]. H=.15,h=.10854503542535461,d0=.0875 and the physical roughness/matching calibration were retained. The exact input functions/parameters were pinned before timing; the derived range audit changes none of them.

| Case | Block0 mean s/step | Block1 | Block2 | Warm linear iterations | Factor setups/block |
|---|---:|---:|---:|---|---:|
| Snu | 0.038118634 | 0.026196329 | 0.021274463 | 3 direct solves/stage | stock native caches |
| Smu | 0.032708908 | 0.027690679 | 0.021306116 | 3 direct solves/stage | stock native caches |
| Nnu | 0.043668492 | 0.039617895 | 0.040009853 | 1.000 mean, [1, 1] | 1 startup |
| Nmu | 0.075884846 | 0.045125761 | 0.051323148 | 1.000 mean, [1, 1] | 1 startup |
| Nspace | 0.052788301 | 0.039175056 | 0.043142183 | 1.000 mean, [1, 1] | 1 startup |
| Ntime | 0.121939577 | 0.073721098 | 0.093749613 | 2.875 mean, [1, 3] | 1 startup |

| Case / repetition | Mean | Median | SD | Min | P05 | P95 | Max | Warm seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Snu / 0 | 0.03811863 | 0.03495092 | 0.00795395 | 0.02375119 | 0.03062542 | 0.04965205 | 0.08687003 | 15.24745 |
| Snu / 1 | 0.02619633 | 0.02575013 | 0.00319093 | 0.02151339 | 0.02348982 | 0.02924536 | 0.06885489 | 10.47853 |
| Snu / 2 | 0.02127446 | 0.02078809 | 0.00318193 | 0.01817499 | 0.01872221 | 0.02449029 | 0.05502461 | 8.50979 |
| Smu / 0 | 0.03270891 | 0.03132232 | 0.00728073 | 0.02290164 | 0.02513032 | 0.04663707 | 0.07039435 | 13.08356 |
| Smu / 1 | 0.02769068 | 0.02627149 | 0.00421362 | 0.02306582 | 0.02435438 | 0.03572921 | 0.06249558 | 11.07627 |
| Smu / 2 | 0.02130612 | 0.02079807 | 0.00303343 | 0.01791273 | 0.01888962 | 0.02454320 | 0.05571537 | 8.52245 |
| Nnu / 0 | 0.04366849 | 0.04267810 | 0.00692272 | 0.03404073 | 0.03499878 | 0.05622964 | 0.07148788 | 17.46740 |
| Nnu / 1 | 0.03961789 | 0.03821908 | 0.00572750 | 0.03304379 | 0.03406271 | 0.05050280 | 0.08680698 | 15.84716 |
| Nnu / 2 | 0.04000985 | 0.03887173 | 0.00503290 | 0.03447113 | 0.03552125 | 0.04669659 | 0.08830010 | 16.00394 |
| Nmu / 0 | 0.07588485 | 0.07732840 | 0.02042913 | 0.04639071 | 0.05055272 | 0.11459214 | 0.21733377 | 30.35394 |
| Nmu / 1 | 0.04512576 | 0.04603856 | 0.00732032 | 0.03454664 | 0.03568815 | 0.05531693 | 0.11979322 | 18.05030 |
| Nmu / 2 | 0.05132315 | 0.04689471 | 0.01097449 | 0.03529422 | 0.03984063 | 0.06665600 | 0.11814450 | 20.52926 |
| Nspace / 0 | 0.05278830 | 0.05237213 | 0.00594414 | 0.04210873 | 0.04564641 | 0.06063051 | 0.09699486 | 21.11532 |
| Nspace / 1 | 0.03917506 | 0.03772699 | 0.00508997 | 0.03378207 | 0.03461542 | 0.05108775 | 0.07710188 | 15.67002 |
| Nspace / 2 | 0.04314218 | 0.03925274 | 0.00900909 | 0.03460424 | 0.03581099 | 0.05994674 | 0.10030343 | 17.25687 |
| Ntime / 0 | 0.12193958 | 0.13232852 | 0.03031261 | 0.04784862 | 0.07189892 | 0.16111860 | 0.21733314 | 48.77583 |
| Ntime / 1 | 0.07372110 | 0.06876690 | 0.01401125 | 0.03541963 | 0.05279266 | 0.09100003 | 0.12961096 | 29.48844 |
| Ntime / 2 | 0.09374961 | 0.08453998 | 0.01897616 | 0.05860513 | 0.06561364 | 0.12598709 | 0.17197425 | 37.49985 |

All candidate block current momentum/continuity residuals pass1e-9. Ntime changes the actual operator over turbulent modulation[.9,1.1], with a retained preconditioner; warm linear iterations average2.875 and range1–3. Nspace/Ntime each use one initial factorization and no warmed refresh in these blocks. Required warmed AWF feedback and pressure-correction minima are:

| Case | min dimensionless C_AWF Pi | min pressure correction/P_ref |
|---|---:|---:|
| Nspace | 1.265368700e-03 | 6.354911991e-05 |
| Ntime | 1.265386093e-03 | 6.354410296e-05 |

Both exceed the preregistered1e-6 nonzero-action/correction criterion. The actual700-cell restart reads/restores all histories at completed step420, then advances to440 with BDF2 history rather than replaying startup.

Across the required Nspace/Ntime warmed work, disjoint attribution is:

| Exclusive category | Fraction |
|---|---:|
| Coefficient/constitutive/source update | 8.106575% |
| Equation assembly/application/solution | 87.816314% |
| Warmed preconditioner setup | 0.000000% |
| True residual/diagnostics | 3.939786% |
| Compact I/O | 0.004512% |
| Controller overhead | 0.132814% |

Attribution closes to 100.000000000% (target>=95%). Inner Krylov convergence actions remain equation-solution work; the final true-residual/sanity checks are separately timed. Timing intervals are non-overlapping. Stock/new mean-of-block-means ratios are1.4406 for Nnu/Snu and2.1092 for Nmu/Smu; no invented stock-ratio gate was applied.

## Engineering screen and output accounting

Ordinary time includes coefficient/source updates, current matrices, BC lifts, convection/history/RHS formation, velocity/pressure solve, true residual/sanity checks and compact logging. Source/operator affine caches are exact for these prescribed functions; their construction is in startup. The largest warmed Nspace/Ntime block mean is used. Host variation is visible in the three block means; the worst block plus25% margin leaves a clear decision, so extending measurements is unnecessary.

```text
N_star = 101005
c_flow = 0.12193957717128796 s/step
C_startup = 2.161009422038 s
C_cold_start_transient_excess = 0.000000000000 s
C_final_output_and_assessment = 0.066757255001 s
C_other_fixed_conservative = 0.718937544618 s
C_actual_restart_read_restore = 0.003523349005 s
C_fixed = 2.950227570662 s
C_recurring = 505 * 0.028191346966 = 14.236630217929 s
C_screen = 1.25 * (C_fixed + 101005*c_flow + C_recurring)
         = 15417.117312468 s = 4.282533 h < 43200 s
```

Recovery occurs every200 accepted steps: floor(101005/200)=505 writes, each conservatively charged at the largest measured required-case write time. Latest-two recovery retention is exercised. Final output is charged once in fixed cost. The5050 compact summaries at cadence20 are already in ordinary means and are not charged again. Fixed cost includes measured startup, final output/assessment, one actual700-cell restart read/restore, and conservative external-process emission/teardown/controller overhead. Any positive cold first20-step excess above20*c_flow is charged once; here that excess is zero because the worst warmed mean already covers it.

The frozen/prescribed kernel screen passes with 27782.883s headroom. T_star and dt_star are historical budget yardsticks. This is not a qualified hill timestep, nonlinear coupled-cost bound or one-day production forecast. The largest remaining frozen-flow cost is equation application/solution (about87.8%); constitutive evaluation dominates the separate small nonlinear-refresh fixture.

## Resources, rejected prototypes and preserved work

The new campaign is CLOSED: 2371.270578 /21600 aggregate Linux operation seconds charged, 19228.729422s unused. The1800s evidence reserve was protected during scientific admission; the final120s publication epilogue is conservatively prepaid. These are operation seconds, not CPU-hours or projected production cost.

| Resource | Measured / bound |
|---|---:|
| Maximum recorded actual total host workers | 42 /72 |
| Largest start admission bound including helpers | 53 /72 |
| Compiler jobs | 4; initial limit8 |
| Peak build tree RSS | 2.457211 /16GiB |
| Peak science/analysis tree RSS | 0.665478 /8GiB |
| New source/build/runtime footprint | 0.833720 /8GiB |
| Generated science footprint | 0.125199 /8GiB |
| Free ext4 space | 70.725994GiB; reserve15GiB |

The isolated healthy ext4 namespace is `/home/jay/build/wj-h26-native-flow-20261117`, authorized because the historical heavy-storage alias was unavailable. No unrelated process was killed. The general Linux checkout remains at f4b70bae with its unrelated developer-guide modification. Both old H25/H25R1 closed-account hashes and4325 pinned old source files are unchanged. The selected4334 source files and320 immutable final input/result/receipt entries were checked, together with exact mapped runtime libraries. No new old-capsule bytecode cache was found. The Mac linked worktree preserves the pre-existing untracked pilot_wj_h25.py.

| Rejected or corrected path | Executed evidence and resolution |
|---|---|
| Immutable journal update | Two attempts failed before child launch; reconciled and conservatively charged. Mutable journal now uses locked atomic updates. |
| XML shorthand quadrature | Native Q3/Q2 grid mismatch guard failed; explicit per-direction basis/points declarations fixed it. |
| MPI argument array | Segmentation during MPI initialization; required argv null terminator added and all final executions pass. |
| Coordinate-keyed constitutive cache | Complete blended-stress weak identity failed; retain per-element quadrature gradients instead of merging C0 limits. |
| First full-step pilot | v3 averaged0.703085529s/step; exact affine caches and gauge-projected sparse preconditioning remove repeated work. Final decision uses the full final-binary matrix, not that short pilot. |
| Incomplete restart identity | H/h were initially omitted; current scientific signature pins them and constants, and the changed-H negative control is executed. |
| Initial collar resolution | n=4->8 rate~2.74 failed2.8 target; n=16 resolves~2.944 without changing the threshold or stress law. |
| Interpolated nonlinear continuity datum | Earlier algebraic g=D U_target was insufficient; current native projected target and final nonlinear problem explicitly use g=0. |

All failed receipts and recovered outcomes remain indexed on Linux; expected changed-H negative-control failures are distinguished from unresolved failures. There are no unknown process outcomes or unfinished scientific operations.

## Reproduction and evidence

Mac source/Git owner: `/Users/jay/Documents/nektar-v2026/.worktrees/wj-earsm-awf`. Linux is authoritative for all build, solver, numerical and timing conclusions. Use `ssh linux-usb`, with the verified fallback `ssh jay@100.105.247.110`. Activate `/home/jay/.local/opt/nektar-v2026-mpi/env.sh`.

Exact executed build/run entry points (the original outputs are immutable and the original journal is now closed):

```bash
R=/home/jay/build/wj-h26-native-flow-20261117
T=$R/source-v8/solvers/IncNavierStokesSolver/Tests/RANS/h26
U=$R/tools-v10/solvers/IncNavierStokesSolver/Tests/RANS/h26
python3 -B "$T/run.py" --root "$R" --label build-v8 --phase build --seconds 300 -- python3 -B "$T/build.py" --root "$R" --version 8
python3 -B "$T/campaign.py" --root "$R" --version 8 --mode verify
python3 -B "$T/campaign.py" --root "$R" --version 8 --mode constant_limit
python3 -B "$T/campaign.py" --root "$R" --version 8 --mode benchmark
python3 -B "$U/analyze.py" --root "$R" --version 8 --stage C --output "$R/publication-c-v8"
```

Fresh execution must use a new isolated parent with basename `wj-h26-native-flow-20261117`, its own fresh accounting journal, source/build/runtime directories, TMPDIR and science storage; do not reopen this closed journal or overwrite old capsules. Copy the frozen source-v8 to the new source-v1, make a new full file-identity manifest with the tested revision, and use version1 in build/campaign commands. Generate fixtures-v2 with fixtures.py and copy only the pinned700-cell profile-v3/case.xml template. The build helper verifies every reused stock object dependency against the copied source and uses four compiler workers. The resource runner enforces the shared-host72-worker/RSS/storage/free-space/finite-allowance limits. Existing successful operations in the original campaign are verified/reused rather than duplicated; each failed attempt needs a new explicit operation label.

Compiler, exact configure/link commands, loaded library hashes, ELF/wrapper hashes, source/dependency manifests, all XML/request identities and raw outputs are supplied by the following artifacts. The final publication manifest verifies every transported compact artifact and all operation request/receipt identities. Large matrices, field histories and recovery binaries remain at the indexed Linux paths.

- [Closure summary and complete numerical/block statistics](wj-h26-linux/checkpoint-c/summary.json)
- [All400 warmed timings/iterations/actions per block](wj-h26-linux/checkpoint-c/raw-timing-samples.json)
- [Linux raw input/result/receipt index](wj-h26-linux/checkpoint-c/linux-index.json)
- [Compiler/source/executable identities](wj-h26-linux/checkpoint-c/identities.json)
- [Loaded-runtime/resource/preservation audit](wj-h26-linux/checkpoint-c/runtime-resource-audit.json)
- [Derived fixed project viscosity range](wj-h26-linux/checkpoint-c/prescribed-model-range.json)
- [Closed finite campaign account](wj-h26-linux/checkpoint-c/accounting-final.json)
- [Final publication manifest](wj-h26-linux/checkpoint-c/publication-manifest.json)
- [Prospective mathematical/design contract](wj-h26-contract.md)
- [Mandatory working-core checkpoint](wj-h26-checkpoint-b.md)
- [Precise H27 coupling contract and H28 boundary](wj-h27-native-coupling-handoff.md)
- [Fresh isolated build/run reproduction recipe](wj-h26-reproduction.md)
