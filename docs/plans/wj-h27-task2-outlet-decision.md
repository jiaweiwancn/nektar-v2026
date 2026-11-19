# Task 2 flat WJ–AWF outlet decision

**Advance B for further flat qualification.** The zero-transverse trace with the existing prescribed normal load passes the corner-directed two-mesh screen and reduces the adjacent raw divergence. The original full-traction policy A remains the default. Its matching M0–M1 comparison misses 16 numerical rows. Both policies still show a separate radial collar sensitivity on M1–M1c, even though both M1c states meet the current equations and physical and inlet-projection gates. [Outlet screen](wj-h27-task2-outlet-cp2.md) · [Isolated collar control](wj-h27-task2-outlet-cp3.md)

## Outlet and experimental region

On the 111×63 mesh, the B near-annulus D RMS is 5.86e−5 compared with 1.44e−3 for A. The annulus spans 11 owners. The full-domain strong divergence RMS is about 1.95e−4 under either law. The B vertex value follows from the prescribed traces. At the outlet, all four matched-mesh B maxima in stress, Pi, shear and gradients pass their inherited numerical limits. They occur in the collar at r=.060–.069. A's corresponding maxima occur at the ground vertex. [Raw outlet and annulus audit](wj-h27-task2-outlet-cp2.md#current-equations-physical-bounds-and-tracings)

Over the fixed experimental window x in [−2.5,8], the matched M1c boundary changes reach 1.264e−5 for U, 3.664e−6 for W, 1.926e−4 for k/Kin and .002903 for Pi/tau0. The summed A and B collar indicators resolve the U, W, Pi and shear RMS effects. Their local maxima remain unresolved. The k and outer epsilon RMS and MAX effects resolve. Both policies satisfy the physical inlet-comparison limits. [Regional effect and both radial indicators](wj-h27-linux/task2-outlet/cp3/outlet-collar-effect.json)

## Second collar refinement

B2 bisects all four M1c children within the original two collar intervals. It preserves every x edge, the wall and terminal widths and the remaining radial edges. The sequence is 111×63, 111×65 and 111×69 with changes confined to the same r interval. The four added B2 radii are .0859375, .0953125, .101602194142255 and .104806582426765. This tests whether the radial differences contract after the localized reduction in outlet sensitivity. [Measured admission and frozen coordinates](wj-h27-linux/task2-outlet/cp4/outlet-second-collar-plan.json)

The first collar bisection reduces the whole-domain strong D RMS by about 37% under either policy, from about 1.95e−4 to 1.24e−4. It leaves A's vertex D near .453. This controls the radial effect separately from the outlet vertex behavior. In the B1c−B1 comparison, the maximum stress change is in Rxx at the upper r=.10 outlet owner, the Pi maximum is at the inlet r=.10, and the gradient maximum is just downstream of the inlet. [Both signed collar comparisons](wj-h27-task2-outlet-cp3.md#radial-sensitivity-under-each-outlet)

The first two B2 attempts stopped at shared-host worker censuses of 95 and 75. They charged 441.785 s and 919.168 s. The third stopped after a bounded suspension with a 174.660 s charge, including 57.124 s suspended. None produced a native iterate. The user then raised the shared-host ceiling from 72 to 80 workers for the remaining outlet work. The wall-time ceiling of 54,000 s and the 3,600 s protected reserve continue to govern the ledger. [Prior worker-capacity readmission](wj-h27-linux/task2-outlet/cp4/outlet-b2-final-admission.json) · [Worker ceiling amendment](wj-h27-linux/task2-outlet/cp4/outlet-worker-amendment.json)

The revised forecast uses the B1c Newton history. Its first two updates used 41 combined linear-check records. The third iteration alone used 274 before the deadline. Shorter native segments retain typed, accepted Newton iterates and use the existing continuation path. The B2 restart after the worker interruptions begins with the sealed B1c guess. Each attempt uses one MPI rank. The watchdog accounts for task-owned stopped workers when projecting the census at resume. [Segmented forecast](wj-h27-linux/task2-outlet/cp4/outlet-b2-worker80-admission.json)

| B2 segment | Charged wall seconds | Peak tree RSS | Peak shared census | Native residual | Outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| Fresh seed from B1c | 761.209 | 8.48 GiB | 81 | 9.2481e−8 | Typed Newton checkpoint at deadline |
| Typed continuation | 595.485 | 8.44 GiB | 74 | 3.4981e−11 | Converged at zero simulated time |

The first of these segments records two watchdog pauses totaling 10.4 s during shared-host pressure. The census briefly sampled 81 while task-owned processes were stopped. The continuation ran without suspension. OpenMPI initially bound both segments to CPU 0, which was busy with other ranks. Only the B2 rank moved to an idle core on the same NUMA node in each segment, using CPU 14 and CPU 15 respectively. [Watchdog receipts](wj-h27-linux/task2-outlet/cp4/operations.json) · [Affinity records](wj-h27-linux/task2-outlet/cp4/b2-affinity.json)

B2 passes all independent current, physical, inlet, boundary, native guard and interface gates. The assembled momentum, continuity, k and outer epsilon residuals are 2.64e−10, 3.13e−14, 3.92e−10 and 3.62e−10. The accounting identity error is 9.52e−15. The ground shear RMS/MAX are .01209/.02198 relative to tau0 and .01021/.01389 relative to Tref(Kin), both scaled by tau0. The whole-domain strong D RMS falls from 1.236e−4 on B1c to 1.194e−4 on B2, and the near-annulus RMS is 5.86e−5 on both. [Independent B2 assembly and raw traces](wj-h27-linux/task2-outlet/cp4/outlet-B2-audit.json) · [Inlet projection](wj-h27-linux/task2-outlet/cp4/outlet-B2-audit-inlet.json)

The second collar bisection reduces misses from 180 to 122 out of 690. All 58 gradient misses clear. The remaining misses are 60 Pi, 58 stress and four shear rows, spanning all 29 stations. The axial and transverse velocities, k and local outer epsilon continue to pass in both comparisons. [All signed rows, extrema and one-sided samples](wj-h27-linux/task2-outlet/cp4/outlet-second-collar-result.json)

| Full-domain field | B1c−B1 RMS / MAX | B2−B1c RMS / MAX | RMS and MAX limit |
| --- | ---: | ---: | ---: |
| R / tau0 | 7.61e−5 / .012257 | 1.62e−5 / .002188 | .001 |
| Pi / tau0 | 1.25e−4 / .009157 | 2.20e−5 / .001353 | .001 |
| Vxr / tau0 | 1.05e−5 / .008530 | 1.35e−6 / .002300 | .001 |
| Velocity gradient in code units | 1.28e−4 / .014687 | 2.76e−5 / .003235 | .01 |

Each norm contracts in 680 of the 690 rows across the two radial differences. The largest remaining stress and shear changes occur at the outlet, on the upper owner of r=.08125 in the final x column. The Pi maximum is in that column at x=14.999467, r=.083773. In the fixed experimental region, the stress and Pi maxima still exceed .001, at .001186 and .001147. The largest interior stress difference is Rrr at x=8 on the upper r=.090625 owner. Resolving the lower half of the collar is the next spatial qualification step for B. [Peak locations and contraction ratios](wj-h27-linux/task2-outlet/cp4/outlet-second-collar-decision.json)

At r=.10 and x=5 the B1c lower and upper owner D values are 3.18e−4 and 4.03e−4. On B2 they are 1.86e−5 and 2.71e−5. At x=15 they decrease from 9.59e−4 to 4.32e−5 on each owner. Removing the planar gradient trace only for constitutive reevaluation at these r=.10 samples reduces the maximum Ryy difference from .02229 to 1.69e−5 under A, from .00754 to 1.69e−5 for the first B bisection, and from .000756 to 6.83e−7 for the second. This isolates a strong trace sensitivity in the normal stress at r=.10. The solver and acceptance tables use the original gradients. [Trace sensitivity with paired raw data](wj-h27-linux/task2-outlet/cp4/outlet-collar-normal-stress-diagnostic.json)

## Outlet-location feasibility

The prospective 118×69 extent test moves the outlet from x=15 to 17.5. It retains all upstream x edges through x=14.375, translates the final grading and retains x=15 as a physical station. Its frozen axes, common observation window and physical RMS/MAX limits are in the feasibility record. The forecast totals 2,464 s including the new build, retained-case parity, native solve, analysis and publication. The forecast peak is 9.92 GiB within the 10 GiB solve cap. The live nonprotected balance at the check was 1,180 s with the evidence reservation in place, a gap of 1,284 s. The extent test is deferred, so outlet-location sensitivity remains open. [Frozen extent test and measured forecast](wj-h27-linux/task2-outlet/cp4/outlet-location-eligibility.json)

## Reproduction and evidence

The registered Linux application uses the Release OpenMPI build on `jay@100.88.70.60`. The numerical binary embeds source `54b2f4268a5903fcf894f032b968bbfcb42a3e999a9c0409d032e999114894a8` and has SHA-256 `3e8db421575122676808adc6c9783f0c897136e9d9d283b9d0cb200c5ec2df45`. The exact-inner table SHA-256 is `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`. The mesh, state, signature, serialized coefficient packet and receipt identities are in the individual audits. The ext4 workspaces, exact session and option files, signed native outputs and stdout are in the Linux archives. [Runtime build and bound sources](wj-h27-linux/task2-outlet/cp1/outlet-runtime.json) · [Final run receipts](wj-h27-linux/task2-outlet/cp4/operations.json)

The final Linux settlement records 49,165.082 cumulative charged seconds, 1,234.918 nonprotected seconds available, 3,600 protected seconds and zero reserved seconds. The independent B2 audit charged 152.960 s and the collar comparison 76.298 s. The archive records all prior stops and the evidence-summary retry. Its SHA-256 is `ff78f310092c58ed682b7a9836acb277f81b5ad7212cb455c5eaa4e439ca29fb`. [Settled account](wj-h27-linux/task2-outlet/cp4/accounting.json) · [Archive and resource receipts](wj-h27-linux/task2-outlet/cp4/operations.json)

```bash
W=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
R="$W/runs/wj-h27-task2r-20261118"
source "$W/env.sh"
D="$R/tools/outlet-python-v12/h27"
S="$R/tools/outlet-control-v6/task2r_run.py"
export PYTHONPATH="$D:$R/tools"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python3 -B "$S" "$R" balance
# B2 uses the separately frozen and sealed plan.
python3 -B "$S" "$R" run outlet-b2-continue-3-inputs analysis 18 -- \
  python3 -B "$D/task2_outlet_case.py" "$R" prepare B2 650 \
    --continuation 3 --guess b1c-after-worker-stop
python3 -B "$S" "$R" run outlet-b2-continue-3-steady outlet-science 765 -- \
  mpirun -np 1 "$R/runtime-outlet-v1/IncNavierStokesSolver" \
    --use-opt-file "$R/science/outlet-b2-continue-3/run.xml.opt" "$R/science/outlet-b2-continue-3/run.xml"
python3 -B "$S" "$R" run outlet-b2-continue-4-inputs analysis 18 -- \
  python3 -B "$D/task2_outlet_case.py" "$R" prepare B2 570 --continuation 4
python3 -B "$S" "$R" run outlet-b2-continue-4-steady outlet-science 675 -- \
  mpirun -np 1 "$R/runtime-outlet-v1/IncNavierStokesSolver" \
    --use-opt-file "$R/science/outlet-b2-continue-4/run.xml.opt" "$R/science/outlet-b2-continue-4/run.xml"
python3 -B "$S" "$R" run outlet-b2-audit analysis 182 -- \
  python3 -B "$D/task2_outlet_analysis.py" "$R" audit "$R/evidence/outlet-B2-audit.json" \
    --case B2 --seconds 168
python3 -B "$S" "$R" run outlet-second-collar-compare analysis 96 -- \
  python3 -B "$D/task2_outlet_analysis.py" "$R" collar "$R/evidence/outlet-second-collar-result.json" \
    --case B2 --seconds 87
```

These commands identify the recorded recipe. The operations JSON gives the exact invocations through the execution lock, including reservation caps and any typed continuation. Reproduction requires fresh output paths and an admitted balance.

| Verified checkpoint | Git SHA |
| --- | --- |
| CP0 baseline and amendment | `a9882692b2eb2efcf4836161d80bb296dafae709` |
| CP1 native outlet verification | `097341dba98bf82cf74ba3d5e3ea421666bbc66d` |
| CP2 matched-mesh evidence | `dacf6a372beb43fa958cd78951474d8ad8708acf` |
| CP3 isolated collar evidence | `e44a2a3db2fbccd1902e4b44c227c2599eec3b35` |

The diagnostic uses zero simulated time throughout. The higher-quadrature and tighter-algebraic controls, including the inner-k gates, remain part of final-recipe qualification. The three-timestep, stock, operational restart and advancing-cost gates follow spatial qualification. Later hill adoption also requires a wake and separation/reattachment sensitivity study. [Task 3 launch hold](wj-h27-linux/task2-outlet/cp4/task3-launch-hold.json)

**TASK2_COMPLETE=NO. TASK3_READY=NO.**
