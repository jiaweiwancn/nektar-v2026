# WJ-4 native finite-strip qualification

Overall qualification is **OBSTRUCTED** for both policies by the L1/Nx2
horizontal divergence budget. Both fine-normal timestep pairs reach the frozen
T4 and pass the engineering field/interface budgets at declared observations.
Their longer results remain **CONDITIONAL finite-strip evidence**.
Stationarity is **NOT_REACHED**.

The producing experiment is retained at
`/media/jay/新加卷/nektar-v2026-runs/wj4-20261110`. Its runtime was preserved
while a separate source-matched observer was built under `readback-runtime`.
The [contract](wj4-contract.md) fixes the equations, engineering budgets,
work limits and comparison definitions. This report distinguishes execution,
accuracy, numerical-policy intervention and stationarity.

The problem is the finite reference-Dirichlet strip **[0,1] x [0,2]**.
Reference differential-field traces are prescribed at both lateral boundaries
and the top. There is no periodic/outflow implementation. The inner reference
retains molecular viscosity, and budget k is nonconstant. Pressure is
pi=p/rho+2k/3 with one global gauge and distinct regional pressure traces.
The selected WJ tensor/root/gradient completion, physical epsilon, independent
scalar diffusivity and uncapped EP/budget sources remain unchanged.

The outer sources remain Sk=Praw-epsilon and
Sepsilon=Ce1*Praw*epsilon/k-Ce2*epsilon^2/k, with scalar diffusivities
nu_d/sigma_k and nu_d/sigma_epsilon, nu_d=Cd*k^2/epsilon. Defaults are
c1=1.8,c2=5/9,Ce1=1.44,Ce2=1.82,sigma_k=1,sigma_epsilon=1.3,
Cd=572/6561,qe=81/(2 sqrt(143)),kappa=.381917863546733. Wall sources are
P_A-p_A,0 for EP and P_A-E_c for budget. There is one outer physical-epsilon
trace, no inner epsilon PDE and no additional independently prescribed flux.

## Implementation and immutable lineage

The reviewed parent is `11733a64fdf9452adc14c45594820c3da2513fd7`.
Mac source/Git work uses the existing `codex/wj-earsm-awf` worktree, preserving
other branches, dirty work, running executables and installed dependencies.
Every local commit was pushed immediately. Both specification originals and
worktree copies retain the exact reviewed identities:

```
1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124  main
0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c  appendix
```

The main eq:integrals-linear literal `quad` is the documented `\quad` erratum;
the intended logarithmic coefficient remains `(d+a*h)/a^2`.

Qualification now requires explicit component results; missing cross runs,
transactions, incomplete targets, unresolved meaningful differences, invalid
statuses, negative/nonfinite norms and nondecreasing meaningful sequences cannot
produce PASS. Re-evaluation of retained WJ-3 evidence passes both policies,
without altering its original receipts. Execution success is separate.

Write-once experiment/request packets bind argv, exact and semantic XML,
mesh/dt/target, initial/restart checkpoint, code, runtime/dependencies and
report/log outputs. Persistent package ledgers reserve work before launch,
retain unknown interrupted reservations and include numerical process-session
wall time, Python algebra, observations, SVD and I/O. Native negative zero is
read losslessly. Arithmetic caches have separate input/output fingerprints.
Measurement-only replay creates a new code/output lineage and verifies original
producer receipts, with zero new governing attempts.

Common physical quadrature splits the union of horizontal element boundaries,
x=.25,.75, normal elements and physical breakpoints. Orders4 and6 independently
check x integration. Full, central-strip and lateral/top trace errors retain
separate measures and maximizing element/face provenance. A single full-domain
pressure-gauge shift is applied before regional extraction. Surface measures
are5 (epsilon3); full areas are2 (epsilon1); central areas are1 (epsilon.5).

The native readback retains spatial residuals at alpha_mass=0, actual saved
callback-stage defects, total continuum error, mesh-initial drift and same-time
differences separately. Both H sides and vertical-face intersections are explicit.
Signed regional weak loads use actual represented fields and callback history;
independent WJ algebra/split quadrature checks source/diffusion and two shared
interface test functions. Their exterior reactions vanish because the tests
vanish at prescribed exterior boundaries. Epsilon uses an outer bubble test;
no additional epsilon flux condition is imposed.

Startup H derivative conditions in RadialProject partly impose initial traction
and k-flux agreement; those initialization identities are labelled separately
from independent evolved load/trace checks. Pressure is an algebraic multiplier:
its discrete consistency adjustment is measured against the unchanged saved
initial pressure, with common-time comparisons in one overall gauge.

## Horizontal result and independent native witness

Both policies confirm startup sizes639,1109,2049,1257,2183,2493,4331 for
L0/Nx1,Nx2,Nx4; L1/Nx1,Nx2; L2/Nx1,Nx2. The4331-unknown cases receive reference
audits only. The3000-unknown guard is explicitly retained in WJ PIMPLE as well
as the older dense solve path.

L1/Nx1 reaches T3=.004096s at dt=3.2e-5,128 steps. Its normalized sampled strong
divergence is about9.63766e-4. L1/Nx2 crosses the fixed1e-3 engineering budget
at step32, t=.001024s:2.12883124e-3 for EP and2.12883176e-3 for budget.
Further Nx2 continuation blocks stop, retaining the accepted state and the
unchanged T3 target. All other accuracy budgets, native hard predicates and
inactive-projection checks pass at that accepted crossing.

The maximum is at x=1,z=h=.1, inner geometry/element29. The central-strip maximum
is also above budget, about1.05847e-3 at x=.5,z=.1 on the explicit right-hand
vertical-face trace. Thus an unchanged global maximum would not characterize
interior sensitivity.

Independent native element-owned BwdTrans/PhysDeriv/StdPhysEvaluate reproduces
the EP raw divergence .001064415620005743 as .0010644156201074153 s^-1;
for budget .0010644158820045482 becomes .001064415882070311 s^-1.
The discrepancy is about1e-13. Both H sides and all face/corner gradient checks
retain their geometry ownership. This is an accuracy obstruction in the
represented field, rather than a radial/x-side observation error or a hard
constitutive/projection failure. Independent evolved loads and regional weak
moments pass at the obstruction. No budget, boundary condition or physics was
changed to continue it.

The horizontal component is **OBSTRUCTED**. The longer fine-normal evidence is
therefore conditional finite-strip evidence and cannot close horizontal
qualification. Evolving L2/Nx2 would require4331 unknowns and is outside this
bounded dense backend.

## Frozen longer design and resource policy

The longer pair uses2493 unknowns, T4=.008192s, T4*u*/H=.004096,
dt=3.2e-5/1.6e-5 and256/512 steps, each with a fresh compatible CP0. Initial
primary/gauge and corresponding physical history values must agree exactly;
dt-specific time labels/signatures remain native. This is timestep sensitivity
against budgets, with no new temporal-order claim.

The measured four-step process costs are24.834s EP and33.903s budget. Before the
pair, frozen estimates were6200.50/7994.78 remaining seconds with1531.25/1918.39
already charged. Both fit the10800-second per-policy cap. Limits remain1500
scheduled attempts per policy, including failed/control/restart work. Numerical
libraries use one thread. Current capacity snapshots precede builds, tests and
the two independent serial policy processes; the total remains below72.

Linux builds use Release GNU C++17, `-ffp-contract=off -fno-fast-math`,
RANS_BUILD_SOLVER=ON and compatible installed dependencies read-only. Initial
compilation used one worker, then four after capacity checks. CTest uses -j1.
The producing native build is0fcf69d, byte-identical in native inputs through
5628d93. The observer build isf93c1c9; subsequent Python-only revisions retain
its native inputs. The exact source-difference check shows only an added
nonvirtual native evaluator and read-only WJMeasurements diagnostics, with
unchanged existing global-field implementations and data-layout declarations.
Both runtimes and their resolved library identities are recorded separately.

The successful build receipts are the clean producing rebuild862.660s, its
final native update92.408s and the fresh observer build875.313s. The executed
commands are `cmake --build "$r/build" --clean-first --parallel 4`,
`cmake --build "$r/build" --parallel 4` and
`cmake --build "$r/readback-runtime/build" --parallel 4`, respectively.
Configuration uses `cmake -S <source>/solvers/IncNavierStokesSolver/Tests/RANS
-B <build> -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release
-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
"-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math"`.

Production changes are concentrated in RANSWallEvolutionDiagnostics,
RANSWallGlobalField, the native resolution probe, assembly-check exports and
the WJ PIMPLE size guard. The shared governing/model equations are preserved.
Python changes repair gates/immutable reuse/accounting, generalize x knots,
measure common physical points and provide bounded native/restart drivers.

## Terminal results

Both policies complete 256/512 steps on L2/Nx1 (2493 unknowns), with native end
times `0.008191999999999991` / `0.008192000000000095` s. Fresh CP0 primary/gauge
and physical history hashes agree between dt choices; native dt signatures and
time labels remain distinct. All observed baseline states have inactive
MassProjection, identical raw/corrected k bytes and zero added k storage load.

The following fine-dt values are normalized by the fixed amplitudes and each
field's own measure. The L2 error and drift budget is 1e-4. Full field/scope
L2 and sampled-maximum tables retain initial representation separately.

| Policy | Field | E2 at T4 | D2 at T4 | Same-T coarse/fine L2 difference |
|---|---|---:|---:|---:|
| EP | velocity | 2.639030e-6 | 1.921809e-6 | 2.144440e-10 |
| EP | k | 8.992589e-8 | 8.992589e-8 | 3.799694e-11 |
| EP | entire outer epsilon | 3.097425e-8 | 2.566816e-8 | 2.019228e-12 |
| EP | reduced pressure | 1.120842e-7 | 1.120842e-7 | 5.973275e-10 |
| EP | mechanical pressure | 2.325414e-7 | 2.325414e-7 | 6.026505e-10 |
| Budget | velocity | 2.639030e-6 | 1.921810e-6 | 2.144412e-10 |
| Budget | k | 4.133290e-7 | 3.026580e-7 | 4.153105e-11 |
| Budget | entire outer epsilon | 3.097523e-8 | 2.567138e-8 | 1.989394e-12 |
| Budget | reduced pressure | 1.120577e-7 | 1.120577e-7 | 5.973003e-10 |
| Budget | mechanical pressure | 9.400936e-7 | 6.926920e-7 | 6.038125e-10 |

The final fine-dt interface quantities below use their declared engineering
normalizations; k/epsilon value-trace defects are absolute physical values.
Pressure-offset departure is measured against -K0/Ne and is not itself the
general-state mechanical coupling residual.

| Quantity | EP | Budget | Budget limit |
|---|---:|---:|---:|
| full traction jump / u*^2 | 2.936982e-5 | 2.936981e-5 | 1e-3 |
| k flux jump / (K0 u*) | 2.424842e-7 | 1.543829e-7 | 1e-4 |
| pressure-offset departure / u*^2 | 1.832727e-6 | 1.832400e-6 | 1e-4 |
| sampled strong divergence H/u* | 8.880164e-5 | 8.880174e-5 | 1e-3 |
| shared k value defect | 0 | 0 | exact trace |
| dependent epsilon trace defect | 2.78e-16 | 1.67e-16 | native trace predicate |

Endpoint K and cutoff remain responsive to evolved traces. At fine T4, EP K
is in [0.8466950351491502,0.8466950354719460] and budget K in
[0.8466950347359299,0.8466950351246831]. Their respective cutoff intervals are
[0.015005743260185285,0.015005743260356210] and
[0.015005743259966491,0.015005743260172332] m. The independently assembled
evolved source/diffusion loads and regional weak moments pass their original
tolerances. Horizontal integration order4->6 changes initial field L2 errors
by at most about1.8e-18 for EP; its declared absolute floor is1e-10.

Every required common-time comparison contributes to final qualification,
including startup and all eight cadence endpoints. Early/middle/end frozen
arithmetic packets use the declared100x reserve and remain finite estimates,
not a nonlinear roundoff theorem. Differences below their estimate receive no
rate. This pair is timestep sensitivity against fixed budgets, with no new
temporal-order claim. The constrained-boundary L2 horizontal comparison also
remains UNRESOLVED under its retained conservative envelope (velocity uncertainty
about1.59e-5 against the1e-5 sensitivity limit), even though prescribed
differential trace differences are zero. It is separate from the measured
divergence obstruction and the unreached Nx2 T3 target.

Stationarity fails independently. Fine velocity window V values are
4.97965e-5,3.56887e-5,2.68921e-5 (both policies to shown precision), exceeding
1e-5 in all three windows. Final normalized spatial u residuals are
2.36518e-5; the actual completed-stage u defects are about1.790e-8.
The small stage defect does not establish a steady state. Fine k/epsilon V
values, all endpoint residual blocks and both continuity regions are retained
in the window/residual tables. No startup steady-residual target was imposed.

Both late [224,256] segments match extra16+16 restarts exactly in all33
retained primary states/times and accepted audit packets, all32 incoming and
callback-history hashes, every raw/corrected k snapshot, final gauge/GLM slots
and final physical samples. Frozen readback preserves saved bytes. Deliberate
step225 rejection restores accepted step224/history/time and keeps its failed
candidate separately. All fully scheduled attempts, including this rejection,
retain one U, one pressure/gauge, one k, one epsilon solve and one logical bound;
additional governing/restoration solves are zero.

Final charged work and publication receipts are recorded in the outer
`publication.json`; the inner receipt preserves its prior immutable snapshot.
Each policy has965 native study attempts, including the4-step cost measurement
and33 extra transaction attempts. Repeated WJ/shared-legacy controls add28
native attempts;28 conservative ordinary legacy-control units are separately
charged, giving1021/1500 cap units per policy. Ordinary control units are not
given WJ one-sweep solve-count claims. Final numerical-study charges, including
117.842s for figure/index/publication work, are6681.529s EP (1.856h) and
9044.226s budget (2.512h), each below10800s. Both ledgers have no active
reservation/session. The minimum declared bound arithmetic reserve on native
study attempts is7.7006464847e-13; measured correction norms are exactly zero.

Across all declared longer observations, the largest normalized full-domain
L2 error/drift is2.639030e-6/1.921810e-6; the corresponding largest sampled
maximum error/drift is1.758951e-5/2.188232e-5, below1e-3. All378 longer
sensitivity leaves per policy (nine common times x three scopes x seven fields
x two norms) pass. The maximum (difference+uncertainty)/sensitivity-limit is
0.001192 EP and0.001109 budget. No order is assigned to below-estimate values.
Independent load, shared weak-load and trace-moment maximum discrepancies are
1.97e-16,3.74e-14 and1.10e-15, respectively, against unchanged declared tolerances.

The overall components for each policy are: mesh inventory PASS; resource
design/actual resource PASS; fresh startup PASS; longer PASS; transactions PASS;
horizontal OBSTRUCTED. Stationarity remains a separate NOT_REACHED outcome.

## Reproduction and immutable evidence bridge

The tested publication source is
`3f063fc48737e5a7e61a0d17100c85a07e7abb6e`, tree
`996ab9fd6d87d34e231a5bf0607de2a15cd893e3`. Its source status is clean.
Native production/observer build identities remain0fcf69d/f93c1c9 as described
above; only Python qualification/publication changed afterward. The separate
publication worktree preserves closed observer source at6666d53. Every relevant
script hash is bound independently of the changing overall Git revision.

The compact [manifest](wj4-linux/manifest.json), [inner receipt](wj4-linux/receipt.json)
and [final publication seal](wj4-linux/publication.json) bind46 transferred
files (5,165,604 bytes), with exact Linux-to-Mac byte equality checked before
adding the local attributes file. Heavy data remain under the experiment root.
The immutable heavy index contains10381 artifacts, including native fields,
expanded options, checkpoints, native/oracle JSON and logs. Its SHA-256 is
`3038587f3b65ac8ece10704e6b215ded4cec96ca58b9f80e8ccab85bf61ca0c4`.
The source-content index SHA-256 is
`da54be300c6eee935e7ea51bcf8f7f20ccf7807b85308d671a4cbe1707c03eff`.
Seventeen executed-test executable/interpreter identities per native build
retain their role-specific resolved library closures. Producer and accepted
readback closures are distinct records.

The primary invocations below expand against the recorded isolated paths;
complete argv/exit/start/time/log hashes are in the copied command receipts.

```bash
r='/media/jay/新加卷/nektar-v2026-runs/wj4-20261110'
producer="$r/source/solvers/IncNavierStokesSolver/Tests/RANS"
observer="$r/readback-runtime/source/solvers/IncNavierStokesSolver/Tests/RANS"
publication="$r/publication-source/solvers/IncNavierStokesSolver/Tests/RANS"
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export WJ4_PACKAGE_ROOT="$r" TMPDIR="$r/tmp"
export LD_LIBRARY_PATH="$r/build/geometry-lib:${LD_LIBRARY_PATH:-}"
# Both original policies used this invocation, replacing EP with Budget.
python3 -B "$producer/run_wj4_native.py" --solver "$r/build/IncNavierStokesSolver" \
  --probe "$r/build/RANSWallResolutionProbe" --policy EP --mode all --output "$r/study"
# EP observer code stayed pinned; Budget used the repaired publication driver.
python3 -B "$observer/run_wj4_native.py" --solver "$r/build/IncNavierStokesSolver" \
  --probe "$r/readback-runtime/build/RANSWallResolutionProbe" --policy EP \
  --mode reprocess --native-study "$r/study" --output "$r/study-readback"
python3 -B "$publication/run_wj4_native.py" --solver "$r/build/IncNavierStokesSolver" \
  --probe "$r/readback-runtime/build/RANSWallResolutionProbe" --policy Budget \
  --mode reprocess --native-study "$r/study" --output "$r/study-readback"
python3 -B "$publication/requalify_wj4.py" --study "$r/study-readback" \
  --policy Budget --output "$r/terminal-gates-Budget.json"
python3 -B "$publication/finalize_wj4.py" --root "$r" \
  --observer-root "$r/readback-runtime" --study "$r/study-readback" --output "$r/compact"
```

The original policies exit0 with honest OBSTRUCTED qualification. Independent
EP/Budget reprocessing exits0 in371.306/352.241s, with zero governing attempts
and unchanged checkpoint/report bytes. Complete-gate certificates add no probe
or governing launch. The final publisher exits0 in110.500s; its log SHA-256 is
`72e9be4b04234f4d50fee70777abb6455762d2fc57c38a57ad7e885deed5347d`.
The inner receipt SHA-256 is
`047e1819824b3e7053fefdc2d8aaa58e61c7c8b83afdb005eb550a6160909245`.

The tables separate [field error/drift](wj4-linux/fields.csv),
[interfaces](wj4-linux/interfaces.csv), [window rates](wj4-linux/windows.csv),
[spatial/stage residuals](wj4-linux/residuals.csv),
[all sensitivity margins](wj4-linux/sensitivity-gates.csv),
[native solve/correction accounting](wj4-linux/run-accounting.csv),
[independent loads](wj4-linux/independent-loads.csv),
[signed weak loads](wj4-linux/weak-loads.csv),
[regional residual norms](wj4-linux/regional-residuals.csv) and
[trace moments](wj4-linux/trace-moments.csv). The
[accepted-stage lineage](wj4-linux/accepted-stage-lineage.json) explicitly marks
startup/frozen-readback callbacks inapplicable and rejected callbacks unaccepted.
The [audit witnesses](wj4-linux/audit-witnesses.json) retain both principal forms
and independently replayed gradient provenance. Exportable PNG/PDF/SVG figures
are in [wj4-linux/figures](wj4-linux/figures).

## Validation and retained development failures

Both producing and independently rebuilt observer runtimes passed all 32
selected tests from 56 registrations: 24 kernel/legacy controls, 3 native/legacy
transactions, 2 read-only-stage tests, 1 metadata test and 2 composite native
observation tests. This includes all 29 WJ-3 selections after shared C++ changes.
The final registered metadata supplement passes15 fixtures, including a failing
early sensitivity check with a passing terminal check, missing nonterminal
coverage, incomplete fields and interrupted-log preservation. Its separate
configuration registers56 tests and selects exactly1 metadata test; it reuses
no native trajectory. The immutable regression-verification packet records every name,
count, command, exit, elapsed time and CTest-log hash.

The executed wrapper is `wj4_validate.py --root <runtime-root> --suite <suite>
--label <label>`, with suites `controls`, `native`, `readonly`, `metadata` and
`composite`. Producing labels are `qualification`; rebuilt-observer labels are
`observer2`, with `alltimes15` in the qualification configuration for the final
metadata supplement. Each wrapper inventories
the exact CTest regex before executing `ctest --test-dir <runtime-root>/build
-j 1 --output-on-failure -R <recorded-regex>`. The compact command receipts retain
the fully expanded absolute argv, not only this invocation template.

The exact reference tangent certificate and both adverse WJ constitutive
witnesses remain in the compiled independent oracle selection. The smooth-budget
reference_nonpositive witness and reason-specific cutoff/scalar/configuration/
signature failures remain in the native selections.

The q=3.5 witness retains N=4.130462147295278 and principal form/nu_eff
-0.00633378103539413 with positive covariance. The commuting 3D witness retains
N=3, R/k=diag(-2/15,-2/15,34/15) and P/epsilon=48/5. These tests pass by
reproducing the adverse values. Finite native covariance and tangent audits are
reported separately; neither is a global realizability or stability proof.

Development receipts retain interrupted preliminary builds, an optional unset
library-path launch error before any governing attempt, and global-ledger label
collisions before any governing attempt. The latter was fixed by namespacing
attempt labels with the immutable experiment path. Failed receipts were not
rewritten or accounting reset. These are implementation/setup failures, distinct
from the independently reproduced horizontal accuracy obstruction.

An accounting reader initially assumed that the standalone legacy stage XML
contained NumSteps. Its failed preflight launched no governing attempt; a
reason-specific fixture now recognizes that declared stage route, while ordinary
missing step targets remain invalid. Append-only supplements retain both original
control receipts and explicitly charge their omitted conservative budget units.
The final all-observation aggregation repair uses a separate publication source
worktree and metadata-only certificates; closed measurement source/files stay
unchanged. No saved signature or history was edited to make a check succeed.

The preliminary build hand-off briefly overlapped two owned builds after a
process-name parser mishandled spaces. Both were stopped (retained -15 receipts),
then a clean four-worker build completed. The concurrent snapshot counted65
threads, below72; the largest retained allocation total is70. The parser was
repaired before subsequent allocation decisions. Build time is recorded
separately from the numerical-study budget.

This study does not establish exact finite-space equilibrium, stationarity,
long-time preservation, global positivity/realizability/stability, native3D/MPI
or terrain validation. The next spatial step must resolve the measured horizontal
accuracy obstruction, with a suitably resolved/scalable backend; separately
qualified homogeneous boundaries remain a different boundary problem.
The bounded next step is to diagnose the velocity/pressure representation at
z=h using the retained face witnesses, and enable a matched L2/Nx2 comparison
(4331 unknowns) through an appropriate scalable frozen-system backend or a
separately justified space choice. Retain the frozen horizontal budgets and
reference-Dirichlet boundary problem for that comparison before extending
duration. The present conditional results do not close horizontal qualification.
