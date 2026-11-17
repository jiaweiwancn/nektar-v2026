# Current quadrature data inverse-cache admission

After source14's passed terminal-output foundation, its retained operator-cache
path still computes every frozen velocity tangent and both scalar inverse
matrices at every Newton base. Only factorization is reused. Admit a separate
opt-in `WJStepPreconditioner=CurrentDataDriftV1` to avoid preparing those inverse
matrices on reused bases. Preserve the fresh control and OperatorDriftV1 paths.
The active approved-v2 equations, reference, inlet, boundaries and guards are
unchanged. V4 remains inactive; no full Task3 run is admitted.

At every Newton base use already current data on the fixed native static
quadrature inventory, retaining each original owner. Let K0=Kin and use the
fixed positive normalization E0(r)=Ew(r;Kin) below H and Ain/r outside. Define
`Dk0=nu+Cd*K0^2/(sigma_k*E0)`, `De0=nu+Cd*K0^2/(sigma_e*E0)` and
`reaction0=1.82*E0/K0`. At each point form the dimensionless ordered tuple

```
(k/K0, E/E0, aux/Dk0,
 Gxx*K0/E0, Gxz*K0/E0, Gzx*K0/E0, Gzz*K0/E0,
 bulk/Dk0, shear/Dk0, Dk/Dk0,
 [De/De0, reactionE/reaction0] for outer points only).
```

Dk/De/reaction and bulk/shear are the actual current coefficients, sampled
at the fixed points before moving-cutoff scalar quadrature is introduced.
The tuple includes every variable input to the frozen-scalar velocity tangent
and every scalar inverse coefficient. Dk0 is a numerical unit for diffusivity;
it does not identify Cd with a momentum coefficient. All scales are strictly
positive fixed reference quantities; use no norm or denominator floor.
Positive k/E entries make every valid full tuple's norm nonzero.

Use `rho=||tuple_current-tuple_prepared||2/||tuple_prepared||2`. This is a
preconditioner refresh heuristic, not an operator-error or tangent bound.
It is distinct from OperatorDriftV1's actual-matrix Frobenius norm and must
be identified separately in native diagnostics. Retain the .02 drift threshold,
eight Newton-base uses, cold/actual-alpha/drift/age refresh order, all-three
factor refresh, saved matching epsilon full operator, current Gamma derivative,
typed same-dt cold retry and single cooperative physical-step deadline.
Every aborted path invalidates both the tuple and factors. Neither tuple nor
inverse factors enters current residuals, FD Jacobian actions, forcing or
accepted-state guards, and neither is serialized in a recovery file.

Also repair source14's zero-NumSteps terminal serialization timing: create
its terminal frame inside the final-export timer. Advancing output scheduling
is unchanged. Test a physical zero-step cold import with terminal output
enabled, exact binary histories and explicit imported-frame provenance;
require frame/index seconds to be contained in final export. Default native
zero residual/guard records are not accepted-step evidence.

Freeze a new source15 capsule after c88e26ac7f9d55c18b4c5139949de628015cc182.
Build the standard main with four jobs. Compare PolynomialScalarsV1/body1
four-step dt=.02 and Historical eight-step dt=.0025 fresh factors against the
new cache, with max uses3 and the same <=1e-8 coefficient/owner/history/Gamma/
cutoff/Ew/action limits. Independently check both current MMS endpoints.
A threshold-zero new-cache run must refresh on every nonzero measured change;
an exactly unchanged tuple may reuse. Independently check its current equations.

On retained700 cells use the existing dt=.005/NL3e-12/linear1e-12/no-Anderson
12-step physical diagnostic. Run fresh, OperatorDriftV1 and CurrentDataDriftV1
with the same eight-use/.02 settings; compare each cached path to fresh using
the fixed<=1e-8 gates. For the new policy qualify prefix3+cold9, retained
recovery9+cold3, terminal-enabled recovery12+zero steps and one retry at step4.
Recovery/science stride3 produces four writes, actual deletion3/6 and retained
9/12. Bind four advancing physical independent current-equation audits and
the physical history checker. Keep complete per-base events, actual matrix/
factor work, all timings and failures. Three short post-startup blocks remain
diagnostics, not representative long-state production forecasts.

Admit materialization60s, four-job build600s, fixtures30s, MMS natives90s,
MMS current audits90s, MMS comparisons60s, zero-drift event audit20s,
physical controls600s, physical comparisons120s, physical history natives580s,
four current physical audits600s, complete history/zero-export audit120s and
evidence90s: <=3060 aggregate additional Linux operation seconds. Allow no
reset of the same43200s campaign account, last published OPEN10630.33897288586s,
and retain3600s protected for evidence. Use one science process/rank, all six
library thread counts1, <=72 actual shared-host workers, science RSS8GiB/
build16GiB, the existing source/build8GiB and science8GiB allocations and
15GiB free reserve. Resolve unknown process outcomes before retries.

All scientific/reference and physical pilot/accuracy gates remain open. A
numerical/cache PASS cannot qualify the failed v2 same-inlet reference, select
a production timestep, finish monitoring/storage or imply one-day affordability.

The first source15 MMS current-oracle command omitted its nested direct-check
time limit; the budget guard rejects its default before equation assessment.
Preserve that charged failure. The corrected command explicitly uses25s total
and20s direct within the admitted30s operation, with the unchanged oracle gates.

All native physical controls, field/history comparisons and four independent
current-equation audits complete. The first complete history assessment then
misclassifies the new scientific index descriptor as a numbered snapshot.
Preserve its27.677273s failed receipt and report. Admit a test-only immutable
validation overlay<=20s and a new read-only complete assessment<=80s, under
the same3060s phase forecast and43200s campaign ceiling. At this admission
the campaign has charged11672.81236364384s, with no open operation. Only the
descriptor classification/binding is repaired; all four required numbered
snapshots and all numerical thresholds remain unchanged. No solver is rerun.

The repaired short history gate passes. Its prefix3/cold4 and recovery9/cold10
endpoints are exactly equal to uninterrupted history. The observed uninterrupted
cache refreshes at4 and10, so those cold imports coincide with normal refresh
bases. Admit one additional offset restart: prefix4, then cold8 steps from its
complete checkpoint to12. At uninterrupted step5 the first base reuses age4
factors; the cold step5 must start unprepared and actually rebuild all three.
Require this distinct cache lineage and the same<=1e-8 full history/owner/
dependent-state comparisons, plus an independent current-equation endpoint
report. No exact endpoint equality is required across the different inverses.

Use a new immutable experiment namespace and a separate read-only checker.
Admit fixture20s, prefix80s, cold150s, current audit150s, history assessment120s
and addon evidence60s, <=580 additional seconds. Actual completed source15
charges plus this remaining forecast remain below the original3060s phase
admission. The43200s campaign account,3600s evidence protection and all worker/
memory/storage limits are unchanged. This addresses an observed coverage gap;
do not repeat the already passed numerical or native matrices.
