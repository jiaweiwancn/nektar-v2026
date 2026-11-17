# Unchanged-v2 Newton inverse-operator reuse admission

Source-v12's final 700-cell startup step spends 4.2660 of 10.9844 native
seconds in factor rebuilds. The current equations and final checks already
pass the independent comparator. This prospective change reduces repeated
inverse preparation; the approved v2 physical model remains unchanged and
the v4 scientific decision remains pending.

## Exact numerical policy

Keep `WJImplicitSolver=CoupledJFNKV1` and default
`WJStepPreconditioner=RebuildEveryNewtonV1`. Add the opt-in
`OperatorDriftV1`. At every Newton base, assemble the current frozen-scalar
velocity tangent K_U and complete current scalar diffusion/reaction matrices
K_k/K_e on their actual signed maps. These are inverse approximations; the
physical residual, FD Jacobian action and final full principal certificate
remain current and independent of cached matrices.

For saved inverse-operator matrices K_i^s define
`rho=max_i(||K_i-K_i^s||_F/||K_i^s||_F)`. Require nonzero finite saved norms;
do not floor a denominator. This is a numerical refresh heuristic, not an
error or tangent bound. Fixed mass/divergence/pressure blocks do not dilute
this coefficient-matrix comparison. Refresh all three factors when cold,
actual BE/BDF2 alpha changes, rho exceeds .02, or the factors have already
served eight Newton bases. Age counts bases, starting at one after a refresh;
it does not count Krylov applications. Compare current matrices even on reuse.

Use the saved epsilon diffusion/reaction-plus-actual-mass matrix for the
Gamma cross action whenever its inverse is reused. Gamma's derivative is
evaluated at the current Newton base and remains fixed within that linear
solve. Cache membership belongs to one immutable geometry/map/reference/source
instance. Numerical policy, threshold and maximum uses bind the history
signature. Derived factors are never serialized.
The normal factory constructs a new Coupling instance before a cold process
restart. An external same-instance history/field transfer must explicitly call
ResetPreconditioners when a cold inverse reconstruction is intended; such a
transfer is not the cold-process evidence in this matrix.

Typed linear failure clears the cache and retries from accepted histories at
the same dt. The second attempt refreshes at every Newton base. The explicit
retry injection first performs a real triple rebuild. Every aborted path
invalidates cached factors. Actual current residual and full principal checks
still precede the only history commit. Preserve the single shared cooperative
step deadline and all source-v12 failure semantics.

Per-step counters and exclusive refresh reasons are always recorded. With
`WJNonlinearTrace=True`, record bounded per-base events containing attempt,
base, alpha, operator drift, age before/use, refresh reason and refresh flag.
Production may disable those detailed events. Their timings and I/O remain
charged when enabled; they do not replace equation checks.

## Bounded qualification matrix

Freeze source-v13 from source-v12 plus the complete current Mac overlay; retain
all earlier code, raw data, failures and accounts. Use a full standard-main
build with four compiler workers. Compare fresh/control and cached Newton
MMS runs with the existing fixed <=1e-8 field/history/trace/cutoff/cache gates.
Use PolynomialScalarsV1/body=1 at dt=.02 for four steps and Historical at
dt=.0025 for eight steps. Set max uses three for those MMS pairs to exercise
age refresh. Also run an operator-drift threshold-zero control, where every
nonzero matrix change must refresh; exact equality may reuse. Independently
check changed-equation/current-state audits at the same fixed tolerances.

On the original retained 700 cells, compare fresh and cached Newton for 12
unforced startup steps, dt=.005, nonlinear 3e-12, linear 1e-12, no Anderson.
Both use threshold .02/max uses eight. This is T=.06, a short numerical/cost
diagnostic; it is not the >=26.1335 physical principal pilot or an accurate
selected production recipe. Report the initial three-step block and three
subsequent three-step blocks, actual variation, all equation/guard checks and
cache work. No one-day or long-duration claim follows from those blocks.

For the cached physical case, test a three-step prefix plus cold nine-step
continuation to step12, recovery9 plus three steps to12, a zero-step cold
import of recovery12, and one retry at step4. Scheduled recovery/science
stride3 produces writes3/6/9/12, actual removal3/6, retained9/12 and retained
scientific outputs3/6/9/12. Verify each retained recovery through its declared
cold process. Independently parse complete binary/source/checksum/history
identities and compare fields, dependent wall state and explicit histories
with the fixed <=1e-8 scales. Zero-step import is not an advancing residual
test. Bind independent current-equation reports for each advancing endpoint.
The new physical history checker must retain pressure level and owners.

Admit source materialization <=60 s, build <=600 s, fixtures <=30 s,
MMS natives <=240 s, MMS audits/comparisons/event control <=360 s, physical native matrix
<=1230 s, physical algorithm comparison <=180 s, physical equation audits
<=900 s, physical history assessment <=300 s and evidence <=90 s:
<=3990 aggregate additional seconds. All failures/overlaps are charged in the
same 43200 s account, last completed 8715.746188972902 s, with3600 protected
for evidence. No allowance is reset. Use one science rank, all six library
thread counts one, <=72 actual shared-host workers, science RSS8GiB/build16GiB
and the existing source/science storage and free-space limits. No full Task3
run is admitted. Final H27-C/D qualification remains open even if this
numerical foundation passes.
