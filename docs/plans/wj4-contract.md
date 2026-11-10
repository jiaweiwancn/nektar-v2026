# WJ-4 qualification contract

The reviewed parent is `11733a64fdf9452adc14c45594820c3da2513fd7` on
`codex/wj-earsm-awf`. Earlier WJ receipts remain immutable. The main and appendix
SHA-256 identities are respectively
`1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124` and
`0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c`.
Original and worktree copies agree. Main eq:integrals-linear retains the
literal `quad` typography erratum for `\quad`; its logarithmic coefficient
is `(d+a*h)/a^2`.

The equations and unsupported routes in WJ-1/2/3 remain unchanged. This is the
finite reference-Dirichlet strip **[0,1] x [0,2]**: reference differential-field
traces are prescribed at **both lateral boundaries and the top**, ground
velocity is zero, and modeled ground k flux is zero. It is not a homogeneous,
periodic or outflow boundary problem. Full WJ covariance, DevStrainRawWorkV1,
physical epsilon, Praw=-R:Graw, independent Cd scalar diffusion, positive wall
nu, Linear full-strain inner stress and uncapped sources remain selected.
EP uses chi=1 and zm=H/2; budget uses zK=zm=H, no physical chi, and transported
nonconstant inner k. The exact constants from WJ-1 remain unchanged.

IMEX1/PimpleSingleSweep performs flow/pressure, k, epsilon once; eta=0, zero
body/scalar forces and gravity, no condensation/search, reuse/verification Off.
WJNativeV1 MassProjection and StopOnBound=On require inactive correction.
Reduced pressure is pi=p/rho+2k/3, with one global gauge, distinct pressure and
gradient traces, shared velocity/k and one dependent outer epsilon trace.

## Qualification and immutable execution

Execution success and component qualification are separate. Required components
use PASS, FAIL, OBSTRUCTED, UNRESOLVED, RESOURCE_LIMIT or NOT_RUN. A missing cross
mesh or transaction, nonfinite evidence, unresolved meaningful temporal
comparison or nondecreasing spatial sequence cannot yield overall PASS.
Requalification of WJ-3 is a new derivative artifact; its receipts are not edited.
Its first-order gate uses resolved differences and p in [0.5,1.5]; spatial u,k,
full epsilon residuals must decrease above the declared 1e-11 absolute floor.
An explicitly identified constant field at arithmetic floor receives no order.

Write-once design and requested-intent packets bind argv, exact/semantic XML,
mesh, dt, target, source code, runtime/dependencies and checkpoint/report/log
identities. Changed intent requires a new root. Measurement-only reprocessing
owns a separate code/output lineage. An unreconciled/crashed attempt is never
silently rerun or overwritten. A persistent package ledger reserves attempts
and wall time before launch; unknown attempts retain their reservation.
WJ-4 caps are 1500 scheduled attempts and 10800 numerical-study wall seconds
per policy, including failures, controls, observations, SVD and I/O. WJ-3's
original 500-attempt evidence is unchanged. Resource exhaustion is a component
RESOURCE_LIMIT. Practical cost estimates are frozen before the longer pair.

## Physical comparisons and predeclared budgets

Cases accept explicit x knots while preserving Nx1 defaults. The principal
horizontal comparison is L1/Nx1 versus L1/Nx2, startup and T3=.004096 s,
dt=3.2e-5, 128 steps. The longer pair is L2/Nx1, T4=.008192 s,
T4*u*/H=.004096, dt=3.2e-5/1.6e-5, 256/512 steps from fresh compatible CP0.
The dense 3000-unknown safeguard remains; L2/Nx2 is not evolved.

Common physical composite quadrature splits the union of horizontal element
boundaries, x=.25,.75, normal elements and physical breakpoints. An independent
horizontal-order increase checks integration. Both H sides and explicit native
geometry ownership at vertical faces/corners retain distinct gradients. The
radial side flag never selects an x side. One overall pressure gauge is applied
before full-domain, boundary or central-strip extraction. Subdomains use their
own measures. Spatial residuals, actual accepted-stage defects, continuum errors,
mesh-initial drift and same-time differences remain separate. The entire outer
epsilon region is included. Measurements hash state/history and saved bytes,
use the actual producing report/callback packet, and perform no governing solve.

Amplitudes are aU=.5, ak=K0=.8466950351321795,
aepsilon=.125/kappa, ap=.25. E2 and D2 divide physical L2 by amplitude times
sqrt(field volume); sampled maxima divide by amplitude. Engineering budgets:

| Quantity | Bound |
|---|---:|
| normalized L2 continuum error and initial drift | 1e-4 |
| normalized sampled maximum error and drift | 1e-3 |
| full traction jump / u*^2 | 1e-3 |
| k flux jump / (K0 u*) | 1e-4 |
| reference pressure-offset departure / u*^2 | 1e-4 |
| sampled strong divergence H/u* | 1e-3 |

These are pilot decisions, not preservation theorems. Pressure-offset departure
from -K0/Ne is distinct from general-state coupling, assessed by full traction
and signed regional weak loads. Same-T timestep/horizontal differences plus
diagnostic uncertainty must fit 10% of corresponding field budgets. A difference
below its uncertainty receives no rate. Passive early/middle/end arithmetic
packets give finite frozen-system estimate envelopes, not a nonlinear roundoff
theorem; spatial residuals/stage defects are excluded from the noise estimate.

Accuracy is checked at most every32 coarse/64 fine steps; a crossing stops
further blocks but retains the accepted checkpoint and frozen target. Native
hard/material-correction exits stop before commit and retain failed candidates.
No unmeasured time is certified. If horizontal accuracy is obstructed, useful
longer finite-strip evidence remains conditional and cannot close that gate.

Stationarity is separate: endpoints160,192,224,256 (fine320,384,448,512)
define three equal final windows. For U,k,epsilon,
V=(H/u*)||delta f||/(delta T*amplitude*sqrt(volume)) must be <=1e-5 in every
window, and normalized spatial u,v,k,full-epsilon residuals <=1e-5 at all window
endpoints, with hard/trace/no-correction conditions. Otherwise NOT_REACHED.
The startup projection residual is not subjected to a steady-state target.

A retained late32-step coarse segment is compared with an additional16+16
restart, plus immutable readback and deliberate post-solve rejection. Exact
same-runtime primary/gauge/time/history, accepted audit/sample equality and
unchanged saved bytes are required. Solve counts remain honest on partial exits.

Linux is authoritative, in a new isolated root using source-matched consumers,
RANS_BUILD_SOLVER=ON, -ffp-contract=off -fno-fast-math. Numerical threads=1,
current capacity checks and <=72 total server workers apply. Shared C++ changes
require all29 WJ-3 selected controls plus focused WJ-4 tests. Heavy data remain
outside Git, indexed by immutable paths/hashes. Conclusions remain finite-strip,
short-duration, serial flat2D observations; no global positivity, realizability,
stability, long-time, native3D/MPI or terrain claim follows.

## Independent gradient replay and observer lineage

A horizontal accuracy crossing is replayed through native element-owned
`BwdTrans`, `PhysDeriv` and `StdPhysEvaluate`. The fast gradient and the
independent gradient retain the same explicit geometry/element at faces and
H/corner intersections. Componentwise differences must fit
`2e-11*(1+max(abs(fast),abs(native)))`; no derivative is replaced in the governing
operator. The complete maximum-divergence witness and comparison defects are
retained. This distinguishes a represented-field accuracy obstruction from a
row/provenance error.

A later observer build remains separate from the producing runtime. All matching
native consumers are rebuilt and the29 prior selections rerun. Its new method
is nonvirtual and adds no field/layout; existing governing methods and model
constants are unchanged. Strict serialized-checkpoint/signature validation is
still required. `--mode reprocess --native-study <producer-study>` verifies the
original argv/XML/input/output/code/runtime receipts, launches no solver and
creates a separate measurement/code lineage. Original summaries are archived
byte-for-byte. Actual accepted callback packets always come from the producing
report. SVD caches have their own immutable input/output receipts. All extra
observer/test work remains charged to the same package-wide ledger.

Final publication has two immutable layers to avoid a hash/accounting cycle.
The inner manifest and receipt bind the measurement-time snapshots. After the
figure/publication commands finish and their reserved wall time is reconciled,
an outer `publication.json` binds the inner receipt, final ledgers and complete
command receipts. Interrupted command logs require a new label. Compact run
accounting includes every scheduled/partial solve count, projection reserve and
raw/corrected hash equality; accepted-stage lineage binds callback history to
its producing report, checkpoint, signature and physical samples.
