# M8-R1 — First failure and joint constitutive feasibility

## Scope and disposition

Bounded diagnosis on `codex/cls-epawf-revised`, authorized by the M8-R1 request.
The initial Mac HEAD exactly matched the reviewed baseline
`ab5c3c6fee19761e8d88692e8e6c331e5889234e`; there were no intervening HEAD
changes to explain before starting. Existing unrelated untracked files and all
M0–M8 source/evidence were preserved. **M8 remains OPEN.** Final bounded
diagnosis, controlled comparisons and kernel feasibility work are complete;
**M8-R1 is complete as a bounded diagnostic result**, with physical
qualification failed/open and the unresolved claims below retained.

The original transition calculation has two distinct problems: an initial
representation that already violates inner covariance, and a later nonlinear
stage whose terminal line search leaves the residual evaluation domain. A
positive source-quadrature k minimum does not describe all the points used by
the dependent epsilon constraint. Improved convergence cannot repair the
retained broad constitutive contradictions.

No terrain calibration, physical-input tuning, field floor, eigenvalue clip,
tolerance relaxation, change of Original/CubicScaled/Regularized identity,
general solver rewrite, or removal of the 3,000-unknown guard is included.
The proposed inner law is an isolated research kernel, not a compiled stress
selection or an adopted weak operator.

## Execution and provenance

Source/Git work used `/Users/jay/Documents/nektar-v2026`. Numerical work used
Tailscale `ssh jay@100.105.247.110`, with large results under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r1` (abbreviated R1).
The unrelated main Linux checkout stayed clean at
`f4b70baec10016dd098005bc5a5507b87e4a57bf`. Separate source/build epochs
preserve the binaries actually used; no historical build was overwritten.

| Purpose | Source revision | External location |
|---|---|---|
| Historical failing flow | `aa854f4c5b05451a81a4e1dd60eaaca98eb547e4` | `m8/source`, `m8/build`, `m8/flow-001` |
| Reviewed startup-corrected preflight/oracles | `ab5c3c6fee19761e8d88692e8e6c331e5889234e` | `m8/initialization`, R1 `archive-source` |
| Original-policy diagnostic replay | `426d48cd44f9f4af017244f6be4c1a6f62d36c01` | R1 `source`, `build`, `reproduction-001` |
| Admissible-difference comparison and guarded controls | `360cbb7ac728b575cdf9af21e6d0963a1d88acb9` | R1 `controlled/source`, `controlled/build` |

The isolated Release builds consume the designated MPI-enabled installation
`/home/jay/.local/opt/nektar-v2026-mpi` and its dependencies. Their single
`geometry-lib` overlay and resolved libraries are fingerprinted. Build commands
use `cmake --build ... --parallel 16`. The allocation ceiling was 19 workers
during builds, well below 72; the later solver/analysis allocation is smaller.
Each native solver uses one MPI rank, explicitly one numerical thread, CPU 0;
focused probes use CPU 1 and independent analyses use CPUs 48–50. All numerical
children set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
`BLIS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, and `NUMEXPR_NUM_THREADS` to 1.
Compilation pins four numerical-library variables; the two additional
analysis settings are an explicit environment difference, not a hidden binary
change. Full before/after manifests and exact commands remain external.

The first build passed. The first controlled build failed under the existing
`-Werror` for misleading indentation in a report-writing statement; its log
is retained. Splitting that statement produced the passing `360cbb7` build
with identical compiler flags. This was a diagnostic formatting correction.
The documented earlier Functions+Capped startup correction was not presumed
to fix transition evolution: the new actual binary replays the old failure.

Six core executed audit helpers are frozen at
`f53bc61112f62be148083fdb5e4b427eeffeac5f` and synchronized to another clean
Linux worktree. [source-binding-receipt.json](cls-epawf-revised-evidence/m8-r1/source-binding-receipt.json)
verifies all 1,801 tracked C/C++/CMake inputs are byte-identical to the
tested `360cbb7` build source and its final runtime still matches the control
pin. New Python/evidence files legitimately change a full input manifest;
no full-manifest identity across that addition is claimed. The exact solver
SHA256 values and all ten build/run manifest epochs are in
[runtime-summary.json](cls-epawf-revised-evidence/m8-r1/runtime-summary.json).
The later control-analysis helper has its own executed-file hash in
`controls-ledger.json`; it adds no compiled input or solver change.

## A. Archived first-failure record

The entire archived `flow-001/qualification/control-coarse/case.evolution.json`,
XML and log were read before new trajectories. Every accepted coefficient
state and complete available history was reconstructed independently using
the clean reviewed source: 324 checks and 1,212 component comparisons pass,
with maximum reported component error `1.448e-11`. See
[archive-ledger.json](cls-epawf-revised-evidence/m8-r1/archive-ledger.json) and
[archive-events.md](cls-epawf-revised-evidence/m8-r1/archive-events.md).
The archived artifact contains no rejected Newton vector and no checkpoint;
these absences are recorded, not reconstructed as historical exports.
An explicitly later archive-only supplement exports **every negative inner
evaluation entry** on the original three constitutive finite sets: 18,180
records over all 12 epochs, with full k/K/nuA/G/S/r_R and eigensystems, plus
the initial field at each identical point. All are on the source-volume set;
the K-sampler and interface sets have none. Its 612 checks match the original
state/history hashes, counts and minima. The 42.24 MB NDJSON remains external,
SHA256 `4259f6dd08247a6a9ffc3136d81751738a2dc71ed927928eca15a85ea3f21e5a`;
[archive-points-index.json](cls-epawf-revised-evidence/m8-r1/archive-points-index.json)
binds the exact file and compact receipts. Enumeration order is not a new
first-minimum event, and this supplement does not rewrite phase-A timing.

Parameters remain Regularized, beta=.4, GammaScale=1, GammaUniform=0,
theta=.2, delta=.1, amplitude=0, eta=0, dt=.005, planned T=.3.
The actual transition weight is a dependent sensor quantity, not inferred
from eta or the requested global chi value.

| Time/event | Recorded observation |
|---|---|
| 0, initial state | Inner covariance already negative: `−.0976947699300747` at `(x,r)=(-.9602898564975362,.5092722755356407)`. Minimum eigenvector `(1,0,1)/sqrt(2)`. Local k=sampled K=`.5105080130598861`, k/K=1; no initial depletion. |
| 0, constitutive membership | Both historical radius .001 and .01 local boxes are already exited. Maximum normalized outer gradient departure `.018332328258716113`; epsilon/epsilon0 reaches `.9833429223226836`. The .01 certificate additionally requires GammaScale=3, which this flow does not have. Membership loss alone does not prove principal failure. |
| .045, first source limiting | Limiter-active physical volume `.001207819098335217`; cap-active volume 0, min lambda=1. |
| .055, last accepted stage | Sampled k and outer epsilon remain positive. Source inner k minimum `.01629424224774017`, inner covariance minimum `−.09983101642016129`. Limiter volume `.09832745735721828`, cap volume 0, lambda=1. Source-ray min b=`.5056987677000203`; trace-ray min b=`.4937805451331741`. |
| .060, failed nonlinear stage | Original stage reports that every admissible Newton trial was rejected. The new diagnostic reproduction identifies the evaluation-domain failure below. |

At the first inner witness,

\[
 G=\begin{pmatrix}0&0&5.191590612064875\\0&0&0\\0&0&0\end{pmatrix},
 \quad S=G+G^T,\quad \nu_A=.08437364924063438,
 \quad r_R=\frac{3\nu_A\lambda_{\max}(S)}{2k}=1.2870516253344717.
\]

The analytical neutral shear at this radius is `1.8959895843783414`.
The initial failure is excessive represented shear with undepleted k.
Horizontal ties occur up to roundoff; the reported x is not a uniquely
identified dynamical trigger. The ledger retains the complete local G/S,
stress eigensystem, sampled/local k, ratios and initial-point comparisons.
There is **no sampled outer covariance loss** in the accepted prefix.

The entire native `robust_source` history was extracted and independently
checked. Lambda is 1 and normalization-cap volume is 0 at every accepted
epoch. Source limiting begins late, but there is **no cap-induced epsilon
trace attenuation** to identify as this failure's cause. Full outer q and
constitutive membership diagnostics are retained; their changing values are
not evidence that an inactive cap caused the nonlinear rejection.
The baseline chi=0 trace itself is lower than incoming epsilon0: preflight
actual-profile ratios span `[.9834528466780685,.9846519510665872]`. For an
initial projected outer interface witness at `(x,r)=(-2.1016667612931865,2)`,
epsilon=`.07142593254506191`, epsilon/epsilon0=`.9833429223226836`,
k=`.5105080130598861`, Gxz=`.45397395788665634`, and
`q_strain=q_rotation=q=3.2447226793354913`. Its full G/S/stress eigensystem
is in `archive-ledger.json#/states/0/sets/interface_outer/worst_point`.
This baseline/trace-projection difference occurs with lambda=1; it is
distinct from cap-induced attenuation and has not been isolated as a cause.

Pressure, continuity, interface and energy are separate observations:

| Quantity | t=0 | Last accepted t=.055 |
|---|---:|---:|
| Weak pressure row norms, inner/outer | 0 / 0 | `5.61e-15 / 1.44e-15` |
| Common gauge residual | 0 | `2.17e-17` |
| Strong divergence L2, inner/outer | 0 / 0 | `.0666061 / .0406064` |
| Full normal traction RMS defect | `.00119693` | `.0558921` |
| Full tangential traction RMS defect | `.000378170` | `.157505` |
| Full traction vector RMS defect | `.00125525` | `.167128` |
| Diffusive k-flux RMS defect | 0 | `.349590` |
| Maximum scaled constraint defect | `5.31e-15` | `1.33e-14` |
| Independent energy accounting remainder | `−1.24e-13` | `1.93e-13` |
| Free energy residual after boundary reaction | `2.58903` | `−3.37e-12` |

The nonzero initial free energy residual is consistent with the substantial
raw unforced residual. At .055, energy work matches the boundary reaction;
eta-anchor work is zero. Weak pressure balance does not imply pointwise
incompressibility or pointwise traction continuity. These interface defects
are reported independently; no controlled result identifies them as the
terminal Newton cause. Mechanical pressures are not required to be equal.

## B. Failure instrumentation and controlled comparison

The original production Newton method uses central-difference columns,
reuses a cached Jacobian when stage alpha permits, and permits only one
stagnation/line-search refresh inside a solve through its `refreshed` flag.
Line-search residual exceptions formerly disappeared into a failed descent
test, and the failure report omitted internal iterates. These observations
motivated instrumentation; none alone establishes a numerical cause.

`RANSWallSolveDiagnostics=Failure` now buffers the solve and writes only its
first problematic stage. It records complete iterate/history/time/alpha,
constraints, field and block residuals, each fraction and evaluation outcome,
exception type/message/context, Jacobian age/rebuild reason, linear residual,
step norms and physical quantities on evaluable trials. Invalid evaluation,
valid non-descent and factor/solve failure have separate classifications.
Default diagnostics are off; detailed mode is serial. Equations, original
central columns, trial fractions, nonlinear tolerances and acceptance gates
are unchanged. Successful stages discard the buffered detailed history.
The reproduced native Newton tolerance is `2e-12`; the later Jacobian-action
reporting threshold is a separate diagnostic quantity.

The optional post-failure Fresh/AdmissibleFresh comparison never advances
native history or returns a substitute production state. The second mode
uses dimensionless variable scales and a bounded central/one-sided stencil
ladder where residual evaluation permits it. This is an explicitly different
diagnostic Jacobian, not a production policy change or a positivity fix.
Its finite-domain label does not make a physically negative-k base state
admissible. A corrected gauge diagnostic scale u*/H=.2 replaces the first
trace's pressure scale .16; that changes diagnostic metadata only.

The first unchanged-policy replay reproduces all 12 archived accepted states
and full available histories, then fails again at .060. Independent analysis
passes 3,656 checks; residual/norm/product discrepancies are at most
`1.04e-15`. All 11 accepted line-search steps descend. There are 131 invalid
evaluations and **no valid residual non-descent trial** in this replay.

The terminal iteration uses a Jacobian nine accepted Newton updates old,
rebuilt within this same .060 stage, with `refreshed=true`. Its linear
residual is about `4.09e-14`; no factorization or linear-solve failure occurs.
This demonstrates an accurate solve of that cached linear system, not that
its matrix equals the current derivative or that stale reuse caused failure.
At that unaccepted base the true scaled velocity-x/velocity-z/k/epsilon
residual norms are `[2.29232184,.44800258,13.63163485,1.14698592]`.
The weak pressure norms are `[5.67e-14,4.61e-14]`, gauge `8.79e-15`, and
maximum scaled constraint defect `.00164974324`. These failed-stage values
are separate from the converged last accepted .055 balance above.

Every one of the 21 terminal fractions, down to `2^-20`, evaluates negative
k on a complete normalization ray used by the **dependent epsilon trace**.
At the base iterate the source-volume GL minimum is `.00226983477`, while
a trace-GLL ray has k=`3.22852626e-9` at
`(x,r)=(1.6046496089512394,4.40175204669948e-7)`. The smallest permitted
fraction gives k=`−1.70055431e-8`, far beyond a `~8.98e-15` roundoff scale.
The unaccepted base iterate also has native-GLL k=`−1.00746589e-6` at the
same x and r=0. External boundary lifting does not produce that change.
The first recorded negative native-GLL k in this .060 solve occurs at
internal iteration 8 (`−1.1713577454207364e-7`). It is an internally accepted
Newton update following iteration 7's fraction `2^-12`, at
`(x,r)=(1.6046496089512394,0)`, never an accepted physical time stage.

Constraints evaluate `DissipationTrace` and its full normalization **before**
the volume residual's exception handler and on a different tangential ray
set. This explains both the formerly empty volume-exception detail and the
misleading positive source-volume minimum. A smaller domain-only fraction
window near `1.52e-7` exists below the original minimum, but it proves neither
descent, convergence nor a physically admissible root. No smaller fraction
was silently accepted.

The strict original fresh central matrix cannot even be assembled at this
base: column 232 minus `2e-6` produces k about `−7.87e-7`, while its plus
evaluation succeeds. That unsuccessful controlled comparison is retained.
It justifies the one additional replay with admissible finite differences,
which again reproduces .055/.060. The latter obtains a fresh matrix but its
same-stage retry does not converge or advance history. Independent action
analysis passes 22,361 bookkeeping, norm, invalid-context and trace checks.
**All 21 fresh-step fractions are also invalid normalization-domain trials**:
zero accepted steps and zero valid non-descent trials. The fresh linear
residual maximum is `4.67e-14`, dimensionless step about `.22920`. At the
smallest permitted fraction the independent trace-ray k is `−1.68871124e-8`,
again far beyond the `8.98e-15` arithmetic scale. Fresh rebuilding therefore
does not rescue this same stage under the unchanged line-search gate. The
refresh flag is not established as its cause; existence of an admissible
stage root remains unresolved.

The fresh matrix uses 293 central columns and 132 lifted constraint
identities. Three columns require smaller differences; 11 invalid coordinate
probes are retained. No one-sided column was ultimately required at the
selected smaller sizes. All 50 directional endpoints in the five-family,
five-level dimensionless ladder evaluate. Their complete field/block error
records include the nonlinear projected dependent-epsilon targets and
constraint derivatives, not a lagged trace.

| Direction | Finest relative fresh-action error |
|---|---:|
| Velocity | `6.27e-9` |
| Free k, independently verified positive physical bubble | `5.97e-9` |
| Epsilon including dependent trace | `1.97e-9` |
| Joint | `1.69e-8` |
| Pressure/gauge | **`2.58e-6`**, fails the fixed `1e-7` diagnostic reporting threshold |

The pressure/gauge discrepancy remains about `1e-6` on coarser differences;
it is an unresolved accuracy observation, not retuned away or confused with
the unchanged native nonlinear tolerance. The coarsest velocity/joint
differences cross five/three outer invariant max branches; finer levels do
not. Limiter/cap patterns and admissible sides are separately checked.
Independent nonlinear trace actions agree to `4.86e-10`. These are finite
derivative observations, not smoothness or conditioning certificates.

At the terminal near-wall normalization point the initial k=K=`.510508`,
Gxz≈11.6694 and r_R≈.02964. At that same point in the last accepted .055
state k≈.0140462, K≈.616178, Gxz≈11.7530 and r_R≈1.19973. This later local
realizability loss tracks severe depletion with little shear increase.
It is separate from the already bad initial represented shear at r≈.5093.
Complete tensors and each invalid context, including formal quantities when
k is negative, are retained with the independently checked failure record.

## C. Initial representation and bounded trajectories

The successful historical p7 zero-step reference used another radial mesh;
it is not a preflight of this flow. The actual fixture was audited at t=0
using two source quadrature levels, keeping GLL points, physical inputs and
each initial coefficient vector fixed. Exactly four zero-step initializations
and 532 independent checks passed. See [preflight-notes.md](cls-epawf-revised-evidence/m8-r1/preflight-notes.md)
and [preflight-summary.json](cls-epawf-revised-evidence/m8-r1/preflight-summary.json).

Both representations use x knots `[-4,-3,-2,-1,1,2]`, p3 velocity/scalars,
p1 pressure and GLL8. Original radials are `[0,2,4]`; the only improved
radials are `[0,.01,.1,.25,.5,1,2,8/3,10/3,4]`. H, samplers, transition,
roughness, beta/gamma/source parameters and acceptance policy remain fixed.

| Preflight, base source quadrature 8/16 | Original | Improved |
|---|---:|---:|
| Cells / actual full unknowns | 10 / 425 | 45 / 1,571 |
| Single dense double array | 1.445 MB | 19.744 MB |
| Indicative FD-volume work / dense factor work | 1 / 1 | 16.63 / 50.51 |
| Raw unanchored velocity-x / velocity-z / k / epsilon physical block norms | 3.71866 / .0272472 / 12.5763 / .132790 | 9.86323 / .025532 / 62.5985 / .111044 |
| Maximum relative inner gradient representation error | 1.7385 | .6539 |
| Maximum outer normalized gradient departure | .01833233 | .00341234 |
| Inner / outer finite covariance minima | −.09773632 / .17433082 | +.07729332 / +.17679689 |
| Outer k/K0 | ≈1 | ≈1 |
| Outer epsilon/epsilon0 | [.98325247,1.00408123] | [.98325247,1.00006559] |
| Maximum finite inner r_R | ≈1.2872 | ≈.7729 |
| Positive finite k/outer epsilon; limiter/cap; min lambda | yes; 0/0; 1 | yes; 0/0; 1 |
| Historical .001 / .01 box membership | fail / fail | fail / fail (old epsilon ±1% also fails) |
| New restricted R1 box membership | fail | fail: gradient departure >.001 |

These are physical mass-dual block norms, not a mixed-unit coefficient
maximum. The improved raw k residual is larger; positive initial covariance
does not make this unforced profile an equilibrium. Pressure/gauge raw
residuals are near roundoff. Native initial coefficients and physics of the
original/base case match the archive exactly. Base/high quadrature pairs have
bitwise-identical initial coefficients. Increasing only independent normal
quadrature changes loads by <`9.44e-15`; raising tangential/outer order as well
changes raw load coefficients by `7.06e-5` / `6.14e-5`. The improved k-block
norm changes `.00215343` on ≈62.6 (mass-dual difference `.00834036`). This is
quadrature sensitivity, not a convergence-order claim. High-quadrature inner
minima are `−.09774725` / `+.07729332`, outer minima `.17432242` / `.17678911`.

The improved initial state is independently finite-admissible in scalars and
covariance, but **not jointly covered by the frozen constitutive certificate**.
The `.0034123386947804146` gradient witness at one-sided r=8/3 remains outside
`.001`; even source GL points exceed it. There was no mesh/box adjustment
after this observation. No polynomial-everywhere admissibility is asserted.

The complete conservative six-slot matrix is bounded as follows. Same-stage
retries are counted as slots even though they do not advance a trajectory.

| Slot | Controlled comparison | Result / interpretation |
|---|---|---|
| 1 | Original mesh, dt=.005, unchanged policy + recording | Reproduces .060 failure, last accepted .055; all accepted histories match. |
| 2 | Same failed iterate, strict fresh original central Jacobian | Assembly leaves finite evaluation domain; no fresh matrix or history advance. |
| 3 | One repeat of original prefix on corrected diagnostic build | Reproduces .060/.055; physical equations/policy unchanged. Required for the justified scaled/one-sided comparison. |
| 4 | Its same-stage admissible fresh-Jacobian retry, max 12 iterations | Matrix obtained; no convergence/history advance. Detailed independent action results retained separately. |
| 5 | Improved mesh, dt=.005, fresh initialization, max T=.06 | First nonlinear stage .005 fails; last accepted time 0. No converged step to qualify. |
| 6 | Same improved mesh, dt=.0025, fresh initialization, max T=.06 | First nonlinear stage .0025 fails; last accepted time 0. Halving dt does not rescue this initialization under the retained solver. |

Slots 5/6 are explicitly unqualified from t=0 because applicability already
exits there. `R1DiagnosticV1` records that first exit and saves/stops the first
audited scalar, covariance or lambda violation before history advance. It
uses source/mass/body/samplers/interfaces, native GLL, inner dense 2p+7 and
outer dense 2p+3 finite sets. Every actual dt change initializes a new IMEX
history; no old multi-step history is reused. There is no long campaign or
reduced-domain replay. The dense 3,000 guard remains in force.

Both improved runs retain their failures, complete internal vectors and
trials. These outcomes do not establish a dt convergence order, transition
passage, absence of an admissible root, or success of spatial refinement for
time evolution. The six allocations are exhausted; no further native
evolution comparison is added.

Independent analysis passes 12,632 implementation/provenance checks. At
dt=.005 the failed first stage contains nine recorded Newton iterates, eight
accepted internal descents and 95 invalid trials. At dt=.0025 these counts
are twelve, eleven and 117. There are **zero valid non-descent trials** in
either run. All 212 exceptions independently reproduce negative k at the
complete normalization points on matching-trace GLL rays. Each final search
rejects all 21 fractions. Source covariance becomes negative by internal
iteration 1; native-GLL k becomes negative by internal iteration 5/7 while
source-Gauss k remains positive. Lambda stays 1 and cap volume zero.
These are unaccepted nonlinear paths, not physical-range guard stops or
accepted physical evolution. Both preserve time/steps=0 and their complete
fresh startup histories; their initial state matches the corresponding
unforced preflight. The strict zero-step guard event/history identities also
match the preflights exactly. See
[controls-ledger.json](cls-epawf-revised-evidence/m8-r1/controls-ledger.json).
The last rejected fractions independently give k=`−1.64477336e-7` /
`−7.45760685e-8`, versus arithmetic scales about `9.1e-15`; sampled K remains
about `.51053`. The offending points are near-wall normal-quadrature points
on the trace-GLL rays, not the values at H. Measured native runtimes are
274.99 / 275.79 seconds; these unsuccessful different-prefix runs are not
a performance comparison. No fresh-J comparison was allocated to either
refined control, so their refresh-policy causality remains unresolved.

The optional finite guard has a retained coverage limitation: its existing
covariance evaluator skips exactly zero inner k, although the old linear
inner tensor is defined for positive sampled K there. No audited exact-zero
case occurs in these controls; negative k is independently flagged. This
prevents claiming a complete exact-zero covariance guarantee from the guard.
No mid-run source change or hidden retry was made to address that edge case.

Focused validation uses the existing `RANSEquilibriumTimeNative`,
`RANSEquilibriumEvolutionNative` and `RANSEquilibriumPatchNative` tests,
sequentially, all passing. These are toy-history, native-field/load and algebra
probes, not additional PDE qualification trajectories. Two zero-step solver
checks also pass: original/DiagnosticV1 stops for covariance at t=0; improved/
LocalV1 stops **only** for `outer_gradient_box_exit`. Each preserves complete
state/history with no accepted audit or history advance. No acceptance
tolerance was changed to obtain these outcomes.

## D. One next-model decision and bounded feasibility

Choose **(a): one explicitly restricted research contract**, with the inner
completion below and the unchanged beta=.4/GammaScale1 outer law only inside
the fixed near-neutral region. This is a concrete next implementation target,
not a currently qualified transition model. The derivation, all dependencies,
pressure/traction conversion, source/trace identities, zero-input cases and
full kernel results are in [candidate-contract.md](cls-epawf-revised-evidence/m8-r1/candidate-contract.md).
The contract and fixed 64-partition plan were frozen before their single
Linux screen: 2,670 checks, zero verification failures, 58.40 seconds.

For the inner candidate set c=2k/3, s=||S||F, a=nuA(K,r)s, rho=4/5,
d=1−rho, and

\[
 t(a,c)=\begin{cases}a,&a\le\rho c,\\
 c-dc\exp[-(a-\rho c)/(dc)],&a>\rho c,\end{cases}
 \quad F=tS/s,\qquad R_A=cI-F.
\]

For k>0 the continuation is C1 at the join, t≥0, t_a≥0, t<c. Its radial
and transverse stress tangents are nonnegative; molecular viscosity supplies
the positive inner momentum margin when turbulent tangents degenerate.
For every symmetric S, covariance ≥c−t>0 at finite a. Trace-free S gives the
stronger rigorous bound `lambda_min(R_A)≥.183503419072273967267571975098*c`.
Full local-k, sampled-K, strain, auxiliary-viscosity, roughness and active
positive-part derivatives are supplied. Uniform turbulent margins require a
declared positive-k/viscosity and bounded a/c regime, rather than unbounded
strain. These claims do not certify a coupled or discrete positive time map.

The complete analytical neutral reference has a/c<`.664850004870988`<rho,
so an open inactive neighborhood retains its stress, profiles, sources,
trace and derivatives. This preserves the analytical neutral reference;
the underresolved coarse projection is outside that neighborhood and cannot
simultaneously retain its negative covariance. At the exact archived first
bad point the candidate kernel gives covariance minimum `+.09997609032554167`.
That is a finite kernel evaluation, not a recomputed native trajectory.

At s=0 use F=0 with the inactive tangent; at nuA=0 use F=0 with the admissible
one-sided viscosity action; at k=0 define F=R_A=0. The joint k/S origin has no
claimed full Fréchet derivative. The same F must supply momentum, production
`P=F:grad(U)−c div(U)`, and every normal/tangential traction component.
The original off-reference analytical AWF elimination becomes a nonlinear
weak stress/source/trace operator when this continuation is active; reusing
its old analytical columns would be inconsistent. That operator is not yet
implemented. Mechanical and reduced pressures remain distinct, with the
derived total normal-traction jump rather than pressure continuity.

Capped pre-source budgets may remain Epre/Dpre, D+=Dpre*w/max(b,delta),
E+=lambda Epre, C+=E+−D+, and the net k source is actual P_F−D+.
The actual dependent outer epsilon trace remains `lambda*Epre(H)`, including
chi=0, with all lambda/profile/constraint derivatives. The source algebra
still cancels redistribution and mean/turbulence exchange; no discrete
positivity theorem is inferred from it.

The single outer applicability region is entrywise normalized gradient
radius .001, local k/K0 in [.99,1.01], epsilon/epsilon0 in [.98,1.02], with
beta=.4, GammaUniform=0, GammaScale=1. The explicit epsilon interval was
selected from the archived actual trace before this fixed proof, and is a
new restricted research box, not a retroactive change to M8's broad claim.
Directed interval/exact-polynomial bounds certify covariance positivity and
the **full incompressible rank-one momentum principal form** ≥
`.53345977306365966796875*nu_t0`. Both limiting max-invariant branches and
their derivative actions are included. k and epsilon are frozen in this
velocity principal certificate; their separate full derivative checks do not
make it a coupled Jacobian coercivity theorem. Positive covariance LDL
pivots in the receipt are not covariance eigenvalue margins. No yx-only test
or sampled scan is presented as this continuum certificate.

All six historical witness families are retained and rescreened, including
the exact necessary-bound conflicts and both invariant-response alternatives.
The unchanged restricted outer law still fails covariance at 37 retained
finite witness directions, outside the restricted region. No unrestricted
gamma/beta tuning is repeated. The broad fixed-gamma family and specified
invariant alternatives remain rejected on their declared box. Common cubic
scaling's qualified transverse repair and capped destruction's earlier
source/trace/energy/native-depletion evidence remain valid within their scope.
For clarity, M8's rejected broad box has gradient-entry radius .25 and both
k/K0 and epsilon/epsilon0 in [.5,1.5]. Its exact beta=.4 necessary bounds
remain GammaScale≥2.44640920758 and GammaScale≤.177525831862; beta=.45 has
the corresponding incompatible 2.27794509320 / .156117116287 bounds. The
new local certificate does not overturn these contradictions or drop their
covariance requirement.

Joint applicability also requires positive sampled K, the **actual** positive
dependent trace, actual outer epsilon bounds and lambda≥1/2. Lambda alone
does not ensure the epsilon box. Lambda approaching zero would drive this
trace toward zero and leave the singular original high-Re outer equations
without an admissible continuation. That model choice is unresolved and
excluded by an explicit applicability stop, not concealed by a floor.

Neither original nor improved flow initialization covers this contract.
Future coverage must be established on actual gradients, k/K, EPS, trace,
branches, samplers and interface/source/physical finite sets before qualified
steps, with independent quadrature/resolution error and a first-exit state.
Finite coverage still would not establish continuum invariance. A scalable
backend may be necessary when the required representation exceeds 3,000.

## Event classification, next action and remaining gates

The [machine-readable first-failure ledger](cls-epawf-revised-evidence/m8-r1/first-failure-ledger.json)
links exact source, executable, XML, full history/checkpoint availability,
result and matrix hashes. Its [event table](cls-epawf-revised-evidence/m8-r1/event-table.md)
separates representation error, inner realizability loss, outer range exit,
scalar inadmissibility, Jacobian/globalization failure, pressure/interface
imbalance and unresolved causes. Successful bookkeeping checks remain
distinct from failed native evolution and the failed pressure/gauge action
accuracy comparison. Failed build/postprocessor attempts are retained.

The recommended next bounded action is a **native weak-operator prototype of
the derived inner completion under the one restricted outer contract**,
including exact/independently checked pressure-gauge and dependent-trace
actions, the full stress/source/traction coupling, and explicit scalar-domain
handling without floors. Its first prerequisite is an actual initial
representation covering the fixed joint region; the current improved fixture
does not. The present data justify this constitutive/admissibility work:
the original same-stage fresh-J comparison leaves its scalar-domain barrier, and
converged residuals cannot repair the initial covariance failure or the broad
outer counterexamples. A native prototype is an implementation experiment,
not a promise of a positive time map. Physical trajectory qualification
requires an independent outer/discrete-admissibility result and coverage.

These later gates remain explicitly open:

- Genuine curved manufactured-solution verification for an adopted law.
- Active capped-source coupling on curved and variable-H patches.
- Outer and fully discrete turbulence admissibility.
- Completed unforced perturbation/refinement and long-time qualification.
- Scalable numerical feasibility where resolution exceeds the dense backend.
- Paired developed floor/hill validation and reserved physical scores.

M8-R1 records a diagnosis and a bounded kernel decision; it does not close M8.
