# WJ-H27 scientific and qualification contract — candidate 1

Status: feasibility investigation, **not qualified for Task 3**. The first
independent reference diagnostic has exposed a negative full constitutive
principal direction. This contract records the candidate selected before that
diagnostic and the prospective gates; it does not select a successful recipe.
Any scientific change needs an explicitly versioned replacement and new tests.

## Authority, identities and scope

The supplied WJ-H27 package is the authority. Mac source/Git work belongs to
`codex/wj-earsm-awf` in `.worktrees/wj-earsm-awf`, reviewed base
`ec0b2f1ba0de65ffabfd41001cc0e0a423ddd9db`. The remote branch was checked at that
exact revision before changes. Preserve the unrelated main checkout and the
pre-existing untracked `Tests/RANS/pilot_wj_h25.py`.

Linux compilation, solver results and timings are authoritative. H27 uses a
new isolated healthy-ext4 namespace,
`/home/jay/build/wj-h27-coupled-20261117`. H25, H25R1 and H26 are closed,
read-only evidence, not H27 initial states or allowances. Source capsules and
operations have create-only manifests/receipts; derived Python bytecode is
disabled. This is a numerical experiment, not a manuscript or a human-reading
attestation. The two independent audit notes document evidence coverage and
the local literature identities.

## Candidate scientific definition

Use the selected full WJ coefficients and physical root, the
`DevStrainRawWorkV1` gradient policy, and raw production `P=-R:grad(U)`.
Molecular viscosity is `6.423173803526447e-5` in both domains; historical pure
WJ calls still use `outerNu=0`. The current positive implicit coefficient
multiplies the **full** symmetric strain. Its complete remainder restores
`V=nu M-R+2k I/3`, with `Pi=p/rho+2k/3`. Positive split diffusion alone does
not establish ellipticity of that complete tensor.

Keep the specified tensor `Rh=P Rw P+(n.Rw.n)nn+bA n+n bA`, where
`bA=-nuA P M n`, and `R=(1-chi)Rh+chi Rw`. In particular, keep its WJ normal
stresses below the collar. Chi is zero for `r/H<=.4`, one for `r/H>=.75`,
and the quintic `s^3(10-15s+6s^2)` between them. The fixed geometry is
`H=.15`, `h=.10854503542535461`, `d0=.0875`, physical length unit `.04 m`,
domain `x in [-10,15]`, and top physical elevation `20`.

Transport Q3 k over the whole domain and Q3 epsilon only outside H; velocity
is Q3 and modified pressure Q2. Ew is dependent. Use conservative scalar
transport, `Dkw=nu+Cd k^2/(sigma_k Ew)` inside and the stated outer k/epsilon
diffusion/source equations. `Cd=572/6561`, `sigma_k=1`, `sigma_e=1.3`,
`Ce1=1.44`, `Ce2=1.82`. The candidate inner source is unamended `P-Ew`.
There is no clipping, floor, EP redistribution, or balance forcing.
Positivity, realizability, physical root and full applicable principal
directions are acceptance conditions, not unconditional claims. Reject a
failed state while preserving accepted history.

The amplitude is `A(K)=K^(3/2)/cl`, with `qe=Cd^(-1/2)`,
`kappa=.3819178635467331`, `cl=kappa qe^(3/2)`. Thus
`Ew=A/max(r,c)` and `Ew(H)=A/H` for `0<c<=H`. This preserves the exact
retained inlet corner `Kin=.008809015145515192`,
`epsilon_in(H)=.002315524054799254`.

The explicitly selected **local counterfactual calibration column** has
constant reference total shear `t=K/qe`, the same finite-nu tensor,
`Q'=E-P`, `k'=Q/Dkw`, `k(H)=K`, and `Q(H)=0`. Determine c by `Q(0)=0`.
Equivalently `A[1+log(H/c)]-integral(P)dr-Q(H)=0`. The independent oracle
solves k and shear at each trial c; no old closed-form cutoff is reused.
The new velocity anchor is `integral(G)dr`, not old `t J0`.
Zero reference endpoint flux is an explicit local calibration choice, **not**
a theorem about the finite-nu outer flow or an extra actual interface law.
Actual storage, transport and interface flux retain the actual imbalance.

Initially bound `K/Kin in [.9,1.1]`, `log(c/H) in [-4,-.2]`. Select only the
increasing planar shear branch `0<=q<=72/13` where chi is positive; at chi=0
the auxiliary shear equation is unique without that q restriction. Reject
out-of-range states; do not clamp. Require independently stable root slope,
split-quadrature budget/weak moments, trace and admissibility before using any
interpolation table. A finite sampled positive principal minimum is not a
continuum certificate; a repeatable negative direction is a rejection.

## Application boundaries and forcing

The complete historical application inlet stays unchanged, including old
inner `tau0 J0_old`, old inner k, full-height outer logarithmic U, constant
outer Kin, and outer `Ain/r`. The new column never replaces these data.
Ground: no slip, zero k diffusive flux. Top: zero normal velocity, prescribed
tangential total traction `tau0=.002601`, zero k diffusion and fixed outward
epsilon flux `Cd Kin^2/(sigma_e rT)+nu Ain/rT^2`, `rT=19.9125`.

Outlet: prescribed complete reference traction, zero k/epsilon diffusion,
and explicit backflow rejection. The separately pinned finite-nu reference
load is `(Vxx_ref-Vzz_ref+Kin/Ne,0,tau0)`, `Ne=81/20`, with the selected
inner column and an outer counterfactual `k=Kin`, `epsilon=Ain/r`, constant
finite-nu total shear. It is not fitted to a developed solution and is not
claimed to satisfy the outer stationary scalar equations. Its construction
must pass reference feasibility before publication as executable input.

Across H require common U and k and conservative mass/traction/k fluxes.
The sole epsilon condition is outer `epsilon(H)=Ew(H)`. There is no old
normal-stress pressure jump or artificial C1 U condition. Determine gauges
from the restricted assembled constant-pressure action under actual
component constraints; the physical open-traction outlet fixes the level.

The flat roughness force is
`b+1_(r<=h)[2 U_H d_x U_H+d_x Pi_H-(2/3)d_x k_H-bx_H] e_x`.
Use the outer H owner, actual transported Q3 k action and complete force work
in mechanical energy. No sampled local time derivative is added.

## Numerical and physical gates fixed before selection

Use Uref=1, Kin for k, tau0 for stress/Pi/traction, positive local inlet
epsilon `Ain/r` outside, and declared physical length/gradient scales. Use
mass/Riesz and independent physical-integral norms; no large time-mass
normalization or dimensional `max(1,rhs)`. Require operator/trace/source and
resolved reference-budget/trace agreement `1e-10`, current momentum/weak
continuity `1e-9`, current k/epsilon `1e-8`, tiny direct-reference fields
`1e-8`. For principal admissibility the necessary prospective lower bound is
strict positivity with a numerical error reserve in units of current
`nu+(1-chi)nuA+chi C k^2/E`; no historical high-Re quality floor is imported.

Verification includes independent evolving coupled continuous MMS, at least
three dt levels with U/k/epsilon order >=1.8, spatial/split-quadrature and
tolerance refinement, independent k/epsilon/mechanical budgets, mixed gauges,
Q3 known-k force/sign, failed/stale/root/invalid-input fixtures, and full
current-data residuals. Shared assembled fluxes cancel to arithmetic tolerance;
physical C0 gradient/traction discrepancies must converge with resolution.

A new same-inlet steady 2D solution needs complete steady residuals and a
finer h/p/q check. Prospective physical development envelope, distinct from
numerical error and experimental uncertainty: U/W scaled RMS<=.01,
MAX<=.02; k RMS<=.05, MAX<=.10; outer epsilon RMS<=.10, MAX<=.20;
ground upward shear RMS<=.05, MAX<=.10. Apply whole-domain and inner/outer
regions and preserve the historical station union in the boundary audit.
Report R, Pi, gradients, Ew and cutoff separately. The envelope is an
engineering criterion; no measured CFD uncertainty supports it.

Run unforced developed control, retained-inlet startup, and finite perturbation
pilots. The principal evolving pilot endpoint is at least
`26.13350407707765`, one bulk transit and at least three interface k/epsilon
clocks. Compare three dt levels at common endpoints for startup and perturbation
and compare a finer h/p/q recipe at common physical, owner-resolved samples.
Numerical RMS/MAX budgets are U `1e-4`, k `1e-3`, stress/Pi/traction `1e-3`,
gradient `1e-2`, local outer epsilon `1e-2` in the scales above.

Pin perturbation peak U `1e-3 Uref`, k and independent outer-epsilon `1%`
of their reference scales; compact support crosses H in
`x in [-2.5,2.5], r in [.075,.45]`. Require resolved support area>=1.5,
half-peak area>=.05 and >=4 cells in each relevant direction. Domain RMS
increments are >=`1e-5 Uref`, `1e-4 Kin`, `1e-4` local epsilon; support RMS
is >=10% of each peak. Interior K changes by 1%; recompute Ew and c.
Dependent epsilon-H changes by `(Knew/Kbase)^1.5`, not an incompatible
independently imposed 1%. Include a second sign/shape and measurable gradient,
k, epsilon, K and cutoff changes. The audit gives the complete support norms.

Cold restart covers all states/time levels/explicit histories/scientific time,
wall/reference branches, model/geometry/parameter identities; rebuilt caches
are checked and charged. Target endpoint scaled differences<=`1e-8`.
Exercise >=3 actual scheduled recovery writes and rotation/deletion, usable
retained recovery files, scientific reduction and terminal retention.

## Resource and cost admission

Aggregate H27 ceiling is 43200 Linux operation seconds, 3600 protected for
evidence. Each operation has a unique request, complete elapsed charge
(including failures), bounded log and resource receipt. Initial access/control
has a conservative 120-second charge. Never reopen the closed older accounts.
Never exceed 72 actual shared-host workers; do not stop unrelated jobs.
Initially one science process/rank and all six numerical thread controls=1;
build at most 4 workers in this implementation (package maximum 8), while
reserving helper processes. Science RSS<=8GiB, build<=16GiB, new source/build
<=8GiB, science<=8GiB, free reserve>=15GiB.

Before expensive implementation/pilots, pass candidate scientific feasibility.
The initial admitted reference diagnostic cap was 300s; its actual charge is
19.4490732899867s. A follow-up feasibility matrix may reserve at most 1800s:
bounded-K/refinement root checks, alternate explicitly labelled inner-D and
endpoint diagnostic, independent planar/native witness and invalid inputs.
This reservation does not admit physical pilots. Their whole matrix must be
forecast from measured complete coupled prefix rates and admitted against the
remaining account before execution: same-inlet steady/finer reference, three
physical classes, both startup/perturbed three-dt windows, finer h/p/q window,
three representative warmed blocks, actual rebuild/retry, cold restart and
three-write retention, plus assessment/publication. Unknown costs remain
unknown; kernel timings cannot close admission or production gates.

Task3 keeps horizon `505.02466521278404` and terminal window
`[447.8483832973273,505.02466521278404]`. For **each** required trajectory,
require `1.25 Ctotal<=86400s`, using the accuracy-qualified coupled dt and
complete initialization, accepted/rejected progress, all recurring work, one
full cold restart and final assessment. Report the sum over the campaign too.
Obtain >=3 representative warmed coupled blocks, actual rebuild and retry
provision; use the largest representative class mean, attribute >=95% without
double counting, and report equation/model work fraction separately (>50%
aim, 70% preferred). No complete cost can be inferred from H26 affine kernels.

## Checkpoints and closure

A records feasibility, this contract, gates and resource plan. B requires the
registered genuinely coupled solver and independent verification. C requires
the developed reference and physical/accuracy pilots. D requires cold restart,
retention, complete cost and executable Task3 manifest. Commit/push and verify
each exact remote checkpoint; never mislabel an unimplemented gate as passed.

`TASK3_READY=YES` requires all five package gates to pass on one pinned recipe.
Do not run the full Task3 horizon. If a genuine model barrier is independently
confirmed, publish its state, direction, convergence/native checks and a
specific proposed remedy for scientific review, close the finite resource
account, mark readiness NO, and preserve the unmet positive-closure goal.
