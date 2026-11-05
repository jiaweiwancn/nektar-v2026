# M8-R10 safeguarded optimization of the destruction model

R10 starts from the R9 **initial** state, the physical R8 iterate-10
endpoint. The reviewed branch head is
3b906f7fd2d2e657526ae291165ab2e4d6b0e901. The remote branch had no
intervening commits when implementation began. Source changes use the
existing branch worktree at /private/tmp/nektar-m8-r8. The separate
Picard checkout and its untracked files are preserved.

The frozen state has SHA256
eb063e15f16b2610f1c5a832280d87894f567e8dab321560202e6dd9af95194c.
Its current physical matrix is the R9 solve1.J0.bin, with SHA256
2194965dee9c5960ba4e1c1604f76e8e85372cabef2d5ced7c8a394f593c8dc6.
The original history has SHA256
b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71.
The loader binds all three to the original residual and constraints
before importing the R9 destruction value and action.

The comparison retains outer6, 2,165 unknowns, p3/p1 spaces and quadrature,
alpha 400, time .0025, eta zero, the original initialized history and zero
forcing. SaturatedR2 with rho .8, beta .4, GammaScale 1, GammaUniform zero,
the original outer response, viscosities, scalar diffusion and projected
trace remain in the physical operator. Capped destruction uses the bound
Linear smoothing selector with theta .2 and delta .1. The total radius is
.125, the scalar and continuous-wall reserves are 5%, and the positive
ground nondepletion restriction remains disabled. Root, action and
original-row gates are 2e-12, 1e-7 and 5e-12.

The numerical selector RANSWallR10Experiment=On enables the comparison with
JointRestorationR5, RefinedR6, ScalarFeasibleR3 and VerifiedR3Action.
Earlier selectors keep their existing paths.

For a total displacement $d$ from the frozen $q$, R10 uses
$$
 M_q(d)=R(q)+J(q)d+D(q+d)-D(q)-J_D(q)d.
$$
The objective is half the sum of squared free differential block norms
in the native consistent-mass dual metric. At an accepted model anchor
$a$, the tangent is
$$
 A_a=J(q)-J_D(q)+J_D(q+a), \qquad
 r_a=M_q(a)-A_a a.
$$
A tangent solve produces a new total displacement $d_\star$.
Backtracking evaluates $a+\eta(d_\star-a)$, using the frozen-base
radius, reserves, outer cuts and full affine algebraic restoration.
Every accepted anchor passes the native and independently accumulated
original-row scans, the actual scalar/source domain and continuous
ground trace. The initial full linear displacement must pass these
checks before it can initialize the incumbent and best candidate.
The scaled accepted native displacement is explicitly checked for its
remaining affine defect and cannot substitute for the restored anchor.

The model Armijo constant is 1e-4. Each accepted anchor must also satisfy
the native slope bound $g_E(d)\le-2\times10^{-4}E(q)$. The slope comes
from physical $J(q)$ and the same mass-dual metric as the native norms.
Active norm ties use every active gradient, including the gauge and
maximum boundary block when needed. A source-valid model-decreasing
trial blocked by that slope can trigger one tangent retry containing
the slope as a general affine inequality. Its original units are
audited separately from spatial observations.

Before line search, R10 checks
$s_\Phi+\|A_ap\|_W^2/2$ against a tolerance of
$10^{-9}\max(1,\Phi,\|A_ap\|_W^2)$, consistent with the retained affine
subproblem KKT check. A focused source-only Richardson check evaluates
$J_D(q+a)p$ at the initial moved anchor. It repeats after an accepted
branch change. Ordered limiter, cap, profile and panel identities
determine whether a centered stencil stays on the selected branch.
An admissible forward stencil can check a matching one-sided branch.

The declared comparison budget is four tangent/QP solves including
the conditional retry, ten model Armijo fractions per tangent, and
50 new source-value evaluations including focused checks. Repeated
displacements reuse the source cache. Model construction stops at 540
seconds, new tangents stop at 500 seconds, and final ray testing starts
only before 565 seconds. The comparison phase has a cooperative
600-second ceiling and a 650-second independent native supervisor
ceiling to preserve final reporting.

Each compared ray has the original R8 selector with at most nine native
trial states, fitting and held-out fractions .125 and .0625, the same
quadratic proposal grid, dyadic fallback and native $10^{-4}f$ decrease
requirement. This selector does not evaluate the full native $f=1$
endpoint. Identical source and control directions share the native
result. Promotion keeps the R9 gain, ground, trace and gain-per-second
criteria.

Authoritative execution is pending access to the preserved Linux volume.
The server responds over Tailscale, but /dev/sda2, labelled 新加卷, is
unmounted. The standard user mount returned NotAuthorizedCanObtain, and
noninteractive sudo returned “a password is required.” The R9 immutable
state, history, physical matrix and source matrices are on that volume,
so the Linux binding check, isolated build and numerical comparison have
not started. The committed implementation has syntax and diff checks
from the source worktree, with native verification pending.

The pinned launch script is queued on Linux as
/tmp/cls-m8-r10-launch.sh, with log /tmp/cls-m8-r10-launch.log. It waits
up to 30 minutes for the archive mount, then creates a fresh isolated
m8-r10/epoch-001, builds with 16 explicit workers, prepares the bound
comparison, runs the native experiment, collects its durable records and
runs the independent terminal audit. A second queued script runs the
three focused native CTests and nine predictor-binding checks after
the runtime fingerprint is published. Every numerical library is
allocated one thread. This queue has performed no native evaluations.

Model monotonicity and tangent consistency are awaiting the Linux
experiment. Native-merit descent is implemented and awaits verification
against the retained control and R9 ascent direction. Practical gain and
cost remain unmeasured. R10 has obtained no coupled root and has accepted
no physical time step. M8-R10 remains incomplete while the archive mount
is unavailable.

M8 remains open. Curved-wall MMS, active-cap and variable-height coupling,
full-element and fully discrete positivity, refinement, long-time
robustness, scalability and physical floor/hill validation remain open.

