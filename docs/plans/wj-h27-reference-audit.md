# H27 independent reference and inlet audit

Status: analytical derivation and an independent numerical oracle prepared;
Linux numerical execution is pending. No numerical validation was performed on
the Mac. Source inspected at `ec0b2f1ba0de65ffabfd41001cc0e0a423ddd9db` in the
branch-owning worktree. This note does not mark H27-A or Task 3 ready.

The scientific choice proposed by the parent is a **local counterfactual
calibration column**, with zero reference endpoint k flux. This is a permissible
definition, provided it is named explicitly, has a qualified admissible branch,
and is never substituted for the new two-dimensional steady application.
There is no derivation that makes the actual interface k flux zero. A global
finite-nu stationary solution can have a nonzero interface k flux, even though
the selected local calibration column has zero flux at both ends.

## Source evidence and the boundary that must change

All line numbers below refer to the inspected revision unless the path is the
new oracle. Literature paths are relative to the main checkout.

| Evidence | Exact location | Consequence |
|---|---|---|
| Physical WJ root, general explicit tensor, raw stress work | `library/SolverUtils/Turbulence/WJEARSM.cpp:54-71,98-152` | Keep the full tensor and raw contraction in the PDE. Simple shear is an independent reference reduction only. |
| High-Re pure kernel rejects nonzero molecular viscosity | `WJEARSM.cpp:78-80` | Add molecular viscosity in the new PDE, not inside this historical constitutive API. |
| Independent scalar Cd and sigma | `WJEARSM.cpp:115,138-142` | Positive momentum splitting does not define scalar diffusion. |
| Existing calibrated tail amplitude and traction | `library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp:116-121,179-183` | A(K) can remain unchanged while the reference is rederived. |
| Old production budget and exponential cutoff | `EquilibriumAWFProfiles.cpp:197-215,469-475` | These are independent of the old cutoff; they cannot be reused for the new nonlinear collar. |
| Original roughness/intercept law | `library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp:92-103,118-125` | Preserve it and the fitted roughness; do not re-fit roughness to the new anchor. |
| Exact high-Re shear reduction and reference constants | `docs/literature/wj_earsm_budget_cutoff_awf_appendix.tex:140-178,240-291` | Defines q_e, Cd, kappa, c_l and the independent planar oracle. |
| Old inner diffusivity and source | same appendix `:312-326,469-510` | Changing inner diffusion to the selected local-k/Cd law is an explicit H27 choice. |
| Old zero endpoint flux follows from its outer constant-k reference | same appendix `:423-440,479-503` | This is not a finite-nu outer theorem. |
| Actual versus reference budget; source positivity counterexample | same appendix sections `ba:actual-budget`, `ba:dynamic-positivity` | Do not impose instantaneous actual integral(P-E)=0 or infer evolving k positivity from a positive reference. |
| Retained inlet formulas, units, viscosity, interface derivative ratio | `solvers/IncNavierStokesSolver/Tests/RANS/wj_h5_neutral.py:18-20,119-138,149-164` | Retain the entire prescribed inlet. |
| Exact sole6mm physical calibration | `docs/plans/wj-h5-linux/calibration.json:513-560` | H=.15, h=.10854503542535461, d0=.0875 in 40mm units; physical H=6mm and matching elevation=9.5mm. |
| Existing physical top/reference context | `solvers/IncNavierStokesSolver/Tests/RANS/WJH5NeutralInlet.md:3-25` | A top shear law is not a prescribed top tangential velocity. |

## Explicit selected column

Use wall/model distance r in [0,H]. Physical elevation is d0+r for this flat
case. Do not confuse H with physical matching elevation d0+H. At fixed sampled
K>0 set

\[
q_e=C_*^{-1/2},\quad C_*={572\over6561},\quad
\kappa^2=\sigma_e(C_{e2}-C_{e1}){C_*^{3/2}\over C_d},\quad
c_l=\kappa q_e^{3/2},\quad
t={K\over q_e},\quad u_*=\sqrt t,\quad A={K^{3/2}\over c_l}.
\]

The selected coefficients have Cd=C_* numerically. Cd remains a distinct
scalar coefficient. The oracle freezes the selected kappa when testing other
scalar choices; it never recalibrates the inlet. Keep the same empirical
intercept b(K,h,nu), z_v=nu b/sqrt(K), and auxiliary viscosity
nu_A(r)=kappa sqrt(K/q_e) max(r-z_v,0). Chi is the specified quintic transition
from zero at .4H to one at .75H, constant outside this interval. These
geometric operations are the declared model, not stress or state clipping.

For trial cutoff 0<c<=H,

\[
 E(r;c)={A\over\max(r,c)},\qquad
 D(r;k,E)=\nu+{C_d k^2\over\sigma_k E}>\nu>0.
\]

This inner diffusivity is an explicit change from the appendix's sampled-K
auxiliary diffusivity. It uses local transported k and the actual dependent E
in the actual PDE; it is identical to the outer k diffusion at a common k/E
trace. It is not the implicit momentum coefficient mu. The proposed healthy
case uses the unamended source P-E. There is no EP or bounded-destruction
redistribution in this definition. Positivity is an admissibility condition
with accepted-history-preserving rejection/continuation, not a theorem of this
sink. If a later amendment is selected, its actual E, source and outer trace
must all change together and this audit must be rerun.

Take a parallel reference U=U(r)e_x. Set G=U'>0 and q=kG/E. The physical planar
WJ root and turbulent shear function are

\[
 N^3-c'N^2-{7\over20}q^2N-c'q^2=0,\quad N\ge c'=9/5,
 \qquad f(q)={3Nq\over5(N^2+q^2)}.
\]

The **complete selected tensor** reduces to

\[
 R_{xz}=-[(1-\chi)\nu_A G+\chi k f(q)],\quad
 R_{xx}=k\left({2\over3}+{3q^2\over5(N^2+q^2)}\right),\quad
 R_{zz}=k\left({2\over3}-{3q^2\over5(N^2+q^2)}\right),\quad
 R_{yy}=2k/3.
\]

The WJ diagonal blocks are retained even where chi=0. Replacing them by
isotropic normal stresses would be a different wall tensor. In the reference,
finite-nu constant total shear determines G implicitly:

\[
 t=\underbrace{[\nu+(1-\chi)\nu_A]}_{\alpha}G+\chi k f(kG/E),\qquad
 P=-R:\nabla U=(t-\nu G)G.
\]

Mean molecular work nu G^2 is excluded from turbulent production. A convenient
new velocity anchor and mechanical pressure are

\[
 U(r)=\int_0^r G(s)\,ds,\quad U_H=\int_0^H G(s)\,ds,\quad
 p(r)/\rho+R_{zz}(r)+gz=\text{constant}.
\]

This anchor is a column diagnostic; it does not overwrite the prescribed
application inlet or prescribe the actual wall traction. The old formula
t J0(H) is no longer the anchor for this nonlinear tensor.

The reference k equations, including a general specified gradient flux Q_H,
are

\[
 Q'=E-P,\qquad k'=Q/D,\qquad k(H)=K,\qquad Q(H)=Q_H,\qquad Q(0)=0.
\]

The cutoff is the extra unknown satisfying the last condition. The physical
upward diffusive flux is -Q, not Q. The equivalent integral equation is

\[
 F(c;K,Q_H)=A[1+\log(H/c)]-\int_0^H P(r;c,k_c)\,dr-Q_H=0.
\]

Thus k(r)=K-integral_r^H Q(s)/D(s) ds is now an **implicit** profile; using the
old explicit profile in this budget would leave stale k and stress. The parent
selects Q_H=0 only for the counterfactual column. If a reference is instead
matched to a specified finite-nu outer stationary solution, its computed
Q_H=D_o k'_o(H) must be included. An outer epsilon flux condition is never
added: the single trace is epsilon_o(H)=A/H.

For every test function phi the independent moment identity is

\[
 -\int_0^H Q\phi'\,dr+\int_0^H(P-E)\phi\,dr
 +Q_H\phi(H)-Q(0)\phi(0)=0.
\]

During actual evolution, reference P is never substituted for actual P and
the actual storage/tangential transport/interface k flux close the actual
energy balance. Q_H=0 is not an extra actual-flow boundary condition.

## Root branch, conditioning and failure conditions

The planar cubic has a unique root N>=c'. Its derivative is
F_N=3N^2-2c'N-(7/20)q^2>0 on the physical branch, and

\[
 N_q={2q[(7/20)N+c']\over F_N},\qquad
 f_q={d\over dq}{3Nq\over5(N^2+q^2)}.
\]

There is no unconditional monotonicity of f for all shear. In fact

\[
 f^2={.35+.65c'/N-(c'/N)^2\over(9/4)^2},\qquad
 f_q\ge0\text{ for }0\le q\le q_{peak}=2c'/.65=72/13.
\]

The first oracle uses this bounded increasing-shear branch wherever chi>0.
For chi=0, G=t/(nu+nu_A) is unique independently of q. On that branch

\[
 m=\partial_G t=\alpha+\chi{k^2\over E}f_q\ge\nu>0.
\]

A trial with no shear root on that interval is rejected. This is not a proof
that no other shear branch exists, and an unexamined alternate branch cannot
be silently selected. Positive m proves only simple-shear invertibility,
not full constitutive ellipticity.

For cutoff conditioning, write lambda=log(c/H) at fixed K and Q_H. The
sensitivities satisfy

\[
 E_\lambda=-E\,\mathbf1_{r<c},\quad
 G_k=-{\chi[f+q f_q]\over m},\quad
 G_E={\chi kq f_q\over Em},\quad
 G_\lambda=G_k k_\lambda+G_E E_\lambda,
\]
\[
 D_\lambda=(D-\nu)(2k_\lambda/k-E_\lambda/E),\quad
 k_\lambda'={Q_\lambda D-QD_\lambda\over D^2},\quad
 Q_\lambda'=E_\lambda-(t-2\nu G)G_\lambda,
 \qquad (k_\lambda,Q_\lambda)(H)=(0,0).
\]

E is continuous at r=c, so moving that split does not generate a delta source.
The cutoff residual derivative is

\[
 F_\lambda=-A-\int_0^H(t-2\nu G)G_\lambda\,dr
          =-Q_\lambda(0).
\]

The extra production sensitivity destroys the old monotonicity argument.
The normalized root condition measure is A/abs(F_lambda). A root with
F_lambda nonzero is locally regular by the implicit function theorem. A
simple sign-changing bracket alone does not prove uniqueness beyond that
bracket. The oracle reports two central-difference slopes and retains every
scanned rejected trial and every sign-changing bracket; continuation and
denser scans must establish the intended bounded K range before caching it.

Prospective local branch screen (to be frozen in the parent contract):

- Start K/K_in in [.9,1.1], extend only after checking the whole enlarged box.
- Search log(c/H) in [-4,-.2], with explicit out-of-range rejection. Broaden
  diagnostic scans to [-12,0] only to identify an obstruction, not to pretend
  the original bounded table was valid there.
- Require positive k and E, positive Reynolds-stress covariance, positive
  full applicable principal directions, a nonzero cutoff slope with stable
  sign under refinement, and budget/trace errors <=1e-10 in their own scales.
- Use k/K_in, E/(A_in/H), covariance/k and the physical tangent scale
  mu_ref=nu+(1-chi)nu_A+chi C(q)k^2/E. Do not normalize by a large time-mass
  coefficient or dimensional max(1,rhs).
- Report the sampled minimum of the symmetric transverse acoustic tensor
  H:(nu(H+H^T)-D_gradR[H]), H=l tensor n, l perpendicular n, divided by mu_ref.
  A negative witness is decisive; a positive finite sampling is not a proof
  of the minimum over all directions. Qualify minima by angular refinement
  and local minimization, with a prospective reserve against derivative and
  interpolation errors. Do not reuse the historical high-Re floor unexamined.

## Inlet/interface corner and the stationary outer obstruction

The retained exact values are K_in=q_e(.051)^2 and
epsilon_in(H)=(.051)^3/(kappa H). Hence the unchanged A(K) gives exactly
E(H;K_in)=epsilon_in(H), for every admissible c<=H. K, A and the outer trace
are not three independent prescriptions. Replacing A by a power of a
finite-nu total traction would generally break this compatibility and is not
permitted without a separate model change.

The retained mean inlet is continuous but has different one-sided derivatives
at H:

\[
 G_{in,-}={u_*^2\over\nu+\nu_A(H)},\qquad
 G_{in,+}={u_*\over\kappa H}.
\]

With the new identical full-WJ collar law this is not a smooth stationary
parallel traction-matched profile. Continuous Dirichlet value data can still
define a weak two-dimensional problem; a derivative kink at the inlet/interface
intersection is not by itself a proof of nonexistence. Its possible corner
regularity loss and one-sided traction mismatch must be retained in spatial
and corner refinement tests. Do not overwrite either prescribed value or
claim artificial C1 velocity continuity. Scalar corner values are compatible;
actual flux matching is a conservative PDE condition away from the boundary
intersection and requires its own convergence evidence.

There is a sharper stationary obstruction to reusing the old constant-k outer
solution. For a stationary parallel state with k=K constant, the k equation
implies P=epsilon and q=q_e. Then turbulent shear is K/q_e constant. Finite-nu
constant total traction forces G constant, hence epsilon constant. The epsilon
equation then reads 0=(C_e1-C_e2)epsilon^2/K, which is impossible for positive
epsilon and the selected unequal coefficients. Thus a new global finite-nu
stationary reference must allow varying k and/or streamwise development. This
does not obstruct the local counterfactual column, which has no outer epsilon
equation, or the requested two-dimensional same-inlet developed reference.

The original sampled inner sink still has a first-zero positivity obstruction.
Neither the selected diffusion nor a positive reference cures it. For k=0 on
an open stationary-flow patch, diffusion/advection vanish and E(K)>0 can give
k_t<0. A healthy bounded branch can be admitted with failure detection; it
cannot be advertised as unconditionally positive or globally robust.

## Independent Linux experiment matrix

Program: `solvers/IncNavierStokesSolver/Tests/RANS/h27/reference_oracle.py`.
It uses an independent polynomial-root evaluation, DOP853 shooting split at
c, .4H, .75H and z_v, independently adapted Gauss–Kronrod source/moment
quadrature, and the full 3D explicit tensor for covariance and principal
direction tests. It does not call native/old reference solvers. The endpoint
flux is an exposed parameter. Every output is create-only; run identities and
all failures belong to the new H27 ledger. Use one process and all numerical
library thread controls=1. No jobs have been launched by this subagent.

| Run | Exact additional arguments | Purpose / limit |
|---|---|---|
| R0 | `--k-factors 1 --scan-points 41 --principal-directions 128` | First finite-cost branch/negative-witness screen; 120s forecast cap. |
| R1 | `--k-factors .9,.95,1,1.05,1.1 --scan-points 81 --principal-directions 512 --rtol 2e-12` | Bounded interface-state branch and tighter numerical audit; admit only using measured R0 cost. |
| R2 | `--k-factors 1 --endpoint-flux-over-A -.05 --scan-points 81 --principal-directions 128` and separate `+.05` output | Detect missing endpoint term; compare root sensitivity and weak moments. This does not change the selected Q_ref=0 model. |
| R3 | `--k-factors .5,1.5 --log-min -12 --log-max 0 --scan-points 121 --principal-directions 128` | Adversarial/out-of-range diagnostics; a success does not enlarge the selected range automatically. |
| R4 | repeat R0 with `--inner-diffusion auxiliary` | Explicit different-model comparison only; detects accidentally using old diffusion. Never mix into selected evidence. |
| R5 | `--k-factors 0,-1` | Invalid positive inputs must be rejected and recorded. |

Prepend `python3 <source>/solvers/IncNavierStokesSolver/Tests/RANS/h27/reference_oracle.py`
and append `--output <fresh-run-dir>/<run>.json`. If R0 is slow, reduce the
first scan to 17 points as a separately identified preliminary search; retain
R0's eventual 41-point qualification. Do not accept an isolated grid scan as
proof that double roots cannot exist between grid points. Verify log-cutoff
slope throughout accepted brackets, add interval/continuation evidence, and
cross-check independent BVP collocation or another quadrature/root method
before claiming branch qualification at <=1e-10.

If no admissible zero-flux root exists, first determine whether the failure is
the bounded search, positivity, covariance, simple-shear inversion or a full
principal negative direction. A compliant next reference candidate is a
counterfactual column matched to a finite-nu outer BVP's Q_H, retaining A(K),
roughness and inlet, with its endpoint term explicitly included. Changing
Q_H is a change of calibration definition that needs a new version and tests;
it must not be chosen from instantaneous actual P merely to zero its budget.
If a full principal obstruction already occurs at the fixed prescribed inlet,
document its state and direction before considering any model amendment.

## Follow-up: the ground principal obstruction is planar

The parent reports the authoritative Linux R0 result c/H=.19696577127569226,
k(0)=.014362547557221588, E(0)=.011755971810747884,
G(0)=22.635162605006414 and q(0) about27.65. The reference budget residual is
about3.7e-12 in its scale, and covariance is positive, but a stable negative
full principal sample occurs at the ground. These values are parent-supplied
Linux results; this subagent has not executed numerical work. The derivation
below independently shows that this is an applicable **in-plane** obstruction,
not a spanwise-only artifact, scalar positivity failure or linear-solver issue.

Hold k, E and sampled K fixed, and let chi=0. At the simple-shear state consider
the unit incompressible rank-one direction

\[
 n=(\cos\theta,0,\sin\theta),\quad
 l=(-\sin\theta,0,\cos\theta),\quad B=l\otimes n,
 \quad\operatorname{tr}B=l\cdot n=0.
\]

The principal quadratic of the full stress is
a=B:[nu(B+B^T)-D_grad R_wall[B]]. The auxiliary off-diagonal block contributes
nu_A cos^2(2theta), but the retained WJ diagonal blocks have a nontrivial,
potentially adverse tangent. With x=cos(2theta), s=sin(2theta), tau=k/E,
D=N^2+q^2, r=c'/N and B_r=.7+2.65r-2r^2, exact planar differentiation gives

\[
 \boxed{a(\theta)=\nu+\nu_Ax^2+\nu_Ws^2
 -\nu_W{q\over N}s(a_0+a_1x)},\qquad
 \nu_W={3N\over5D}{k^2\over E}={k\over G}f(q),
\]
\[
 a_0={.7+.65r\over B_r},\qquad
 a_1={.7-.05r-2r^2\over B_r}.
\]

This follows from the exact planar WJ formula
a_W=beta1 S+beta4(S Omega-Omega S), with
beta1=-6N/(5D), beta4=-6/(5D), and the physical-root differential

\[
 \delta N={q\tau[2.7N x+2(N-c')]\over
                 3N^2-2c'N-.35q^2},\qquad
 \delta I_S=q\tau x,\quad\delta I_\Omega=q\tau.
\]

The raw and deviatoric strain conventions agree on this direction because its
trace is zero. It is an admissible velocity polarization and wave normal for
the supported two-dimensional x-z application. Pressure cannot remove it:
the polarization is transverse to the wave normal.

There is a simple analytic certificate needing only one fixed angle, avoiding
angular-sampling uncertainty. The exact planar shear function has f(q)<=3/10.
For q>=18sqrt(2), the cubic implies N>=18 and r<=1/10. On 0<=r<=1/10,

\[
 a_0\ge17/21,\quad a_1\ge5/7,\quad {q^2\over D}\ge2/3.
\]

The first two inequalities follow from monotonicity: the derivative numerator
of a0 is -1.4+2.8r+1.3r^2<0 on this interval, and that of a1 is
-1.89-5.4r^2<0. At theta=30 degrees, x=1/2 and s=sqrt(3)/2, hence

\[
 \boxed{a(\pi/6)\le\nu+\nu_A/4+
 {k\over G}\left({9\over40}-{7\sqrt3\over30}\right).}
\]

The reported reference state satisfies the q condition. A negative value of
this upper bound certifies a negative in-plane principal direction even if
the reference's scalar diffusivity, covariance and positive split mu pass.
The independent verifier prints the bound and cross-checks the exact planar
quadratic against three finite-difference steps of the separate full 3D
tensor implementation.

### Why changing only positive inner diffusivity cannot fix this family

This conclusion can be made stronger for the **retained inlet/startup
profiles**, without assuming the reference's own ground k remains unchanged.
Let d=nu+nu_A(0)>0, a=kappa sqrt(K/q_e), and define

\[
 J_P(r)={1\over A}\int_0^r P_A\,ds
 =\log(1+ar/d)+{\nu\over d+ar}-{\nu\over d}.
\]

On chi=0 the reference G=t/(d+ar) and P_A are independent of local k, E and
the choice of positive scalar diffusion. Suppose c<.4H. In the tail collar,
E=A/r. Starting from k(H)=K,Q(H)=0, the physical increasing-shear branch
implies G<E/t whenever k>=K: at G=E/t the auxiliary and WJ shear contributions
are each at least t, and molecular traction is strictly positive. If the
comparison q exceeds the bounded branch, the branch endpoint itself gives
G<E/t. Therefore P=(t-nu G)G<E, Q<0 on backward integration and k>=K. This
barrier closes consistently for any D_kw>0. Thus

\[
 {1\over A}\int_0^H P\,dr
 <J_P(.4H)+\log(2.5)=B_{\max},\qquad
 {c\over H}\ge\exp(1-B_{\max})=:\gamma_{\min}.
\]

If c>=.4H the same lower bound remains true whenever gamma_min<=.4. It
therefore covers the entire c<=H family with the stated hypotheses. This is
not the obsolete claim that general new production is cutoff-independent;
it uses the actual chi=0 lower segment and the proved collar inequality.

The retained inlet has a fixed ground k_in(0), reconstructed independently by
the historical integral profile, and fixed G_in(0)=t/d. Consequently its new
dependent dissipation gives

\[
 q_{in}(0)={k_{in}(0)G_{in}(0)c\over A}
 \ge {k_{in}(0)G_{in}(0)H\gamma_{\min}\over A}.
\]

If this lower bound exceeds 18sqrt(2) and the boxed principal upper bound is
negative for k_in(0),G_in(0), **every ordinary positive inner-diffusion choice
in the zero-reference-flux family has an inadmissible retained inlet-derived
initial state**. No fitted change of D can repair it. The verifier checks
these inequalities with the fixed calibrated inputs and reports all values.
At the actual initial profiles, the mean flow is parallel; a developed
boundary solution could have additional normal-to-inlet derivatives, so this
argument specifically certifies the required inlet-derived startup field,
not every conceivable developed boundary gradient.

The existing named bounded-destruction amendment is inactive where local
k>=theta K with theta<1, including the retained inlet and this healthy
reference. When it reduces actual dependent E by lambda<=1, it increases q
at this frozen state and cannot invalidate the large-q certificate. It
cannot be used as an undeclared increase of E, source replacement or stress
regularization. Redefining the source or stress to manufacture a pass would
violate task-package lines105–107 and249.

### Bounded verification and honest next step

After synchronizing the updated oracle, execute under the parent's Linux
runner (one process; all library threads1):

```sh
python3 <source>/solvers/IncNavierStokesSolver/Tests/RANS/h27/reference_oracle.py \
  --principal-audit-input <R0>/reference-kin-v1.json \
  --output <fresh-run>/planar-principal-obstruction-v1.json
```

This performs no new column solve. It reconstructs the retained inlet ground
k independently, prints analytic and full-tensor FD witnesses at six planar
angles, minimizes the analytic planar form, and evaluates the uniform
diffusivity bound. A negative certificate is a scientific failure of the
selected tensor/reference/required-startup combination. Larger mu, smaller
dt, more accurate assembly, extra implicit iterations or stress-remainder
freezing cannot change this continuous principal symbol.

A potentially compliant remaining **reference** investigation is an independent
finite-nu outer BVP supplying a nonzero Q_ref(H), with the endpoint budget term
retained and A(K), inlet and roughness unchanged. It must supply the flux on
scientific boundary grounds; tuning Q_ref from instantaneous actual production
would violate the package. Positive Q_ref can algebraically allow a smaller
cutoff and larger plateau E, but existence, sign, magnitude and admissibility
of such an outer reference have not been established. This alternative must
be versioned and independently tested, not assumed to repair the issue.

If no compatible reference branch repairs the required startup, the minimal
honest closure is H27_MODEL_AND_REFERENCE=NOT_MET and TASK3_READY=NO, with the
negative planar certificate preserved. A change to the lower wall tensor
(for example replacing its adverse WJ normal-stress response) would be a
physical-model amendment beyond the specified tensor and must be presented
for separate scientific review; it is not authorized as a silent numerical
fix. No such change is implemented here.
