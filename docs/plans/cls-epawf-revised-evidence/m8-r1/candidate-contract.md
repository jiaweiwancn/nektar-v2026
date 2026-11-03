# M8-R1 phase D — one bounded next-model contract

## Material Passport

- Authority: the bounded M8-R1 user request; baseline `ab5c3c6fee19761e8d88692e8e6c331e5889234e`.
- Scope: kernel feasibility and an explicit applicability decision only.
- Production adoption: **none**. No C++ closure, solver or physical input is changed.
- Numerical evaluation: Linux through Tailscale; one worker on CPU 50; numerical threads one.
- Status at derivation: formal candidate specified **before** its numerical screen.
- Original M8 broad box and all failed witnesses remain authoritative failures.

## Inner completion specified before evaluation

Use the existing physical auxiliary viscosity `nuA=nuA(K,r)>=0`, the local
transported `k>=0`, and the full symmetric strain `S=grad(U)+grad(U)^T`.
Let `s=||S||F`, `c=2k/3`, `a=nuA*s`, `rho=4/5`, and `d=1-rho=1/5`.
For positive k define

\[
 t(a,c)=\begin{cases}
 a,&a\leq\rho c,\\
 c-dc\exp[-(a-\rho c)/(dc)],&a>\rho c,
 \end{cases}\qquad
 F_k(S)=\begin{cases}t(a,c)S/s,&s>0,\\0,&s=0.\end{cases}
\]

The candidate Reynolds and physical mean stresses are

\[
 R_w=cI-F_k(S),\qquad
 T_w=-p_wI/\rho_f+\nu S-R_w
     =-\pi_wI+\nu S+F_k(S),\quad
 \pi_w=p_w/\rho_f+c.
\]

This is a constitutive replacement for the inner full-strain completion. It
is not eigenvalue clipping, a field floor, an iterative viscosity split or a
change to the scalar AWF reference profiles. Negative k remains invalid.

### Complete differential and monotonicity

On the active branch set `z=(a-rho*c)/(d*c)>0` and `q=exp(-z)`. Direct
differentiation gives

\[
 t_a=q,\qquad t_c=1-q(1+dz).
\]

On the inactive branch `t_a=1,t_c=0`. Both derivatives agree at `z=0`, so
the continuation is C1 there. It is generally not C2. Since
`exp(-z)*(1+d*z)` decreases from one to zero for `z>=0`, one has
`0<t_a<=1`, `0<=t_c<1` at positive k and finite a. Also `0<=t<=a` and
`t<c`; the active equality gap is `c-t=d*c*exp(-z)`.

For `s>0`, `h=t/s`, and arbitrary admissible changes of all inputs,

\[
 \begin{aligned}
 \delta c&=\tfrac23\delta k,&
 \delta s&=(S:\delta S)/s,&
 \delta a&=s\,\delta\nu_A+\nu_A\delta s,\\
 \delta F&=h\delta S+
 (\nu_A t_a-h)\frac{S:\delta S}{s^2}S
 +t_a\,\delta\nu_A S+\frac{t_c}{s}\delta c\,S,\\
 \delta R_w&=\delta c\,I-\delta F.
 \end{aligned}
\]

The fixed-k/fixed-nuA derivative on symmetric tensors has radial eigenvalue
`nuA*t_a` and tangential eigenvalue `h`. Both are positive if `nuA>0`;
`nuA*t_a<=h<=nuA`. Consequently F is the gradient of the convex radial
potential `integral_0^s t(nuA*r,c) dr` at fixed positive k and nuA. The full
inner mean stress includes the actual molecular term. For every
incompressible rank-one gradient `H=u tensor xi`,

\[
 H:D_A[\nu S+F][H]\geq
 (\nu+\nu_A t_a)|u|^2|\xi|^2\geq\nu|u|^2|\xi|^2.
\]

This is an inner frozen-state principal result, not positivity of the whole
coupled Newton matrix or a fully discrete turbulence theorem. Cross-field
and sampled-K derivatives above must remain in the actual coupled action.

### Covariance and bounded margins

For any symmetric S, `||F||2<=||F||F=t`, hence

\[
 \lambda_{\min}(R_w)\geq c-t>0\quad(k>0,\ a<\infty).
\]

For the physical trace-free 3-by-3 strain there is a stronger uniform bound.
If its largest eigenvalue is l, the two other eigenvalues sum to -l and
their squares sum to at least `l^2/2`. Thus `l<=sqrt(2/3)*s`, and

\[
 \lambda_{\min}(R_w)\geq c-\sqrt{2/3}\,t
 \geq c[1-\sqrt{2/3}]>0\quad(k>0).
\]

This covariance margin is uniform in strain even though the turbulent
radial tangent can tend to zero. For a controlled weak-divergence ratio
`abs(trace(S))/s<=eta<1`, the corresponding coefficient is
`eta/3+sqrt(2/3)*sqrt(1-eta^2/3)<1`. No such bound is silently inferred
from a weak pressure residual; the actual trace/strain ratio must be measured.

If `a/c<=L` and `k>=kmin>0`, a uniform sufficient covariance margin is
`d*(2*kmin/3)*exp[-max(0,(L-rho)/d)]`. A turbulent radial tangent margin
also requires a positive lower bound for nuA; without that bound the
physical molecular viscosity still supplies the inner mean-stress margin.
The uniform margins and any finite acceptance threshold must be stated
explicitly for the bounded kernel suite. No uniform positive turbulent
principal margin under unbounded strain, or strict covariance margin as k
approaches zero, is asserted. The trace-free covariance bound above remains
uniform in strain when k has a positive lower bound.

For exactly incompressible S the trace is `trace(R_w)=2k`. With weak discrete
incompressibility the full S is retained and its nonzero trace remains a
reported constraint/trace defect. No undeclared deviatoric projection is used.
The covariance norm bound itself remains valid for any symmetric S.

### Sampled K and auxiliary-viscosity derivatives

K remains the selected physical sampler of the represented k field, so
`delta K` is that same sampler applied to `delta k_h`. It is independent of
the local delta k at a generic quadrature point. The original AWF viscosity is

\[
 \nu_A=b_A\sqrt K\,[r-z_v(K)]_+,
 \quad b_A=c_{\mu,A}c_{l,A},\quad
 z_v=\nu b(K)/\sqrt K.
\]

On the positive-viscosity branch,

\[
 \partial_K\nu_A=b_A\{r/(2\sqrt K)-\nu b'(K)\},\qquad
 \partial_r\nu_A=b_A\sqrt K.
\]

For a smooth wall `b(K)=10.7`, so `b'=0`. For roughness h>0 set
`x=h*sqrt(K)/(70*nu)`,
`p=max(.5-.4*x^.7,1-.79*x^(-.28))`, and `b=10.7*(1-x^p)`. Then

\[
 x_K=x/(2K),\quad
 p_{1,x}=-.28x^{-.3},\quad p_{2,x}=.2212x^{-1.28},\quad
 b_K=-10.7x^p\{p_K\log x+p/(2K)\}.
\]

Use the actual active or one-sided max branch. At `r=z_v`, the positive-part
directional derivative is `b_A*sqrt(K)*max(delta r-z_v'(K)*delta K,0)`;
there is no unique two-sided Jacobian. Geometry changes and map derivatives
also change r and the physical strain and must be included if tested.

### Zero inputs and the joint origin

- At `S=0,k>0`, an open inactive neighborhood gives `F=nuA*S` and
  `D_S F=nuA*I`; the k and nuA partial derivatives of F vanish at S=0.
- At `nuA=0,k>0`, `F=0`; the admissible one-sided nuA action is
  `delta nuA*S`. The mean-stress molecular response is retained.
- At `k=0`, define `F=0,R_w=0` for all finite S and nuA. For fixed nonzero S
  and strictly positive nuA, the one-sided k action tends to
  `(2/3)*delta k*S/s` in F. At fixed nuA=0 this k action is zero. There is no
  strictly positive covariance margin at k=0.
- At the joint origin `(k,S)=(0,0)` with nuA positive, the directional
  derivative depends nonlinearly on the ratio of the incoming k and strain
  directions. No full Fréchet Jacobian is claimed there. The analogous joint
  k/nuA origin also requires one-sided treatment.
- The high-Re outer equations remain singular if their transported k or
  epsilon vanishes. This inner extension does not supply an outer continuation.

## Exact open preservation of the flat neutral reference

At the neutral reference `k=K0=Qstar*u_star^2` and
`U0'=u_star^2/(nu+nuA)`. Therefore

\[
 \frac{a_0}{c_0}
 =\frac{3\sqrt2}{2Q_*}\frac{\nu_A}{\nu+\nu_A}
 <\frac{3\sqrt2}{2Q_*}<\rho=0.8.
\]

This strict bound is uniform over the finite wall layer and is independent
of the roughness branch. Continuity supplies an open all-ray inactive
neighborhood. There the complete inner stress, production and every
existing directional derivative agree with the original completion,
including auxiliary-viscosity one-sided branches. The neutral AWF profiles,
cutoff, matching value and reference pressure offset remain unchanged.

## Production, capped budgets, actual trace and chi=0

Use precisely the same Rw in momentum, production and total traction:

\[
 P_w=-R_w:\nabla U=ts/2-c\nabla\cdot U.
\]

The divergence term must remain in a weakly incompressible implementation.
With exact incompressibility `P_w=ts/2>=0`; at k=0 the actual production
vanishes. The mean-flow turbulent work and the k-production term cancel
algebraically in total energy using this same stress. Physical molecular
strain dissipation remains `nu*S:S/2`; prescribed roughness-force work is
still reported separately and is not assumed nonpositive.

The pre-repair reference profiles may remain
`Epre=(1-chi)*EA+chi*EC`, `Dpre=(1-chi)*EA+chi*pA0`, with equal metric ray
budgets B. Every counterfactual flat reference is inactive by the same
strict ratio above, so pA0 and its reference-profile derivatives remain the
original prescribed quantities. Apply the existing capped normalization
to the actual local k/K, then

\[
 D_+=D_{pre}w/m,\quad E_+=\lambda E_{pre},\quad
 C_+=E_+-D_+,\quad S_{k,w}=P_w-D_+.
\]

Thus `integral(J*C+) = 0`, `integral(J*D+)=integral(J*E+)=lambda*B`, and
the dependent outer trace is the actual `epsilon(H)=lambda*Epre(H)`.
The finite existing quadrature/budget mismatch must still be measured;
it is not silently renormalized away. In the conservative form the
redistribution flux has zero endpoints; use it or the source representation,
not both. All K, local k, chi, cutoff, normalization, sensor and projected
trace derivatives remain coupled to the stage.

At chi=0 the same new inner completion and source repair define the adapted
baseline. They do not recover the unrestricted original scalar AWF response.
Once saturation is active, parallel shear has nonlinear traction
`tau=nu*U'+sign(U')*t(sqrt(2)*nuA*abs(U'),c)/sqrt(2)` with positive slope
`nu+nuA*t_a`. A nonlinear inner solve replaces the original off-reference
analytical AWF elimination. Keeping the old closed columns as an exact
active-branch reconstruction would be inconsistent.

## One outer applicability decision fixed before the kernel screen

The next decision is option (a): combine the inner completion above
with the existing beta=.4, GammaScale=1 outer law **only on an explicitly
stated near-reference region**. This is a new restricted research contract,
not repair of the rejected M8 broad box. Original, CubicScaled, Regularized
and both failed invariant-response alternatives and their exact witnesses
remain unchanged in the archive.
In particular, the historical broad box still means gradient-entry radius
`.25` with both k/K0 and epsilon/epsilon0 in `[.5,1.5]`; its incompatible
necessary fixed-Gamma bounds remain a failure. The restricted box below is
a separate applicability proposal, not a revised definition of that domain.

The completed original t0 ledger gives `epsilon/epsilon0` in
`[.9833429223226836,1.003851648978953]`, k/K0 equal to one to roundoff,
and maximum normalized gradient-entry departure `.018332328258716113`.
Consequently the original initial gradient is outside the region below;
the scalar range also exceeds the old +/-1% certificate. This is measured
initial representation coverage, not a claim of a constitutive failure.

The single R1 outer region is now fixed as

- beta=.4, GammaUniform=0, GammaScale=1, original C(q), no outer molecular term;
- entrywise `abs(A*-e_x tensor e_z)<=.001` in the incoming orthonormal frame;
- `k/K0 in [.99,1.01]` and `epsilon/epsilon0 in [.98,1.02]` separately;
- both limiting max-invariant branches and every incompressible rank-one direction;
- target momentum principal margin `.1*nu_t0` and positive covariance;
- fixed proof algorithm: 64 scalar k/epsilon-ratio partitions, exact rational
  polynomial cancellation, directed coefficient intervals and null-form LDL.

The expanded epsilon interval is selected from the actual archived scalar
range before evaluation. No gamma/beta retuning or open-ended box search is
permitted. If this proof is inconclusive, record that failed proof attempt
and leave the contract unqualified. The improved phase C representation
must independently fit the fixed region; the region is not moved to fit it.
A passing interval proof is separate from finite representation coverage
and from forward invariance of the PDE state set. Every retained exact
historical witness is evaluated and explicitly classified; witnesses outside
this restricted region are retained, not reclassified as passing broad tests.

Applicability requires `lambda>=1/2`, positive sampled K (with the bounded
kernel controls using `K/K0 in [.5,1.5]`), the actual positive dependent trace,
and actual outer epsilon/volume/transfer values in the region above. The
lambda bound is a declared stop condition, not a replacement field value.
It covers the archived prefix's measured lambda=1; no causal epsilon
attenuation mechanism is inferred from source limiting alone. The contract
exits before lambda approaches zero or the original high-Re outer equations
become singular. No epsilon floor, clipping or undeclared continuation is
allowed; what happens beyond exit remains unresolved.

No kernel numerical screen or candidate adoption is implied by this formal
derivation. The eventual fixed witness suite, failures, exact source, command,
results and final one-region decision are appended after Linux evaluation.

## One completed Linux kernel screen

The pre-screen contract is retained verbatim in
[kernel-contract-preregistered.md](kernel-contract-preregistered.md), SHA256
`082a791fb086a08bd83285e68ae5e05923fdb688c7919410263bb894a6f22313`.
The source driver and the 64-partition plan were frozen before this one run.
Root executed the suite through Tailscale on CPU 50 with all six numerical
thread settings equal to one. It ran no solver or build. The complete
[validation receipt](kernel-validation.json), [runner](kernel-runner.json),
[log](kernel-run.log) and compact [summary](kernel-summary.json) record the
command and exact input hashes. All mathematical imports came from the clean,
unchanged `ab5c3c6fee19761e8d88692e8e6c331e5889234e` archive; the new driver is
bound separately by SHA256, rather than attributed to the historical binary.

The run passed 2670 verification checks with zero failures in 58.40 seconds.
The **one fixed outer region** passed its rigorous interval certificate:
the full incompressible momentum principal form is at least
`0.53345977306365966796875*nu_t0`, exceeding the requested `.1*nu_t0`.
Both complete max-invariant derivative branches are enclosed, so this includes
active branches and limiting tie combinations. Interval LDL proves strict
positive covariance throughout this same finite state box. Its positive
pivots are factorization certificates, not covariance eigenvalue bounds.
The complete 64-leaf bound/LDL records and small covariance intervals are in
[kernel-outer-proof.json](kernel-outer-proof.json); the full interval matrices
remain externally stored and SHA-bound by [kernel-curation-index.json](kernel-curation-index.json).
No molecular outer viscosity was added to obtain this principal margin.
The certified principal derivative freezes the positive k and epsilon while
differentiating the entire stress with respect to all nine gradient entries,
including C'(q), the quadratic/cubic factors, and both formal invariant
branches. Separate analytic and finite-difference checks include the k and
epsilon coefficient actions. Those cross-field checks are not a coercivity
certificate for the coupled momentum/k/epsilon Jacobian.

For the outer law the exact same-stress production relation is

\[
 P_{reg}=-R_{reg}:\nabla U
 =P_\beta+\Gamma(S-S_{ref}):\nabla U.
\]

Its opposite is the turbulent stress work in the mean-energy equation, so
the local exchange cancels algebraically when the actual same stress and
total traction are used in both equations. Gamma and Sref remain prescribed
spatial data in the state derivative. In particular, write `Rreg=cI+Breg`
(Breg denotes the remaining tensor, without assuming it is trace-free at
nonzero discrete divergence). Then

\[
 P_{reg}=-c\nabla\cdot U-B_{reg}:\nabla U,
 \qquad T_o=-\pi_oI+\nu S-B_{reg},\quad
 \pi_o=p_o/\rho_f+c.
\]

Thus a reduced-pressure implementation must retain the isotropic
`-(2k/3)*div(U)` production term and its matching pressure/stress work.
Mechanical pressure is `p_o/rho_f=pi_o-c` for total-traction bookkeeping;
these are not independent extra forces. Likewise the inner production is
`F:grad(U)-c*div(U)`, rather than `F:grad(U)` away from exact incompressibility.
No weak pressure residual is silently replaced by pointwise zero divergence.
The unchanged outer production contraction is also covered by the earlier M8
[indexed-kernel receipt](../m8/proof-kernel-fresh.json) and the
[archived-stage accounting](archive-ledger.json) for the native
beta=.4/GammaScale1 law; this reuse is not a new native test of the inner
completion. This identity does not imply global
nonnegative production, a discrete energy inequality or positivity of the
outer high-Re equations. The completed numerical outer screen tests
covariance and full state derivatives/principal forms; it does not supply an
additional all-state production-sign certificate. That distinction is retained
instead of reading energy or scalar admissibility out of a momentum margin.

The inner continuum proof has two explicit dimensionless margins. Directed
interval evaluation bounds the neutral `a0/c0` above by
`0.664850004870987993954772952780164787607499586529604475224952`, leaving
more than `.135149995129012` below rho=.8. For every trace-free strain and
positive local k, its covariance is at least
`(2k/3)*0.183503419072273967267571975098036202678017506447776623855768`.
Its physical momentum principal margin is at least the actual molecular nu.
These are frozen-state continuum constitutive statements. The checks of all
local-k, sampled-K, r, strain, roughness, limiter and actual-trace actions
confirm the specified differentiated formulas; they do not turn the full
coupled Newton matrix or a discrete time step into a positive map.

At the **exact first archived bad inner point**, applying the proposed law
only as a kernel evaluation gives minimum covariance
`+0.09997609032554167`; the original stress at that same point has
`-0.0976947699300747`. The candidate is active there (`a/c=1.8201658640223455`).
This does not change the archived trajectory or demonstrate that a native
nonlinear wall solve using the proposed law would converge.

The fixed witness suite retained all six preregistered families: 17 M1
constitutive cases, 72 M8 native kernel cases, 16 exact fixed-Gamma witnesses,
48 invariant-response witnesses, eight invariant analytic controls, and two
outside-envelope controls. The unchanged beta=.4/GammaScale1 target still
has 37 finite certified negative covariance directions in these retained
inputs; 94 tested inputs lie outside the proposed restricted box. Those
failures are retained in [kernel-retained-witnesses.json](kernel-retained-witnesses.json).
Their original models,
counterexamples and broad-box rejection are not overwritten by the local
interval certificate. A finite witness list is not an envelope proof.

## Coverage and the one next-model decision

The recommendation is option (a), **a restricted research contract**: use
this inner completion together with the unchanged beta=.4/GammaScale1 outer
law solely inside the fixed region and actual EPS/trace/positive-K conditions
above. The mathematical kernel has feasibility; production adoption is not
part of this result. There is no claim that the original broad-box model is
repaired, nor a claim of full discrete or global nonlinear admissibility.
This is a recommendation for the next research/implementation contract,
not a currently usable qualified transition model. Its near-reference
physical relevance and kernel margins do not establish that a native initial
representation satisfying the joint contract is presently available.

Neither initial native transition representation presently satisfies this
joint contract. The independent [preflight summary](preflight-summary.json)
and [refined higher-quadrature audit](preflight-resolved-high.json) retain
the actual coefficient and input hashes. The original gradient-entry departure is
`.018332328258716113`. The sole improved representation has
`.0034123386947804146` at both base and higher quadrature, exceeding `.001`.
Its higher-quadrature outer epsilon ratio
`[.9832524731447201,1.0000658783917953]` fits the expanded scalar box, and its
finite inner/outer covariance minima are positive; this does not repair the
failed gradient coverage. Its chi=0 trace ratios
`[.9834528466780685,.9846519510665872]` explain why the old +/-1% scalar
certificate was insufficient. No mesh or state-box adaptation was made after
these measurements.

The coverage plan is therefore an explicit gate, not a claim of completed
qualification: require actual native gradient, local k, sampled K, outer EPS,
dependent EPS trace and both active max branches on every declared physical,
source, transfer, mass and interface set to satisfy the fixed bounds before
a qualified step. Save the state and exit qualification at the first failure;
any continuation is diagnostic. A continuum coverage or forward-invariance
proof remains open. The available refined initial state already fails this
gate, so current trajectory runs remain unqualified from t=0.

An implementation of the inner candidate would have to replace the original
active-branch analytical wall elimination with a nonlinear wall solve,
differentiate its full residual and dependent traces, and use the same F in
momentum, production and **all** traction components. In particular, normal
traction continuity requires

\[
 (p_w-p_o)/\rho_f
 =\nu\,n\cdot(S_w-S_o)n-c+n\cdot F n+n\cdot R_o n.
\]

At fixed geometry the reduced-pressure inner weak stress terms become
`integral[(nu*S+F):grad(v)]-integral[pi_w*div(v)]`, with the actual total
traction at boundaries and interfaces. Their directional action therefore
contains `nu*delta S+delta F` and `delta pi_w`, including all sampled-K and
local-k terms in delta F. The scalar source uses the same actual P and D+;
its dependent trace is differentiated from lambda*Epre(H). These nonlinear
stress/source/trace and interface terms replace the old active-branch
analytical columns; they have not been assembled or verified in a new native
weak operator. Curved metrics, normals, variable H, and any corresponding
geometry derivatives need their own future implementation and MMS audit.

This is a physical pressure jump, not pressure continuity. Reusing the old
off-reference analytical pressure or traction columns as exact after
saturation would be inconsistent. Genuine curved MMS for the adopted law,
active-cap/variable-H coupling, full discrete admissibility, unforced
refinement and long-time evolution, the scalable backend and physical terrain
validation remain open. The lambda>=1/2 guard and actual outer EPS box are
both required; no lambda->0 singular high-Re continuation has been supplied.
M8 remains open.
