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
explicitly for the bounded kernel suite; no uniform margin for unbounded
strain/depletion is asserted.

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
