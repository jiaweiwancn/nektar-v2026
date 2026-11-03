# M8-R2 restricted native prototype contract, version 2

This contract authorizes a research implementation and a bounded implicit-stage
test. M8 remains open. It separates constitutive feasibility, native weak
operator correctness, convergence of the unchanged coupled residual, and
physical applicability of a converged stage. None implies the next.

The reviewed source is `5186d63bd825aadd198488d93cff02fc9496a232`.
The R1 [candidate contract](../m8-r1/candidate-contract.md) is retained byte for
byte, SHA256
`c635154db92a07be5bf27b5b50817340e9d6b5301d778c6a1099ac1d220c8879`.
Its pre-screen contract has SHA256
`082a791fb086a08bd83285e68ae5e05923fdb688c7919410263bb894a6f22313`.
This version corrects the molecular-viscosity and dimensional-gamma notation
in R1 and specifies the native prototype and solver domains. Historical
failures, certificate bounds and producing source/builds are not revised.

## Physical stresses and the matching condition

Let `G=grad(U)`, `S=G+G^T`, `c=2k/3`, and use mechanical pressures `p_w,p_o`.
The selected molecular viscosities are **nu_w=nu and nu_o=0**. The high-Re
outer model acquires no molecular diffusion. In particular,

\[
 T_w=-p_w I/\rho_f+\nu_w S_w-R_w,\qquad
 T_o=-p_o I/\rho_f+\nu_o S_o-R_o.
\]

Writing `R_w=c_w I-F` and `R_o=c_o I+B_o`, the reduced pressures are
`pi_i=p_i/rho_f+c_i`. Thus

\[
 T_w=-\pi_w I+\nu_w S_w+F,\qquad
 T_o=-\pi_o I+\nu_o S_o-B_o.
\]

All components of `(T_w-T_o)n=0` are required at the physical interface.
The mechanical normal-pressure relation is

\[
 \frac{p_w-p_o}{\rho_f}
 =n\cdot(\nu_w S_w-\nu_o S_o-R_w+R_o)n
 =n\cdot(\nu_w S_w-\nu_o S_o+F+B_o)n-(c_w-c_o).
\]

The common continuous k field gives `c_w=c_o` at matching points; it does
not imply mechanical-pressure continuity. No `nu*S_o` term is inserted by
using reduced pressure. Normal and tangential traction, pressure work and
the isotropic production term use these same conventions.

## Fixed outer law and applicability region

Retain beta=.4, the original CLS `C(q)`, GammaUniform=0 and
GammaScale=1. The dimensional coefficient is

\[
 \gamma(x)=\mathrm{GammaScale}\,\nu_{t0}(x),\qquad
 \nu_{t0}=C_{\mu,0}K_0^2/\epsilon_0=\kappa u_*r,
 \quad K_0=Q_*u_*^2,\quad\epsilon_0=u_*^3/(\kappa r).
\]

GammaScale is dimensionless; gamma has viscosity units. The selected outer
stress is `R_o=R_CLS^beta-gamma(x)*(S-S_ref(x))`, with the prescribed
incoming Cartesian flat-reference extension `S_ref`. It is independent of
the current state and current curved-wall normal. At fixed geometry the
state derivative has `delta gamma=delta S_ref=0`. In dimensional production,

\[
 P_o=-R_o:G=P_\beta+\gamma(S-S_{ref}):G
     =-c\,\mathrm{div}(U)-B_o:G.
\]

Mean turbulent work cancels this same-stress production algebraically.
This identity is not a production-sign or discrete-energy theorem.

The R1 restricted box is unchanged. In the incoming orthonormal basis,
`A*=Q^T G Q/G0`, `G0=u*/(kappa*r)`, require

- every entry of `A*-e_x tensor e_z` in `[-.001,.001]`;
- local `k/K0` in `[.99,1.01]` and actual lifted `epsilon/epsilon0` in
  `[.98,1.02]`, including the dependent matching trace;
- both limiting max-invariant branches, positive sampled K and
  `lambda>=.5` for the capped source.

The existing R1 interval proof certifies this frozen-state outer region:
momentum principal form at least
`.53345977306365966796875*nu_t0` and strict positive covariance. It does not
certify state invariance, a complete coupled Jacobian, or a time map. The
inner bounded kernel controls use `K/K0 in [.5,1.5]`; a finite-set check of
that range is explicitly reported and is not inferred from local k.
The rejected broad boxes and invariant alternatives retain their failures.

## Selected inner stress

The new explicit selector is `RANSWallInnerStress=Linear|SaturatedR2`.
Linear remains the default and retains `F=nuA*S`. SaturatedR2 uses the same
physical auxiliary viscosity `nuA(K,r)>=0`, local transported k, and full
symmetric strain. Fix rho=.8, d=.2, `s=||S||F`, `a=nuA*s`, `c=2k/3`:

\[
 t(a,c)=\begin{cases}a,&a\leq\rho c,\\
 c-dc\exp[-(a-\rho c)/(dc)],&a>\rho c,
 \end{cases}\qquad F=tS/s,\quad R_w=cI-F.
\]

At `s=0` set F=0. At `nuA=0` set F=0. At `k=0` set F=0 and Rw=0;
negative k is invalid. No field floor, eigenvalue clipping, post-step
repair, or undeclared deviatoric projection is introduced. Covariance is
evaluated at exact inner k=0 for either selector instead of being skipped.
The Linear law may then fail covariance; zero is not silently accepted as
its covariance matrix.

For positive c on the active branch, `z=(a-rho*c)/(d*c)` and `q=exp(-z)`:
`t_a=q`, `t_c=1-q*(1+d*z)`. Inactive values are `t_a=1,t_c=0`.
For s>0, `h=t/s`,

\[
 \delta F=h\delta S+(\nu_A t_a-h)
 \frac{S:\delta S}{s^2}S+t_a\delta\nu_A S
 +\frac{t_c}{s}\frac23\delta k\,S.
\]

The sampled `delta K`, local `delta k` and the actual one-sided roughness
and positive-part branches in `delta nuA` remain distinct. At fixed
nonzero S and positive nuA, the one-sided k derivative at k=0 is
`(2/3)*delta k*S/s`; it differs at nuA=0. The joint k/strain origin is
not assigned a nonexistent full Frechet derivative. The zero-strain action
at k>0 is `nuA*delta S`. Branch-aware verification must state any narrower
evaluation domain required by the implementation.

At fixed k and nuA this radial map is monotone. The full physical inner
principal form is bounded below by the actual molecular nu. For any
symmetric strain, `lambda_min(R_w)>=c-t`; for trace-free strain its stronger
bound is `c*(1-sqrt(2/3))`. Neither statement certifies the coupled Newton
matrix or a fully discrete positive map. No strict margin is claimed as
k tends to zero or as the turbulent radial tangent tends to zero.

Use this same selected F in every actual inner momentum load, full normal
and tangential traction, pressure bookkeeping, production, mean/turbulence
energy, covariance, coupled directional action and condensation. With weak
discrete divergence, the actual production is

\[
 P_w=-R_w:G=F:G-c\,\mathrm{div}(U)=ts/2-c\,\mathrm{div}(U).
\]

The isotropic divergence term is retained. The existing native patch already
solves the represented velocity/pressure/k system nonlinearly. Its actual
callers must be audited before any architectural change; active stress
cannot be reconstructed by treating old analytical linear columns as exact.
The scalar diffusion coefficients retain their current molecular/nuA/nu_t
definitions. In particular, F/s is not substituted for nuA in scalar
diffusion. Unsupported entry points must reject the selector explicitly.

## Reference, capped source and dependent trace

The analytic neutral reference is inactive because
`a0/c0 < 3*sqrt(2)/(2*Qstar) < .664850004871 < rho`.
Its roughness-dependent auxiliary/reference profiles remain unchanged.
The projected discrete initial field may activate saturation: its a/c,
active volume and stress change are measured and are not labeled an exact
discrete equilibrium.

Retain the existing preprofiles
`Epre=(1-chi)*EA+chi*EC`, `Dpre=(1-chi)*EA+chi*pA0` and their metric ray
budgets. The same capped normalization gives
`D+=Dpre*w/m`, `E+=lambda*Epre`, `C+=E+-D+` and `Sk=P_actual-D+`.
Theta=.2 and delta=.1 remain fixed. The integral budgets and actual trace
`epsilon(H)=lambda*Epre(H)` are checked, including chi=0. Finite quadrature
mismatch is reported, not silently renormalized. A conservative flux and
its equivalent redistribution source are never applied simultaneously.
Local k, sampled K, chi, moving cutoff, normalization, sensor and projected
trace dependencies remain in the actual coupled actions. F has no direct
chi dependence; pressure acts through the existing sensor/profile/source/
trace/constraint path as well as its direct affine terms.

## Scalar-domain policy and physical acceptance

The new opt-in `RANSWallNonlinearPolicy=ScalarFeasibleR2` keeps Legacy as the
default. Both evaluate the same true coupled residual and use the original
Newton acceptance tolerance. The scalar-feasible policy starts from an
admissible initial state, rather than a terminal R1 Newton iterate.

The common scalar-domain inventory is declared in each report. It includes
actual source, mass and prescribed-load quadrature; physical samplers,
interfaces and transfers; native GLL8 and declared dense reconstruction;
and **complete normalization rays for both volume sources and dependent
matching constraints**, including their baseline/current profiles and
actual trial-dependent cutoffs. Finite-set positivity is not positivity
everywhere inside a polynomial element.

Mathematically the inner extension permits k>=0, while the outer high-Re
fields require k>0 and epsilon>0 and all physical samples require K>0.
Actual lifted epsilon, not an unlifted auxiliary field, is checked. Any
stricter numerical boundary rule must be explicit in the policy receipt.
Trial-dependent points, normalization and dependent trace are rebuilt and
checked before residual evaluation and before internal acceptance. Check
the represented k/K values before a dependent-trace exception can obscure
the limiting point. A fraction-to-boundary or equally explicit domain
globalization changes the step, not the field or residual. Use admissible
one-sided actions at exact domain boundaries.

Report Jacobian rebuild policy, independent actions, step limits, limiting
locations/sets, smallest scalar margins, branch changes and all true block
norms. Scalar-feasible internal Newton states are not physical trajectory
samples. A line-search collapse, stagnation, or nonzero residual minimizer
is a failed stage search, not an accepted stage or proof of no root.
Pressure actions retain the full nonlinear sensor/trace path while checking
direct affine pressure, continuity and gauge contributions exactly. The
existing independent `1e-7` action threshold is not relaxed.

`R2LocalV1` denotes the proposed strict converged-stage physical guard:
the unchanged R1 outer box, positive sampled K, actual dependent EPS,
lambda>=.5 and covariance at all declared sets, including inner k=0.
ScalarDomain is the internal evaluation-domain check, not that physical
guard. A comparison stage probe collects post-convergence physical failures
without discarding the convergence evidence; such a state is reported as
converged but unqualified. A native qualified stage may not advance history
through a physical exit. The selector, physical law version and applicable
checkpoint signature are recorded; cross-law restart is rejected. The
nonlinear policy/version is recorded separately from the physical law.

## Frozen representation and experiment bounds

All physical inputs remain those of the R1 improved control: H=2, roughness
.1, sample/tail height1, physical samplers2.5/3.5, original transition,
ustar, viscosity, direction, beta/gamma/source parameters and eta=0.
Keep X=`[-4,-3,-2,-1,1,2]`, all inner knots
`[0,.01,.1,.25,.5,1,2]`, p3 velocity/k/epsilon, p1 pressure, native GLL8,
and the 3000-unknown guard. At most two targeted changes are permitted:

| Order | Outer r knots | Field sizes u,v,k,epsilon,piw,pio | Total with gauge |
| --- | --- | --- | --- |
| 1 | 2+j/3, j=0..6 | 592,592,592,304,42,42 | 2165 |
| 2, only if needed | 2+2j/9, j=0..9 | 736,736,736,448,42,60 | 2759 |

Base source quadrature is tangential8/normal16; high is12/24. GLL8 stays
fixed. Reconstruct both one-sided gradients at internal outer interfaces,
including the original r=8/3 witness. The expected cubic-derivative error
factors1/8 and1/27 are planning estimates, not measured bounds. Relative
to R1's1571 unknowns, dense-factor estimates are2.62 and5.42, while FD
volume-work estimates are1.84 and2.93. No open mesh search, inner refinement,
box enlargement, backend rewrite or dense-guard change is authorized.

Each preflight reports actual DOFs, quadratures, raw **unanchored** mass-dual
field norms, complete scalar/covariance/outer/dependent-trace margins, inner
a/c and saturation physical volumes/locations, and same-state independent
quadrature sensitivity. Positive finite coverage is not a continuum proof.
Failure after these two choices is an explicit coverage/capacity blocker.

After operator/action verification and preflight, freeze two matched
first-stage problems: Linear and SaturatedR2, both ScalarFeasibleR2, same
represented initial state and true constraints, eta=0, no forcing/anchor,
actual startup dt=.0025 and fixed target time t*=.0025. Let

\[
 R_\alpha(x)=\alpha M(x-h_0)+L(x,t_*).
\]

Here h0 is the identical actual initial history for every level and both
laws. The only permitted continuation is alpha=`[3200,1600,800,400]`, at
most four levels per law. Alpha400 is the original target stage. Previous
level solutions may initialize the next guess but never replace h0 or
advance physical time/history. These are auxiliary stage-mass solves,
not time-step refinement. The frozen execution protocol sets an explicit
Newton limit no greater than16 per level before execution; no failed-level
retry or extra level is implicit. Verify the final original target residual.

Only if the new-law target converges and passes every scalar, covariance,
outer, trace, source/energy/continuity/traction/constraint gate may the two
short independent dt/dt2 native controls run, horizon at most four coarse
steps and first physical-exit stop. If both target roots are unresolved,
replace that short-run allocation with at most one frozen-velocity/pressure
scalar-stage localization, clearly labeled diagnostic. Bounded failure is
not proof of root absence. Genuine curved MMS, active-cap curved/variable-H
verification, full discrete positivity, long-time evolution, scalability
and paired terrain remain open.

Every numerical receipt binds the producing source, binary/runtime,
selected tool imports, XML, command/environment, results and histories.
Execution is Linux through Tailscale with explicit worker allocations and
numerical threads one. Contract preparation itself invokes no solver.
