# WJ-1 theory-to-code contract

The reviewed parent is `ba0f3068c334bdbe43ce0cadf5ad3ad45340baae` on
`origin/codex/cls-epawf-pimple-one-sweep`. The fetched remote head on 10 November
2026 was exactly that commit. Work uses `codex/wj-earsm-awf` in the isolated Mac
worktree `.worktrees/wj-earsm-awf`. The original Mac checkout and the Linux main
checkout, branches, uncommitted files, installation and running campaigns are
preserved. Linux execution evidence is in `wj1-linux-evidence.md`.

## Reviewed specifications and erratum

Both previously untracked specifications are committed as exact copies; their
original files were not edited. Actual SHA-256 values match the requested review:

- `docs/literature/wj_earsm_neutral_abl_ep_awf_note.tex`:
  `1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124`.
- `docs/literature/wj_earsm_budget_cutoff_awf_appendix.tex`:
  `0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c`.

Main-note `eq:integrals-linear` contains the literal `quad`; this is a typographic
erratum for `\quad`. The intended formula is
`Delta J1 = (z2-z1)/a - (d+a*h)/a^2 * log(y2/y1)`, for `z1>=h`.
`AnalyticalWallFunction::Integrals` already implements its stable, segment-local
form using positive endpoint viscosities, `log1p`, and small-ratio moments.
The specification bytes deliberately retain the erratum.

## Outer constitutive evaluator and shared action

`WJEARSM.{h,cpp}` implements main-note `eq:wj-normalization`, `eq:wj-invariants`,
`eq:wj-cubic`, `eq:wj-basis`, `eq:wj-coefficients`, `eq:wj-ceff`, the outer scalar
equations and `eq:wjstab-root-derivative`. This is the full 3D explicit tensor
with the 2D cubic surrogate on every 3D state, not the planar two-term tensor or
the implicit 3D ARSM. With `tau=k/epsilon`,

```
a = R/k - 2I/3
S = tau sym(G), Omega = tau skew(G)
s = tr(S^2), w = tr(Omega^2) <= 0
(N-9/5)(N^2-2w) = (27/10) s N, N >= 9/5
Ceff = 3N/[5(N^2-2w)]
nu_eff = Ceff k^2/epsilon, nu_d = Cd k^2/epsilon
P = -R:G
Sk = P-epsilon
Se = (Ce1 P-Ce2 epsilon) epsilon/k
Dk = nu_d/sigma_k, De = nu_d/sigma_epsilon
```

The physical root is bracketed using monotonicity of `F/N`; a pure rotation or
zero gradient gives `N=9/5`. Long-double arithmetic precedes checked conversion.
No covariance, production, root or diffusivity clipping is used. Invalid or
unrepresentable values throw. `RANSClosureValues` returns full covariance,
anisotropy, generic deviatoric stress, generic stress remainder, actual raw work,
deviatoric work, scalar sources, separate viscosities and root/invariant/beta
diagnostics. WJ leaves CLS `quadraticStress` and `cubicStress` zero.

`ApplyWJClosureDirection` is the shared action for later native integration.
First-order directional arithmetic differentiates all tensor products, every
coefficient (including the zero-at-reference I4 coefficient), production, sources
and diffusivities. Only the root value is iterated. Its derivative is explicitly

```
FN = (9/5)(N^2-2w)/N + 2N(N-9/5) > 0
dN = [(27/10) N ds + 2(N-9/5) dw]/FN
dnu_d = 2 Cd k/epsilon dk - Cd k^2/epsilon^2 d(epsilon)
```

`WJParameters::Defaults()` is distinct from `RANSParameters`. It uses c1=1.8,
c2=5/9, Ce1=1.44, Ce2=1.82, sigma_k=1, sigma_epsilon=1.3 and
Cd=572/6561. Only the selected c1/c2 law is supported; changes reject. Cd may be
an explicitly supplied positive **constant coefficient**, and scalar transport
coefficients may be supplied explicitly. Neutral calibration then uses their
actual values and requires Ce2>Ce1. The WJ parameter type exposes no Yap,
damping, max-invariant, timescale floor, clipping or CLS-variant knobs. All
nonzero outer molecular-diffusion requests reject. Physical wall nu stays
strictly positive and separate.

The omitted-parameters `EvaluateRANSClosure(WJHighRe,...)` dispatches through
`RANSClosureKernel::WJ()` using WJ defaults. Passing legacy `RANSParameters` to
WJ rejects, even if values resemble WJ coefficients. Both bool branches of
`EvaluateRANSClosureLinearization` reject WJ; the bool is never ignored. Legacy
robust variant entry points likewise reject WJ. Standard and CLS defaults retain
Ce2=1.92 and their existing conventions.

## Off-constraint choice and kernel factory

The strict evaluator rejects a nonzero trace (with a relative roundoff allowance
of 64 double epsilons) in the gradient or its direction. The explicitly named
`DevStrainRawWorkV1` factory completion instead uses

```
theta = tr(Graw)
S = tau [sym(Graw)-theta I/3], Omega = tau skew(Graw)
Pdev = -k a:dev(Graw)
Praw = -R:Graw = Pdev - 2 k theta/3
```

Sources use Praw. This completion preserves tr(R)=2k and actual raw stress work;
it is a choice away from the incompressible manifold and agrees with the notes
on that manifold. Adding lambda I leaves N, a, R and Ceff unchanged and changes
Praw by -2k lambda. Nonnegative-production assertions apply only under the
incompressible hypotheses. This choice does not alter CLS or the AWF full-strain
inner stress, which still uses local transported k and G+G^T.

`RANSClosureKernel::WJ` provides actual compiled factory/adapter dispatch for both
values and directional actions. No session/PDE WJ selection is added. The native
session guard still accepts only StandardKEpsilon and CLSHighRe.

## Neutral reference, embedding and pressure

`RANSNeutralReference::WJ` implements `eq:equilibrium-constants`,
`eq:outer-reference-family`, `eq:kappa-compatibility` and `eq:outer-covariance`:
Ne=81/20, C*=572/6561, qe=81/(2 sqrt(143)),
`kappa^2=sigma_epsilon (Ce2-Ce1) C*^(3/2)/Cd`. Default kappa is
0.381917863546733. `AtHeight` anchors U to the supplied (H,UH), gives constant
k=qe ustar^2, epsilon=ustar^3/(kappa z), and returns full covariance. The epsilon
gradient flux is `-(Cd/C*) ustar^4/(sigma_epsilon z)`; Cd is not implicitly
replaced by stress Ceff. No 0.4/0.09/independently rounded 0.087 enters defaults.

Planar physical velocity (U,W) embeds into Cartesian (x,y,z) as (U,0,W);
gradient indices are physical x and z (0 and 2), retaining all three covariance
rows, especially Ryy=2k/3 at neutral shear. A native mesh's coordinate mapping
must make this embedding explicit.

`WJOuterTraction` returns the full `(-p/rho I-R)n` for a unit physical normal.
`WJInterfacePressureJump` prepares normal matching using
`D_A n.(G_w+G_w^T).n - 2k_w/3 + n.R_o.n` with the actual tensors.
Its flat reference is -K/Ne. Tangential matching must also use full traction.
`WJReferenceMechanicalPressure` on flat WJ profiles returns
`p_o(H)/rho-K/Ne+g(H-z)+2(K-k_ref(z))/3`, so
`d[p_w/rho+gz+2k_ref/3]/dz=0`. Curved profiles reject this flat pressure helper.
`WJMechanicalPressure` explicitly recovers p/rho from stored
`pi=p/rho+2k/3+gravitationalPotential` before sampled forcing is differentiated.
This is a specified storage convention, not an assumption about every native
pressure variable.

## Calibrated finite-nu walls

`WJCalibratedAWFParameters` separates calibration from source policy.
`EquilibriumAWFProfiles::WJ` requires explicit WJCalibratedEPAWF or
WJCalibratedBudgetAWF, supplying c_mu_wall=C* and
c_length_wall=kappa qe^(3/2), alpha=kappa/sqrt(qe). Direct construction with
legacy wall coefficients rejects. Both select endpoint K=k(H). The existing
`AnalyticalWallFunction` supplies unchanged signed vector momentum reconstruction,
roughness correlation, negative zv, smooth h=0 limit and stable piecewise
finite-nu J0/J1/J2/positive production moments. The baseline cutoff is
`2 c_length_wall nu/sqrt(K)`, never the old literal 5.1.

The wall viscosity is `nu_t,A=alpha sqrt(K) (z-zv)+` and
`D_A=nu+nu_t,A`, `D_k,w=nu+nu_t,A/sigma_k`. Production is turbulent only:
`nu_t,A |U'_A|^2`. Interface traction is `t_A+F(H-h)`, preserving reversal and
vector cross terms. The budget uses counterfactual zero-forcing production at
current sampled K; actual matching velocity or forcing never changes the cutoff.
On flat geometry `Pbudget=t0^2(J0-nu J2)` is evaluated through the stable positive
production moment, `Lambda=Pbudget/A`, `zc=H exp(1-Lambda)` and the applicable
`zc<=zm` gate is enforced without clipping. On curved local metrics the analytic
metric dissipation integral is matched to independently integrated J*pA0 by a
strictly monotone logarithmic root and its metric gate. This is a normal-column
kernel, not a complete curved momentum/pressure or terrain-equilibrium solver.

EP uses Epre=(1-chi)EA+chi Ec and Dpre=(1-chi)EA+chi pA0, with the simplified
source PA-Dpre. At chi=1 and matching reference data the arithmetic cancellation
is pointwise, hence valid for every scalar test mode. At chi=0 it constructs the
calibrated baseline directly, without evaluating Ec, its gate or its root. Active
EP fails under RequirePreservation; the explicitly selected UseOriginal policy
returns direct baseline data on a failed cutoff/regime gate. EP's declared
minimumReynolds is a separate applicability input, not outer damping. An
unrepresentable numerical cutoff rejects rather than clipping.

Budget has no chi blend, EP redistribution, Reynolds switch or fallback. It
requires zK=zm=H and its budget gate. The options object's default
minimumReynolds=100 is an unused legacy placeholder; an explicit nondefault value
is rejected for WJ budget. Q0 and k_ref implement appendix `ba:Q`/`ba:k-profile`:
`Q0=int_0^z(Ec-pA0)`, `k_ref=K-int_z^H Q0/D_k,w`. Both Q0 endpoints and all
nonconstant-mode diffusion/source balances are tested. No k(0) condition, local-k
coefficient replacement, reference clipping or evolving-field reset is introduced.
There is one outer epsilon value trace, dependent inner epsilon, no inner epsilon
PDE and no independently prescribed epsilon flux.

## Independent verification and limitations

`WJKernelProbe` runs compiled API/rejection checks and a stream protocol using
production factory/adapter/profile/load APIs. `validate_wj1.py` compares these
with independent NumPy polynomial roots, complex-step differentiation, centered
differences over four steps, and SciPy quadrature split at physical breakpoints.
`wj1_reference_certificate.py` contains the note's complete exact tangent and
rational null-minor/LDL certificate. It proves the all-direction incompressible
principal bound nu_eff/250 with k,epsilon frozen. All 8 trace-free tangent columns
are compared with the compiled action; the note's arbitrary extension outside
that subspace is not mistaken for the named raw-gradient completion.

The required q=3.5 adverse direction remains negative despite positive covariance
and production. The commuting 3D example remains diag(-2/15,-2/15,34/15), N=3,
P/epsilon=48/5. The rough canonical budget profile is positive. Independent and
compiled evaluation also expose a smooth-wall budget reference with
k_ref(0)=-0.330127502951175 for the same K,nu,H and h=0. Its cutoff gate passes;
reference positivity is a separate qualification and no clipping conceals its
failure. The API returns that unqualified reference as documented. General
sampled-K sinks can also destroy local positivity in an evolving depleted field;
WJ-1 selects no positivity amendment or outer regularization.

Kernel verification is not native coupled equilibrium preservation, a global
realizability/stability proof, a dynamic positivity guarantee, or terrain-flow
validation. No long-time, terrain, optimization or stabilization campaign is run.

## Native integration audit and bounded WJ-2

The following source paths were inspected and deliberately remain unavailable
for WJ PDE dispatch:

- `RANSEquilibriumWallEvolution`: model selection, variables, initialization,
  source/boundary function names and output assume eps_tilde; calibrated profile
  and PIMPLE guards are CLS-specific. Epsilon must become model-aware.
- `RANSWallEvolutionPatch`: OuterClosure and wall/profile creation must dispatch
  one model configuration; preserve the physical 2D-to-3D embedding and complete
  full traction/relative pressure modes before enabling WJ.
- `RANSWallEvolutionNonpressure`: duplicated bool-cls/strain-branch derivatives
  must use the shared WJ action; scalar diffusion must use dnu_d rather than its
  current dnu_t/sigma. The two generally differ off reference.
- `RANSWallEvolutionResolution`: reconstructs CLS stress from implicit viscosity
  minus quadratic/cubic groups; use actual full WJ stress and generic remainder.
- `RANSWallEvolutionPicard`, `Robust`, `Feasible`: residual, frozen solve,
  retraction/admissibility and numerical derivative paths must dispatch WJ
  consistently. CubicScaled, beta=0.4, CLS regularization qualification and
  SaturatedR2 are not WJ defaults.
- `RANSWallEvolutionDiagnostics`: stress decomposition and scalar signatures
  contain CLS identities; report WJ invariants, direct epsilon, scalar viscosity,
  actual raw stress work and full covariance instead.
- `RANSWallEvolutionContact`: preserve sampled-K/epsilon trace sensitivities and
  scalar diffusivity in contact/constraint rows; do not use a CLS-only tangent.
- `RANSWallPhysicsSignature` and `RANSWallEvolutionPatch::Signature`: replace CLS
  identity strings with explicit model, tensor/root/completion, parameters,
  wall source policy, calibration, pressure convention and epsilon ownership.

WJ-2 is bounded to this sequence:

1. Complete native residual/action/diagnostic dispatch.
2. Full traction and pressure jump, and one outer epsilon trace.
3. Model-aware signatures and restart compatibility.
4. Zero-step **uncancelled** reference residual, including all constraints and
   interface rows, without a stored defect correction hiding a baseline defect.
5. Two-step versus 1+1 restart and rejection rollback.
6. Then spatial/reference refinement and short temporal tests at fixed physical
   H, zm and sampling. Long-time/terrain qualification remains later work.
