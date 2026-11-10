# WJ-2 native serial flat 2D contract

This package continues reviewed WJ-1 HEAD
`7a06e1c1f7e4e81d42f53f5c5d94bea33668bb85`. It uses the existing full 3D
WJ tensor, cubic-N surrogate, implicit directional action and calibrated
finite-viscosity profiles. The [WJ-1 contract](wj1-theory-to-code-contract.md)
continues to specify their equations and independent exact tangent certificate.
This package enables a bounded native IMEX1/PimpleSingleSweep route and tests
its transactions. A frozen matrix is not the full nonlinear R3 Jacobian.

The specification SHA-256 identities remain:

```
1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124  wj_earsm_neutral_abl_ep_awf_note.tex
0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c  wj_earsm_budget_cutoff_awf_appendix.tex
```

The originals and worktree copies are unchanged. The main-note
`eq:integrals-linear` typography has literal `quad` where `\quad` was intended:
`Delta J1=(z2-z1)/a-(d+a h)/a^2 log(y2/y1)` for `z1>=h`.

## Owned configuration and field semantics

`RANSWallPatchParameters::wj` owns `WJParameters::Defaults()`.
`LoadWJNativeConfiguration` validates the reviewed constants and creates the
WJ neutral reference. WJ code never obtains constitutive parameters from the
legacy `RANSNeutralReference::Parameters()` facade. Native calls explicitly
select `DevStrainRawWorkV1`; the strict incompressible kernel remains available
for theory tests.

| Property | WJ-2 choice |
|---|---|
| Model | WJHighRe, full 3D tensor with 2D cubic-N surrogate |
| Constants | c1=1.8, c2=5/9, Ce1=1.44, Ce2=1.82, sigma_k=1, sigma_epsilon=1.3, Cd=572/6561 |
| Reference | qe=81/(2 sqrt(143)), kappa derived from actual Cd (0.381917863546733) |
| Gradient | S=(k/epsilon) dev(sym G_raw), Omega=(k/epsilon) skew G_raw |
| Scalar transport | direct physical epsilon, nu_d=Cd k^2/epsilon; zero outer molecular diffusion |
| Stress | Original full WJ outer covariance; full-strain Linear inner covariance |
| Wall EP | WJCalibratedEPAWF, uniform chi=1, InterfaceTrace K=k(H), zm=H/2 |
| Wall budget | WJCalibratedBudgetAWF, InterfaceTrace K=k(H), zm=zK=H, no physical chi |
| Wall viscosity | physical nu>0, WJ calibration c_mu=C*, c_length=kappa qe^(3/2) |
| Integration | serial flat 2D tensor-cubic basis, IMEX1/PimpleSingleSweep, eta=0 |
| Numerical bound | explicit MassProjection with WJNativeV1 inventory/signature, inactive on tested startup and baseline steps |
| Forcing/reuse | zero prescribed body/scalar forces and gravity; evaluation reuse Off |

In native 2D, velocity components `(u,v)` and mesh coordinates `(x,y)` embed
in physical `(x,z)`. The physical y velocity and its derivatives are zero.
All constitutive tensors, including R_yy, remain 3 by 3. Canonical session,
output and checkpoint scalar identity is `epsilon`; legacy `eps_tilde`
sessions keep their original semantics. Neutral and physical function-reference
initialization use the same typed configuration and direct epsilon name.
A resting function-defined state and its immutable restart test canonical
function/session/checkpoint semantics separately from the neutral reference.
WJ rejects legacy epsilon variable lists and checkpoint identities.

WJ transitions, physical destruction capping, CubicScaled, gamma regularization,
SaturatedR2, nonlinear Picard/R3 searches, condensation, curved maps, 3D and MPI
are rejected. General kernel parameter/action support does not imply native
support for changing Cd or the gradient completion. These native knobs reject.
Legacy capabilities and defaults remain available on their existing routes.

## Full stress, raw work and frozen systems

The shared legacy robust modifier now refreshes every generic returned stress
field after changing R. For legacy raw-strain CLS,
`deviatoricStress=R-2kI/3` may have nonzero trace off constraint; it is not
projected. Always `anisotropy=deviatoricStress/k` and
`productionDev=P+(2/3)k tr(G_raw)=-deviatoricStress:G_raw`.
Original, CubicScaled and Regularized tests preserve the original physical R
and production. They check a beta!=1, gamma>0 state.

WJ's trace-free deviatoric stress instead splits naturally against
`dev(sym G_raw)`. Its returned `stressRemainder` uses that split. An implicit
raw-symmetric-gradient matrix requires the additional isotropic remainder
`(2/3)nu_eff div(U) I`. Native assembly defers the difference between its
implicit operator and the complete returned stress; it never inserts the WJ
remainder into a legacy raw-strain reconstruction. A focused compiled test
checks both reconstructions and simultaneous k/epsilon/G actions at divergent
states with nondefault constant Cd and sigma values.

The physical residual and resolution/energy loads use full R,
`Praw=-R:G_raw=Pdev-(2/3)k div(U)`, independent `nu_d` diffusion and direct
epsilon sources. No general production is inferred from N. Frozen flow, k and
epsilon blocks use the same values, state, history, time, constraints and
coefficients. Pressure columns use `AffinePressureAction`; the WJ wrapper
does not call the unsupported complete R3 action. `AssembleNonpressureR3`
explicitly rejects WJ.

After vector velocity and pressure/gauge, WJ solves **k then epsilon**. This
first-order single-sweep ordering gives the sole epsilon solve the final
endpoint K and retains its trace mass columns. The physical epsilon trace
therefore remains consistent after the k update. There is no post-solve Lift,
trace reset, nonlinear restoration or additional equation solve. The scalar
order is bound in the WJ signature. Legacy epsilon-then-k ordering is preserved.
The uncapped k block freezes the incoming outer epsilon, while the epsilon
block freezes the newly transported k; each matrix is checked against the
physical residual at exactly its own assembly state.

Uncapped EP and budget carry the **complete** `source.kSource` in the frozen
k load with zero local destruction reaction. Only a capped legacy selection
uses the declared production/implicit-rate split. At resting mean velocity,
the EP sink is `-p_A,0` and the budget sink is `-E_c`. Both physical and frozen
loads retain these sinks in nonconstant modes. Normalization B/m/lambda is
inapplicable and cannot affect WJ sources or epsilon_H. Actual
`CappedDestruction()` controls normalization in PIMPLE, Contact, Source,
Feasible and physical sample reporting. Legacy transition/forensic diagnostics
reject WJ rather than applying CLS qualification identities.

## Interface, pressure and reference qualification

The fitted interface shares velocity and k values. Each region retains its
one-sided gradients and pressure trace. Weak assembly uses full mean traction
with opposing normals. Eddy viscosity and mechanical pressure need not be
continuous. There is one dependent outer epsilon value trace and no inner
epsilon PDE or extra prescribed epsilon flux. The epsilon transport matrix
retains the dependent trace columns and their target.

Native stored pressure is `pi=p/rho+2k/3`, without gravitational potential.
Sampled forcing recovers the mechanical pressure gradient as
`grad(pi)-2 grad(k)/3` before evaluating wall forcing. With zero gravity,
budget inner pi is constant while mechanical pressure varies with the
transported k. The continuum interface pressure offset is `-K/Ne`, Ne=81/20;
the pressure representation retains two relative traces and one gauge.

The pilot uses nu=1.5e-5, u*=0.5, H=1, h=.1, K=0.8466950351321795,
top=2 and one horizontal element. The graded normal knots are recorded in
`wj2_cases.py`, including .0025,.005,.01,.02,.04,.08,.1 and the fitted H.
The counterfactual cutoff uses endpoint K, independent of actual forcing.
Evolved k is transported. EP's reference k is constant; budget's reference is
`K-integral_z^H Q0/D_k,w`, with neither ground k=0 nor a local-k coefficient
substitution.

The appendix's sufficient rough-reference conditions are **zv<=0, d>=2nu,
Lambda>=1**, where d=nu-a zv. They imply Q0<0 in the interior, k_ref>K and
positive reference covariance. Independent quadrature records these actual
hypotheses and both Q0 endpoints. This statement is about the continuum
reference. The represented native polynomial is audited separately over the
complete native consumer/interface/reconstruction/ground inventory.

At the same K,nu,H but h=0, the valid cutoff is about
.000580427216979312 and k_ref(0)=-.3301275029511697. Neutral initialization
rejects this with `reference_nonpositive` and an exact location/value witness
before projection or flooring. The low-level kernel keeps the true negative
reference. Cutoff failure and invalid scalars have separate tests/reasons.

## Audit and independent verification

`WJCompleteAudit` uses NativeElementSideV2 geometry, element, region, side,
gradient and scalar provenance, including both interface sides. It reports
scalar/trace positivity and profile gates separately from covariance spectra,
root branch/residual/F_N, WJ s/w/I4 ranges, nu_eff/nu_d, divergence and
Praw/Pdev/work. CLS max-invariant branch, beta=.4 counterfactual, R1/R2 box
and cap-lambda qualification are marked inapplicable.

The momentum tangent diagnostic uses the shared WJ action at frozen positive
k and epsilon on incompressible rank-one increments. In-plane directions
receive a 32-angle bracket plus bounded minimization; transverse polarization
uses its symmetric 2 by 2 quadratic form. This is finite spatial sampling and
a numerical x-z propagation diagnostic, not a global 3D stability certificate.
WJ-1 separately retains the exact all-direction **reference** bound nu_eff/250.
Both adverse WJ-1 witnesses remain explicit: positive covariance with negative
principal form at q=3.5, and the N=3 state with two negative variances. They
pass by reproducing the adverse values, with no clipping.

`RANSReturnedFieldsChecks` verifies the common returned-field repair.
`WJNativeEP` and `WJNativeBudget` run the real solver and strict immutable
checkpoint probe. Independent NumPy/complex-step WJ algebra and separate
split SciPy quadrature check extracted divergent states/actions and actual
represented source/diffusion loads. New higher-order native field evaluations
export basis values, coefficient IDs and gradients; expected stresses and
loads come from independent code. Endpoint-K perturbations check cutoff and
epsilon-trace responses. Pressure-only differences include inner, outer,
relative, gauge and a nonconstant outer mode producing roughness forcing.

Root/tensor identities retain about 1e-12 scaled tolerance. Actions use
complex-step comparisons at roughly 3e-10 relative, 5e-11 absolute; integrated
loads use 1e-10 relative and 2e-10 absolute. Frozen residual comparison uses
1e-8 absolute, pressure columns 1e-10, centered pressure differences 1e-7.
These checks record actual defects; tolerances are not convergence claims.

## Zero step and transactions

Startup reports the uncancelled spatial residual with alpha_mass=0, eta=0.
`raw_reference_norms` is ordered `(u,v,k,outer epsilon,inner continuity,
outer continuity,gauge,boundary/trace)`. Each unconstrained regional load uses
the native consistent-mass dual norm, divided by its physical rate scale and
sqrt(region volume). Rate scales are u*^2/H, K0 u*/H, epsilon0 u*/H, and
u*/H for continuity. Gauge scales by u*^2; constrained rows use their field
amplitudes. The complete outer epsilon region is included. Interface traction,
k flux, pressure jump and strong divergence are separately reported. No eta
anchor or equilibrium source cancels representation error. One normal-mesh
refinement and a source quadrature-order check keep all physical inputs fixed.

The declared dt is 1e-6 s, dt*u*/H=5e-7. Fresh immutable WJ checkpoints supply
two uninterrupted steps and one step/restart/one step with the same runtime.
Primary coefficients, gauge, time, step count, every native GLM/history vector
and time label, last implicit history, accepted audit packet and physical
samples must match exactly. An accepted checkpoint has a zero-step native
readback and a strict probe readback. WJ restart restores the saved reference
directly without constructing/projecting/lifting another reference; strict
diagnostics perform no Lift, governing solve, state change or history advance.

End-of-step deliberate hard rejection occurs after the scheduled equations
and audit. It must restore the last accepted state/history/time and packets,
retain the failed candidate separately, and retain one vector-U, one pressure/
gauge, one epsilon, one k solve and one logical bound operation. Genuine early
failures retain partial counts. Raw/corrected k candidates, transport matrices,
correction budgets and stage hashes remain available; baseline corrections
must be inactive. MassProjection is a numerical operation, not physical
destruction capping.

Positive changed wall viscosity and valid changed EP/budget choices must
reach and fail checkpoint **signature mismatch** validation. Unsupported Cd,
gradient completion, cap, sampler, transition, epsilon semantics, map and
reuse tests are configuration rejections, separately identified. No saved
signature is edited. The existing short Regularized CLS budget transaction
is a distinct registered shared-path control.

The [Linux evidence report](wj2-linux-evidence.md) binds the tested source,
compiler policy, commands, inventory, native/oracle JSON and checkpoints.
Successful short native steps and transactions do not establish long-time
preservation, full spatial/temporal convergence, global positivity,
realizability/stability, native 3D/MPI or terrain validity. The bounded WJ-3
recommendation in that report follows the actual residual and profile results.
