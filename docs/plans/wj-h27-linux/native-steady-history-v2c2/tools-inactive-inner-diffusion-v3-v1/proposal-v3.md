# WJ-H27 inner k diffusion amendment — inactive v3 proposal

Status: **INACTIVE; NOT APPROVED; NOT IMPLEMENTED; NOT QUALIFIED.** This
document proposes one scientific change for review. It does not amend the
active [v2 contract](wj-h27-amended-candidate-v2.md), authorize a production
run, or activate a different equation through an option or reference table.
The existing v2 sources, executable identities, tables and results remain
unchanged. Implementation requires a separately reviewed, explicitly approved
version. No claim is made that the proposal will pass H27.

## Evidence and its limits

The source-v11 physical v2 steady calculations on 700 and 2800 cells have
converged and independently met their current discrete steady equations. The
reported retained-inlet comparisons are:

| Quantity, scaled by Kin | 700 cells | 2800 cells | Frozen limit |
| --- | ---: | ---: | ---: |
| Inner k departure RMS | 0.132331 | 0.132795 | 0.05 |
| Inner k departure MAX | 0.459519 | 0.440730 | 0.10 |
| Coarse/fine inner k difference RMS | 0.013824 | Same comparison | 0.001 |
| Coarse/fine inner k difference MAX | 0.129202 | Same comparison | 0.001 |

These are the parent audit's reported values; the immutable coarse/fine
assessment paths and hashes must be attached before this proposal is submitted
for scientific approval. The existing
[prospective refinement admission](wj-h27-physical-reference-refinement-v1.md)
records why both physical departure and numerical error must be assessed.
The reported U/W, outer-epsilon and ground-shear physical envelopes pass.

The persistent inner-k adjustment is a material scientific concern. However,
the coarse/fine k differences fail the numerical budget, especially in MAX;
these two solutions do not establish a mesh-resolved limiting profile or a
conclusive continuum-model barrier. Algebraic convergence cannot replace
spatial, quadrature and tolerance convergence. Preserve the v2 failures and
perform the bounded resolution/owner-aware station checks justified by those
data. Do not retune an envelope, discard the inlet region, or declare v3 a
successful remedy from this comparison.

The proposed mechanism is to use the same roughness auxiliary viscosity for
inner k diffusion that enters the retained historical wall construction.
This is a hypothesis about the transport law, not an attribution of the
observed departure to diffusion alone. The current full-WJ collar, finite-nu
reference and developed two-dimensional flow differ from the historical
parallel construction; this change cannot make that inlet an exact steady
solution by assertion.

## Exact proposed equation change

Let r=z-d0, K(x,t)=k(x,H,t), and retain the current roughness function
nuA(K,r). For all 0<=r<H, propose exactly

```
Dkw_v3(K,r) = nu + nuA(K,r).
```

There is no adjustable multiplier and no blending by chi. This coefficient
is distinct from the momentum split coefficient
`mu=nu+(1-chi)*nuA+chi*nuW`. In particular Dkw_v3 remains `nu+nuA` in
the full-WJ part of the inner collar. It is not `mu`, and it is not
`nu+Cd*k^2/(sigma_k*Ew)`. With the retained sigma_k=1, the formula above
is the exact proposal; varying sigma_k or inserting a new turbulent
Prandtl/Schmidt factor is outside its scope.

On the existing bounded negative-intercept roughness branch,

```
qe = Cd^(-1/2),  a(K) = kappa*sqrt(K/qe),
X = h*sqrt(K)/(70*nu),
m = max(0.5-0.4*X^0.7, 1-0.79*X^(-0.28)),
b = -10.7*expm1(m*log(X)),  zv = nu*b/sqrt(K) < 0,
nuA(K,r) = a(K)*(r-zv(K)).
```

This is the existing roughness correlation, with its existing branch and
derivative checks; it is not a new fitted law. For positive nu and the
admitted branch, Dkw_v3 is positive without a floor. At fixed geometry its
directional change is `delta Dkw_v3=nuA_K*delta K`; the local transported
k and dependent Ew do not add separate direct derivative terms. This changes
the scalar nonlinear coupling and must be represented consistently in MMS,
residuals, any Jacobian/preconditioner and independent verification.

Retain conservative scalar transport exactly:

```
inner:  dt(k) + div(U*k) = div(Dkw_v3*grad(k)) + Pactual - Ew,
outer:  dt(k) + div(U*k) = div(Dko*grad(k)) + Pactual - epsilon,
        dt(epsilon) + div(U*epsilon)
          = div(De*grad(epsilon))
            + Ce1*Pactual*epsilon/k - Ce2*epsilon^2/k,
Dko = nu + Cd*k^2/(sigma_k*epsilon),
De  = nu + Cd*k^2/(sigma_e*epsilon).
```

Keep independent Cd=572/6561, sigma_k=1, sigma_e=1.3, Ce1=1.44,
Ce2=1.82 and physical nu=6.423173803526447e-5. Keep the approved v2
isotropic lower-wall/full-WJ-collar tensor, geometry-only quintic chi,
DevStrainRawWorkV1, selected physical WJ root, raw
`Pactual=-R:grad(U)`, all outer equations and all momentum/AWF terms.
No EP, source adjustment, actual P-E balance constraint, clipping, floor or
extra forcing is proposed. Manufactured loads remain verification-only and
must be rederived from the amended continuous equations.

Keep A(K)=K^(3/2)/cl, `cl=kappa*qe^(3/2)`, and
`Ew=A(K)/max(r,c(K))`. Keep the entire historical application inlet, its
old inner k/U construction and its full-height outer profiles. Keep Kin,
roughness h=.10854503542535461, H=.15, d0=.0875, domain, physical units,
top flux including molecular viscosity, all other boundary laws, and the
frozen physical/numerical limits. The new reference column is never an inlet
replacement.

## H interface and epsilon inventory

Continue k and velocity across H, with common conservative mass, traction
and k flux. A diffusion jump is allowed:

```
Dkw_v3(H-) = nu + nuA(K,H),
Dko(H+)   = nu + Cd*K^2/(sigma_k*(A(K)/H)).
```

These expressions need not coincide. For a smooth solution on either side,
using the same upward normal, require
`Dkw_v3*partial_n k|H- = Dko*partial_n k|H+`. Equivalently the two
outward diffusive contributions cancel. Common U and k also give the common
advective flux. Do not impose C1 k or average the coefficients to remove the
jump. Actual interface flux may be nonzero; the reference endpoint condition
below does not prescribe it. Native weak cancellation and independently
sampled one-sided physical flux convergence are separate checks.

Transport epsilon only on the outer domain. Its sole H condition remains
the dependent Dirichlet trace `epsilon(H)=A(K)/H`, with its exact trace
variation `A'(K)*delta K/H`. There are no inner epsilon degrees of freedom,
no inner epsilon transport equation and no additional epsilon flux-matching
condition. Keep the existing common velocity/traction formulation and
pressure-level determination; add no pressure jump or artificial C1 velocity
condition.

## Rederive the counterfactual reference and outlet data

Retain the explicit local reference choice `t=K/qe`, `k_ref(H)=K`,
`Q_ref(H)=Q_ref(0)=0`, with `Q_ref=Dkw_v3*k_ref'` and
`Q_ref'=Ew-P_ref`. Solve the same finite-nu v2 shear equation and the
amended scalar column on the declared bounded shear/root branch: retain the
increasing branch `0<=q<=72/13` where chi>0, with q=k_ref*G/Ew. The
general endpoint-flux budget is

```
A*(1+log(H/c)) - integral_0^H P_ref dr - (Q_ref(H)-Q_ref(0)) = 0.
```

This equation calibrates a counterfactual reference only. The actual flow
retains its full storage, transport, production, destruction and interface
flux imbalance. Do not constrain its instantaneous integrated P-E to zero.

The conditional fast derivation can be checked anew. If the newly obtained
branch satisfies `0<c<.4H`, then Ew=A/r on [.4H,H]. Integrating the new
column backward from H through that interval is independent of c. Below
.4H, chi=0, `G=t/(nu+nuA)` and `P_ref=nuA*G^2` are independent of
the local k and Ew. Put `d=nu-a*zv`; the lower production integral is

```
(1/A)*integral_0^(.4H) P_ref dr
  = log1p(a*.4H/d) + nu/(d+a*.4H) - nu/d.
```

With `B_P=integral_0^H P_ref dr`, the selected zero-endpoint reference
then gives `log(c/H)=1-B_P/A`; for general Q_ref(0)=0 it is
`1-B_P/A-Q_ref(H)/A`. Under these branch hypotheses, the normalized
ground-flux residual Q_ref(0)/A has derivative +1 with respect to log(c/H);
the displayed integral-budget residual divided by A has derivative -1.
Recompute B_P using the proposed
diffusivity; v2 coefficients cannot be reused. Reject a branch exit or a
nonpositive k/shear/root margin, rather than extrapolating the formula.
Integrate the complete lower column as well: a formal cutoff alone does not
establish positive k, a ground zero flux, or an admissible full column.

Before an executable v3 table exists, independently compare full trial-IVP
roots, split quadrature budgets and velocity anchors against this conditional
construction over K/Kin in [.9,1.1]. Retain all failed trials. Qualify the
cutoff and its K derivative, endpoint traces, roughness branch and state
margins with refinement and independent held-out points. Preserve the
existing reference budget/flux tolerance 1e-10, derivative error check 1e-8 in
`Kin*dlog(c/H)/dK`, and traction tolerance 1e-10 in tau0 units. Reject
out-of-range runtime K and invalid cutoffs; do not clamp or extrapolate.

Regenerate and independently certify the fixed outlet reference traction
using the new Kin column:
`(Kin/Ne-2*chi*k_ref*d_WJ(q_ref),0,tau0)`, Ne=81/20, where
`d_WJ(q)=3*q^2/[5*(N(q)^2+q^2)]` uses the unchanged selected WJ root.
Below .4H its x component
remains Kin/Ne. Keep the outer counterfactual `k=Kin`, `epsilon=Ain/r`
and finite-nu constant-total-shear branch. Check all segment endpoints and
H continuity. Preserve the physical top epsilon flux and the inlet/H corner
binding. Issue a new whole reference identity even where the old and new
outer profiles agree; never mix a v2 cutoff/traction segment with a v3
calibration silently.

## Versioning and qualification required after approval

Use a distinct full-model identifier for the new scalar law while retaining
the explicit v2 tensor identifier, plus new source, reference-table, manifest,
MMS and checkpoint identities. A cold restart must reject a v2 identity;
no v2 time levels or explicit histories are valid v3 histories. An expressly
labelled v2 field may be considered only as a fresh steady initial guess,
with recomputed dependent trace/reference/advection data, step/time zero and
no inherited convergence or qualification. Preserve all v2 evidence.

The complete [H27 gates](wj-h27-contract.md) remain binding:

1. Reopen scientific feasibility for the amended column, table and diffusion
   jump before application use. Check positive scalar diffusion, current
   covariance, root and applicable full planar principal margins, including
   the collar; unchanged tensor code does not certify changed solution states.
2. H27-B requires the registered coupled solver with the exact new equation,
   independent continuous evolving MMS and directional/current residual
   checks, dependent Gamma elimination, signed ownership, conservative scalar
   and mechanical budgets, direct-reference checks, invalid/failure/history
   controls, temporal order and spatial/quadrature/tolerance refinement. Verify
   both sides of H with a flux-compatible manufactured case; an imposed flux
   jump belongs only to an explicitly derived verification load.
3. H27-C requires a newly resolved same-inlet steady 2D reference and the
   unforced developed, retained-inlet startup and finite-perturbation pilots.
   Keep the station union, whole/inner/outer envelopes, owner-aware profiles,
   perturbation support/amplitudes, common endpoint, three dt levels and all
   numerical budgets. Keep inner k RMS/MAX limits .05/.10 and numerical
   RMS/MAX budget .001 in Kin units. A less-adjusted but unresolved profile
   still fails qualification.
4. H27-D requires new-model cold restart, recovery/rotation/retention and
   complete measured coupled cost on the same accuracy-qualified recipe,
   followed by the executable manifest. Preserve the original horizon,
   terminal window and per-trajectory `1.25*Ctotal<=86400s` gate. Do not
   run the full Task3 horizon during this qualification.

No resource allowance is reset by this proposal. All prior H27 charges,
failures and the protected evidence reserve carry forward. Forecast and admit
the full remaining matrix from complete coupled measurements before expensive
work, within the existing Linux ledger and worker/storage/RSS limits.
Readiness remains **TASK3_READY=NO**. A future failure must be reported with
its witness and uncertainty; neither this proposal nor positive reference
samples close H27-B/C/D.
