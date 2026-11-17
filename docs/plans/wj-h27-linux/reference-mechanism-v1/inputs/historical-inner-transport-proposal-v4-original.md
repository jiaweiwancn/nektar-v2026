# H27 historical inner transport — inactive v4 candidate

Status: **NOT APPROVED; NOT IMPLEMENTED; NOT QUALIFIED.** Approval would
authorize a separately versioned scientific experiment and its full H27
qualification. It would not accept an application result, loosen a gate or
authorize full Task3. The active solver remains the user-approved v2 tensor
with its original scalar diffusivity and counterfactual cutoff calibration.
The task package section 8 requires separate review before an additional
physical-model change. The prior tensor approval does not cover these two
new scalar/reference changes.

## Evidence and decision

The [source-v11 checkpoint](wj-h27-steady-history-checkpoint-v2.md), commit
`1d2df25fb15381eaff1a3dbbfc17fdf7d024b5b1`, preserves independently converged
700/2800-cell v2 equations and failed physical and numerical reference gates.
Inner-k physical RMS departures are .132331/.132795 in Kin units, against
.05; MAX departures are .459519/.440730 against .10. The inter-mesh inner-k
RMS/MAX differences .013824/.129202 exceed .001, so this is not a resolved
continuum barrier proof. An unchanged-physics x refinement is in progress.

The [v3 diffusion-only proposal](wj-h27-inner-diffusion-amendment-proposal-v3.md)
does not have affirmative mechanism evidence. Its inactive zero-endpoint-flux
column differs from the historical inlet more than v2. A separate independent
factorial column control now distinguishes the two proposed changes:

| Inner diffusion | Cutoff rule | Inner k RMS / MAX departure in Kin | Counterfactual QH/A |
| --- | --- | --- | ---: |
| v2 scalar law | full-WJ zero-endpoint-flux calibration | .099384 / .326425 | 0 |
| v2 scalar law | historical cutoff | .116530 / .486285 | .04451675 |
| nu+nuA | full-WJ zero-endpoint-flux calibration | .188467 / .424999 | 0 |
| nu+nuA | historical cutoff, proposed v4 | .018939 / .021467 | .04487878 |

All four full Kin-column controls pass ground-flux, independent integrated
budget, positive-k, declared shear-root and traction checks. The v4 control's
inner velocity MAX departure is .006323 Uref. The last row is within the
existing inner k and velocity envelopes **for that counterfactual column
only**. It is not a same-inlet 2D application reference, a full principal
certificate, an all-K table, or a physical pilot. The result supports testing
a concrete hypothesis; it does not predict success of the 2D solver.

This proposal selects two already defined historical transport ingredients,
not coefficients fitted to the observed 2D errors. Their joint effect is
tested explicitly above. Full-WJ collar stress, physical viscosity and the
application inlet remain binding. If v4 fails, preserve and report the failure.

## Exact proposed equations

Let r=z-d0, K=k(x,H,t), qe=Cd^(-1/2), cl=kappa*qe^(3/2),
A(K)=K^(3/2)/cl, and retain the existing negative-intercept roughness branch:

```
a(K) = kappa*sqrt(K/qe),
X(K) = h*sqrt(K)/(70*nu),
m(K) = max(.5-.4*X^.7, 1-.79*X^(-.28)),
zv(K) = -10.7*nu*expm1(m*log(X))/sqrt(K),
nuA(K,r) = a(K)*(r-zv(K)),
d(K) = nu-a(K)*zv(K),
J(K,r) = log1p(a(K)*r/d(K))
           + nu/(d(K)+a(K)*r) - nu/d(K).
```

Change exactly the inner diffusivity and dependent cutoff to

```
for 0 <= r < H:
    Dkw_v4(K,r) = nu+nuA(K,r),
    c_v4(K) = H*exp(1-J(K,H)),
    Ew_v4(K,r) = A(K)/max(r,c_v4(K)).
```

Retain sigma_k=1; no new multiplier, adjustable sigma, momentum-mu substitution
or blend by chi is allowed. The proposed Dkw remains nu+nuA throughout the
inner collar, including its full-WJ portion. The cutoff is the historical
auxiliary parallel-column calibration at the current K. It is no longer
selected by imposing two zero fluxes on a full-WJ reference column. Reject an
invalid branch, nonpositive coefficient or cutoff outside its qualified
0<c<.4H interval; never clamp or extrapolate. Retain K/Kin in [.9,1.1].

Actual conservative equations remain

```
dt(k) + div(U*k) = div(Dk*grad(k)) + Pactual - Eactual,
Pactual = -R:grad(U),
Eactual = Ew_v4 (inner), epsilon (outer),
Dk = Dkw_v4 (inner), nu+Cd*k^2/(sigma_k*epsilon) (outer),
dt(epsilon) + div(U*epsilon)
    = div((nu+Cd*k^2/(sigma_e*epsilon))*grad(epsilon))
      + Ce1*Pactual*epsilon/k - Ce2*epsilon^2/k  (outer only).
```

Keep Cd=572/6561, sigma_e=1.3, Ce1=1.44, Ce2=1.82,
nu=6.423173803526447e-5, Kin=.008809015145515192, H=.15,
h=.10854503542535461, d0=.0875 and all physical units/domain/inlet parameters.
Keep A(K), the entire historical inlet, all outer equations and coefficients,
all velocity/pressure/AWF terms, conservative EX2 advection and BDF2 history.
No P-E balance is imposed on the actual flow. No EP, floor, clipping,
equilibrium source or new application forcing is introduced.

Keep the approved v2 tensor exactly:

```
P = I-n*n^T, M = G+G^T, bA = -nuA*P*M*n,
Riso = 2*k*I/3 + bA*n^T + n*bA^T,
Rwall = (1-chi)*Riso + chi*RW(k,Ew_v4,G),
V = nu*M-Rwall+2*k*I/3.
```

The existing geometry-only quintic chi is zero below .4H and one above .75H.
RW uses the complete full-3D DevStrainRawWorkV1 WJ tensor and unchanged selected
physical root. New Ew changes its inputs, so unchanged tensor code does not
certify changed states. Preserve current covariance and full supported planar
principal checks with their existing reserves. Keep pure-WJ outer nu=0 API
convention and explicit physical molecular viscosity in both PDE regions.

## Interface, derivatives and counterfactual reference

Keep common k, velocity and conservative interface mass/traction/k flux.
Dkw_v4(H-) generally differs from Dko(H+); allow the physical gradient jump
while checking common diffusive flux. Do not require C1 k or average its two
owners. Epsilon remains outer-only, with the sole dependent condition
epsilon(H)=A(K)/H and derivative A'(K)*deltaK/H. No inner epsilon degrees of
freedom or extra epsilon-flux matching condition is introduced.

At fixed geometry the exact variations are

```
delta Dkw_v4 = nuA_K*deltaK,
delta log(c_v4/H) = -J_K(K,H)*deltaK,
delta Ew_v4 = (A'/r)*deltaK                 when r>c_v4,
delta Ew_v4 = Ew_v4*(A'/A-c_v4'/c_v4)*deltaK when r<c_v4.
```

Retain the moving-cutoff split treatment and its one-sided derivative handling.
Derive and independently check J_K across every admitted roughness branch;
reject a non-differentiable branch crossing instead of inventing a derivative.
No finite-difference/model Jacobian may silently keep the v2 scalar law.

For outlet calibration/verification only, keep t=K/qe, k_ref(H)=K,
Q_ref(0)=0 and the exact v4 Dkw/cutoff. Solve for **Q_ref(H)** on the bounded
declared reference branch; do not also force Q_ref(H)=0. Require

```
Q_ref = Dkw_v4*k_ref',
Q_ref' = Ew_v4-P_ref,
A*(1+log(H/c_v4)) - integral(P_ref,0,H)
    - (Q_ref(H)-Q_ref(0)) = 0.
```

The Kin diagnostic selects Q_ref(H)/A=.04487877859763451 in [0,.4]. Certify
existence, uniqueness/conditioning and derivatives throughout [.9,1.1] Kin,
including held-out roots and quadrature refinement. This is a counterfactual
inner column, not a whole-domain stationary solution: its nonzero H flux is
not compatible with an outer constant-k parallel column. It must never be
substituted for the independent coupled 2D physical reference or applied as
an actual H boundary flux. That reference must solve common H flux itself.

Regenerate the complete fixed outlet traction from the new Kin reference
column, including its x-normal stress and continuous H endpoints. Keep the
outer finite-nu reference branch, physical top shear/epsilon flux, no-slip
ground, zero ground k flux, scalar outlet conditions, no backflow and actual
pressure-level determination. New reference coefficients require a new whole
identity even where an outer segment is unchanged. Certify endpoint/budget
errors <=1e-10, derivative errors <=1e-8 in the existing normalized units and
traction errors <=1e-10 tau0. No executable v4 table exists at proposal time.

## Implementation and full qualification after approval

Give v4 its own complete model, table, source, MMS, checkpoint and manifest
identities while retaining v2 as the tensor identifier. Reject v2 histories
and restarts. Start fresh from the retained inlet-derived fields and rebuild
all dependent traces and explicit histories; no v2 acceptance carries forward.

Reopen H27-A for the new reference/root/cutoff and interface feasibility.
Re-derive independent continuous coupled MMS and all changed scalar, trace and
boundary evaluations. Verify directional and signed current actions, root/trace
failure controls, temporal order, spatial/quadrature/tolerance separation,
conservative/mechanical and two-owner interface accounts. Existing v2 evidence
is a foundation, not v4 qualification.

H27-C requires a newly resolved, same-inlet steady 2D reference passing every
frozen physical envelope and numerical budget, then the actual developed,
startup and both finite-perturbation pilots. Retain amplitudes/support,
26.13350407707765 minimum principal evolving duration, common-time three-dt
checks and refined h/p/q owner-aware comparisons. Do not change the inner-k
.05/.10 physical or .001/.001 numerical RMS/MAX limits.

H27-D requires physical cold restart, actual recovery/retention/terminal and
monitoring/storage qualification, complete warm blocks and rebuild/retry cost,
>=95% disjoint attribution, and a pinned executable launch manifest. Keep
T=505.02466521278404, terminal window [447.8483832973273,505.02466521278404]
and per-trajectory 1.25*Ctotal<=86400 s, plus total campaign disclosure. Do not
execute the full Task3 horizon during H27.

No allowance resets. The same 43200 aggregate Linux-operation seconds,
protected 3600 evidence seconds, <=72 actual shared-host workers, one science
rank/thread policy and existing RSS/storage/free-space caps remain binding.
Forecast/admit each complete remaining matrix from measured coupled costs.
`TASK3_READY=NO` until every gate passes on one executable recipe.
