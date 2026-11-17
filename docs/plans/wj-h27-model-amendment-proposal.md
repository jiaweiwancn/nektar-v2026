# H27 model amendment proposal: isotropic wall normal stress and full-WJ collar

## Status and decision boundary

**PROPOSAL ONLY. No active PDE/model change, no production qualification, and no Task3 launch authorization.** This is a concrete candidate for separate user review after publication of the specified-model obstruction. It does not convert H27's failed model gate into PASS. The candidate tensor must not enter the active solver until the user authorizes this scientific change.

This proposal follows the independently checked negative principal action of the specified H27 wall tensor. The coordinator reports an authoritative Linux native/oracle comparison over 23 states with maximum error 1.68e-14, while the model fails its principal gate. On the retained inlet the normalized ground minimum is about -0.9350801, and positive-width witnesses at r=.0015,.003,.006 give approximately -.8220,-.6901,-.4240. The separately derived zero-endpoint reference family bound c/H>=.1619259 and inlet-direction upper bound -.518046 rule out curing that family by arbitrary positive inner diffusivity. The finite-nu outer reference gives Q_H/A=-.03958397595 and c/H=.21367837, which worsens this mechanism. These numbers are coordinator-supplied summaries of Linux evidence; the closure report must bind their exact receipts. No numerical run was performed on the Mac for this proposal.

The source coverage limitation is documented in `wj-h27-boundary-audit.md`, section 8: historical inner points were checked for covariance but skipped the outer full-WJ principal search. H26's tensor/work/collar tests were not a certification of the revised wall normal-stress derivative.

## Exact candidate

Retain physical nu>0, the old imposed application inlet, H=.15 (6 mm), d0=.0875 (3.5 mm), hs=.10854503542535461, the scalar coefficient convention, and the geometry-only C2 transition:

\[
 \chi=0\quad(r/H\le .4),\qquad
 \chi=1\quad(r/H\ge .75),\qquad
 \chi=s^3(10-15s+6s^2),\quad s=(r/H-.4)/.35
\]

in the intermediate interval. Let n be the physical wall normal, P=I-nn^T, M=grad U+(grad U)^T, and use the existing positive auxiliary viscosity nu_A. Define

\[
 b_A=-\nu_A PMn,\qquad
 R_{iso}=\frac23kI+b_A\otimes n+n\otimes b_A,
\]
\[
 \boxed{R_{proposal}=(1-\chi)R_{iso}+\chi R_W(k,E_w,\nabla U).}
\]

Above the wall layer, use the same full WJ tensor with transported epsilon. Local k remains the actual transported scalar, Ew the actual dependent destruction. This replaces the full WJ tangential block and normal-normal component below the collar with an isotropic normal-stress block; it preserves the signed auxiliary tangent-normal shear. It is a model amendment, not numerical stress clipping or an added artificial viscosity. It is also distinct from the historical full-strain linear inner model, since R_iso supplies only tangent-normal turbulent shear.

Trace R_iso=2k because b_A.n=0, and the construction is covariant under simultaneous rotation of the physical tensor and normal. Full WJ is recovered identically in the open collar (.75H,H), with the same dependent epsilon trace and common limiting interface state. The candidate does not require C1 finite-element velocity continuity.

## What the candidate proves at chi=0

Freeze scalars, geometry and nu_A in the local momentum principal block. For a unit incompressible rank-one velocity gradient B=a tensor xi with |a|=|xi|=1 and a.xi=0,

\[
 \delta R_{iso}[B]=-\nu_A\{P(B+B^T)n\otimes n
             +n\otimes P(B+B^T)n\},
\]
\[
 \boxed{B:\{\nu(B+B^T)-\delta R_{iso}[B]\}
       =\nu+\nu_A|P(B+B^T)n|^2\ge\nu>0.}
\]

Thus the candidate removes the identified chi=0 normal-stress-derivative mechanism without a fitted floor. The result covers all propagation/polarization directions, not only the planar shear derivative. It is a local frozen-scalar principal statement, not a global nonlinear well-posedness or positivity theorem. Below the collar, modes without tangent-normal shear have only molecular coercivity; conditioning and attainable numerical resolution require measurement.

The three covariance eigenvalues of R_iso are `2k/3-|b_A|`, `2k/3`, `2k/3+|b_A|`. Consequently,

\[
 \boxed{R_{iso}\succeq0\quad\Longleftrightarrow\quad 2k/3\ge|b_A|.}
\]

This condition must be tested on actual fields; positivity of k alone does not ensure it. In the blend, covariance follows from convexity if both constituent tensors are realizable, but the full actual tensor and margins must still be reported. For fixed geometry the principal quadratic is the convex combination of the candidate chi=0 quadratic and the full-WJ quadratic including molecular viscosity. The latter need not be positive over all states: **collar admissibility remains unproved and must be tested**. Geometry derivatives and derivatives of chi contribute lower-order terms and remain in the complete stress divergence/work accounts.

## Reference identities retained, physical equations changed

For a parallel column U=(U(r),0,0), the original specified H27 tensor and this candidate have identical tangent-normal shear at every chi. Their contracted production is therefore identical:

\[
 P_{parallel}=\{(1-\chi)\nu_A+\chi C(q)k^2/E_w\}G^2,
 \qquad G=U',\quad q=kG/E_w.
\]

If scalar diffusivity, tail amplitude, endpoint flux, shear branch and all other reference inputs remain fixed, the parallel k/reference-budget/cutoff equations and finite-nu velocity anchor are unchanged. Existing scalar reference-root results can be used as controls, then must be independently rechecked with the candidate tensor. This identity does not permit reuse of old evolving states or a change to the imposed historical inlet.

For general flow, raw production changes to

\[
 P_{proposal}=(1-\chi)\left[\nu_A|PMn|^2-\frac23k\nabla\cdot U\right]
                 +\chi[-R_W:\nabla U].
\]

The `-2k div(U)/3` term remains for weakly incompressible fields. Do not preserve the old actual production while changing momentum stress. The k and epsilon equations, mechanical-energy budget and dependent-destruction policy must use the same actual stress consistently. General tangential/tangential turbulent stresses, normal traction, modified pressure and their feedback into sampled mechanical-pressure forcing also change.

With `V_ref=nu M_ref-R_ref+2k_ref I/3`, keep the prescribed-reference outlet principle, but rederive the candidate's fixed loads:

\[
 \Pi_{ref}=V_{zz,ref}-K_{in}/N_e,\qquad
 t_{x,ref}=V_{xx,ref}-V_{zz,ref}+K_{in}/N_e,\quad t_{z,ref}=u_\tau^2.
\]

For planar pure shear, writing the WJ anisotropy as `R_W,xx/k=2/3+d`, `R_W,zz/k=2/3-d`, the candidate gives `Vxx=-chi*k*d`, `Vzz=chi*k*d`, hence `Pi_ref=chi*k*d-K_in/Ne` and `t_x,ref=K_in/Ne-2chi*k*d`. It recovers +K_in/Ne below the transition and the full-WJ reference continuously in the collar. These differ from the specified H27 tensor's full WJ diagonals at chi=0. Fixed loads must be regenerated and hashed, never evaluated from the current evolving state. Ground traction remains an output and the top remains normal impermeability plus prescribed tangential traction; the previously derived fixed finite-nu epsilon top flux is unchanged if its calibration inputs are unchanged.

## Concrete review and requalification package

The separately named `h27/proposal_principal.py` is **proposal_only**. It loads the exact retained inlet read-only, evaluates an independent full-3D WJ expression and candidate tensor, and screens physical r points, planar directions and 256 Fibonacci sphere normals with two transverse modes. It reports finite-difference step sensitivity, covariance and the exact chi=0 coercivity identity. Its finite positive coverage is not a continuum certificate or a coupled-flow result. The source and any output must retain their proposal-only identity; output is exclusive-create. It is prepared for an admitted Linux analysis run, never Mac execution.

Before any active implementation, the user must decide whether to adopt this exact scientific tensor change. After authorization, qualification must include:

1. A versioned amended model contract, regenerated reference normal pressures/outlet loads, strict old/new signature separation, and fresh initial states/history. Retain the exact imposed inlet and original physical geometry.
2. Independent analytic/native tensor-direction agreement; planar and full-3D covariance, trace, rotation, production and current-data tests; resolved collar principal-direction searches and adversarial states. The proof at chi=0 cannot replace these gates.
3. Rechecked reference roots, conditioning, endpoint fluxes, epsilon-H compatibility, same-inlet corners and complete reference/scalar budgets, followed by a separately resolved 2D steady/developed reference with finer h/p/q confirmation.
4. Complete coupled BDF2/scalar integration, all current residuals, independent continuous manufactured convergence, mixed boundaries/gauge, Q3-k pressure forcing and signed mass/scalar/mechanical-energy balances.
5. Genuine startup, developed-control and finite admissible perturbation pilots under the unchanged prospective amplitudes/support/physical envelope, plus matched dt and spatial/quadrature accuracy evidence. No amended tensor is qualified by quiet prescribed-coefficient loops.
6. Cold restart, actual recovery rotation/scientific retention, current-state monitoring, storage and full coupled one-case/campaign cost qualification. Keep T=505.02466521278404 and its original terminal window; do not launch that full horizon in this review package.

Until those gates pass on one executable recipe, the amended model remains unqualified and `TASK3_READY=NO`. If the collar screen fails, publish its exact witness and reassess this proposal rather than changing chi, coefficients, inlet, roughness or floors without a separately reviewed scientific basis.

## Observed proposal-only Linux screen

The frozen diagnostic was subsequently executed under the H27 account in
4.830688584013842 seconds. It reports `POSITIVE_FINITE_COVERAGE_ONLY`:
minimum sampled planar principal .004435791335539807, full-3D two-mode
principal .00512317614873607, and covariance/k .2817664872069159. Both FD
levels pass the declared step-sensitivity screen. The
[raw output](wj-h27-linux/publication-feasibility-v1/raw/proposal-screen-v1.json)
retains every point, direction witness and proposal-only flag; the
[closure report](wj-h27-feasibility-closure.md) binds its receipt and identity.
The active PDE/model remains unchanged and user scientific approval remains
pending. These finite positive results do not close any coupled H27 gate.
