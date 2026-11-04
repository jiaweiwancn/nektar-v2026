# M8-R9 — nonlinear inner destruction in the coupled direction model

The fixed comparison starts from the retained R8 iterate 10, an in-box
numerical guess with $E=2.911641510023811$, a continuous ground minimum
of $6.733016805\times10^{-6}$, and a scaled retained/projected trace
discrepancy of $4.242968946\times10^{-4}$. The startup history is the
original initialized $h_0$. The preceding matrix belongs to iterate 9.

The experiment uses the same 2,165-unknown outer6 representation and
quadrature, $\alpha=400$, target time .0025, $\eta=0$, zero forcing, and
SaturatedR2 with $\rho=.8$, $\beta=.4$, GammaScale 1, GammaUniform 0,
and the original outer $C(q)$. The viscosity conventions, scalar
diffusion and projected trace are unchanged. The bound XML explicitly
selects capped destruction with **Linear** smoothing, $\theta=.2$ and
$\delta=.1$. The native residual, the $2\times10^{-12}$ root gate, the
$10^{-7}$ action gate, and the $5\times10^{-12}$ original-row gate
apply unchanged. The numerical opt-in is `RANSWallR9Experiment=On`
with `DirectionExperimentR5`, `ScalarFeasibleR3`,
`VerifiedR3Action`, `JointRestorationR5` and `RefinedR6`.

## Direction model and comparison budget

The extracted $D$ assembles the native positive destruction load on free
$k$ rows after external lifting. It shares the native profile,
normalization and point source. Its action differentiates the local
limiter, sampled $K$, transition sensors, profile and moving quadrature,
including the test-basis evaluation and complete-ray $m$. The
normalization response is active on uncapped rays when $m=b$.
The source correction vanishes on algebraic and overwritten rows.

At the frozen base $q$, the model is
$$
M_q(d)=R(q)+J(q)d+D(q+d)-D(q)-J_D(q)d.
$$
The same $J(q)$ defines the matched linear direction and the native
linear prediction. Each source refinement uses
$A_j=J(q)-J_D(q)+J_D(q+d_j)$ and
$r_j=M_q(d_j)-A_jd_j$, with refined forward and transpose solves
for $A_j$. The constrained solve uses the existing consistent-mass
metric and restores $C(q)+J_C(q)d=0$. Radius .125, the 5% scalar and
continuous-wall reserves, and the outer cuts apply to the total
displacement from $q$ throughout. Strictly positive ground values
can move outward subject to the reserve and the native checks.
An invalid full tangent seed falls back to the accepted control
displacement. A refinement that fails to reduce the model objective can
supply a checked midpoint for the next tangent. Only full reconstructed
directions compete for the best model candidate, with at most three
refinements.

Both rays use native fitting and held-out samples at .125 and .0625,
the same R8 quadratic proposal and dyadic fallback, and the same
continuation decrease requirement $10^{-4}f$. At most nine native trial
states are available per ray. Source-only evaluations and component
actions have separate counts. Each actual candidate must pass the
finite scalar inventory, continuous ground trace, box and native
physical audits. The matching trace is evaluated on the actual trial
and contributes to $E$.

A source candidate qualifies for promotion only with at least 25% more
native residual reduction than the linear control, an absolute gain of
at least $10^{-4}E(q)$, a ground minimum at least as high as the
control, a scaled boundary/trace block no more than 10% larger, and
better residual gain per second after charging the model construction
and trial cost. These comparison thresholds select a numerical policy.
The native root and physical qualification criteria remain separate.

## Endpoint and destruction verification

The fresh complete physical J10 was assembled on the Linux endpoint.
The initial snapshot reproduced the R8 state, history, residual and
constraints byte for byte. The physical matrix and the destruction
action have separate hashes in the decision ledger.

The extracted value matched the native residual decomposition exactly
at the base and the .125 control trial. The independent quadrature
reconstruction of the base component differed by $2.19\times10^{-12}$
in the scaled $k$ block. Excluded-row support and $M_q(0)-R(q)$ were
exactly zero. At the base, the added value derivative and the subtracted
$J_D(q)$ cancel on the checked branch, leaving the physical $J(q)$.

| Destruction direction | Richardson step | Scaled absolute error |
|---|---:|---:|
| Coupled linear control | $2^{-9}$ | $1.79\times10^{-9}$ |
| Pressure part of control | $2^{-9}$ | $7.18\times10^{-11}$ |
| Selected source-model direction | $2^{-11}$ | $1.60\times10^{-8}$ |

The action checks used centered value evaluations at both step sizes
inside each Richardson pair with matching ordered branch and panel
inventories. At the larger step on the source direction the unextrapolated
pair differed by $2.53\times10^{-6}$. At $2^{-11}$ the difference was
$1.64\times10^{-7}$. The table uses the smaller check for that direction.

At the selected linear trial, retaining nonlinear destruction reduced
the $k$-block vector error by about 500 times. The remaining maximum
block error was in streamwise momentum. This check used one unique
native trial state. The progress stream logged it twice because the
selected control fraction coincided with the fitting fraction. The
reported source-evaluation cost includes the duplicate.

| Scaled vector error at the linear trial, $f=.125$ | Linear $R+fJd$ | Nonlinear destruction $M_q(fd)$ |
|---|---:|---:|
| Streamwise momentum | 3.80218e-5 | 3.80218e-5 |
| Normal momentum | 3.61164e-7 | 3.61164e-7 |
| $k$ | 6.15990e-3 | 1.23205e-5 |
| Outer epsilon | 9.41692e-8 | 9.41692e-8 |
| Essential and matching constraints | 5.98423e-6 | 5.98423e-6 |

At the base, $b$ and $m$ ranged from .8272669574 to 1. At .125 their
minima fell to .8272643807, with $\lambda=1$ in both evaluations. The
ordered moving volume inventory recorded one limiter switch and no
cap, profile-category or panel-identity switch. Maximum profile weight
changed from .9999461927 to .9999715631. Uncapped normalization and
profile changes therefore enter the finite-step source calculation.

## Bounded direction result

All three source tangents restored the full affine algebraic equations.
The largest algebraic block was below $2.22\times10^{-19}$ in their
respective tangent models. Both complete original-row scans passed for
each solve. The worst violation was $4.635\times10^{-12}$ across the
linear control and the three refinements. In both tested directions
the scaled coefficient step reached .125. The wall reserve was active
near $x=1.5639$, and both directions also had active retained and
projected outer constraints. The normalized QR multiplier of the wall
reserve was 2.96 for the control and 3.05 for the selected source-model
direction. The ledger gives the governing radius and outer rows with
their gaps and multipliers.

| Source refinement | Affine tangent $E$ | Full source-model $E$ | Full model differential objective |
|---|---:|---:|---:|
| 1 | 2.903479 | 28.775805 | 414.054649 |
| 2 | 2.903329 | 6.233355 | 19.458583 |
| 3 | 2.903305 | 6.557977 | 21.534675 |

Refinement 2 supplied the lowest objective of the three reconstructed
directions. Its full-model objective still exceeded 4.321036 at the
initial source-model anchor. The bounded construction resolved the
affine constraints and found a feasible direction for native testing.
It did not obtain a useful minimum of the nonlinear model.

| Ray and fraction | Linear $E$ | Source-model $E$ | Actual $E$ | Actual trace error / $E_0$ | Ground minimum |
|---|---:|---:|---:|---:|---:|
| Linear, .125, selected | 2.909883843 | 2.910252782 | 2.910253518 | 3.772440e-4 | 5.933475e-6 |
| Source, .125 | 2.974136634 | 2.976893288 | 2.976894318 | 4.206138e-4 | 5.933679e-6 |
| Source, .0009765625 | 2.911874171 | 2.911874231 | 2.911874231 | 4.238863e-4 | 6.726772e-6 |

The source predictions in the last two rows independently reconstruct
the destruction component at already tested states. They reduce the
$k$-block vector error from .013586 to $1.79\times10^{-5}$ at .125 and
from $8.54\times10^{-7}$ to $1.10\times10^{-9}$ at the smallest fraction.
All eight source-ray trials showed native non-descent. At the smallest
fraction the predicted and measured increases per unit fraction are both
about .2383, so damping that direction did not address its local ascent.
The scalar, wall and box observations passed at the listed source
trials. The smallest rejected fraction kept 99.907% of the base ground
minimum at $x=1.563841069$. The .125 linear trial kept 88.125%, reduced
$E$ by .04767%, and was selected and saved.

The linear R8 selector proposed the cached .0625 trial after its held-out
check, while the already evaluated .125 trial gave the lowest actual
$E$ and was saved. The source selector found no predicted decreasing
proposal and used the bounded dyadic fallback. Each ray had the same
nine-trial ceiling and selector. The linear ray used two evaluations and
the source ray used eight. No continuation was launched.

The native experiment took 155.65 seconds and exited normally. The
measured direction and trial cost was 21.80 seconds for the control
and 98.68 seconds for the source policy including its baseline
construction. The source work included 25 source-only evaluations
(10.03 seconds), four source actions (2.65 seconds) and three model
factorizations (.523 seconds). The independent finite-step comparison
used three offline term reconstructions and zero additional native
trials.

## Selected state and next decision

| Qualification | Independent result on the saved linear trial |
|---|---|
| Coupled root | Reconstructed $E=2.9102535179430$, above $2\times10^{-12}$ |
| Actual retained/projected trace | Scaled discrepancy 3.7724401502e-4 |
| Scalar and continuous wall | Finite inventory admissible. Cubic ground minimum 5.9334752793e-6 at $x=1.5638416679$, checked at 80 and 110 digits |
| Original box | In box. Max gradient departure .000999984943. Retained and projected epsilon minima .980000032867 and .980000032940 at $(1.75,0,2)$ |
| Covariance and lambda | Measured minima .000119099406 and .176877001374. Minimum $\lambda=1$ |
| Same-state accounting | Source budget $E_{\rm pre}-D-C=-1.94\times10^{-16}$. Integrated redistribution $1.94\times10^{-16}$. Weak-minus-direct energy $-7.82\times10^{-14}$ |
| Physical history | Native final report matches the selected state and residual. Original history slots verified. Zero accepted physical steps |

These are finite-inventory and wall-trace observations. The source and
energy identities concern the same retained non-root state. The ledger
contains the complete block norms, matching traction and flux
observations, selected-state audit, and immutable Linux payload
references. The isolated authoritative Linux build used 16 compiler
workers. Native runs used one rank and one numerical thread per library.
The three focused native CTests and nine predictor binding tests passed.

Keep the verified destruction correction and test a descent-safeguarded
model solve at this frozen endpoint before any new continuation. Require
monotone reduction of the evaluated source model from a feasible anchor
and verify descent of the complete native merit at the base before
spending trial evaluations. Preserve the frozen-base coupled constraints
and reserves. The current evidence supports the destruction curvature
model and identifies an unresolved direction construction.

M8 remains open. Curved-wall MMS, active-cap and variable-height
coupling, full-element and fully discrete positivity, refinement,
long-time robustness, scalability, and physical floor/hill validation
remain outstanding.
