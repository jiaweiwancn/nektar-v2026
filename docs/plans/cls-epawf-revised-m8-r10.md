# M8-R10 safeguarded optimization of the destruction model

R10 completes the bounded comparison at the R9 initial state, the physical
R8 iterate-10 endpoint. Model globalization is monotone and the resulting
direction satisfies native-merit descent. The new policy fails the
unchanged practical promotion gates. The linear control is retained,
and the conditional continuation is omitted. M8 remains open.

The reviewed HEAD was 3b906f7fd2d2e657526ae291165ab2e4d6b0e901.
The remote branch had no intervening commits before implementation.
The producing Linux source is
6d6a92bb0f9b4e7102d9c3bdefd5066d0a1bba07. Later report and analysis
changes leave that producing source and build preserved. Development
used the existing branch worktree at /private/tmp/nektar-m8-r8. The
separate Picard checkout and the shared Linux checkout's existing
modifications were preserved. The Linux data volume was restored
before the isolated build.

The compact [decision ledger](cls-epawf-revised-evidence/m8-r10/decision-ledger.json)
pins the executable, XML, runtime, commands, every accepted model anchor,
native proposals, retained state and independent audits. Large payloads
remain below
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r10/epoch-001.

## Frozen base and compatibility

| Bound quantity | SHA256 |
|---|---|
| R9 initial state, physical R8 iterate 10 | eb063e15f16b2610f1c5a832280d87894f567e8dab321560202e6dd9af95194c |
| Current physical matrix, R9 solve1.J0.bin | 2194965dee9c5960ba4e1c1604f76e8e85372cabef2d5ced7c8a394f593c8dc6 |
| Original initialized history | b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71 |

The matrix remains associated with the initial state. The R9 saved linear
trial has a different state hash. R10 reuses the bound physical matrix,
D(q) and J_D(q). The initial residual and constraints reproduce their
R9 bytes. The new executable also reproduces the R9 linear direction,
.125 trial state and trial residual byte for byte.

Compatibility checks retain identical residual, Norms and destruction
value function bodies, plus the complete nonpressure-action, transition,
robust-source and residual-resolution source files. The private geometry
library was rebuilt in its isolated epoch and has a different binary
hash. Its 80 SpatialDomains source inputs match R9. All other changed
runtime artifacts are the rebuilt solver and two probes. Current and
historical runtime fingerprints are retained.

The target uses 2,165 outer6 p3/p1 unknowns and the existing quadrature,
alpha 400, target time .0025, eta zero, original h0 and zero forcing.
SaturatedR2 retains rho .8, beta .4, GammaScale 1, GammaUniform zero and
the original outer response. Viscosities, scalar diffusion, transition,
samplers, profiles, projected trace and capped source stay fixed.
Destruction smoothing is Linear, theta .2 and delta .1. Total radius
.125, 5% scalar and continuous-wall reserves, the disabled positive-value
inward rule, original outer box and 3,000-unknown guard are retained.
Root, action and original-row gates are 2e-12, 1e-7 and 5e-12.

## Opt-in direction policy

RANSWallR10Experiment=On selects the comparison with
JointRestorationR5, RefinedR6, ScalarFeasibleR3 and VerifiedR3Action.
Earlier selectors retain their paths.

For total displacement d from q, the verified model is
$$
M_q(d)=R(q)+J(q)d+D(q+d)-D(q)-J_D(q)d.
$$
Its objective is half the sum of squared native free differential
consistent-mass dual norms. At accepted anchor a, the affine tangent is
$$
A_a=J(q)-J_D(q)+J_D(q+a), \qquad r_a=M_q(a)-A_a a.
$$
Refined forward and transpose solves apply to this actual A_a. Every
tangent and anchor retains full frozen-base algebraic restoration,
including the nonzero dependent-trace equation.

The full linear direction initializes both incumbent and best only after
the original native and extended row scans, affine restoration, slope,
actual source-domain and continuous-ground checks pass. All constraints,
radius and reserves continue to refer to q. Backtracking evaluates
a+eta(d_star-a) with up to ten dyadic fractions. Actual source-model
Armijo decrease uses c=1e-4. Rejected trials preserve the incumbent.
Every accepted intermediate point competes for best retention.

The native slope is computed using physical J(q) and the native
mass-dual scaling. The helper handles active block ties conservatively,
including boundary and gauge blocks. The bound is
g_E(d) <= -2e-4 E(q) = -.0005823283020047622.
A source-valid model-Armijo trial blocked by this bound can trigger
one retry with the slope included as a general affine QP inequality.
That row has a null spatial witness and is scanned in original units.

The R9 ascent regression gives +.236151665288554 independently and is
rejected before native backtracking. The current control gives
-.014934753521388. Independent Cholesky reconstruction confirms all
accepted-anchor slopes. A .125-scaled native displacement leaves affine
defect .000371259782791366, exactly (1-.125) times the base trace block,
and fails restoration as a model anchor.

## Model verification and accepted anchors

The tangent identity r_a+A_a a=M_q(a) has largest scaled difference
9.31e-18. Algebraic correction support remains zero. Across the four
tangents, s_Phi+||A_a p||_W^2/2 is -.0747625, -.0726830, -.0705890
and -.0685321. All satisfy the declared numerical consistency test.
Both complete original-row scans pass for the control, all tangents
and accepted anchors, with maximum violation 4.635e-12.

At the initial moved anchor, the source-only action check has error
6.19e-9. An accepted limiter-pattern change triggers a second check,
with error 1.66e-9. Both use matching ordered limiter, cap, profile and
panel inventories on centered Richardson stencils and meet 1e-7.
The final accepted anchor changes a branch but is not subsequently
linearized in this bounded comparison.

| Accepted model anchor | Accepted eta | Phi | Native-base slope |
|---:|---:|---:|---:|
| Initial verified linear direction | — | 4.3210357075 | -.01493475352 |
| 1 | .03125 | 4.3190921048 | -.008630059813 |
| 2 | .03125 | 4.3171253688 | -.003575275636 |
| 3 | .015625 | 4.3154992235 | -.001636016409 |
| 4, slope-constrained retry | 1 | 4.2992543536 | -.0005823283020 |

The first full candidate gives Phi=414.054649 and its midpoint gives
19.604809. Both are rejected and the verified initial incumbent remains
in place until a decreasing fraction passes every gate. Runtime records
and compile-time regressions cover incumbent retention, rejected
midpoints, nonzero-base restoration along feasible segments, scaled
native displacement rejection and the R9 ascent rejection.

The third tangent has a source-valid Armijo point at eta=.03125 with
Phi=4.3150060181, but its native slope is +.0003032428173. This supplies
the declared trigger for the one slope-constrained retry. The retained
smaller step is feasible, and the fourth tangent includes the slope row.

For the recovered general row, the native, extended coefficient and
transformed gaps are respectively 4.95e-16, 5.86e-16 and 4.44e-16 in
original units. Its normalized multiplier is .05894. The continuous-wall
reserve, projected outer epsilon and coefficient radius are also active.
The wall reserve multiplier is 2.9684, with the governing point near
x=1.563842. These are direction-model restrictions with their existing
physical and numerical roles.

The comparison stops at its four-tangent budget after continuing
monotone model progress. It does not establish constrained stationarity.

## Matched native trials and promotion

Both rays use the existing R8 selector, designated fitting and held-out
fractions .125 and .0625, quadratic proposal grid, dyadic fallback and
native 1e-4*f sufficient decrease. The selector does not test f=1.
Native trials retain the actual trace residual, complete scalar
inventory, recomputed continuous ground extrema, original box,
covariance and lambda checks.

| Ray | Fraction | Actual E | Outcome |
|---|---:|---:|---|
| Linear fitting sample | .125 | 2.910253517943 | Selected control |
| Linear held-out sample | .0625 | 2.910834135384 | Eligible, weaker |
| Source fitting sample | .125 | 2.911671604342 | Residual increase |
| Source held-out sample | .0625 | 2.911629893491 | Decrease below the required sufficient-decrease margin |
| Source model proposal | .04541015625 | 2.911628012134 | Accepted proposal |

The source held-out sample decreases E, but its reduction is smaller
than 1e-4*f*E(q). The actual .04541015625 proposal passes that rule.
It is distinct from the linear state selected and saved.

| Compared selected ray state | Linear | Source proposal |
|---|---:|---:|
| Actual reduction in E | .001387992081 | .0000134978894 |
| Ground minimum | 5.933475279e-6 | 6.442557036e-6 |
| Ground retained relative to base | 88.1251% | 95.6860% |
| Actual trace block | .0003772440150 | .0004052718817 |
| Construction and trial cost | 31.18 s | 163.80 s |

The source proposal has 102.83 times less residual gain and 5.25 times
the cost. It passes the ground and trace comparison gates. It fails
the 25%-greater-gain, absolute-gain and gain-per-second gates. No
threshold was changed after observing these results. The conditional
continuation is omitted.

The comparison phase took 168.66 s. It used four tangent/QP solves,
20 source Armijo fractions, 30 new source values and nine cache hits.
There were four source actions and five total native trials, two for
the control and three for the source ray. Source values took 12.78 s,
actions 2.72 s, tangent factorization and constrained solves 106.93 s,
and source native trials 6.27 s. The native process exited normally
after 205.20 s including initialization and final reporting. The 600 s
cooperative comparison and 650 s supervisor limits were available.

The isolated Linux build used 16 compiler workers. Numerical runs
used one rank and one thread each for OpenMP, OpenBLAS, MKL, BLIS,
vecLib and NumExpr. The three focused native CTests and nine predictor
binding tests passed. The independent retained-state audit took
248.34 s. Read-only source-proposal accounting took 179.08 s and
created no additional native search trials or model-construction values.

## Independent retained state and conclusions

The selected linear state and residual reproduce R9 byte for byte.
Its independently reconstructed E is 2.91025351794298, above 2e-12.
The scaled retained/projected trace discrepancy is .000377244015022.
Weak continuity and gauge blocks remain at roundoff.

The original finite scalar inventory is admissible with zero unresolved
signs. Its minimum inner k is .0006027608924. The separate continuous
ground polynomial has minimum 5.9334752793e-6 at x=1.5638416679,
confirmed at 80 and 110 digits. These observations retain their finite
and ground-trace scopes.

The fixed outer box passes. Maximum normalized gradient departure is
.0009999849430. The retained and projected epsilon minima are
.9800000328669 and .9800000329396 at (1.75,0,2). Measured covariance
minima are .0001190994061 and .1768770013736, and minimum lambda is one.

At the retained state, the source accounting defect Epre-D-C is
-1.94e-16, integrated redistribution is 1.94e-16, and same-state
weak-minus-direct energy work is -7.82e-14. At the accepted source
proposal those quantities are -2.13e-15, 3.52e-16 and -2.31e-14.
These identities are evaluated on their respective actual retained
coefficients. They describe accounting at non-roots.

| Conclusion | Result |
|---|---|
| A. Model globalization | Verified monotone, with consistent tangents, checked moved actions and retained feasible anchors |
| B. Native-merit descent | Verified slope safeguard and an actual sufficient-decrease source proposal |
| C. Practical progress | Source policy fails promotion against the matched linear control |
| D. Root and physical advancement | No coupled root. Original startup slots verified. Zero accepted physical steps |

The narrow observed limitation is inadequate native gain and cost after
the four-tangent budget, with the native-slope safeguard active at the
final model anchor. Further model progress remains unresolved. These
results do not establish physical instability or absence of an
admissible root.

For the next bounded test, retain R10's safeguards and include the
verified native-slope row from the first tangent at this same physical
base. The single constrained retry accepted a full model step and
delivered most of the objective decrease, while the first three
tangents required substantial damping. Judge that test using the same
native gain, depletion, trace and cost gates before any continuation.

M8 remains open. Equilibrium preservation, constitutive applicability,
nonlinear root convergence, full-element positivity and long-time
physical robustness retain separate status. Curved-wall MMS, active-cap
and variable-height coupling, full-element and fully discrete
positivity, refinement, long-time robustness, scalability and physical
floor/hill validation remain outstanding.
