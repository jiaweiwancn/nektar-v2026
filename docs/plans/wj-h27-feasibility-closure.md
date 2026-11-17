# WJ-H27 feasibility result: selected reference/startup rejected

**Positive H27 closure was not achieved. TASK3_READY=NO. The requested positive-closure goal remains open.**

The specified wall tensor, selected finite-nu reference and required
inlet-derived startup fail a necessary constitutive principal-admissibility
test. Independent algebra, full-tensor finite differences and the actual
native C++ adapter/directional kernel agree. A scientifically derived outer
reference flux was investigated and does not repair this candidate. This is
a scientific prerequisite failure; coupled solver, physical-pilot, restart
and production-cost gates have not been represented as passed.

| Required gate | Result | Evidence / outstanding work |
|---|---|---|
| H27_MODEL_AND_REFERENCE | NOT_MET | Negative applicable 2D principal directions in the selected reference and retained-profile startup |
| H27_COUPLED_NUMERICS | NOT_MET | Registered transported-scalar coupled solver and its verification not implemented |
| H27_PHYSICAL_PILOTS_AND_ACCURACY | NOT_MET | No same-inlet developed 2D reference, physical advancement, or qualified h/p/q/dt |
| H27_RESTART_MONITORING_STORAGE | NOT_MET | No coupled cold restart, recovery rotation or terminal-retention qualification |
| H27_COMPLETE_COST | TARGET_NOT_MET | Complete coupled cost is unmeasured; affordability has not been established |
| TASK3_READY | NO | Every prerequisite must pass on one executable recipe |

## Checkpoints and tested source

| Checkpoint | Exact commit | Remote verification |
|---|---|---|
| H27-A: scientific contract, prospective gates, reference/BC audits and native feasibility tools | `529eb9337bde9eb1be3606160f0ee1617ca4b54a` | Exact branch ref verified before continuing |
| Explicit scientific amendment proposal and proposal-only screen | `86b8f0d2b3493a147b93ebeee475e7e313eca63e` | Exact branch ref verified |
| Final raw evidence, account and this handoff | Commit containing this report | Exact final remote SHA is reported with delivery |

H27-B/C/D were not reached. Their names are not assigned to diagnostic or
proposal checkpoints. Reviewed parent `ec0b2f1ba0de65ffabfd41001cc0e0a423ddd9db`
was preserved. Every committed change through the proposal checkpoint adds
new files; historical executable/model paths and defaults were not modified.

Authoritative execution was on `jay-Super-Server`, in
`/home/jay/build/wj-h27-coupled-20261117`. Native source capsule `source-v3`
was built with four compiler workers and byte-bound afterward to the H27-A
source. The audit verifies 1452 inherited current native inputs and 20 exact
Mac-to-executed source/calibration bindings. The unchanged H26 capsule's 4334
files and all three older closed accounts were rehashed successfully.

Executed native ELF SHA256:
`85860a03bff3e42d077dfc853fc6f4c3e3c3c411a246591b16acd06a11d54b59`.
This executable is **WJH27PrincipalProbe**, a constitutive diagnostic, not a
coupled EquationSystem or a Task3 executable. Compiler/configure/link,
source, loaded-library and input identities are retained in the evidence.
Analysis used Python 3, NumPy 1.21.5 and SciPy 1.8.0. Every science/analysis
process was serial; all six numerical-library thread controls were one.

## Exact candidate and independent reference results

The [contract](wj-h27-contract.md) preserves finite molecular viscosity in
both domains, physical roughness/height and the entire historical inlet.
Use full selected WJ, physical root, `DevStrainRawWorkV1`, raw `P=-R:grad U`,
independent `Cd=572/6561`, complete symmetric-strain diffusion/remainder and
the specified wall tensor retaining WJ normal blocks. Chi is zero below .4H
and full WJ above .75H with the H26 C2 quintic between them.

The candidate selects positive inner scalar diffusion
`nu+Cd*k^2/(sigma_k*Ew)`, unamended `P-Ew`, amplitude
`A(K)=K^(3/2)/cl`, and `Ew=A/max(r,c)`. No EP, clipping, floor, fitted
roughness change or actual-flow equilibrium forcing was used. Its explicitly
local counterfactual column has reference traction `K/qe`, `k(H)=K`,
`Q(H)=0`, and cutoff determined by `Q(0)=0`. It is not a stationary outer
solution or a substitute for the requested 2D developed application.

| Quantity at Kin | Independent Linux result |
|---|---:|
| c/H, original diagnostic | .19696577127569226 |
| c/H, tighter independent shooting/quadrature | .19696577127512963 |
| Finite-nu column velocity anchor U(H) | .42150691440063115 |
| Reference ground k | .014362547557221588 |
| Reference ground Ew | .011755971810747884 |
| Original independently integrated budget error / A | -3.7233798792697865e-12 |
| Refined independently integrated budget error / A | -8.955715018129626e-13 |
| Original analytic normalized cutoff slope | .9999999999879299 |
| Dependent epsilon(H) | .0023155240547992535 |

The retained inlet/interface scalar values are compatible with this unchanged
tail amplitude. At K/Kin=.9,1,1.1 the refined cutoffs are respectively
.20395270192101528, .19696577127512963, .19107351175357618. Budgets and weak
moments pass 1e-10, but the roots fail full principal admissibility. A finite
scan is not claimed to prove global root uniqueness.

An explicitly different auxiliary-inner-diffusion control gives
c/H=.1971067395536903 and a negative principal witness. Endpoint-flux
controls -0.05A and +0.05A give c/H=.2179270634901486 and
.17329993793722342, with the endpoint term included and budget errors below
1e-10. Neither control repairs the reference witness. Invalid K=0 and K<0
are rejected and retained in raw output.

## Principal failure and its scope

For unit in-plane wave normal xi and polarization a perpendicular to xi,
let B=a tensor xi. The necessary frozen-scalar momentum quadratic is

`B : [nu*(B+B^T) - D_grad(Rwall)[B]]`.

The scale is the current positive split coefficient, not time mass. At the
ground it is `mu=1.1490971129249645e-4`. The negative contribution comes from
the WJ normal-stress derivative retained below the collar. Positive mu and
positive covariance do not establish the sign of the complete tangent.

| State / positive wall distance | Native principal / mu | Minimum covariance / k |
|---|---:|---:|
| Selected reference, r=0 | -.5455586175395419 | .25577096566172025 |
| Selected reference, r=.00234375 | -.41704589380911133 | .26235515532882986 |
| Selected reference, r=.0046875 | -.28815828340789645 | .2697083516786046 |
| Retained inlet startup, r=0 | -.9350801177065614 | .25213809271992316 |
| Retained inlet startup, r=.0015 | -.8220007398190949 | See complete native table |
| Retained inlet startup, r=.003 | -.6901006539084942 | See complete native table |
| Retained inlet startup, r=.006 | -.4240283824677582 | See complete native table |

The native inlet ground quadratic is `-1.0744978636101457e-4`, at an
in-plane angle about .409856 radians. Across 23 negative and positive control
states, native values/directions agree with independent algebra to a maximum
scaled error of `1.683176636693416e-14`; actual-adapter finite differences
also agree. These are frozen constitutive checks, not coupled PDE residuals.
Strict negative values at positive r exclude a boundary-only artifact.
Pressure does not eliminate a divergence-free transverse velocity mode;
positive scalar diffusion and lower-order source/AWF coupling cannot repair
this second-order velocity symbol.

The [independent derivation](wj-h27-reference-audit.md) gives a fixed 30-degree
negative certificate without angular sampling. For every positive inner k
diffusivity in the declared zero-endpoint-flux family, it bounds c/H below by
.16192588845890177. The unchanged inlet then has q(0)>=27.28589600825159,
above the certificate threshold 18*sqrt(2). Its principal upper bound is
`-5.95285360749929e-5`, or `-.518046172124358 mu`. Adjusting only positive
inner diffusivity cannot cure this family. The existing bounded-destruction
amendment is inactive on these healthy k>=K profiles and does not supply a
permitted increase of Ew.

This is not an all-cutoff or all-possible-reference nonexistence theorem. The
51-cutoff diagnostic includes positive small-cutoff controls; as c tends to
zero, the ground tangent tends to a positive molecular/auxiliary form. No
cutoff-independent negative tail witness was found. Arbitrarily prescribing
such a small cutoff would discard the independently required reference budget.

## Scientifically derived nonzero-flux reference investigation

An independent stationary parallel **outer** BVP retains k(H)=Kin,
epsilon(H)=A/H, constant finite-nu total shear, zero top k diffusion and the
fixed finite-nu reference epsilon top flux. It supplies Q(H) instead of tuning
it from instantaneous actual production. Viscosity continuation starts from
the exact high-Re column; there are no state clamps or fixed retry counts.

| Initial log-radius nodes / tolerance | Final nodes | Nonlinear iterations at each positive-nu stage | Q(H)/A | Independent k budget error |
|---|---:|---|---:|---:|
| 101 / 2e-7 | 160 | 3 | -.039583975820884514 | -3.658813529927585e-11 |
| 201 / 2e-9 | 591 | 4 | -.03958397595157432 | -9.210451845653722e-12 |

Both use six viscosity stages, including the one-iteration zero-nu control;
all attempt records are retained. The flux refinement difference is
-1.3068980619923565e-10. The refined maximum sampled ODE defect is
2.24905948920974e-9. This is useful reference evidence, not the missing 2D
physical reference.

The resulting inner root is c/H=.21367837350095728, with independent budget
error -4.627680183499905e-12. The negative endpoint flux increases the cutoff;
its generalized lower bound is .16846410962722858. It therefore worsens the
same inlet mechanism instead of repairing it.

## Physical pilots, complete cost, restart and Task3

No accepted physical step was advanced. No new same-inlet developed 2D field
or physical adjustment to the inlet is available; counterfactual column
anchors above must not be relabelled as developed-flow measurements. No
manufactured source was used to conceal the obstruction.

There is no accuracy-qualified coupled dt, current coupled iteration table,
three-block complete-step timing, cold coupled restart, recovery write or
rotation result. Coupled recovery writes=0 and rotations=0. The feasibility
JSON/TSV retention is not scientific terminal-output qualification.
Per-trajectory and total Task3 campaign forecasts are **unavailable**, not a
claimed cost failure measured on a working solver. H26's .121939577 s/step
and kernel-only forecast are not extended to this task.

The required Task3 horizon remains 505.02466521278404 and terminal window
[447.8483832973273,505.02466521278404]. The
[launch hold manifest](wj-h27-linux/publication-feasibility-v1/task3-launch-hold.json)
contains null coupled executable/state/recipe/forecast/commands and explicit
pre-launch abort reasons. **No valid Task3 launch command exists at this
checkpoint.** No full Task3 horizon was launched.

## Proposed remedy and decision

The [separately named proposal](wj-h27-model-amendment-proposal.md) changes
the lower-wall tensor to `Riso=2k I/3+bA tensor n+n tensor bA`, blended with
full WJ using the same chi. At chi=0 its exact incompressible rank-one form
is `nu+nuA*|P(B+B^T)n|^2>=nu>0`. Its covariance still needs
`2k/3>=|bA|`. Parallel shear/production and the selected scalar reference
equations are identical; general production, pressure and normal tractions
change and require consistent rederivation.

A proposal-only Linux screen over the retained profiles, two FD levels,
planar directions and 256 sphere normals/two polarizations found positive
finite coverage: minimum sampled planar value .004435791335539807 and full
3D two-mode value .00512317614873607; minimum covariance/k .2817664872069159.
This is not a continuum certificate or coupled qualification. The active
model was not changed. User scientific approval is pending because the task
specifies a different tensor and forbids silently implementing different
physics. If approved, adopt a versioned new model and execute every remaining
B/C/D gate, including full collar/adversarial checks; this proposal does not
itself qualify them.

## Closed resource account and complete operation table

The H27 feasibility account is **CLOSED** at 703.6786409680499 /43200 seconds,
with 42496.32135903195 seconds unused. This includes a conservative initial
120-second control charge and a prepaid 300-second control/publication
provision; the actual sealing/publication script took 1.1236091719474643s.
Independent overlapping operations are charged separately; idle parent waits
are not counted as another scientific operation. Failed operations are kept.

| Operation | Result | Charged seconds |
|---|---|---:|
| Bootstrap/control | FULLY_CHARGED | 120.000 |
| First independent reference | PASS (execution) | 19.449 |
| First source materialization | FAIL, transfer not complete | .561 |
| Verified delta source-v2 | PASS | 1.817 |
| Native build-v2 | FAIL, print-helper name collided with Nektar::Matrix | 26.720 |
| Analytic planar/uniform-D bound | PASS (rejection certificate) | .631 |
| Exact retained-inlet scan | PASS (diagnostic execution) | 8.883 |
| Refined K=.9,1,1.1 reference matrix | PASS (execution; model rejected) | 115.730 |
| Auxiliary scalar-diffusion control | PASS (execution; model rejected) | 19.058 |
| Negative endpoint control | PASS (execution; model rejected) | 20.248 |
| Positive endpoint control | PASS (execution; model rejected) | 20.243 |
| Invalid-input fixtures | PASS (invalid inputs rejected) | .624 |
| Repaired source-v3 materialization | PASS | 1.343 |
| Native build-v3 | PASS | 26.261 |
| Initial independent native comparison | PASS | 1.371 |
| Independent outer BVP and inner follow-up | PASS (execution; remedy rejected) | 14.512 |
| Native comparison with positive-width witnesses | PASS | 1.397 |
| Proposal-only finite-direction screen | POSITIVE_FINITE_COVERAGE_ONLY | 4.831 |
| Final control/publication provision | CONSERVATIVELY_PREPAID | 300.000 |

PASS in an operation receipt means the diagnostic executed successfully;
scientific gates use the separate matrix at the start of this report. The
first materialization ran before its source transfer completed; that failed
receipt is retained. The slow transfer was stopped, and a small verified
delta replaced it. The native compile collision was fixed by renaming a print
helper; the failed capsule/log remains intact. No numerical target or physical
equation was changed to repair either operational failure.

Measured maximum actual shared-host workers=38/72; maximum admission bound
including helper reservation=39/72. Peak build RSS=1694711808 bytes (<16GiB),
analysis RSS=65740800 (<8GiB), source/build/runtime=213291374 (<8GiB), raw
science=1045335 (<8GiB), free ext4 space=75730812928 bytes (>15GiB). No
unrelated job was stopped. General Linux source remains f4b70baec10016dd098005bc5a5507b87e4a57bf
with the original developer-guide modification. The Mac branch worktree
preserves its unrelated untracked pilot file.

Closed account SHA256:
`f82e7b89d4d2e68ab177a5e9821703fad888be917cfecb5df3b06ef3565ae0c6`.
The 176059-byte transport archive SHA256 is
`48c9e3732705eaddbaeaed22197640d030db3d3c57ae9cf55d31fe8fe0a28611`.
All 73 transported artifacts and the publication/account manifests were
rehash-verified on the Mac. No operation outcome is unresolved.

## Evidence and continuation

- [Complete machine-readable summary and native/outer iteration tables](wj-h27-linux/publication-feasibility-v1/summary.json)
- [Source/runtime/preservation/resource audit](wj-h27-linux/publication-feasibility-v1/integrity-resource-audit.json)
- [Closed account](wj-h27-linux/publication-feasibility-v1/accounting-final.json)
- [All raw diagnostic results](wj-h27-linux/publication-feasibility-v1/raw/)
- [Every operation request, receipt and log](wj-h27-linux/publication-feasibility-v1/operations/)
- [Transport identities](wj-h27-linux/transport.json)
- [Reproduction commands and Task3 hold](wj-h27-reproduction.md)

Do not reopen the closed account or import historical flow histories as new
states. An authorized amended-model continuation can use a separate ledger
carrying forward the 703.6786409680499-second charge and original aggregate
43200-second ceiling. The preserved unused allowance is not an extra budget.
Positive closure still requires the full requested solver and qualification.
