# H27 CP1: source19 coupled numerical qualification

H27_COUPLED_NUMERICS=PASS for the pinned source19 numerical implementation.
The independently assembled current equations, exact A→B→A return, inverse
reuse, typed retry and atomic history failure tests complete the outstanding
numerical-foundation requirements. The source19 physical reference, pilots,
production operations and complete cost remain separate gates.

The diagnostic uses the same public C++ Coupling implementation and static
libraries as the registered solver. The application ELF remains
`37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b`
and its numerical source identity is
`56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf`.
The probe ELF is
`fa2bf63e186fce3dd47293ef11b58be6439967c34d4d62beaf8fc4d9268584ca`.
The approved v2 model and reference table are pinned by CP0. The probe calls
Coupling directly. The separate source19 normal-main MMS and physical-equation
reports exercise the registered entry point. [CP0](wj-h27-continuation-cp0.md)

## Current actions, live inverse reuse and rejected histories

A and B are saved accepted native histories at steps 10 and 18 of the
eight-cell Historical MMS crossing fixture, with dt=.02. Their two-sided
cutoff intersection moves from x=.504795249975577 in owners 1/3 to
x=.3720752367119697 in owners 0/2. Their independently sampled changes in
transported k/Kin, K/Kin, cutoff/H, local Ew and local outer epsilon are
.00293978, .00227524, .000144400, .00416360 and .00381244, respectively.
The original-owner velocity, gradient and essential boundary data also
change by the admitted amounts. [Admission](wj-h27-coupled-return-admission-v1.md)

The warm instance reconstructs the saved histories and evaluates the
complete current zero-time-mass equations at each exact input and clock.
The largest difference from independent 20×24/24×32 reconstructions is
4.694194908508949e-14 in fixed scales. Every identical input has an
identical native residual norm in the cold evaluator. The direct and
assembled native all-row differences in flow, k and epsilon are at most
4.95e-19, 1.70e-17 and 8.78e-18. Boundary and dependent trace gates pass.
The frozen-clock residual is evaluated as a diagnostic at a time-dependent
state. The actual replays use the full BDF2/EX2 equations.

The first B Newton base reuses the inverse at age 3 and measured drift
.0015778957872658265. A returns with exactly the saved four-history payload,
reuses at age 6 with zero drift, and later refreshes at the age limit.
Across the four matched replay pairs, all seven coefficient/action history
roles, original-owner fields, lifts, derived K/Gamma/c/Ew and weighted
all-row advection duals meet 1e-8. The largest reference-scaled coefficient
difference is 6.328355421304852e-14. No pressure alignment is applied.
All six independent complete current-step audits pass with maximum
momentum 3.20e-12, continuity 2.02e-17, k 1.61e-13 and epsilon 1.68e-13.
The predeclared limits are 1e-9 for momentum/divergence and 1e-8 for scalars.

After B, a correctly checksummed A binary carrying B's explicit action
fails the native history reconstruction guard. The original input and
accepted histories remain byte-bound. The frozen prior-scalar/source mutation
is detected on all six accepted replays. The injected typed linear failure
triggers one retry with a new factor at every retry base. Its incomplete
engine count flag is preserved. A one-Newton-iteration failure receives
zero retries and writes an identical input and output history.
[Independent report](wj-h27-linux/continuation-cp1-source19/science/coupled-return-v19r1/independent.json)

The transient Newton data memo is local to one attempt at one time.
Its key is the complete free flow/pressure/k/epsilon vector with exact
equality. Essential velocity and k and dependent Gamma are rebuilt by
`FromFree` from that vector, the scientific clock and fixed topology/table.
Finite-difference points cannot overwrite the validated current-state memo.
MMS loads are evaluated at the requested clock. Physical lift/load reuse is
confined to the instance's fixed physical boundary contract. The advancing
residual adds the supplied full history and computes current actions.
The separate live inverse uses drift/age/alpha, and the final accepted
residual and full principal guard use freshly evaluated current data.

## Numerical requirements and verified evidence

| Requirement | Coverage |
| --- | --- |
| Independent small actions, mixed boundaries, true pressure gauge and Q3 k sign | Native signed action tests against the independent assembled mixed matrices and mutation controls. Current path and tested dependencies matched through the source reuse audit. CP1 also compares all current native rows and independent boundary actions. |
| Evolving positive coupled continuous MMS and tiny direct solution | Historical and PolynomialScalarsV1 source19 normal-main cases retain separate exact/current cutoff roots. Q8 tiny direct nonlinear comparisons pass 1e-8. Both independent bulk and sampled-H body contributions are tested. |
| Temporal, spatial, quadrature and algebraic convergence | Source19 three-level temporal U/k/epsilon orders exceed 1.8, with spatial, q and tolerance controls. The steady nested study resolves Q3 field and gradient rates. Source19 retains mild pressure errors separately. |
| k, epsilon and mechanical accounts, common traces and owner flux | Independent budgets distinguish residual/reaction identities from physical gradient flux. Steady spatial normalized absolute k jumps fall 9.11e-5 → 1.57e-5 → 3.27e-6 and epsilon jumps 5.92e-3 → 2.13e-3 → 6.71e-4. Traction-power jumps fall 3.83e-8 → 4.03e-9 → 3.43e-10. |
| Current cutoff, admissibility and adversarial inputs | Source19 cutoff panel rows, moving owners, inverse certificate and no-root/multiple-root/stale/invalid controls. The unchanged WJ/planar implementation's independent value, covariance, principal and derivative negative controls are byte-bound. The final native guards cover the supported planar application. |
| Complete current coupling, cold reconstruction, cache return and failed step | CP1 six independent equation checks, exact A/B audit, warm inverse across two owner columns, stale history rejection, complete typed retry and atomic iteration-limit rejection. Normal-main source19 physical current equations and cold retained continuation independently pass on their diagnostic recipes. |

The [machine-readable coverage](wj-h27-linux/continuation-cp1-source19/numerical-coverage.json)
binds these reports and the exact source policy files. On any subsequent
numerical source change, review and requalify the affected paths before
transferring this decision to the selected production recipe.

## Timing and resource account

| Case | Native advancing solve s | Newton updates | Krylov iterations | Block factor builds |
| --- | ---: | ---: | ---: | ---: |
| First A | .349992 | 3 | 5 | 3 |
| B with reuse | .346703 | 3 | 5 | 0 |
| Return A | .348394 | 3 | 5 | 3 |
| Injected retry A | .385801 | 3 | 5 | 12 |
| Cold A | .349724 | 3 | 5 | 3 |
| Cold B | .349132 | 3 | 5 | 3 |
| Rejected A | .135117 | 1 | 1 | 3 |

The full native operation is 5.250963 s including setup, replays, exact
audits, binary I/O and export. The eleven qualification operations and
publication total 124.91539241705323 s. No per-step subtimer is added twice.
The earlier controller import failure occurred before the operation runner.
Its command failure, stale/rejected-input results and all new receipts are
preserved. Peak observed workers were 12, process tree RSS 1,012,813,824 bytes.
One science rank and thread limits of one applied, with four compiler workers.

The archive has 248,372 bytes and SHA256
`ef793337cb318774d09668cd5ff5ad1da9d2542441d191a48ff84d12d5743173`.
All 66 inventory files were verified. Exact binaries and omitted raw native
reports are bound to the Linux archive and original run paths. The replay
driver and frozen command receipts support reproduction in a new namespace.
[Inventory](wj-h27-linux/continuation-cp1-source19/inventory.json)

The same aggregate account remains OPEN at 17116.471991412807 seconds,
with 22483.528008587193 seconds available outside the 3600-second evidence
reserve and zero reservations. The physical reference still exceeds the
inner-k envelope and its spatial convergence is open. TASK3_READY=NO.
[Final receipt and account](wj-h27-linux/continuation-cp1-source19/post-publication.json)
[Launch hold](wj-h27-task3-launch-hold-v2.json)
