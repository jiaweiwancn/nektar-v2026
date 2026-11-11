# WJ-6 finite dense capacity and matched qualification

Reviewed parent: `6ee8f356d0cb31d9b99f86127d676ffa408d34e1`, fetched without
newer remote commits. The Mac WJ worktree is clean at entry; other checkouts
and live campaigns are preserved. WJ-4/WJ-5 outcomes, receipts and closed
ledgers remain immutable. WJ-5 remains OBSTRUCTED; local M1/M2 and the half-dt
control do not qualify the missing uniform L2/Nx2 cell.

Original/worktree specification SHA-256 values agree:

```
1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124 main
0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c appendix
```

Main eq:integrals-linear retains its literal `quad` / `\quad` typography
erratum and the correct `(d+a*h)/a^2` logarithmic coefficient. No specification
or model change is authorized. Exact source/tree/content, runtime/dependencies,
contract and XML identities are frozen in write-once run designs before use.
An observer repair requires a separate code/output lineage.

## Fixed problem and execution option

The complete 3D WJ tensor and prescribed cubic-N surrogate, Original stress,
DevStrainRawWorkV1, direct physical epsilon, Praw=-R:Graw and independent Cd
scalar diffusion remain unchanged. Constants are c1=1.8,c2=5/9,Ce1=1.44,
Ce2=1.82,sigma_k=1,sigma_epsilon=1.3,Cd=572/6561,qe=81/(2 sqrt(143)),
kappa=.381917863546733. Molecular viscosity is omitted only in the outer
high-Re model. Positive nu=1.5e-5 remains throughout the AWFs; inner stress
is full-strain Linear and destruction is uncapped.

EP uses chi=1,zm=H/2,source P_A-p_A,0; budget uses zK=zm=H,source P_A-E_c,
no physical chi and transported nonconstant inner k. K=k(H) is an endpoint
sample. There is no inner epsilon PDE or extra epsilon flux. Canonical
u*=.5,H=1,top=2,h=.1,K0=.8466950351321795 are fixed.

The finite strip [0,1]x[0,2] retains reference differential Dirichlet traces
at both laterals/top, no-slip ground and zero modeled ground k flux.
Q3 velocity/scalars, regional continuous Q1 pressure, one overall gauge,
pi=p/rho+2k/3, shared U/k, separate pressure/gradient traces and one dependent
outer epsilon trace remain fixed. IMEX1/PimpleSingleSweep retains its one
U/pressure/k/epsilon solve and logical WJNativeV1 MassProjection per completed
attempt, eta=0, uncondensed, zero body/scalar forces/gravity, reuse/verification
Off, StopOnBound=On. Qualified runs require identical raw/corrected k and no
material correction. All native hard predicates remain active.

`RANSWallDenseEvolutionCapacity` is a SolverInfo execution option. Its default
is3000; the only extended range is finite exact integers3001..4331. Parsing
must reject fractional/nonfinite/out-of-range input before casting. Larger
capacity requires serial flat2D WJ PimpleSingleSweep, eta=0 and uncondensed
evolution. One shared effective-capacity policy guards both Patch::Solve and
the WJ one-sweep entry. Other routes keep3000. There is no constructor size
guard: large reference audits and strict readbacks retain their legacy use.

The opt-in binds a64-hex `RANSWallDenseEvolutionResourceDesign` SHA-256 to an
external finite resource-design file. Both controls are reported, included in
exact/semantic XML and requested intent, and excluded from physical checkpoint
identity by using execution-only SolverInfo fields. No checkpoint is edited.
Physical off/on signatures, accepted states/gauge/histories/audits and solve
counts must agree on a small otherwise identical case.

No pressure-space/continuity change, stabilization, stress/source clipping,
coefficient change, roughness smoothing, new boundary, governing correction,
general sparse backend, native3D/MPI or stationarity qualification is included.

## Mesh and observation identities

x knots are [0,1] for Nx1 and [0,.5,1] for Nx2. All normal knots are the exact
existing uniform L2 vector (two bisections of the original graded vector):

```json
[0.0, 0.000625, 0.00125, 0.001875, 0.0025, 0.003125, 0.00375, 0.004375, 0.005, 0.00625, 0.0075, 0.00875, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.07, 0.08, 0.08499999999999999, 0.09, 0.095, 0.1, 0.125, 0.15000000000000002, 0.17500000000000002, 0.2, 0.25, 0.30000000000000004, 0.35000000000000003, 0.4, 0.47500000000000003, 0.55, 0.625, 0.7, 0.7749999999999999, 0.85, 0.925, 1.0, 1.0375, 1.075, 1.1124999999999998, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4125, 1.475, 1.5375, 1.6, 1.7000000000000002, 1.8, 1.9, 2.0]
```

Expected native dimensions (verified by preflight and actual assembly):

| Mesh | State | Free velocity | Pressure+gauge | Assembled flow | Free k | Epsilon solve |
|---|---:|---:|---:|---:|---:|---:|
| L2/Nx2 |4331|1790|187|1977|900|242|
| M2/Nx2 control |2975|1250|133|1383|630|122|
| L2/Nx1 |2493|716|125|841|360|98|

Common physical observations use normal order8 and independent order10,
horizontal4/6, actual compared mesh unions, zc/zv/h/zm/H and x=.25,.75.
All actual normal/x face owners and intersection owners remain explicit;
duplicate maximum/face observations have zero integration weight. Quadrature
norm differences must fit the existing absolute1e-10 floor. Every field/scope,
the entire outer epsilon domain, one full-domain pressure alignment, sampled
K/cutoff and independent evolved source/diffusion/interface weak loads remain.
Pressure-test projection/complement, element modes, Gxx/Gzz and extrema are
retained. Any extra maximum-neighborhood view is passive and element-owned.

Spatial residual alpha_mass=0, actual completed callback defect, continuum
error, initial-representation drift and same-time differences remain separate.
Measurements use the actual checkpoint producer/callback, hash state/history
before/after, and never solve, Lift, reset, project U or advance history.

## Conditional run phases and gates

Each policy owns fresh native CP0 at each dt. Coarse dt=3.2e-5 and fine
dt=1.6e-5 use distinct signatures/time labels with exactly equal initial
physical coefficients, gauge and corresponding physical history values.
A fresh same-runtime L2/Nx1 coarse baseline is the conservative default;
old WJ-4/WJ-5 states remain context, not an unproved runtime-equivalence claim.

1. Capacity controls and35 existing unique WJ-5 regression selections, plus
   focused parsing/route/preflight/signature/transaction/metadata checks.
2. Per-policy actual4331 pilot: coarse four accepted steps, complete passive
   observation and frozen-system arithmetic/SVD. Fine pilot is eight steps to
   the same time. Pilot states count in the trajectory; no early Nx1 callback
   is fabricated. The fresh Nx1 baseline also executes its own first4 steps
   and genuine observation/SVD for measured mesh costs and a common pilot
   comparison; these4 are included in its32-step target. Native startup
   observations are required on every mesh/dt.
3. Freeze conservative measured cost projection before primary continuation.
   Primary target is Tx=.001024: coarse32 TOTAL, fine64 TOTAL, coarse Nx1
   baseline32 TOTAL. Retain CP24 and accepted24->32 exactly. Hard/material
   correction exits stop immediately. An engineering failure permits only this
   frozen Tx pair as explicitly adverse diagnostic, never T3 continuation.
4. One actual-capacity transaction per policy, conditional on hard-admissible
   uncorrected24->32 and resources: private CP24, rejection at absolute25 after
   scheduled solves, verify rollback and immutable original producer bytes,
   frozen readback, restart24->28 then28->32. Cost is9 extra attempts. Compare
   all9 accepted states, all8 incoming/callback histories, audits/samples,
   raw/corrected k, gauge/time/counts with an explicitly indexed slice of the
   full immutable producer. An abort need not publish a success report or CP25.
5. Only if BOTH policies pass every primary startup/pilot/cadence accuracy,
   horizontal/time sensitivity, quadrature/arithmetic/transaction gate and the
   remaining cost fits, extend SAME Nx2 histories to T3=.004096 (coarse128,
   fine256 TOTAL) and coarse Nx1 to128. No Tx step is repeated. No T4 is
   authorized. Stop optional continuation on any failed gate. Primary PASS
   with unaffordable T3 is reported separately from extension RESOURCE_LIMIT.

Accuracy cadence is no coarser than32 coarse/64 fine, including startup and
the actual4/8 pilot pair. Primary arithmetic is early4/8, middle16/32 and
end32/64. Optional extension adds required64/128 and128/256 packets, retaining
prior ones. A conservative coefficient envelope covers EVERY actual observation
row; each trajectory's uncertainty contributes. Preserve separate k/epsilon
units, dependent trace maps, unconstrained pressure/gauge, prescribed-data/lift/
basis/geometry/accumulation errors and factor100. These are finite frozen-system
estimates, not a nonlinear stability theorem. Missing/stale/nonfinite evidence
cannot pass; a late passing endpoint cannot erase an earlier failure.

Amplitudes remain aU=.5,ak=K0,aepsilon=.125/kappa,ap=.25. E2 and D2 divide
physical L2 by amplitude*sqrt(own field/scope measure). Unchanged limits:
E2/D2<=1e-4; normalized sampled maxima<=1e-3; full traction jump/.25<=1e-3;
k-flux jump/(K0*.5)<=1e-4; reference pressure-offset departure/.25<=1e-4;
sampled |div U|*H/u*<=1e-3. Horizontal/time difference PLUS both finite
uncertainties must fit10% of the corresponding field budget. No rate is
inferred from this two-grid/two-dt design or from below-floor differences.

## Finite resource design and terminal outcomes

New isolated root: `/media/jay/新加卷/nektar-v2026-runs/wj6-20261111`.
Per-policy persistent limits:800 charged native/conservative control units
and10800 numerical-study seconds. Accepted/rejected/development/control work,
observations, SVD, I/O and transactions count; unknown reservations stay charged.
Build time is separately recorded. Core with fresh Nx1 is128+9 at Tx or512+9
TOTAL at T3 before controls. A finite closure/transaction reserve is maintained
before optional extension, with at least25% conservative projection headroom.

Initial per-process pilot ceilings:600 s native four/eight-step invocation,
480 s full observation,600 s arithmetic/SVD,8 GiB peak aggregate process-tree
RSS,8 GiB invocation output and32 GiB cumulative policy output. Resource
monitoring may terminate only its own process group. Current host memory,
storage and thread inventory are captured before launch. Pilot timings/RSS/
volume and actual block sizes determine a write-once next-phase projection;
later native blocks use finite measured bounds within the same package caps.
No bound or target is silently raised. A ceiling yields RESOURCE_LIMIT.

Linux is authoritative: RANS_BUILD_SOLVER=ON, matching consumers, Release,
-ffp-contract=off -fno-fast-math, serial rank and numerical-library threads1.
Start with one compiler worker; any later allocation requires a current load
check and total scientific/compiler/MPI-helper workers<=72. Other campaigns
are preserved. Exact commands/exits/times/logs, source/runtime/mesh/XML/CP0/
producer/stage/history lineage and independent byte verifications are sealed.

Terminal component statuses are PASS, OBSTRUCTED, UNRESOLVED, RESOURCE_LIMIT,
NOT_RUN (implementation/metadata failure uses FAIL until repaired under a new
lineage). No overall qualification PASS with a missing/failed required gate.
WJ-4/WJ-5 negatives and stationarity NOT_REACHED remain unchanged. The scope is
a finite serial flat2D reference-Dirichlet strip, not exact discrete equilibrium,
general nonlinear stability/realizability, scalable production or long-time
homogeneous ABL validation.
