# WJ-H7: bounded neutral transport, T=.04 only

Fetched origin and clean WJ worktree both start at reviewed
`a2642861740cdf1c9bf4f8adce44d15dd4b4a70d`. Other worktrees and the main
checkout's unrelated files are preserved. This contract precedes numerical work.
Use a new root `wj-h7-neutral-transport-20261113` on UUID5AB8A669B8A642FD,
verified at `/media/jay/新加卷1/nektar-v2026-runs`. Never modify H6 ledgers,
reports, checkpoint bytes or signatures. No C++ change/build is planned.

Freeze the three H6 references and their actual Function/mesh bytes. H_w is
5/5.5/6mm (.125/.1375/.15 code), h_s is respectively
4.381467178582103/4.359758967365407/4.341801417014184mm. Each has 16 inner,
34 outer cells, Nx1 on [0,8], Q3/Q1, 2328 states, q12/normal16 and the
3000-state guard. Keep WJNeutralBudgetReferenceV1, Original full3D WJ,
DevStrainRawWorkV1, physical epsilon, raw work and independent Cd diffusion;
Ne=81/20,qe=81/(2sqrt143),Cd=572/6561,kappa=sqrt(1.3*.38*sqrt(Cd)),
u_tau=.051,nu=.051*.005/3.97,d0=.0875,top20,K=.008809015145515192,
z0_e=.0056762703212168 (reuse exact stored H6 values). Positive wall nu,
nonconstant inner k, uncapped Budget/current-K cutoff and dependent A(K)/H
epsilon are unchanged. No inner epsilon or new flux/pressure constraint.
Reference Functions/physical coordinates, CompositeTunnelV1/OpenTunnelV1,
complete outlet and ReferenceTractionV1 top data, IMEX1/PimpleSingleSweep,
eta0,ExactStepV1,PimpleVerifyOff,StopOnBound and IO_CheckSteps1 remain fixed.

Specifications retain SHA256 main
1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124 and
appendix0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c,
including the quad/\\quad typography erratum and (d+a*h)/a^2 coefficient.

The ONLY matrix is 5.5mm pair first, then6mm, then5mm. Per height reuse H6
.01 CP0 unchanged and run one normal four-step block. Create exactly one
fresh .005 CP0 (NumSteps0, no restart), verify state and fixed-reference
binary64 arrays bitwise against H6, then one normal eight-step block.
Max36 attempts, three fresh CP0 reservations, no retry, extra level/mesh,
segmentation or evolution beyond .04. Each XML clones its H6 parent; only
TimeStep (.005 branch only), NumSteps, restart and output locations may change.
Do not relabel checkpoint dt/history. Require zero Project/Lift on restart.
Inputs are write-once, hashed before launch; no cached-intent substitution.

Shared comparisons at integer-aligned times0,.01,.02,.03,.04 use128eps
max(1,|t|) time allowance, retaining actual timestamps. Assess every accepted
fine odd step and every available rejected candidate. Bind each state to its
genuine producing report, audit packet, callback/history hashes, time/alpha,
fixed reference and checkpoint. Never manufacture history from endpoints.
Native hard failure ends that branch; other predeclared branches may run.
An expanded-inventory missed negative or principal sign discrepancy within
2e-11 stops all unlaunched blocks. Preserve native decisions and distinguish
retrospective qualification from precommit admission. A quality-only loss
(margin<half exact) is reported QUALITY_LOST, not a changed physical guard.

Start each height's inventory with the complete H6 common physical union,
including old H5 faces, both owners, x0/4/8, all Fig.3/Fig.12 heights, H,
exterior/interface points and1e-6 physical neighborhoods. Add every new native
witness and its same neighborhoods, retaining all legitimate owners. Scan
the current union after each block, then replay all earlier states on the
final union for that height. Cache only immutable geometry/basis rows; every
state's U/G/k/epsilon/K/stress/cutoff is freshly evaluated. Batch independent
full3D algebra/action must agree with the existing angular calculation and
compiled action within2e-11 at representative/decisive witnesses.

Publish reference error, CP0 drift and same-point/time dt differences separately.
For outer E_G=max|Gxz/Gxz_exact-1|,E_epsilon=max|epsilon/epsilon_exact-1|,
E_U=max|U-U_exact|,E_k=max|k/K-1|. Inner U/k/R use their own finite-nu
reference; no inner epsilon. Report full physical G, transverse U, R,
reduced/mechanical pressure with fixed scales1,K,A/H,u_tau^2 and absolute
errors at exact-zero components. Physical L2 norms use common elementwise
Gauss8 quadrature and full physical Jacobian, split at actual element/physical
breakpoints, with outer-only epsilon. There is no added field-preservation
percentage gate or convergence-order claim. Cross-height outer comparisons
(if reported) use only r>maxH and are not mesh-error estimates.

Keep native admissibility, half-reference quality, and numerical evidence
statuses separate: PASS/QUALITY_LOST/OBSTRUCTED/UNRESOLVED/NOT_RUN,
RESOURCE_LIMIT for denied work. Exact in-plane margin is independently
computed (~.004137527846908); require half at every assessed accepted outer
owner. Transverse/scalar/covariance/root/cutoff/work checks remain distinct.
For complete attempts expect one U, pressure, k,epsilon solve and one bound
CALL, zero ACTIVATIONS/material corrections; preserve partial counts/failure.
Retain raw/corrected k snapshots and all actual step/callback histories.

Terminal accepted/predecessor strict open_checks uses the actual producing
report and checkpoint. A failed producer needs allow_failed_producer=true;
its rejected candidate remains ineligible. No accepted step means reuse a
qualified CP0 view, no synthetic stage. Report native frozen-block residuals,
genuine completed-stage recombined defect and alpha0 uncancelled spatial
residuals separately, with native consistent-mass dual definitions. Keep
signed weak components, actual versus counterfactual production, K/cutoff,
epsilon/traction/k-flux traces, weak continuity and strong divergence.
Independent q12/24/32 weak checks use actual projected exterior test traces;
unchanged component2e-9+2e-8|value| and regional4e-9+2e-8|sum| limits.
No adaptive q, governing diagnostic solves, projection or restoration.

Reuse verified governing producer c344dc2209f030e426e2ee5501d5e5b72a6d2c46
and observer e62b5f4d265f068dd5af3f090bfb0e385c3ae53e, source/ELF/actual
loaded-library checks and -ffp-contract=off -fno-fast-math. One rank;
OMP/OpenBLAS/MKL/BLIS/NumExpr/VECLIB threads1; inspect current workers and
reserve conservatively below72. A fresh persistent ledger caps1800s,
8GiB process-tree RSS and8GiB new output;60s setup and240s closure are
reserved. Operation caps:60s fresh CP0,120s native block/endpoint view,
180s scan batch. Failed/unknown reservations stay charged. Before subsequent
matrix blocks, use actual center-pair governing/every-state/endpoint costs
with25% remaining-work headroom and full closure reserve. Do not equate CP0
time with timestep cost; a nonfitting projection gives RESOURCE_LIMIT/NOT_RUN.

Rerun15 H5/H6 methods plus focused H7 matrix/dt-lineage/time/odd-step,
accepted/rejected callback, zero-safe metrics, bound accounting and status
fixtures. Keep heavy scalar inventories/checkpoints/history outside Git.
Only compact per-step results, witnesses, input/runtime/receipt identities and
plots are published. No historical all-run replay/new provenance framework.

The horizon is only .005 strip flow-throughs. The full-height log remains an
idealized control, not finite-depth inflow. The7mm observation remains inside
AWF and the whole-paper outside interval empty. No physical H is selected.
Old Blackadar/20mm, H4 backend/capacity, EP, hill and longer work stay paused.
H7 may recommend a separately bounded next decision; it authorizes none beyond
this matrix and never overwrites earlier denials or scientific outcomes.
