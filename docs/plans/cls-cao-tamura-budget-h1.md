# Cao–Tamura budget AWF H1: finite reference, two-region geometry and native compatibility

H1 implementation and focused Linux verification are complete. A fresh affine
rough-flat finite-profile case accepted **8 steps to T=0.00008 s** with no
positivity correction. Split/uninterrupted restart and deliberate rejection
rollback agree bitwise. The thin-strip/deep-outer geometry is a verified native
kernel; mapped hill governing evolution remains guarded.

Source: `1124e15677a1efa2fc8b07dae73ec0aca31f983e`, based on reviewed
`cb85f8a90fe4955cb620697bece49070cbf04b58`, on the unchanged
`codex/cls-epawf-pimple-one-sweep` branch. The remote had no intervening commits
when fetched. All 22 changed source/header/test files match the isolated Linux
source export by SHA256. Mac source and Git work did not modify the shared
Linux checkout, installation, other worktrees, or historical runtimes.

[Native evidence](cls-cao-tamura-budget-h1-evidence/native.json),
[source identities](cls-cao-tamura-budget-h1-evidence/source-revisions.json),
[case manifest](cls-cao-tamura-budget-h1-evidence/case-manifest.json),
[kernel/preflight results](cls-cao-tamura-budget-h1-evidence/kernels-preflight.json),
[strict readbacks](cls-cao-tamura-budget-h1-evidence/readbacks.json) and the
[artifact index](cls-cao-tamura-budget-h1-evidence/artifact-index.json) retain the
underlying numbers. Heavy fields, all native histories, checkpoints, stage
packets, frozen transport matrices and logs remain under:

```text
/media/jay/新加卷/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111
```

## Preserved baseline and experiment data

The completed refined P1 steady candidate still qualifies through T=4; its
actual retained checkpoint is step 6432/T=4.02. Its completed small inner-k
response test has amplification .591594 at elapsed .5, with zero correction.
Those long trajectories were reused. A one-step native regression against the
frozen prior runtime confirms identical committed/post-flow/pre-k/raw/corrected
vectors, callback history, reference, defect, audit values, accepted state,
complete native histories and time labels (12 checks).

The separate P2 campaign remains under its existing first-three-passing-windows
or T=6 stopping policy. Its dated live snapshot is in
[latest-status.json](cls-cao-tamura-budget-h1-evidence/latest-status.json);
supervisor, driver, native child, observer and executable identities were checked.
No new P2 continuation, runtime replacement or cap extension was made. H1 does
not depend on its completion.

At **2026-11-11 22:04:43 Asia/Shanghai**, the verified P2 frontier was step 5619,
T=3.511875, running toward step 6400/T=4. Its supervisor 626365, driver 626371,
native child 752066 and observer 689309 matched their recorded start identities.
The observer had assessed endpoints 4000/4800/5600 and remained active. The
immutable copied checkpoint has SHA256
`0116e918ecbe8012b7e5e668f61cfbcf7f631a78bd9389c053ac2e790917eaa2`.
This is a dated live frontier, not a completed P2 qualification claim.

The frozen `cao-tamura-2006-v1.1` dataset passed the Linux auditor: 98 series,
1,151 points, 39 scalar facts, 1,133 visible points and 18 zero-weight occluded
points. All nine snapshot file hashes match. The local 497,768-byte paper has
the required SHA256
`adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df`.
No digitization was repeated. Visibility, original normalized height weights,
coordinate conversions, extraction envelopes, source ambiguities and holdouts
are unchanged. See the [data receipt](cls-cao-tamura-budget-h1-evidence/data-audit.json).

Only rough Fig.3a/b/c enter this adapter. Fig.3d and every hill quantity remain
reserved. Paper vertical y becomes solver z; the frozen Fig.3c conversion already
includes assumed displacement d0. The displacement is not subtracted again from
physical wall distance. The crest-speed-up height .01 m is vertical above the
crest; reattachment near x/Hh=6.5 is a ground-extrapolated held-out observation.

## New adapter and explicit completion assumptions

`cao_tamura_budget_h1.py` uses accepted-budget-AWF case/topology utilities. It
does not import the old EP experiment generator, its P/N sampler, collar,
fallbacks, cMu=.09 or cLength=2.55.

| Input | H1 convention |
|---|---:|
| Hill height Hh; half footprint L | .04 m; .10 m |
| Cube height b; provisional effective AWF h_eff | .005 m; .000878311859293 m |
| Aerodynamic z0; assumed displacement d0 | .0002 m; .0035 m |
| Nominal incoming delta; computational roof | .25 m; .8 m |
| Selected normal AWF thickness Hw | .015 m |
| Velocity scale U(Hh); Uinf | 1 m/s; 1.456234076945 m/s |
| u_star; molecular nu | .0742679379242 m/s; 3.80952380952e-6 m²/s |
| Calibrated wall cMu; cLength | .0982278953282; 2.510862475220 |

U(Hh)=1 is an explicit dimensional realization of the frozen Hh/U(Hh)
convention, not a measured absolute speed. Re_H=10500 uses U(Hh). Linear
interpolation of the sealed CSV gives U(Hh)/Uinf=.686702787575; the original
source-audit precision differs by only 5.33e-11. The resulting roughness Reynolds
number is 3.899066741, versus reported 3.97. This source discrepancy is retained.

The immutable inlet table uses C1 Hermite U with its exact derivative, positive
log-linear k/epsilon, and a finite constant continuation after its last
zero-slope knot at z=.8. U is already constant above .4, retaining the last
measured U/Uinf=1.00010352. An interpolated Hh knot fixes the declared velocity
scale without altering any measured mean value.

For the unmeasured turbulence completion, Rxx is inferred from Fig.3c,
Ryy=.5 Rxx and Rzz=.25 Rxx are assumed, hence k=.875 Rxx. Spanwise variance,
epsilon and these ratios are not measured or uniquely determined. Below the
first measurement the inferred k is held constant for the raw inlet table.
A weak isotropic freestream k=1.5(.003 Uinf)² is reached by z=.4. Epsilon knot
values are `cMu_wall^.75*k^1.5/ell`, with
`ell=min(kappa_CLS*max(z,z0),.025 m)`. The .025 m length cap and .3% freestream
intensity are explicit provisional choices, not calibrated hill parameters.

Fresh native initialization separately uses AWF no-slip U scaled to finite
U(Hw), a nonconstant budget-reference inner k anchored at K(Hw), natural ground
k flux, and an outer-epsilon projection constrained to exact E_C(Hw). Thus the
transported inner seed differs from the raw RMS completion below Hw. Three
constrained radial projections and one setup Lift/one iteration are recorded
separately from physical governing solves. Actual native h0, fixed reference
and all histories are saved; h0/reference match bitwise in this pilot.

The raw finite epsilon completion at Hw is .105482244385 m²/s³; the accepted
AWF trace is .105045058782 m²/s³. Log interpolation does not enforce their
equality. The archived input manifest's contrary prose is corrected by an
explicit erratum, preserving its original bytes/hash. Solved epsilon is never
replaced after k. Legacy neutral-log range diagnostics in the first pilot's raw
reports are likewise labelled explicitly in the final implementation; finite
profile errors come from the strict finite reader.

The completed profile fits measured U knots by interpolation. Its selected-law
weighted RMSEs are .000907143 m²/s² for -Rxz and .00844453 m/s for sigma_u,
well above the extraction-only envelopes .0000402849 and .000870463. These are
real closure/completion discrepancies. No statistical measurement-error claim,
near-wall-isotropy fit or exact RANS-equilibrium claim is made.

## Immutable Regularized reference and native consistency

`RANSWallRegularizedReference=FiniteInflowV1` selects the new const provider;
`NeutralLogV1` remains the unchanged default. The finite initializer is separately
selected by `RANSWallReference=FiniteInflowV1`. File bytes require a verified
`RANSWallReferenceSHA256`.

For prescribed Cartesian e_in, n_in, X_in:

```text
zeta = (x - X_in) dot n_in
S_ref = U_ref'(zeta) (e_in tensor n_in + n_in tensor e_in)
gamma = Gamma * nu_t_OriginalCLS(U_ref', k_ref, epsilon_ref)
Gamma=1, gamma0=0, cubic beta=.4
```

The experimental frame is e_in=(1,0,0), n_in=(0,0,1), X_in=(-.2,0,0).
Incoming height is physical z above the upstream floor, independent of terrain
normal distance and outer eta. Original CLS viscosity is evaluated only on the
fixed completion, then frozen. A positive dummy distance is used solely for that
distance-independent viscosity calculation; no dummy Yap/source output is used.

One shared evaluation supplies `OuterClosure`, Picard's positive implicit
viscosity, cancellation-stable residual components and constitutive diagnostics.
The legacy independent R3 nonpressure action explicitly rejects finite mode.
Neutral initializers, analytical neutral weak loads and neutral error reporters
remain explicitly neutral-specific. Cache/configuration identity includes actual
profile content hash, version, units, interpolation/extrapolation, effective
origin/frame and gamma interpretation; defaults also participate. Ordinary
changed-profile, digest, origin, gamma, dt and affine-mesh restarts are rejected.

Native kernel tests recover the neutral viscosity at exact neutral input knots
to 2.22e-16 relative error. This does not claim that a coarse Hermite table is an
exact logarithm between knots. S_ref is symmetric/trace-free; gamma's stress
correction is zero when S=S_ref. Uniform weak turbulence with U'=0 has exactly
zero shear stress and production, while retaining its physical scalar sinks.
The reference is unchanged by evolving fields. Fixed-branch directional tests
respect the strain/rotation switch; coarse/fine action disagreement is 1.22e-13.

At all native outer points provider/assembly reference values agree exactly.
Actual selected stress-work production agrees to roundoff, as does epsilon's
production term. Strict native versus cancellation-stable spatial residual
components differ by at most 4.52e-16. Production includes the complete stress,
including isotropic, nonlinear and Regularized terms.

## Two-region geometry and independent checks

`RANSTwoRegionHill` implements `CosineNormalVerticalV1` for the C1 cosine hill,
with explicit feet and one-sided curvature. It keeps four queries distinct:
physical map, inner normal r/metric, nearest physical rough-wall distance, and
prescribed incoming-frame height.

Only 0<=r<=Hw uses the certified normal map. Its upper graph Gamma_H has
`x_H' >= 1-Hw*C > 0`, where C=Hh*pi²/(2L²)=19.7392088022/m. The outer map is
`(x_H(s), z_H(s)+eta*(Z-z_H(s)))`; its Jacobian is
`x_H'(s)*(Z-z_H(s))`. For the selected Hw and Z=.8, a conservative global lower
bound is .52441434 m, while the sampled minimum is .66666403 m. At each fixed r
the horizontal map is monotone with fixed side endpoints; at fixed physical x
height increases strictly with r. This establishes injectivity of the inner
strip, unique outer vertical columns, and no gap/overlap except their common
interface. It does not rely solely on sampled positive Jacobians.

Nearest-wall distance uses bounded interval subdivision of the complete
cosine/flat floor and returns a distance error bound and minimizing footpoints.
A fluid-side hill/flat tied minimum is tested; no normal averaging is used.
The nominal roof is not included as a no-slip Yap wall in H1. Facility top/roof
modelling remains an explicit H2 boundary decision.

Native tests verify flat recovery, height/footprint/foot-side metrics, invalid
normal reach, map/inverse, physical gradients, interface points/measures/opposite
area vectors, independent graph volume and divergence flux, physical observation
lookup and geometry identity rejection. For the [-.2,.2] test domain:

| Geometry check | Error/result |
|---|---:|
| Map/inverse | 1.25e-16 m |
| Independent finite-difference Jacobian | 6.31e-10 m |
| Physical coordinate-gradient action | 2.09e-9 |
| Interface point/opposite area-vector disagreement | 0 / 0 |
| Native volume / independent graph volume | .316 / .31600000000000006 m² |
| GL32–GL64 volume / divergence-flux difference | 1.11e-16 / 1.11e-16 m² |
| Normal-distance recovery / certified interval bound | 1.51e-16 / 1e-9 m |

The same kernel supports the planned [-1,1.2] m extent and correctly locates crest
and wake observations there. Geometry identities are strict; native affine
cross-mesh restart rejection is also tested. Mapped solver restart/evolution is
not claimed. Curved metrics use independent quadrature checks; flat GL4 polynomial
exactness is not assumed for terrain.

## Joint zero-step preflight and provisional selection

h_eff is provisionally fitted only to U_ref(.015) using permitted flat data,
reported u_star/Uinf and completed K(.015) in the calibrated AWF J0. It is held
fixed across the height screen. It is neither cube height b nor aerodynamic z0.
The accepted budget's counterfactual friction scale at selected K is .0885432,
versus incoming .0742679; its derived flat z0 is .0000657265 m, versus nominal
.0002 m. These differences are retained, not removed by retuning.

| Hw (m) | C Hw | K_ref(Hw) (m²/s²) | Flat cutoff (m) | E_C(Hw) (m²/s³) | Minimum budget margin* |
|---:|---:|---:|---:|---:|---:|
| .01 | .197392 | .0238329014 | .000727688 | .146535364 | 2.146589 |
| **.015** | **.296088** | **.0250145796** | **.000714123** | **.105045059** | **2.465087** |
| .02 | .394784 | .0250234083 | .000715859 | .078825507 | 2.641739 |

*Margin is productionBudget/minimumBudget minus 1, minimized across five signed
curvatures [-C,-C/2,0,C/2,C] and K factors [.5,1,2] at each height. All 45 cases
pass. Independent physical J-weighted integration, split at roughness, viscosity
intercept and cutoff, closes the budget within4.4e-16 relative error. Hw=.01 is
slightly below the first mean observation; .015/.02 lie inside the measured
height span. The fixed-profile covariance screen is positive; the complete
native audit below is stronger than this knot-only screen.

The cutoff always uses the counterfactual flat parallel zero-forcing p_A0(K),
not finite-reference or terrain production. Actual inner source remains
P_actual-E_C, with no prescribed redistribution or artificial E_C/k reaction.
The only epsilon interface datum is E_C(Hw). Failure remains an explicit
applicability result; there is no EP fallback or clipped cutoff.

Hw=.015 is the provisional H2 candidate: it has measurement coverage, a safe
normal strip and a resolved affordable compatibility mesh. This is not a
calibrated developed rough-flat solution.

For each screened height, a comparable flat four-panel tangential mesh with
normal knots `[0,.00025,.001,.003,Hw,.03,.06,.1,.2,.25,.3,.4,.5,.8]` has four
inner and nine outer panels and an estimated 2000 coupled unknowns. Only the
selected .015 case received native initialization/evolution. The roughly
.0005–.00104 m screened cutoffs and h_eff lie in the first resolved inner panels;
source quadrature still splits at their exact locations. The full hill mesh
would require additional tangential/outer resolution and a different backend.

## Executed native compatibility, representation and transaction tests

The affine physical pilot uses x=[-.2,-.1,0,.1,.2] m,
z=[0,.00025,.001,.003,.015,.03,.06,.1,.2,.25,.3,.4,.5,.8] m, cubic
velocity/scalars, linear pressure in both regions and one common gauge.
There are 2000 coupled unknowns, below the unchanged 3000 limit. Source quadrature
is unchanged (16 normal/8 tangential), retaining cutoff/roughness panels.
Molecular diffusion remains inside the strip and is omitted from the chosen
high-Re outer equations. ExactStepV1, NativeElementSideV2, unit relaxation,
physical floor/reserve and 2000/20000/2000/900 bound controls are retained.

Both streamwise boundaries and the top are prescribed by `ExternalConstraints`.
Ground U is no-slip; ground k retains its native natural condition. This is a
finite-profile compatibility pilot, not an experimental precursor. Endpoint
K=k(Hw) and z_K=z_m=Hw hold throughout.

An earlier 1678-unknown setup without outer knots .25/.4 had epsilon
`-1.51555746e-7` at `(x,z)=(-.1980144928,.30397101435)` in the native projection
of a positive completed profile. It was rejected at T=0, before any governing
solve or native history advancement. The exact seed and read-only witness are
[retained](cls-cao-tamura-budget-h1-evidence/initialization-rejection.json).
Adding those two outer representation knots and freshly initializing resolved
the undershoot; no source, profile, floor or acceptance threshold was changed.

| Completed primary run | Result |
|---|---:|
| dt / accepted steps / final time | 1e-5 s / 8 / 8e-5 s |
| Outer/U/p/epsilon/k/final-bound calls | 8 / 8 / 8 / 8 / 8 / 8 |
| Additional governing/restoration solves | 0 / 0 |
| Bound activations / signed, positive, negative, absolute supply | 0 / all 0 |
| Worst corrected regional continuity/common gauge | 1.10e-17 |
| Worst raw k / epsilon frozen linear residual | 1.36e-19 / 9.49e-20 |
| Lowest audited inner/outer covariance over steps | .01004107 / 1.90595e-5 m²/s² |
| Final refreshed epsilon matching diagnostic | 3.21298e-6 |
| Last momentum splitting defect / identity-check error | 5.94e-10 / 3.02e-19 |
| Mean recorded native attempt cost | 2.246 s |
| Whole8 wall / user / system / peak RSS | 39.30 / 34.33 / 2.49 s / 836420 KiB |

Pressure mass-inverse actions are 76 per step, counted separately. Bounding
performs no active-set/KKT/mass-factor work in this trajectory. All eight saved
packets verify solved epsilon bits unchanged from pre-k to final, and inactive
correction returns the exact raw vector. The imposed frozen trace closes to
roundoff; the refreshed trace remains a one-sweep lag diagnostic.

The full harness took 143.01 s, including separate initialization, split4+restart4,
startup checks and rejection tests. Split/uninterrupted state/history/time bits
agree; forcing hard rejection after all four solves and the bound restores the
accepted seed and all histories exactly. Failed candidate/history packets are
preserved. Altered profile bytes/digest, frame, gamma interpretation, timestep
and a still-affine same-size mesh all reject ordinary restart. Four corrupted
strict-readback inputs (checkpoint digest, state, native history and signature)
also reject. No readback calls Lift, advances history or performs a PDE solve.

| Strict finite readback | Native h0 | T=8e-5 |
|---|---:|---:|
| Permitted-data weighted U RMSE (m/s) | .005181394 | .005179620 |
| Permitted-data weighted -Rxz RMSE (m²/s²) | .000844626 | .000843161 |
| Permitted-data weighted sigma_u RMSE (m/s) | .009940074 | .009942235 |
| U physical area-L2 error vs declared composite completion | .00021460075 | .00021395052 |
| dU/dz physical area-L2 error | .53054570 | .51122365 |
| k / outer-epsilon physical area-L2 error | 2.72464e-5 / 3.37736e-5 | 2.73048e-5 / 3.37186e-5 |

All 45 permitted native data queries are valid; no weight is dropped or
renormalized away. The composite inner target is the declared AWF completion,
not an assertion that the measured inlet is exactly that wall-law curve.
Errors separate native representation from subsequent drift. Actual area-L2
drifts from native h0 are 5.00932e-6(U),2.44653e-6(k),7.62057e-7(outer epsilon).
Integration splits every native interface, fixed-profile knot and physical AWF
breakpoint; GL16/32 results agree within 6e-14 relative. Whole-element unsplit
quadrature had a visible gradient-norm error and was replaced only in this
read-only diagnostic. Production quadrature and states were untouched.

These eight steps cover only .002 Hh/U(Hh) and do not qualify temporal accuracy,
streamwise development, long-time stability or experimental accuracy. The finite
profile is not an exact equilibrium of the complete selected equations.

## Build, runtime and reproduction

Authoritative Linux Release build: GCC 11.4, `-fno-fast-math -ffp-contract=off
-D_GLIBCXX_ASSERTIONS`, four build workers (two for final focused rebuilds).
All pilots use one rank and one computational-library thread; OMP, OpenBLAS,
MKL, BLIS, VECLIB and NUMEXPR are explicitly 1. Concurrent worker/process load
was inspected and kept below the aggregate 72 limit. Shared installs were not
rebuilt. The three focused CTests pass in .55 s; kernel/preflight execution has
209 passing assertions.

An early incremental development build exposed stale objects because copied Mac
mtimes preceded Linux object mtimes. Authoritative `build-final` was configured
in an empty directory; actual object rebuilds and final source hashes were
verified. Earlier setup/build/harness failures remain indexed as development
evidence. A final Python path-type error was resolved by validating its already
completed native reports read-only, without repeating physical solves.

| Immutable artifact | SHA256 |
|---|---|
| Eight-step solver | `b2bdfc9aeffd423798197928c41f9d826097b0076409c300d5958605bf1549ac` |
| Final finite-diagnostic-label solver | `52b77c587f61cc2a1f4263cb68ece6a395d3ee0ad3fbeadce67c2ef9e5e38e67` |
| Final reference/geometry kernel | `5e4ee65ad1ea8affc2e1258c9f79e59f101cc29f5776a5894de68c5d6b675510` |
| Fixed finite profile | `9375601441e2c70fb87cf8f31de305afd4b67fd2f82a750f9274f6bb68eed697` |
| Native h0 checkpoint | `6d6e3c1757eb0eb542c68051171354a6935eaf8b5463d161a6954ed330f911f7` |
| Accepted step8 checkpoint | `6c0f3b91b5fed45fc7bb071e82ea932f9061ad270c053fa2b99aa05820d6044b` |

The final source's finite-only audit naming change was validated with an isolated
one-step old/new comparison: every governing stage, history and numerical audit
value agrees bitwise after renaming. The eight-step frozen runtime is preserved.
Complete linked-library hashes and readback runtime identities are in the compact
evidence; checkpoint signatures/history bytes were never edited.

Example reproduction on Linux, using new unused output/build directories:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
R='/media/jay/新加卷/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111'
S="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
cmake -S "$S" -B "$R/reproduce-build" -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_BUILD_SOLVER=ON \
  -DRANS_PREBUILT_SPATIAL_LIBRARY="$R/compatibility-resolved/runtime-solver/lib/libSpatialDomains.so.5.9.0" \
  -DCMAKE_CXX_FLAGS='-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$R/reproduce-build" --parallel 4 --target \
  IncNavierStokesSolver RANSCaoTamuraBudgetH1Checks RANSWallResolutionProbe \
  RANSAcceptedAWFProbe RANSHistoryRejectionChecks
ctest --test-dir "$R/reproduce-build" --parallel 1 --output-on-failure \
  -R 'RANSCaoTamuraBudgetH1|RANSAcceptedAWFKernel|RANSHistoryRejection'
python3 -B "$S/run_cao_tamura_budget_h1.py" \
  --solver "$R/reproduce-build/IncNavierStokesSolver" \
  --probe "$R/reproduce-build/RANSCaoTamuraBudgetH1Checks" \
  --output "$R/reproduce-compatibility"
python3 -B "$S/read_cao_tamura_budget_h1.py" \
  --probe "$R/reproduce-build/RANSWallResolutionProbe" \
  --campaign "$R/reproduce-compatibility" --output "$R/reproduce-readbacks"
python3 -B "$S/check_cao_tamura_h1_contracts.py" \
  --solver "$R/reproduce-build/IncNavierStokesSolver" \
  --campaign "$R/reproduce-compatibility" --output "$R/reproduce-contracts"
```

`check_cao_tamura_h1_neutral.py` records the separate one-step frozen neutral
regression. The indexed `collect_cao_h1_evidence.py` preserves the executed
collation. No historical campaign is repurposed by these commands.

## H2 handoff and scaling

The next priority is a genuine rough-flat precursor with explicit boundary roles,
followed by a bounded mapped hill test. Remaining implementation work is concrete:

1. Wire the two-region map into U/p/epsilon/k weak assembly: physical metric
   gradients, masses, advection, diffusion, source weights and all boundary fluxes.
   Retain normal AWF distance only in the strip; use nearest physical wall distance
   for outer CLS/Yap and incoming-frame height solely for the fixed reference.
2. Implement/qualify full terrain traction and tangential roughness-pressure
   coupling, complete G, both continuity blocks/common gauge and a metric-consistent
   pressure-correction mass inverse. A geometry-only positive Jacobian does not
   validate this coupling or its conservation.
3. Extend shared k positivity/constitutive inventories to physical mapped rows,
   endpoint functionals, moving cutoff rays, ground extrema and one-sided radial
   and tangential native provenance. Preserve the original raw-centered consistent
   physical mass objective and actual accepted epsilon state.
4. Replace the compatibility boundary policy with qualified inlet/top/outlet weak
   fluxes, outlet/backflow handling and pressure/gauge constraints. Decide the
   physical roof/freestream role and wall-distance convention explicitly. Develop
   the rough-flat flow and assess reserved Fig.3d before treating h_eff or the
   missing-turbulence completion as calibrated.
5. Add terrain observations and physical norm panels in the frozen global/vertical
   coordinates, independent of normal r/outer eta, with curved-metric quadrature
   verification. Keep hill profiles, crest speed-up and reattachment held out.
6. Replace the dense pilot's scaling bottleneck with a qualified blocked/sparse
   or matrix-free solve and reusable geometry/distance catalogs before a useful
   resolved precursor/hill domain. Raising the 3000 limit alone is insufficient.

For illustration,20 tangential by 16 vertical panels (8 inner/8 outer), cubic
U/scalars and linear pressures require about 10871 unknowns: one dense double
matrix alone is .945 GB, before other matrices, factors and caches. A 40 by 24
example requires about 33495 unknowns and 8.98 GB per matrix. The current 2000-unknown
pilot already peaks near .8 GiB and takes about5 s/step including startup/output.
These are size/cost estimates, not benchmarks of a scalable hill backend.

The current findings support proceeding to the boundary/metric/backend H2 work,
with Hw=.015 as a provisional choice. They also identify a substantial flat
stress/completion mismatch to assess there. They do not support fitting hill
holdouts, changing the accepted sink, or claiming completed terrain evolution.
EP-AWF, further neutral optimization, Original-CLS rehabilitation and optional
long neutral campaigns remain deferred.
