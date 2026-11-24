# C1 finite-viscosity WJ–EARSM on a normal-coordinate hill

The flat implicit baseline and curved governing checks passed on Linux. The full-height physical pilot reached .02 of the 1.48 target before the principal-response guard rejected further advancement. The inputs, source and binary identities, independent results, per-attempt diagnostics and raw Linux locations are recorded with the results. [Machine-readable evidence](wj-h27-c1-evidence.json)

## Geometry and implementation

Lengths use $L_0=40$ mm. The physical floor is $Z_b(a)=H_h\cos^2(\pi a/(2L_h))$ for $|a|<L_h$, and is flat on either side, with $H_h=1$ and $L_h=2.5$ in code units. Define $\mathbf X_b=(a,Z_b)$, its upward unit normal $\mathbf n$, and $\mathbf X_{\mathrm{eff}}=\mathbf X_b+d_0\mathbf n$ with $d_0=.0875$. The wall layer is $\mathbf X=\mathbf X_{\mathrm{eff}}+r\mathbf n$ up to $H=.15$. The flat matching level is 9.5 mm above the original floor. The fixed finite-viscosity table supplies $h=.11541061001735275$ (4.6164 mm), $K_{\rm in}=.008809015145515192$, $c(K_{\rm in})=.0243960321782528$ and $z_{0e}=.005676270321216798$.

With $\ell=\sqrt{1+Z_b'^2}$ and $\kappa_b=-Z_b''/\ell^3$, the effective surface has $ds/da=\ell(1+d_0\kappa_b)$ and $\kappa_g=\kappa_b/(1+d_0\kappa_b)$. The inner determinant in $(a,r)$ is $\ell(1+(d_0+r)\kappa_b)=(ds/da)J$. Each foot is a mesh face. The surface and tangent agree there, and volume and face calculations use the curvature from the owning patch.

The global bound $M=H_h\pi^2/(2L_h^2)$ on $|Z_b''|$ gives $\mathbf n(a)\cdot(\mathbf X_b(b)-\mathbf X_b(a))\le M|\mathbf X_b(b)-\mathbf X_b(a)|^2/2$. A point on a physical-floor normal at distance $d_0+r$ has squared distance to any other foot larger by at least $(1-(d_0+r)M)|\mathbf X_b(b)-\mathbf X_b(a)|^2$. Thus the complete displaced tube has unique normal feet.

Above matching, the separate ruled mesh joins $\mathbf X_H(a)$ to $(a,Z_T)$, with $Z_T=20$ for the physical case. Its horizontal derivative is positive. A lower bound for its determinant in the outer chart is $(1-(d_0+H)M)(Z_T-H_h-(d_0+H)\sqrt{1+\max|Z_b'|^2})/(Z_T-d_0-H)$. At a fixed physical horizontal coordinate, the vertical coordinate increases along this map, giving global injectivity with vertical side boundaries outside the two feet. At $H$, both owners integrate $|\partial_a\mathbf X_H|\,da=J_H\,ds$. The mapped analysis uses the exact piecewise analytic chart on every native quadrature point. The remaining geometry integration error is measured with quadrature refinement.

An experimental sample at fixed laboratory $x$ and vertical height $y'$ above the physical hill is placed at $(x/L_0,Z_b(x/L_0)+y'/L_0)$. The inverse of the appropriate patch then determines its chart coordinates. The outer coordinate labels the ruled mesh. The converter and geometry audit use the laboratory convention in [Cao and Tamura](https://doi.org/10.1016/j.jweia.2005.10.001).

| Model operation | Implementation |
| --- | --- |
| Chart, certification, and invariant mass and continuity | `Geometry::At`, `Geometry::ValidateNormalMesh` and `Flow::BuildInvariant` in `WJNativeFlow.cpp`, using `RANSWallNormalMap` |
| Cartesian gradients, physical face cofactors and cutoff panels | `WJNativeRules.cpp` |
| Shared matching ownership and dependent outer epsilon | `Coupling::BuildTopology`, `DirichletState`, and `DependentTraceDirection` in `WJNativeCoupling.cpp` |
| Sampled pressure, k and convection | `Flow::NormalSampledAction`, `MixedFlow::SetTransportedKSpace` and `Coupling::MomentumAdvection` use the tangential matching derivative. The velocity derivative includes the derivative of the frame. |
| WJ stress, work and admissibility | `AdaptAmended`, `CheckPlanarState`, and `EvaluatePlanarVelocityTangent` receive the physical wall normal. Gradients and momentum remain Cartesian. |
| Current nonlinear equations | `AdvanceNewton` in `WJNativeTimeNewton.cpp` assembles the new-time BE or BDF2 residual, with current conservative scalar transport and both sampled momentum components. |
| Initialization and recovery | `InitialVelocity`, `ProjectContinuity`, and `WriteHistory` in `WJNativeCoupling.cpp`. Full state and advection histories bind to the reference, chart, topology and numerical signature. |
| Verification and sampling | `WJNativeHillMMS.cpp` derives continuous sources in a smooth virtual tube. `c1_geometry.py`, `c1_mms.py`, and `c1_oracle.py` reconstruct physical locations, derivatives, fluxes and budgets independently from chart-basis exports. |

The native XML and FLD files use affine chart elements for storage and connectivity. The result packet includes the physical chart parameters and signed element coefficient maps. A point sampled in physical space must be mapped to a chart coordinate before evaluating the stored polynomial.

The matching spaces share $\mathbf U$, $k$ and $\Pi$. The velocity and k residuals impose traction and total k-flux balance weakly. Outer epsilon has the dependent Dirichlet trace $A(K)/H$ and a one-sided reaction. The full local momentum advection and stress divergence remain active alongside the sampled load

$$
\frac{\mathbf t\,\mathbf 1_{0<r<h}}{J_H}\,\partial_s\left(\Pi_H-\frac{2K}{3}+u_{t,H}^2\right).
$$

The mixed pressure column is $G_h-C_\Pi$ in the direct solve and its pressure Schur action. Current advection and both scalar transports see the refreshed fields and lifts. The cutoff uses the qualified reference family evaluated at the whole matching trace.

## Reference and reproduction

The approved exact-inner input has SHA-256 `3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1`. The governing TeX note reviewed in the main checkout has SHA-256 `61e5b616daf368583884df9798fdc2cadabfbe105151fc0ee466c1a843f545a9` and the Cao–Tamura PDF has `adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df`. The live GitHub tip at the start was `57845c53a29ba32cca49c31a89e0f4ec4da4f300`. The historical Newton pilot used explicit advection through step 16. C1 uses `CoupledJFNKV1`, `SparseLU` and `WJNewtonAdvection=Implicit` from the first step, with full refreshed residual tolerance at most $10^{-10}$. [Reference table](wj-h27-linux/task2r-r1/tables-exact-inner-v1.json)

The authoritative Linux machine was `jay-Super-Server` via `jay@100.105.247.110`. The prescribed NTFS mount was absent and an ordinary mount required local authentication. The controls and bounded guard diagnostics used the isolated ext4 build and run root `/home/jay/build/wj-c1-20261124.7zDGkMC3`. The Release executable uses GCC 11.4, CMake 3.22 and the existing MPI-enabled Nektar++ 5.9 libraries with Open MPI 4.1.2. The build uses six workers. Each solver and audit uses one rank or process and one thread per numerical library. Set `SRC` to a clean checkout of the tested revision and `RUN` to a new absolute experiment directory.

The final numerical sources were compiled at `3274486e345e0b90df5225918aa40806bbf5330e`. Linux validated the fixture and independent audit code at `780db051f93dfdaa6bee1d464860bc33abd79c35`. The standard-main ELF has SHA-256 `061b577d3222d53fac57c35eee5e5272b6653e707cb5f747df9d59de9d951594` and its embedded numerical source identity is `f8b43960b7b5a9c7c8161898f4bcc751bbddb18b42dee5fc6cd2480693ec3e9c`. The three-step flat and curved coefficient vectors, input histories and residuals agree exactly with the earlier independently audited replays. The final revision's same-signature cold and uninterrupted checkpoint hashes are `97fd497c3f3b7ca7175a46c52586461bcad24633c90d35eaf2cd823de6d8de17` (flat) and `10fe512fcf1272b9458f179372cce5fdbb751529fd23d49eb890f4954e81eade` (curved manufactured).

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$RUN/tmp"
export TMPDIR="$RUN/tmp"
tests="$SRC/solvers/IncNavierStokesSolver/Tests/RANS/h27"
table="$SRC/docs/plans/wj-h27-linux/task2r-r1/tables-exact-inner-v1.json"
cmake -S "$tests" -B "$RUN/build" -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$RUN/build" --parallel 6 --target IncNavierStokesSolver WJH27HillGeometryProbe
solver="$RUN/build/IncNavierStokesSolver"
run_case() {
  (cd "$1" && /usr/bin/time -v mpirun -np 1 "$solver" --use-opt-file run.xml.opt run.xml > run.log 2>&1)
}

python3 "$tests/c1_fixtures.py" flat "$RUN/flat" --steps 3 --science 1 --step-seconds 180
run_case "$RUN/flat"
python3 "$tests/c1_fixtures.py" flat "$RUN/flat-cold" --steps 1 --step-seconds 180 \
  --restart "$RUN/flat/state.recovery-2.bin"
run_case "$RUN/flat-cold"
cmp "$RUN/flat/state.final.bin" "$RUN/flat-cold/state.final.bin"
python3 "$tests/physical_transient_oracle.py" --table "$table" --seconds 180 \
  --native "$RUN/flat/state.science-1.json" --out "$RUN/flat-be-audit.json"
python3 "$tests/physical_transient_oracle.py" --table "$table" --seconds 180 \
  --native "$RUN/flat/state.result.json" --out "$RUN/flat-bdf2-audit.json"

python3 "$tests/c1_fixtures.py" manufactured "$RUN/curve" --steps 3 --science 1 \
  --step-seconds 180 --trace-derivative-check
run_case "$RUN/curve"
python3 "$tests/c1_fixtures.py" manufactured "$RUN/curve-cold" --steps 1 \
  --step-seconds 180 --trace-derivative-check --restart "$RUN/curve/state.recovery-2.bin"
run_case "$RUN/curve-cold"
cmp "$RUN/curve/state.final.bin" "$RUN/curve-cold/state.final.bin"
for spec in "be state.science-1.json" "bdf2 state.result.json"; do
  read -r name frame <<< "$spec"
  python3 "$tests/c1_oracle.py" --table "$table" --native "$RUN/curve/$frame" \
    --order 20 --out "$RUN/curve-$name-audit.json"
done
for height in 0 1; do
  d="$RUN/probe-$height"
  python3 "$tests/c1_fixtures.py" hill "$d" --height "$height" --steps 0
  "$RUN/build/WJH27HillGeometryProbe" "$d/run.xml" "$table" "$d/geometry-probe.json"
done
python3 "$tests/c1_assess.py" --table "$table" \
  --geometry "$RUN"/probe-{0,1}/geometry-probe.json --out "$RUN/geometry-audit.json"

for level in 0 1 2 3; do
  d="$RUN/mms-spatial-$level"
  python3 "$tests/c1_fixtures.py" manufactured "$d" --level "$level" --steady-seconds 180
  run_case "$d"
  python3 "$tests/c1_oracle.py" --table "$table" --native "$d/state.result.json" \
    --order 20 --seconds 300 --out "$d/audit.json"
done
for spec in ".08 4 2" ".04 8 4" ".02 16 8" ".01 32 16"; do
  read -r dt steps stride <<< "$spec"
  d="$RUN/mms-time-$dt"
  python3 "$tests/c1_fixtures.py" manufactured "$d" --level 2 --dt "$dt" \
    --steps "$steps" --science "$stride" --recovery "$stride"
  run_case "$d"
done
python3 "$tests/c1_assess.py" --table "$table" --spatial "$RUN"/mms-spatial-*/audit.json \
  --temporal "$RUN"/mms-time-{.08,.04,.02,.01} --out "$RUN/convergence.json"
```

For the steady quadrature witness, generate level 3 again with `--quadrature 20 24`. For the iterative witnesses, generate level 3 with `--tolerance 1e-13` and the .01 transient with `--tolerance 1e-12`. Compare with `c1_assess.field_norms` on the unchanged spaces. Use `c1_oracle.py --order 24` to repeat each independent audit with the finer integration rule. The solver outputs and independent audit orders for the measured cases are recorded individually.

For the developed flat control, run `c1_fixtures.py flat "$RUN/flat-steady" --steady-seconds 180 --step-seconds 180`, then generate the three-step flat input above under a new path. Set `WJInitialFields` in that input's `SOLVERINFO` to the absolute `flat-steady/state.result.json` path. Compare the three coefficient vectors of the typed steady and transient result packets directly.

The prospective smooth $L^2$ rates were 4 for Q3 velocity and both scalars and 3 for Q2 pressure, with 3 for the Q3 velocity gradient. Curved forcing uses a divergence-free streamfunction, a varying matching $K$, a nonzero matching $k$ diffusion flux and all three sampled forces. The manufactured pressure error removes its additive gauge. The time comparison uses .16 and .32 as common positive times. Comparing fields on an identical mesh also isolates time increments from the spatial error floor.

## Measured governing checks

The independent physical flat audits with 20- and 24-point rules measure the following residuals. The evolving state completed three steps at $\Delta t=.02$. The restart begins from the step-2 checkpoint and includes its earlier BDF2 state. The uninterrupted and cold step-3 checkpoints match byte for byte within each tested revision. Changing the target advection scheme to `Explicit` rejects the implicit checkpoint with a scientific identity error.

The retained outer inlet has physical finite-viscosity momentum and epsilon defects $\nu U_{\log}''$ and $2\nu A_{\rm in}/r^3$. Those terms are present in the evolving equations. The audit measures accepted-step algebraic error and integration differences after reconstructing the current residual. The manufactured refinements independently measure spatial representation and time error. Outlet development and long quiet preservation require a developed state of the same discrete equations.

A separate small flat steady solve with the same retained inlet reaches refreshed residuals `[5.77e-15, 8.04e-17, 7.99e-13, 1.23e-14]`. The typed steady state has exactly the same flow, k and epsilon coefficients after three implicit steps to .06. This quiet control uses the developed discrete solution. Its short duration leaves the long-flat gate open.

| Flat control | Momentum | Continuity | k | Outer epsilon |
| --- | ---: | ---: | ---: | ---: |
| BE at .02 | 5.97e-14 | 2.69e-16 | 5.54e-12 | 1.02e-13 |
| BDF2 at .06 | 1.63e-13 | 2.82e-16 | 6.29e-12 | 1.69e-12 |
| Limit | 1e-9 | 1e-9 | 1e-8 | 1e-8 |

Independent reconstruction of the native curved basis gives maximum gradient and divergence defects of $3.94\times10^{-13}$ and $1.61\times10^{-13}$ on the full-height geometry probe. Position, inverse, determinant, normal and matching measure defects are at most $3.56\times10^{-15}$. A round trip through laboratory coordinates at the sides, roof, floor, matching surface and hill feet has maximum chart error $4.90\times10^{-13}$ against $2\times10^{-11}$. The full-height lower bounds for the physical-floor tube stretch, effective $J$, and outer determinant are approximately .812, .760 and .770. The rotation control covers the lower tensor, collar and outer WJ states, including nonzero discrete divergence. Its principal minimum before the orientation-dependent arithmetic reserve changes by at most $5.00\times10^{-16}$. Rotation of stress, work, traction and covariance agrees within $3.34\times10^{-16}$. The maximum relative difference for the represented cutoff, its slope and the dependent dissipation derivatives at 39 held-out trace values is $6.44\times10^{-16}$ against $10^{-8}$.

| Curved spatial mesh cells | Velocity RMS | Pressure RMS | k/Kin RMS | Epsilon relative RMS | Raw traction jump RMS | Raw k flux jump RMS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 14 | 1.561e-6 | 2.461e-7 | 2.367e-5 | 3.308e-5 | 1.807e-7 | 3.126e-10 |
| 56 | 9.924e-8 | 1.473e-8 | 1.166e-6 | 1.838e-6 | 1.614e-8 | 1.476e-11 |
| 224 | 6.221e-9 | 8.893e-10 | 4.501e-8 | 7.943e-8 | 1.469e-9 | 7.697e-13 |
| 896 | 3.908e-10 | 5.69e-11 | 1.529e-9 | 2.806e-9 | 1.363e-10 | 4.315e-14 |

The finest pair gives measured spatial rates 3.99, 3.97, 4.88 and 4.82. The pressure and scalar rates are specific to this smooth fixture. Changing the 896-cell native panel counts from 12×16 to 20×24 changes the four field norms by $2.13\times10^{-18}$, $1.29\times10^{-18}$, $3.05\times10^{-16}$ and $1.36\times10^{-16}$ respectively. These are norms of the field difference.

| Common time | Same-space rate for velocity | Pressure | k | Outer epsilon |
| ---: | ---: | ---: | ---: | ---: |
| .16 | 2.029 | 2.061 | 1.990 | 1.991 |
| .32 | 2.041 | 2.024 | 2.015 | 2.017 |

Those rates compare the .04 to .02 field differences with the .02 to .01 differences. At .32 the finest differences are $1.22\times10^{-8}$, $2.32\times10^{-8}$, $1.19\times10^{-6}$ and $2.06\times10^{-6}$ respectively. The absolute velocity error approaches the measured spatial floor at the .01 step.

The independent curved BE residuals are `[3.16e-15, 3.22e-17, 1.75e-13, 1.36e-13]` for momentum, continuity, k and epsilon. The corresponding cold-tested BDF2 frame has `[8.71e-16, 2.62e-17, 4.32e-14, 3.09e-14]`. The independently integrated interior epsilon lift defect is at most $2.33\times10^{-16}$, against $10^{-10}$. Its pointwise polynomial approximation is reported separately and is $1.24\times10^{-12}$ in the curved BDF2 frame. Native refreshed lift defects are zero at the reported precision. The relative directional difference for the represented dependent lift is $3.04\times10^{-11}$. The cutoff crosses the radial element face at .02435, moving from chart abscissa .4671463702 at BE to .4316143358 at BDF2.

The maximum shared-trace discrepancy over the manufactured refinement is $3.47\times10^{-18}$. The greatest owner measure discrepancy is $4.45\times10^{-16}$. The destruction-integral identity agrees within $6.29\times10^{-18}$ and the layer mass coordinate identity within $4.78\times10^{-19}$. The raw layer normal-mass defect drops from $1.23\times10^{-6}$ to $4.63\times10^{-10}$. The maximum inner weak k budget defect is $1.03\times10^{-16}$ and the maximum power accounting discrepancy is $1.14\times10^{-21}$. The power accounting uses the discrete $U^{n+1}\cdot D_t U$, the weak divergence term with $|U|^2/2+\Pi-2k/3$, physical boundary work, separate sampled work, molecular dissipation, production and the manufactured load. The curved cold restart and uninterrupted run give byte-identical checkpoints at .06 within each replay revision.

The direct flat-map operator comparison gives relative defects of $1.81\times10^{-16}$ for velocity mass, $1.79\times10^{-16}$ for pressure mass, $9.81\times10^{-17}$ for divergence, $2.07\times10^{-14}$ for sampled pressure, and $3.16\times10^{-14}$ for sampled k. The matched one-step applications have maximum advection-coefficient differences below $8.11\times10^{-14}$. The two nonlinear solutions have a maximum velocity/pressure coefficient difference of $1.39\times10^{-10}$ with the prescribed solver tolerance and separate initialization projection.

The refined raw column k budget RMS drops from $3.44\times10^{-9}$ at 14 cells to $1.45\times10^{-12}$ at 896 cells. The absolute weak-divergence work drops from $1.75\times10^{-13}$ to $4.69\times10^{-18}$ over the same levels. At time .32 the absolute BDF kinetic-work correction decreases from $4.07\times10^{-10}$ at .08 to $7.78\times10^{-12}$ at .01. The equation identity includes both corrections directly.

At 896 cells the raw matching epsilon flux is the outer Dirichlet owner's $3.31\times10^{-7}$. The largest independent free matching row for momentum or k is $1.39\times10^{-19}$. At .32 with $\Delta t=.01$ the largest such row is $5.21\times10^{-17}$, the weak layer k defect is $1.24\times10^{-16}$ and the power accounting defect is $1.63\times10^{-21}$. Repeating the independent time audit at 20 and 24 quadrature points changes the four scaled residuals by at most $4.74\times10^{-16}$ and the four absolute field error estimates by at most $6.78\times10^{-17}$.

Tightening the spatial nonlinear target from $10^{-10}$ to $10^{-13}$ returns the same 896-cell field. Both solves reach a maximum refreshed residual of $3.31\times10^{-14}$. At .16, tightening the .01 timestep target from $10^{-10}$ to $10^{-12}$ changes the fields by $3.88\times10^{-16}$, $5.78\times10^{-16}$, $6.87\times10^{-14}$ and $9.00\times10^{-14}$ for velocity, pressure, k and outer epsilon. The corresponding loose and tight maximum residuals are $2.21\times10^{-12}$ and $1.54\times10^{-14}$. At .32 the changes are $5.15\times10^{-16}$, $5.31\times10^{-16}$, $1.18\times10^{-13}$ and $1.58\times10^{-13}$. The loose and tight maximum residuals there are $2.36\times10^{-12}$ and $1.54\times10^{-14}$.

## Physical hill protocol

The planned physical endpoint is $t=1.48$ in code units, 74 steps of .02, using one hill-height travel time at the retained approach velocity of approximately .678 at original elevation 40 mm. The top is 800 mm and uses zero normal velocity, the retained tangential traction, zero k flux and the retained outer epsilon flux. The `ZeroTransverseV1` outlet is at $x/H_h=15$, beyond the experiment's last wake-profile station at 7.5. The native solver rejects advective backflow on the complete outlet polynomial. The physical inputs set a constant matching height and retain the full finite-viscosity inlet and fixed reference table. The manufactured loads are confined to the declared verification mode.

The physical seed is a projected streamfunction. The constrained Sobolev projection preserves the essential velocity and reduces the represented continuity defect. The optional slowdown and outer epsilon perturbation apply at initialization. The latter vanishes at the inlet and matching lift. The advancing equations use the same stress, diffusion, wall dissipation and source terms at every step.

### Shallow physical control

At hill height .002 (0.08 mm), three physical steps completed to .06 with BE and BDF2. The continuity residual of the represented seed dropped from $8.26\times10^{-8}$ to $3.02\times10^{-13}$. A cold start from the full .04 history gives the same .06 checkpoint as the uninterrupted run, SHA-256 `173f684a9a9d5b8f5863daff1f53ad83da9165acc541962767960985602108bd` on `9fb0382656165926aee383f89992cd922b31e541`.

| Independently audited shallow state at .04 | Measured | Limit |
| --- | ---: | ---: |
| Momentum | 1.18e-12 | 1e-9 |
| Continuity | 1.21e-15 | 1e-9 |
| k | 8.86e-12 | 1e-8 |
| Outer epsilon | 9.44e-12 | 1e-8 |
| Dependent epsilon trace projection | 1.53e-16 | 1e-10 |

At .06 the accepted minimum $k$ is .008808647 and minimum outer epsilon is $1.744274\times10^{-5}$. The covariance and principal margins are .27929 and .0032913, against $10^{-10}$ and $10^{-8}$. Whole-face $K/K_{\rm in}$ is $[.9999598640,1.0000258094]$, giving margins of .0999598640 and .0999741906 to the certified edges. The cutoff range is $[.0243958784,.0243962714]$.

At .04 the independent signed matching mass flux toward the outer region was $1.21\times10^{-8}$ on both sides. The raw k fluxes toward the outer region were $-5.19\times10^{-9}$ and $-1.96\times10^{-8}$ from the inner and outer polynomials. The maximum weak matching row was $4.94\times10^{-15}$. The three wall-force powers from pressure, k and convection were $-7.40\times10^{-6}$, $1.04\times10^{-9}$ and $-2.25\times10^{-7}$. The mean-power accounting discrepancy was $9.17\times10^{-16}$, and the independent outlet minimum normal velocity was positive at $2.35\times10^{-6}$. The evidence records the matching momentum fluxes, one-sided outer epsilon flux, and full layer balance.

### Full-height failure diagnosis

| Targeted change | Observed result |
| --- | --- |
| Represent the normal streamfunction on the coarse mesh | The first full-height state fails the principal guard at $a=-2.2621$, $r=.13953$. The margin is $-3.85\times10^{-4}$ with positive covariance. |
| Slow the initial streamfunction and project weak continuity in the mass norm | The projection creates a large very-near-wall shear correction. The guard finds a negative covariance eigenvalue. |
| Use an $H^1$ continuity projection | The near-wall covariance failure clears. The first sampled failure shifts to the outer transition. |
| Recover horizontal far-approach streamlines and bisect the streamwise mesh | With the stronger supported initial epsilon perturbation, the first failed principal margin improves to $-6.57\times10^{-5}$ at $a=-9.5765$, $r=.79058$. |
| Give the epsilon perturbation its own entrance ramp | The refined state clears the early entrance failure and reaches the upslope outer transition at $a=-2.1981$, $r=1.3501$. The independent margin is $-1.76905\times10^{-2}$. The coarse representation fails between-node epsilon positivity. |
| Spread velocity recovery higher in the outer region | The refined .65 slowdown reaches $a=-1.5731$, $r=3.1549$ before the first failed guard. Independently differentiating the continuous preprojection streamfunction confirms a negative principal response at that point. |
| Shorten the outer extension of the matching tangent | Two refined starts fail earlier, at $r=1.3501$. The narrower transition was reverted. |
| Add a smooth initial outer epsilon perturbation at the broad velocity transition | The first meshes fail the between-node positivity certificate. |
| Resolve the radial perturbation transitions and bisect the abscissa | The refined full-height initial state passes all guards and BE accepts the first physical step at .02. The second step triggers its final principal guard. The resolved radial mesh with coarse abscissae fails an upstream initial principal check. |

The broad .65 seed's first failing point is identified by `hill-broad-refined-65/run.log` at laboratory $(x,z)/L_0=(-1.673733,3.471990)$. Its native principal margin is $-0.0276623$ against the required $10^{-8}$, with covariance over $k$ of .232952. Independent directional differences of the same state agree with the native sign and magnitude. Differentiating the continuous velocity seed at the same physical point gives $-0.02490125$ at both tested spacings. The corresponding coarse continuous point has positive principal response despite a negative projected response there, so spatial representation can also cause initial guard failure. The refined startup problem persists in the continuous seed. The computed negative margin applies to this state. The existence of another admissible full-height trajectory remains open.

Reproduce that refined full-height guard case with the executable for `ce57a9c250a3d8cc78b63c23187dcd64b6e765d6` and the build recipe above. Exit status 1 is the observed guard result. All native inputs, the full guard record, the raw script for independent continuous-shape evaluation and the Linux executable have content hashes in the evidence.

```bash
python3 "$tests/c1_fixtures.py" hill "$RUN/hill-guard" --steps 0 --step-seconds 180 \
  --bisect-streamwise --initial-slowdown .65 --initial-epsilon-bump 6
status=0
run_case "$RUN/hill-guard" || status=$?
test "$status" -eq 1
python3 "$tests/c1_assess.py" --table "$table" --guard "$RUN/hill-guard/run.log" \
  --out "$RUN/hill-guard-audit.json"
```

At the broad-seed preflight, the whole matching trace remains at $K/K_{\rm in}=1$ with margins .1 and .1 to the qualified interval boundaries. The corresponding cutoff is .0243960321782528. The projected broad .65 seed passes weak continuity at $8.86\times10^{-13}$, and the constitutive failure occurs before the first BE iterate.

### Resolved full-height mesh and first accepted step

The resolved 40 mm mesh has $30\times117$ Q3/Q2 quadrilaterals. Its initialization reaches weak continuity $2.54\times10^{-13}$ and zero native dependent-lift defect. It has minimum $k=.0088090147$, minimum transported epsilon $1.74427\times10^{-5}$ and covariance and principal margins .275642 and .00372008. Its whole-face table range is $[1,1]$ in $K/K_{\rm in}$ and the margins to both certified edges are .1. The initial steady residual is $5.50$ in outer epsilon, consistent with an evolving initial profile. The first BE step is accepted at .02 with maximum refreshed residual $9.87\times10^{-12}$.

| 40 mm state at .02 | Momentum | Continuity | k | Outer epsilon |
| --- | ---: | ---: | ---: | ---: |
| Fresh native residual | 3.20e-13 | 1.30e-15 | 2.64e-12 | 9.87e-12 |
| Independent 6×8 volume rule | 1.46e-11 | 1.29e-13 | 1.91e-11 | 3.14e-10 |
| Independent 20×20 volume rule | 1.14e-5 | 3.37e-15 | 5.39e-6 | 1.62e-2 |
| Scaled audited limits | 1e-9 | 1e-9 | 1e-8 | 1e-8 |

The 6×8 comparison independently reproduces the discrete equations. The 20×20 result exposes underintegration of this physically evolving field, especially in outer epsilon. Its dependent trace projection still agrees within $3.34\times10^{-16}$. The pointwise polynomial approximation error is $1.07\times10^{-10}$. At .02 the whole-face $K/K_{\rm in}$ is $[.9946692853,1]$ with edge margins .0946692853 and .1. The cutoff ranges from .0243960322 to .0244279405. The accepted minima are $k=.00863812$ and epsilon $1.74427\times10^{-5}$, with covariance and principal margins .276146 and .00287745.

Independently integrated matching mass flow is $-1.07415\times10^{-5}$ toward the outer region for both owners. The raw k fluxes are $5.59834\times10^{-5}$ and $5.59997\times10^{-5}$. The inner weak k balance defect is $1.62\times10^{-15}$ on the 20-point audit. The difference between the independently evaluated power balance and the residual work is $3.15\times10^{-16}$. At this integration order the residual work is $-3.71\times10^{-5}$. The three sampled forces do work of $-1.88364\times10^{-3}$, $-3.97\times10^{-8}$ and $5.94384\times10^{-5}$ in pressure, k and convection order. The independent outlet minimum of $2.35\times10^{-6}$ indicates outward flow in the audited frame.

At the .04 trial the nonlinear equation solver reached $3.08\times10^{-12}$ before the final principal guard rejected the candidate at chart $(a,r)=(-7.5844131,.11535282)$, in the upstream flat wall layer near fixed $h$. The native principal margin is $-0.03418757024$ and the independent directional finite-difference values round to the same number. Covariance remains positive at .261095. The preserved accepted checkpoint exactly matches the earlier .02 recovery file, with SHA-256 `48f292398bd5f80488e23f2968df3b3ad83a2d9a68a3a17a3596c5359ac2887e`.

The first parallel attempt at the timestep and cold comparisons exhausted cooperative deadlines with overlapping MPI core bindings. Repeating with `mpirun --bind-to none -np 1` recovered the scientific comparison. The .02 cold continuation reproduces the failed point and all printed gradient and constitutive values exactly. Each smaller timestep passes through .02, then rejects a converged candidate at the same coordinate. This rules out the .02 timestep, the original cold process and the nonlinear stop as sole explanations of the observed loss.

| Full-height physical timestep | Last accepted time | First rejected candidate time | Principal margin at the rejected point |
| ---: | ---: | ---: | ---: |
| .02 | .02 | .04 | -3.41876e-2 |
| .01 | .02 | .03 | -1.75590e-2 |
| .005 | .02 | .025 | -7.54935e-3 |

The .005 accepted BDF2 state passes the matched-order independent audit with maximum scaled residual $3.17\times10^{-10}$. Lowering the seed slowdown from .65 to .50 still loses positive principal response at the same upstream location at .04. Its margin is $-1.80505\times10^{-2}$. Raising the native quadrature to 12×16 reduces the higher-order mismatch. The independent 20×20 residuals at .02 become `[1.67e-7, 3.46e-15, 1.15e-7, 8.59e-4]`. The same-signature high-quadrature .04 continuation still fails the principal guard at chart $(a,r)=(-7.523049,.1153952)$, with native margin $-0.0359244$.

Reducing the initial slowdown further to .35 moves the first .04 guard failure to the rising slope at $(a,r)=(-1.573052,.11535282)$. Its principal margin is $-8.56\times10^{-6}$ and covariance over k is .265785. This separates sensitivity to the initial velocity profile from the fixed model and illustrates why the tested loss is local to the attempted evolution.

The .005 trajectory also restarts from its genuine BDF2 history at .015 and accepts .02. The original and cold runs have the same numerical signature and maximum refreshed residuals $8.37\times10^{-12}$ and $8.38\times10^{-12}$. The maximum coefficient differences in flow, k and epsilon are $6.94\times10^{-13}$, $1.56\times10^{-16}$ and $4.48\times10^{-17}$. The rejected .04 cold case preserves its .02 predecessor byte for byte.

Across the 11 distinct accepted full-height positive-time states, the whole-face $K/K_{\rm in}$ stays in $[.9946692852,1]$ with minimum lower and upper margins .0946692852 and .1. The observed cutoff range is $[.0243960322,.0244279405]$. The smallest k, epsilon, normalized covariance and normalized principal response in those states are .00863533, $1.74427\times10^{-5}$, .256923 and .00280212.

The separate steady preparation route exhausted its bounded solves. The unshifted run reached a residual of 5.47585 in 377.75 s of solve time. One repeat used an unrecognized option name and ran unshifted. The corrected run used `WJSteadyPreconditionerMass=50`, confirmed by the executable, and reached 5.49103 in 360.04 s. This mass belongs to the inverse preconditioner. Both main comparisons preserved the original state and their unconverged candidate evidence. The steady-plus-perturbation route requires a converged startup state.

Generate this physical case and its preflight with the final source and build recipe above.

```bash
profile=(--step-seconds 180 --bisect-streamwise --resolve-outer-seed \
  --initial-slowdown .65 --initial-epsilon-bump 6)
python3 "$tests/c1_fixtures.py" hill "$RUN/hill-preflight" --steps 0 "${profile[@]}"
run_case "$RUN/hill-preflight"
python3 "$tests/c1_fixtures.py" hill "$RUN/hill-dt02" --dt .02 --steps 2 --science 1 "${profile[@]}"
status=0
run_case "$RUN/hill-dt02" || status=$?
test "$status" -eq 1
cmp "$RUN/hill-dt02/state.recovery-1.bin" "$RUN/hill-dt02/state.failed-accepted.bin"
python3 "$tests/c1_oracle.py" --table "$table" --native "$RUN/hill-dt02/state.science-1.json" \
  --order 20 --seconds 500 --out "$RUN/hill-be-audit.json"
python3 "$tests/c1_oracle.py" --table "$table" --native "$RUN/hill-dt02/state.science-1.json" \
  --order 6 --radial-order 8 --seconds 300 --out "$RUN/hill-be-matched-volume-audit.json"
```

For the .01 and .005 repeats use `--dt .01 --steps 4` and `--dt .005 --steps 8` with the same `profile` options. Run each under `/usr/bin/time -v mpirun --bind-to none -np 1 "$solver" --use-opt-file run.xml.opt run.xml` in its own case directory. For the accepted restart, generate another .005 case with `--steps 1 --restart "$RUN/hill-dt005/state.recovery-3.bin"` and compare its .02 packet to `hill-dt005/state.science-4.json`. The .02 higher-quadrature pair uses `--quadrature 12 16 --step-seconds 360`, one fresh step and one additional cold step from the fresh `state.recovery-1.bin`.

For the steady comparison, generate the same resolved seed using `--steady-seconds 360`. For the shifted run, set `WJSteadyPreconditionerMass = 50` in the session parameter block and confirm `preconditioner_mass_only` in the steady evidence packet.

## Measured reference cost

The flat three-step replay took 13.22 s with peak resident memory 278,140 KiB. Its measured construction, initialization and complete advancement took 1.40, .65 and 10.20 s. Recovery took .081 s in total, scientific frames .386 s, and final export .040 s. The three-step curved manufactured replay took 10.03 s and 86,016 KiB. Construction, initialization and complete advancement took .61, .67 and 7.62 s. Recovery took .066 s, scientific frames .075 s and final export .014 s. The 896-cell steady manufactured check took 123.80 s and 1,580,636 KiB.

| Disjoint final-step work in seconds | Flat replay | Curved manufactured replay |
| --- | ---: | ---: |
| Constitutive evaluation | .5332 | .1107 |
| Wall evaluation | .1257 | 1.3286 |
| Quadrature and assembly | .9292 | .2797 |
| Factor rebuilds | .3909 | .00433 |
| Mixed momentum and pressure application | .0390 | .00024 |
| Scalar application | .00726 | .00013 |
| Guards and refreshed residual | 1.0006 | .2639 |
| Manufactured forcing | 0 | .2419 |

The mixed direct solve handles momentum and pressure together. The separate pressure-Schur timing counters are zero in these `SparseLU` runs. Factor rebuilds include scalar factors, which are also available as a subcategory in the evidence. Timing samples include the observed Linux host load and measure a single rank.

The three-step shallow physical run took 58.48 s and 1,030,624 KiB. The one-step cold continuation, including load, signature verification and history reconstruction, took 36.48 s and 880,628 KiB. These small physical timings inform the pilot budget. The production cost requires a measured full-height trajectory.

The shallow final BDF2 step spent 1.081 s in constitutive evaluation, .347 s in wall work, 2.816 s in assembly, 1.776 s in rebuilds, .148 s in mixed solves, .029 s in scalar applications and 1.982 s in guards and residuals. Its complete logged step times, including scoped output and recovery, were 20.03, 8.53 and 8.50 s. This gives an order-of-magnitude scale for the 74-step pilot once an actual-height first step supplies its own warm rate. The refined full-height broad-seed preflight used 78.47 s and 1,899,964 KiB and produced no advancing cost sample.

The 40 mm preflight took 49.14 s with 2,545,764 KiB peak RSS. The first accepted step plus the rejected .04 candidate, their startup, and recovery took 159.36 s with 3,129,380 KiB peak RSS. The accepted BE solve spent 7.714 s in constitutive evaluation, 1.752 s in wall evaluation, 17.962 s in assembly, 26.147 s rebuilding factors, 2.409 s in mixed applications, .478 s in scalar applications and 8.349 s checking guards and residuals. The rejected BDF2 solve spent 57.76 s. The 12×16 first step including setup and output took 208.24 s with 3,813,300 KiB peak RSS. These costs bound the short diagnostics. A warm-step rate for a qualified actual-height pilot remains open.

## Gate decisions and next package

| Flag | Status |
| --- | --- |
| `CLEAN_IMPLICIT_BASELINE` | PASS |
| `CURRENT_CURVED_GOVERNING_VERIFIED` | PASS |
| `REAL_HILL_PILOT_PASSED` | BLOCKED |
| `C2_READY` | NO |
| `LONG_FLAT_QUALIFIED` | OPEN |
| `PRODUCTION_AFFORDABLE` | OPEN |
| `EXPERIMENT_VALIDATION_COMPLETE` | OPEN |

The full-height accepted trajectories end at .02. The measured principal-response loss survives both timestep reductions, the finer native quadrature, and different starting velocity tapers. The lower-amplitude taper shifts its physical location. The higher-order audits identify a separate spatial integration error. Together these results establish a blocker for the tested starting states and meshes. The model's admissible branch over other full-height states remains to be investigated.

Complete C1 with a geometry-continuation initialization. Begin from an admissible developed shallow state. Remap the volume flux as the chart changes, impose the fixed inlet and dependent trace, and project continuity at each height. Validate constitutive admissibility and outer quadrature at each stage. Start a fresh BE/BDF2 history on the final mesh and advance to one hill-height travel time with a bound and checkpoints sized from the measured accepted-step cost. This proposed change concerns preparation of the state. The finite-viscosity equations and reference family govern each stage.

1. In C2, compare spatial and temporal witnesses at common times for the physical hill and reevaluate the top and outlet truncation effects. Assess the developed whole-face $K$ range and validate any needed fixed-$h$ extension against independently solved reference columns.
2. Qualify the developed, fresh-start, and signed perturbation flat classes through at least code time 26.13350407707765. Keep the established production horizon at 505.02466521278404 and the terminal window beginning at 447.8483832973273.
3. Reduce factorization, assembly and recovery costs using measured workload attribution. Verify $1.25C_{\mathrm{total}}\le86400$ seconds per production trajectory with setup, rejected work, recovery, output and assessment included.
4. Use the qualified hill trajectory in H1 for comparisons at the original laboratory stations, with the vertical datum mapping above.
