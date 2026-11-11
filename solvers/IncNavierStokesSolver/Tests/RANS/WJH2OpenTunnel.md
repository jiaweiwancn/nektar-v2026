# WJ Budget physical flat-tunnel sessions

Generate a fresh, portable production session with the measured-profile
completion and its separately recorded assumptions:

```bash
python3 -B wj_h2_cases.py --output /path/to/fresh-case \
  --kind physical --steps 4 --dt .02 --reuse ExactStepV1
IncNavierStokesSolver /path/to/fresh-case/case.xml
```

Use the source-built library environment recorded in `wj-h2-linux-evidence.md`.
The generator refuses an existing output directory. It writes the mesh,
Functions and `case-design.json`, including actual normal knots, conditional
roughness fit and epsilon completion. A refinement adds `--refined`; the H2
dt=.02 refined case is a retained hard-failure witness, not a qualified baseline.
The prepared longer continuation was denied before execution and requires
separate approval. These examples do not authorize bypassing that stop.

The selected production entrypoint and boundary roles are:

```xml
<I PROPERTY="SolverType" VALUE="RANSEquilibriumWallEvolution" />
<I PROPERTY="TurbulenceModel" VALUE="WJHighRe" />
<I PROPERTY="RANSWJGradientPolicy" VALUE="DevStrainRawWorkV1" />
<I PROPERTY="RANSWallFormulation" VALUE="WJCalibratedBudgetAWF" />
<I PROPERTY="RANSWallNativeGeometry" VALUE="BudgetTerrainV1" />
<I PROPERTY="RANSWallMapContinuation" VALUE="CompositeTunnelV1" />
<I PROPERTY="RANSWallBoundaryMode" VALUE="OpenTunnelV1" />
<I PROPERTY="RANSWallOpenTop" VALUE="SlipV1" />
<I PROPERTY="RANSWallTemporalCoupling" VALUE="PimpleSingleSweep" />
```

`BoundaryValuesInner` supplies u/v/k at the inlet; `BoundaryValuesOuter` supplies
u/v/k/epsilon. Functions use x and **physical vertical y** in the 2D expression
interface; reported stress axes are physical x, spanwise y, vertical z. The
outer computational rho is not the physical height. Ground no-slip takes
precedence at the inlet corner. The shared H velocity/k values must agree and
epsilon(H) must equal A(K)/H; that epsilon degree belongs to the dependent trace.

`OutletTractionInner` and `OutletTractionOuter` each require u/v components of
the complete physical mean traction. The physical pilot uses zero components.
That natural condition fixes the pressure level through stress; it is not a
Dirichlet condition pi=0. No additional p function or global mean-pressure row
is allowed. Outlet k/epsilon diffusive flux is zero. Sampled material backflow
stops explicitly; no backflow data/extrapolation or positive-velocity clamp is
implemented in this restricted route.

`SlipV1` fixes top normal velocity, zero tangent traction and zero scalar
diffusive flux. Do not add conflicting top Functions. The ceiling is at physical
z=20. The ground has zero modeled k flux. There is no inner epsilon equation
or additional interface epsilon flux.

For the separate analytic neutral control:

```bash
python3 -B wj_h2_cases.py --output /path/to/fresh-control \
  --kind analytic --steps 2 --dt .0001 --reuse Off
IncNavierStokesSolver /path/to/fresh-control/case.xml
```

This selects `ReferenceTractionV1`, supplies full reference outlet traction,
`TopTractionOuter/u`, and outward `TopDiffusiveFluxOuter/k,epsilon`. Normal top
traction remains a reaction. The required transverse outlet traction and top
shear are +u_ref²; the epsilon diffusive flux is nonzero. Its K=q_e*u_ref²,
transported inner budget reference, physical-r outer log/tail and pressure
constants are distinct from the measured-k physical completion. No experimental
z0 is imposed independently on that analytic control.

`test_wj_h2_native.py` runs the compiled boundary, pressure-column, weak-load,
reuse, restart/rollback, signature, corner and backflow controls. Passive
`open_checks` requires the exact original producing report/checkpoint and binds
their state/history hashes. A zero-step restart audit cannot replace the
preceding completed callback. The flat open route, one serial rank, p=3 and the
3000-state evolution limit are deliberate restrictions; curved physical open
boundaries and a developed hill are not qualified by these sessions.
