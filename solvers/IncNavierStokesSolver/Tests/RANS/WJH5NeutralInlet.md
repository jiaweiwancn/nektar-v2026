# Neutral WJ/Budget reference, zero-step use

`wj_h5_neutral.NeutralReference` is **WJNeutralBudgetReferenceV1**. It does not
construct the measured-rms `Completion` or use its Blackadar tail. Lengths in
the code use H_h=40mm; velocities use U_infinity. The common mean fit uses the
three Fig.3(a) centers at physical z<=25mm, fixed u_tau=.051 and d0=3.5mm.
Viscosity stays `.051*.005/3.97`, independently of the fitted intercept.

```bash
python3 wj_h5_cases.py --height-mm 5.5 --output /absolute/new/case-directory
IncNavierStokesSolver /absolute/new/case-directory/case.xml
```

The generator always writes NumSteps=0 and has no restart/step option. Supported
organized candidates are5,5.5,6mm. This is an idealized full-height neutral
equilibrium control, with compatible complete outlet traction, top shear and
epsilon flux. Top normal traction is an impermeability reaction. It is not the
physical finite-depth, slip-ceiling precursor. Functions use physical z and
r=z-d0, including the composite outer-map inverse. No built-in Neutral selector
is used on CompositeTunnelV1.

The inner k is the nonconstant integrated budget reference. Both stress and
scalar diffusion retain molecular viscosity. The only outer epsilon value trace
is A(K)/H. Inner/outer total tractions match in the continuum with pi_w=-K/Ne,
pi_o=0; turbulent shear and velocity gradients need not be continuous.

`test_wj_h5_neutral.py` is registered as WJH5NeutralInlet. It checks the fit,
dimensional conversions, exact constants, roughness inverse/empty intervals,
split quadrature/profile positivity, full traction, low-H mesh ownership and
projected native exterior test traces. `study_wj_h5.py` provides continuum,
generation and passive modes. `run_wj_h5.py` enforces the H5-only zero-step
allocation. `publish_wj_h5.py` preserves raw and corrected passive lineages.

Keep the Fig.12 z=7mm rms target in the observation table. Its r=3.5mm lies
inside every candidate AWF. Requiring it outside gives no feasible H when
h_s<H is retained. The Fig.3-only reference interval is not a whole-experiment
matching-height selection.

The H5 Linux initializations used1174 unknowns each. The5mm CP0 passed native
admission;5.5/6mm stopped on principal-negative projected face states at T=0.
No governing timestep, mesh tuning or extra CP0 was used. See the H5 evidence
for the continuum/representation distinction before using these configurations.
