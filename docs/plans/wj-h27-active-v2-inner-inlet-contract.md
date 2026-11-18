# Active-v2 inner inlet fixture

Reviewed head ec2c2a1496730ed6530782721e5d6ad2ac40d2b3 includes the
failed row-D test with the full historical inlet. This fixture changes
only the permitted inner inlet U and k. The historical outer U, k and
epsilon for r≥H, the top, outlet and ground data, physical geometry,
displacement and solver spaces are fixed. The closure and external loads
come from qualified table
b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d.
Use active ScalarCdV2, the amended wall stress and its full current work,
finite molecular viscosity, current v2 cutoff, auxiliary viscosity and
roughness, and transported outer epsilon with the dependent H trace.
Keep the ordinary historical-inlet default and TASK3_READY=NO.
[Approved v2 model](wj-h27-amended-candidate-v2.md)
[Original boundary contract](wj-h27-boundary-audit.md)
[Old diagnostic](wj-h27-inner-k-fixture-final.md)

At fixed Kin, solve $U'=G$, $k'=Q/D_{k,w}$ and
$Q'=E_w-(T_{\rm in}-\nu G)G$ with constant total shear
$T_{\rm in}=\nu G-R_{xr}$. Impose $U(0)=Q(0)=0$, $k(H)=K_{in}$ and
$U(H)=U_{\rm outer,old}(H)$. Treat $T_{\rm in}$ and $Q(H)$ as outputs.
The viscosity calibration remains $a=\kappa\sqrt{K_{in}/q_e}$.
Integrate directly, splitting at the fixed cutoff and model/material
breaks. Continue from the canonical v2 shear column with the joint
endpoint Jacobian in $(T_{\rm in}/\tau_0,Q(H)/A)$. Check the increasing
shear branch and preserve rejected trials. Refine the IVP and validate
with a forward integration, independent full-tensor shear, and
split Gauss budgets. Apply 1e-10 scales to endpoint, shear and
$\int(P-E)dr+Q(H)-Q(0)$. Check sampled positivity, covariance and
supported planar principal response. Compare one-sided derivatives,
total traction and k flux against the finite-nu outer inlet. The
stationary composite is a separate diagnostic.

Freeze values and derivatives as panel polynomials bound to an independent
SHA-256 before the native runs. Target 1e-8 artifact and independent
profile/gradient agreement. Keep one-sided H derivatives and common
U, k and epsilon values. Bind inlet, closure, external loads and full
boundary identities in the native diagnostic and its history. Verify
historical v2 parity, profile actions, loads, outer traces, inlet corners,
rejection of opposite histories and stale artifacts. On the existing
28×100 grid require the new inner trace's independently measured
representation RMS/MAX to be within 10% of the respective numerical
budget for U, k and their relevant gradient actions. Report interpolation
and finite-element projection errors separately.

First solve the retained 28×100 Q3/Q2 steady physical problem with
matched inlet. Use a fresh steady guess, zero elapsed physical time
and no change to the closure or interface or exterior loads. Independently
check current momentum and continuity at 1e-9, k and epsilon at 1e-8,
the trace, roots, positivity, covariance, principal response and weak
accounting. Keep separate physical owner gradients and fluxes. Compare
inner U and k to the frozen inlet and to the historical profile under
distinct labels. Use the original outer references, zero W and tau0
for ground shear. Retain all 29 frozen stations, panels, owners and
corners. Apply the inherited physical U/W .01/.02, k .05/.10, outer
epsilon .10/.20 and ground shear .05/.10 RMS/MAX envelopes. Use constant
Kin scale for every k comparison.

Choose directional or targeted refinement from the measured 2D budgets
and error locations. Compare common owners and fixed stations, with
independent split quadrature for inlet departures and exact squared
polynomial integration where available. Keep the original field
limits, including inner-k .001/.001. Run separate quadrature and
algebraic controls at .0001/.0001 for inner k and report all inherited
coupled spatial rows. Refine endpoint and owner samples. Physical
closure requires adequate spatial resolution. A locally steady inlet,
coupled steady state and transient qualification are separate gates.

Use ext4 root /home/jay/build/wj-h27-coupling-r1-20261117. At opening
the verified live account is 22,283.290719605 of 43,200 seconds, with
17,316.709280395 non-protected and 3,600 protected, and no reservation.
The old package's 2,510.914170112 seconds remain charged. Admit at most
6,000 additional seconds for this package, with 300 reserved for final
evidence. Initially allocate 250 seconds for column work and about
650 for the isolated build and binding checks. Reallocate unused time
under the same limit and freeze each run cap before launch. Use one
science process and six thread limits at one, four build workers, at
most 72 shared workers, 8/16 GiB science/build RSS stops, two 8 GiB
storage limits and at least 15 GiB free. Prior 28×100 and 48×100 steady
measurements, 466 and 627 seconds at 2.874 and 5.506 GB, are forecasts
only. Keep the Linux primary and deployed source19 untouched. Retain
each attempt in a create-only operation directory with exact input,
output and resource receipts.
