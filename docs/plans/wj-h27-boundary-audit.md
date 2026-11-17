# H27 retained inlet, finite-viscosity boundaries and prospective physical gates

This independent audit preserves its preparation-time scope. Subsequent Linux
outcomes and current readiness are in
[the H27 feasibility report](wj-h27-feasibility-closure.md).

## Material Passport

- Artifact: independent boundary/inlet audit for H27-A; source review and algebraic derivation, not solver validation.
- Audited historical revision: `ec0b2f1ba0de65ffabfd41001cc0e0a423ddd9db` in the branch-owning Mac worktree. No historical input, source PDF, solver, runner, model contract or Git history was changed by this audit.
- Authority: the supplied WJ-H27 package, including retained application inlet, finite outer molecular viscosity, common full-WJ collar, physical mixed boundaries and prospective physical gates. The coordinator explicitly selected the H27 reference-traction/reference-flux extensions below before pilots.
- Local source coverage: H5/H10 inlet contracts and calibration; H13 boundary, geometry, clocks and observation protocol; H25/H25R1 contracts/closure; historical boundary/forcing/inlet sources; both WJ TeX notes; relevant local paper pages and associated Markdown extracts. No Linux command or Mac solver/numerical validation was executed.
- PDF integrity: `cao2006` SHA256 `adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df`, 19/19 pages; `richards2011` SHA256 `2a8d769dde2d85a6034f2cc64a5da70da31feb9dfdf2c236484cf38520541780`, 10/10 pages. Structural preflight PASS with bundled pypdf. Earlier default-Python UNAVAILABLE sidecars are preserved separately. Sidecars/extracts/rendered checks are in `/private/tmp/wj-h27-boundary/`; structural PASS does not claim whole-paper semantic review or human reading.
- Direct image checks: Cao-Tamura PDF/printed pp.5-7, Richards-Norris PDF p.3/printed p.259. Relevant text also checked on Cao pp.4,13,15 and Richards PDF p.5/printed p.261. Extraction has mathematical-font errors; equation signs and the probe-calibration percentage were checked on rendered pages.
- TeX identities unchanged: main `1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124`; appendix `0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c`.

## 1. Units, geometry and the exact retained application inlet

Use length unit `L0=H_h=0.04 m`, velocity unit `U_ref=U_infinity`, time unit `L0/U_ref`, k/stress/modified-pressure unit `U_ref^2`, dissipation unit `U_ref^3/L0`, and viscosity unit `U_ref L0`. `U_ref=1` in code is a nondimensional reference, not an asserted experimental speed of 1 m/s. The sources do not establish a dimensional U_infinity value to convert code time to seconds.

Use physical elevation `z`, displaced wall distance `r=z-d0`, and a separate mesh/chart normal coordinate `rho`. In the flat composite chart, `z=d0+rho` inside and `z=d0+H+(rho-H)*(Z-d0-H)/(Z-H)` outside. Do not confuse the outer chart coordinate with physical wall distance. The numerical ground is physical `z=d0`; the original tunnel floor is physical `z=0`.

| Fixed input | Nondimensional value | Physical meaning |
|---|---:|---|
| x domain | [-10,15] | [-400,600] mm, length 1000 mm |
| top elevation Z | 20 | 800 mm above original floor |
| d0 | .0875 | 3.5 mm, paper's assumed 0.7 times 5 mm cube height |
| modeled wall distance H | .15 | 6 mm; matching plane at physical 9.5 mm |
| top wall distance rT | 19.9125 | 796.5 mm |
| fitted aerodynamic intercept z0e | .005676270321216798 | .2270508128486719 mm |
| fixed sand-grain roughness hs | .10854503542535461 | 4.341801417014184 mm; neither cube height nor z0e |
| u_tau | .051 | .051 U_ref |
| nu | 6.423173803526447e-5 | .051*.005/3.97, calibrated using reported .2 mm, not fitted z0e |
| N_e, q_e | 81/20, 81/(2 sqrt(143)) | equilibrium WJ constants |
| C_d | 572/6561 | independent scalar coefficient, .08718183203779911 |
| kappa | .3819178635467331 | sqrt(1.3*(1.82-1.44)*sqrt(C_d)) |
| K_in | .008809015145515192 | q_e*u_tau^2 |
| A_in | .0003473286082198882 | u_tau^3/kappa |
| tau0 | .002601 | u_tau^2, kinematic stress |

These exact stored calibration bytes are in `docs/plans/wj-h5-linux/calibration.json:514-542`; `wj_h13_cases.py:26-33` requires their binary64 identity. Domain knots and physical geometry are in `wj_h13_domain.py:15-16` and `wj-h13-fullheight-paired-domain-evidence.md:65-74`. Unit provenance is `wj_h1_experiment.json:10-11`; the chart formula is in `wj-h10-flat-inlet-initialization-contract.md:16-19`.

The application inlet at x=-10 is unchanged for its entire height:

\[
 U_{in}(r)=\begin{cases}\tau_0J_{0,old}(r),&0\le r<H,\\
 (u_\tau/\kappa)\log(r/z_{0e}),&H\le r\le r_T,\end{cases}
 \quad W_{in}=0,
\]
\[
 k_{in}(r)=\begin{cases}k_{w,old}(r),&r<H,\\K_{in},&r\ge H,\end{cases}
 \qquad \epsilon_{in}(r)=A_{in}/r\quad(r\ge H).
\]

Here `old` denotes the stored high-Re-derived application calibration, including its original finite-nu inner profile. It is not the new H27 reference problem. Precisely,

\[
 a=\kappa u_\tau,
 \quad z_v=-.0026018310348560384,
 \quad D_{old}(r)=\nu+a\max(r-z_v,0),
 \quad J_{0,old}(r)=\int_0^rD_{old}^{-1}\,ds,
\]
\[
 P_{old}=\tau_0^2\nu_A/D_{old}^2,
 \quad E_{old}=A_{in}/\max(r,r_{c,old}),
 \quad r_{c,old}=.026419898309769336,
\]
\[
 Q_{old}(r)=\int_0^r(E_{old}-P_{old})\,ds,
 \qquad k_{w,old}(r)=K_{in}-\int_r^H Q_{old}(s)/D_{old}(s)\,ds.
\]

Do not recompute this imposed inlet using the new collar/reference cutoff or the new inner diffusivity. `wj_h5_neutral.py:119-138`, `wj_h2_cases.py:69-94` and the appendix `:476-510` identify these functions. The fitted z0e is the original fixed-slope, mean-only three-point fit at physical heights <=25 mm; no rms-derived k, finite-depth taper or free-stream cap is permitted (`wj-h5-neutral-inlet-evidence.md:34-50`, `wj-h13-fullheight-paired-domain-evidence.md:33-37`). The old full-height logarithm exceeds U_ref near the 800 mm ceiling; that intentional idealization remains.

## 2. Source boundary law versus selected H27 finite-nu extension

All stresses below are per unit density. For the CURRENT H27 fields define

\[
 M=\nabla U+\nabla U^T,\quad
 \Pi=p/\rho+2k/3,\quad
 V=\nu M-R+2kI/3,\quad
 T=V-\Pi I=\nu M-R-(p/\rho)I,
\]
\[
 j_k=-D_k\nabla k,\qquad j_\epsilon=-D_\epsilon\nabla\epsilon.
\]

Use the revised wall tensor below H and full WJ above H; the wall collar must limit to the same WJ tensor. Outer diffusivities are `Dk=nu+Cd*k^2/(sigma_k*epsilon)` and `De=nu+Cd*k^2/(sigma_e*epsilon)`. The coordinator's selected inner law is `Dkw=nu+Cd*k^2/(sigma_k*Ew)`, a positive H27 choice distinct from the old auxiliary-viscosity diffusion. This changes the new reference/transport equations, not the imposed historical k inlet.

The historical implementation uses complete prescribed reference-traction Functions, not zero velocity gradients. `wj_h13_cases.py:35-46` gives top tangent traction tau0, zero normal velocity, fixed outward epsilon diffusion `B/rT` with `B=Cd*K_in^2/sigma_e`, zero top k flux, complete outlet tractions `(K_in/Ne,0,tau0)` inner and `(-K_in/Ne,0,tau0)` outer, and zero outlet scalar diffusion. `RANSWallEvolutionBoundary.cpp:57-74,99-115` applies those prescribed values off-reference; it does not recompute a current-state outlet normal-stress condition. Its pressure seed is not an inlet pressure condition.

H27 explicitly retains the physical reference-load principle and re-derives its data with finite viscosity and the full normal-stress tensor. The reference data and changes below must be published at A and pinned before developing the actual 2D solution. This is not a claim that the old numerical top flux or outlet constants are unchanged.

### Ground and inlet

At r=0: `U=0`, `j_k.n=0`; no imposed k=0, no imposed ground shear, and no inner epsilon PDE. The computed ground traction is `T(-e_z)`; its tangential force on the fluid has the opposite sign from upward shear. For flat no slip and the below-collar tensor, the tangential shear reduces to `(nu+nu_A)*(partial_r U_t+grad_t W)`. Its value is an output, not `u_tau^2` prescribed again. The scalar ground law uses current Dkw, not momentum mu. Ground pressure/normal stress is a reaction. No slip does not require k=0 in this modeled rough layer (appendix `:505-510`).

At x=-10 impose the retained U,W,k and outer epsilon values above. Do not impose Pi or p, nor duplicate the dependent epsilon-H row. Record inlet scalar diffusion and traction as reactions when forming budgets. Conservative scalar boundary flux is `U_n*q+j_q.n`; zero prescribed diffusion at another boundary does not eliminate its advective scalar flux.

### Top: selected mixed velocity and fixed finite-nu reference flux

At r=rT, n=e_z and Pn=I-nn^T:

\[
 U\cdot n=0,\qquad P_nTn=\tau_0 e_x,
 \quad j_k\cdot n=0,
\]
\[
 j_\epsilon\cdot n=J_{\epsilon,T}^{ref}
 =\left[\nu+\frac{C_dK_{in}^2}{\sigma_\epsilon(A_{in}/r_T)}\right]
 \frac{A_{in}}{r_T^2}
 =\frac{B}{r_T}+\frac{\nu A_{in}}{r_T^2}.
\]

Thus `nu*(partial_r U+partial_x W)-R_xz=tau0`. There is no top tangential Dirichlet velocity, no extra normal traction and no extra pressure condition. `n.T.n` is the impermeability reaction. In the solve, the flux is the fixed number above: `-De(current)*partial_n epsilon=Jref`; it is not a fixed gradient or a coefficient-updated numerical load. Relative to H25, the exact added outward scalar load is `nu*A_in/rT^2`; its ratio to the old load is `nu*sigma_e/(kappa*u_tau*rT)`. Record both terms explicitly in the manifest and the epsilon budget.

The supporting literature is Richards-Norris (2011), PDF p.3/printed p.259, Eq.(23): driving upper shear, zero k diffusion, and nonzero outward epsilon diffusion are needed for the homogeneous surface-layer reference. PDF p.5/printed p.261, Eq.(48) retains the independent scalar-diffusion coefficient in the stress-transport analogue; Sec.3 tests inlet/outlet profiles in an empty domain. These are high-Re results; the extra molecular term is H27's derivation, not a quotation from the paper. Do not import that paper's k-e constants or optional volumetric top buffer sources.

### Outlet: selected PRESCRIBED new-reference complete traction

At x=15, n=e_x, prescribe `Tn=t_ref(r)`, `j_k.n=0`, and outer `j_epsilon.n=0`. Preserve explicit rejection of sustained outlet backflow; do not silently clamp the convective flux or invent a backflow scalar state. The historical diagnostic limit is `U_n<-1e-12` (`RANSWallEvolutionBoundary.cpp:155-169`); any H27 scaled tolerance must be explicitly declared.

Construct a distinct, pinned finite-nu/collar counterfactual calibration column with zero forcing, reference total shear tau0 and the retained physical normal-stress level:

\[
 V_{xz}^{ref}=\tau_0,\qquad
 \Pi^{ref}=V_{zz}^{ref}-K_{in}/N_e,
 \qquad T_{zz}^{ref}=K_{in}/N_e,
\]
\[
 t_{out}^{ref}(r)
 =\left(V_{xx}^{ref}(r)-V_{zz}^{ref}(r)+K_{in}/N_e,\ 0,\ \tau_0\right).
\]

Evaluate V_ref with the full finite-nu strain, revised wall/collar stress, reference local k and actual reference Ew; outside use the declared counterfactual outer k/epsilon/gradient branch. Solve/check the reference shear relation `nu*G_ref-R_xz,ref=tau0`; do not insert the old high-Re gradient and silently assert its total shear is tau0. The reference column, root branch and its one-sided interface limit require the independent checks assigned to the model/reference implementation. It is boundary calibration only, not a replacement for the same-inlet 2D developed solution and not an unproved stationary outer scalar solution.

The derivation recovers the old loads: old isotropic inner `Vxx=Vzz=0` gives `+K_in/Ne`, while old planar outer `Vxx=-K_in/Ne,Vzz=+K_in/Ne` gives `-K_in/Ne`. Their jump was a consequence of the obsolete inner tensor, not an independent physical pressure jump to preserve. A common reference tensor in the H27 collar produces a common limiting outlet traction. Publish the new profile, its units, interpolation/assembly error and hash; then keep it fixed throughout startup, physical reference and perturbation cases.

For current fields the imposed normal equation is

\[
 \Pi=2\nu\partial_xU-R_{xx}+2k/3-t^{ref}_{out,x},
\]

with vertical traction `nu*(partial_x W+partial_r U)-R_zx=tau0`. These are total-traction conditions, not simultaneous prescribed Pi and stress gradients. The scalar diffusion equations use their own current coefficients. Advected momentum/scalar and traction loads must all appear in the boundary budgets.

Rejected alternatives: (a) retaining the old piecewise constants with the new common collar creates an avoidable outlet/H traction-data discontinuity; (b) using `Vxx_current-Vzz_current+K_in/Ne` as the prescribed t_x changes the law to `Pi-Vzz=-K_in/Ne`, a current-state normal-stress-pressure condition. It is not the historical fixed complete-traction law and is not selected. No boundary load may be fitted to the eventual developed profile.

### Pressure nullspace and interface

A constant shift of Pi leaves interior gradients and the AWF pressure-gradient sampler unchanged, but changes full outlet traction by `-c n`. Since outlet normal velocity test traces are free and total outlet traction is prescribed, this physical problem fixes the pressure level; no additional gauge is imposed. By contrast, if every boundary normal velocity were essential and only tangent tractions remained free, a constant pressure could remain a gauge. Therefore inspect the restricted assembled pressure action on the constant field, including boundary and AWF terms, with the actual per-component constraints. The existence of a first-component Neumann edge alone is not a nullspace criterion.

Across H use common U,k and common oriented mass, total-traction and advective-diffusive k fluxes. Use the sole outer epsilon trace `epsilon_o(H)=Ew(H)`. With the retained tail amplitude and cutoff <=H,

\[
 E_w(H;K)=\frac{(K/q_e)^{3/2}}{\kappa H},\quad
 E_w(H;K_{in})=A_{in}/H.
\]

Do not prescribe a second epsilon flux, old `Pi_w-Pi_o=-K/Ne` row, or artificial C1 U continuity. The old jump is derived from old normal stresses in appendix `:738-765`; H27 replaces those stresses. The conservative general interface laws are appendix `:1018-1056`. Discrete one-sided C0 gradients/tractions need convergence and owner-resolved reporting, not pointwise forced equality.

## 3. Corner compatibility and complete sampled roughness forcing

The inlet/H value conditions are exactly `U_inner(H)=U_outer(H)`, `W=0`, `k_inner(H)=K_in`, and `epsilon_in(H)=Ew(H;K_in)`. Reject a new reference/amplitude branch that breaks the last equality. The inlet's known gradient kink is `U_inner'(H)/U_outer'(H)=1-z0e/H`; it cannot be removed by changing the prescribed inlet. It is a regularity issue for the new common-tensor interior, not permission to overwrite either trace. Verify weak convergence with explicit inlet/H neighborhood support; do not claim a globally smooth manufactured-order rate for the physical corner.

The inlet/ground corner has U=0 on both sides and `k_in'(0)=0`, consistent with zero k flux. The inlet/top corner has W=0 and k derivative zero. The selected finite-nu epsilon flux exactly matches the retained epsilon derivative there because `-De_ref*(-A_in/rT^2)=Jref`. The retained parallel inlet velocity does NOT exactly satisfy the selected total top shear: its shear is `tau0+nu*u_tau/(kappa*rT)`. A smooth parallel initial field therefore fails natural traction compatibility. This is expected startup/development data, not an exact stationary solution. Keep the essential inlet, assemble top traction weakly, and quantify inlet/top corner convergence; do not set the top tangent velocity or erase molecular stress. The derivative mismatch is finite and explicitly known.

The newly derived outlet-reference traction removes the old tensor-induced load jump, but does not establish existence or regularity of the complete same-inlet 2D problem. No corner gate is closed solely by value checking; the physical solve must also pass its current budgets, convergence and admissibility tests.

For flat 2D, the historical sampled convective term is **conservative tangential acceleration**, not merely `U_H*partial_x U_H`:

\[
 F_{A,x}=\left[\partial_x(U_t^2)+\partial_x\Pi
       -\tfrac23\partial_x k-b_x\right]_{H^+}
 =\left[2U\partial_xU+\partial_x\Pi
       -\tfrac23\partial_x k-b_x\right]_{H^+},
\]
\[
 f=b+\mathbf1_{r\le h_s}F_A.
\]

This is appendix `:927-950` and native `RANSWallEvolutionPatch.cpp:771-790`: `ut*div(ut)+grad(ut)*ut`, plus mechanical-pressure gradient and minus physical body acceleration, then tangent projection and constant ray extension. If another advective convention is chosen, its difference must be declared and justified. No sampled local time derivative is included. Use the actual Q3 k gradient for `-2/3 grad k`, with outer element ownership at H. Gravity/hydrostatic contributions cancel when represented consistently. Charge `integral U.f_A` in mechanical energy and do not count the roughness force twice in pressure/traction work.

## 4. Fixed horizon and observation contract

The original domain clocks remain relevant because inlet, domain, units, K_in and A_in are unchanged. Source: `wj_h13_cases.py:48-63` and `wj-h13-fullheight-paired-domain-evidence.md:267-284`.

| Clock in L0/U_ref | Historical value |
|---|---:|
| retained inlet bulk U | .9566264028836502 |
| T_bulk=25/U_bulk | 26.13350407707765 |
| slowest outer transit=25/U_in(H) | 57.17628191545675 |
| top transit | 22.935076496133906 |
| old conservative inner diffusion clock | 195.80590488759879 |
| interface K/epsilon=K_in*H/A_in | 3.8043289281565653 |
| top K/epsilon=K_in*rT/A_in | 505.02466521278404 |

Historical horizon is `max(3*T_bulk,2*T_slowest_outer,T_inner_diffusion,T_top_k/epsilon)`, hence `T=505.02466521278404`; terminal window `[T-T_slowest_outer,T]=[447.8483832973273,505.02466521278404]`. Retain those exact values. Recompute additional H27 diffusion/developed-state clocks for interpretation; a slower new clock cannot justify shortening the original requirement. Physical seconds are these numbers times `0.04/U_infinity[m/s]`, not the wall-clock computational runtime. Steady-solver iterations contribute no physical spinup time.

For the principal entry pilot require at least `max(T_bulk,3*K_in*H/A_in)=26.13350407707765` completed code-time units. Declare identical physical endpoint times for dt comparisons and retain startup and perturbed matched refinement windows. This entry pilot is much shorter than the 505-unit Task3 horizon; neither proves stationarity alone.

Preserve all H13 stations: the union of

`[-10,-7.5,-5,-2.5,-1.25,0,1.25,2.5,4,5.5,6.5,7.5,10,12.5,15]`,

`[-2.5,0,2.5,1.25,3.75,5,6.25,7.5]`, and

`[3.75+0.25*i for i=0,...,17]`.

At every station retain physical height and element-side ownership, both H traces and all face/cutoff/roughness sides. Record whole-profile and inner/outer norms. Retain the frozen Fig.3 physical heights and the physical 7 mm near-floor observation (`r=3.5 mm<H=6 mm`), rather than reclassifying it as outer data. Fig.3's flat reference is at the virtual crest x=0, not x=-10. The source is `wj_h13_cases.py:65-86`; Cao pp.6-7 establish the reference location and distinct figure axes; p.13 specifies the rough-case 7 mm measurements.

Observables: U,W,k,outer epsilon, Pi, all six independent Rij, physical velocity gradient, full total shear and ground traction; current interior K, Ew, cutoff and branch; actual P=-R:grad U versus counterfactual production; epsilon-H identity and epsilon reaction flux; common mass/k/interface traction flux; signed scalar/mechanical storage, advection, diffusion, source, roughness work and boundary terms; weak continuity and strong divergence separately; positivity/root/realizability/principal margins and rejected attempts. Pressure keeps the outlet-determined level.

## 5. Prospective physical development envelope, distinct from numerics

H13 explicitly left `application_error_budgets=None` (`wj_h13_cases.py:73-80`); H25 states that its comparison budgets are numerical, not wind-tunnel uncertainty (`wj-h25-production-strategy-contract.md:92-103`). No established same-inlet physical-development acceptance tolerance was found in these contracts. Freeze the following **engineering criterion**, selected before H27 pilot/reference results:

| Physical adjustment to retained inlet / forcing scale | RMS limit | sampled MAX limit |
|---|---:|---:|
| U and W, each divided by U_ref | .01 | .02 |
| k difference divided by K_in | .05 | .10 |
| outer epsilon difference divided pointwise by A_in/r | .10 | .20 |
| upward ground tangential shear minus tau0, divided by tau0 | .05 | .10 |

Apply U/k/epsilon limits at every station, separately in inner and outer regions, and in whole-domain and regional physical norms. For W the inlet comparator is zero. Apply wall-shear norms along the ground, with all selected wall stations included in MAX. An RMS is `sqrt(integral |delta q/S|^2 / integral 1)` using physical measures and split quadrature; a station uses its physical line measure. MAX means the prospectively pinned sampling inventory, not a theorem about the continuum supremum. Use absolute positive fixed scales, never divide by zero wall velocity. Add corner neighborhoods to the numerical audit but do not silently drop unfavorable samples from the physical envelope.

This allows percent-scale physical adjustment while requiring much smaller independently assessed discretization errors. It is motivated by a bounded engineering definition of acceptable approach-flow change; it is not deduced statistically from the paper and is not evidence that this new model will pass. The selected tolerances preserve a meaningful comparison near the matching height, where finite-nu effects and the revised collar matter most. Failure must be reported as unacceptable physical development (or an unresolved numerical error if refinement does not pass), without retrospective tolerance tuning.

R, Pi, gradients, inner dependent Ew, cutoff and flux/source changes must be reported with scales `u_tau^2`, the declared reference gradient, and their own positive reference quantities. Pi is not an imposed inlet datum; do not fabricate a pressure-inlet preservation requirement using obsolete old pressure seeds. Stress evaluated with the revised law on retained profiles can be a diagnostic comparator, explicitly distinct from the old stress field. The same-inlet developed H27 solution is the numerical preservation reference; compare it separately with the imposed profiles above.

The H27 numerical budgets remain U RMS/MAX 1e-4 U_ref, k 1e-3 K_in, R/Pi 1e-3 u_tau^2, gradients 1e-2 of declared scale, outer epsilon 1e-2 of local positive reference epsilon, and traction 1e-3 u_tau^2. They cannot be replaced with this looser physical envelope. Finer h/p/q verification and time-refinement comparisons use common physical points and both element-side owners.

Cao p.5 reports sensor-voltage pitch sensitivity below +/-1.5% over 1-8 m/s, not a uniform uncertainty bound on mean velocity, Reynolds stress, dissipation or CFD predictions. H5's figure-reading bounds are deterministic digitization effects: z0e joint bound [.2171964,.2370378] mm, not confidence limits; d0/delta/u_tau uncertainties are unquantified (`wj-h5-neutral-inlet-evidence.md:43-50`). Do not call this engineering envelope experimental validation or claim a dissipation uncertainty absent from the experiment.

## 6. Finite perturbation amplitude, support and energy contract

The historical seed imposes velocity .001 U_ref and 1% k/local-epsilon increments before dependent-trace lifting (`seed_wj_h25.py:137-151`); H25R1 did not qualify that finite transient (`wj-h25r1-closure.md:95-109`). Its old peak-point normalization alone gives no mesh-independent energy lower bound. Freeze the following initial-data rules before H27 evolving pilots:

1. Use a compact smooth streamfunction perturbation, `delta U=partial_r psi`, `delta W=-partial_x psi`, or its native weak-divergence-compatible counterpart. It vanishes at exterior essential boundaries; construct fresh coupled histories after admission. Do not perturb pressure just to satisfy a desired norm.
2. Initial candidate support is `x in [-2.5,2.5]`, `r in [.075,.45]` (physical x +/-100 mm, elevation 6.5-21.5 mm). It crosses H=.15 and includes interior K. Require at least four streamwise element intervals and four resolved normal intervals spanning the interface; enlarge/support-adjust prospectively if the mesh cannot resolve that shape. The essential inlet, top and ground remain unchanged.
3. The nonzero support must have physical code area >=1.5 (2400 mm^2 per cross-section) and a half-peak superlevel set of area >=.05 for velocity and k; outer epsilon requires >=.05 in the outer region. Freeze the support geometry and source shape before comparing meshes/dt. Native representation must independently meet these measures; a few nonzero nodal coefficients are insufficient.
4. Require `MAX |delta velocity|/U_ref = .001` to representation tolerance, `RMS_domain |delta velocity|/U_ref >=1e-5`, and `RMS_support >=.1*MAX`. Thus the incremental kinetic energy norm `0.5*integral |delta velocity|^2/(|Omega| U_ref^2)` is at least `5e-11`. Also record component norms and gradient change. A streamfunction shape may be rescaled once as initial-data construction, before histories, to meet the fixed peak; not after observing pilot stability.
5. Require `MAX |delta k|/K_in=.01`, `RMS_domain >=1e-4`, `RMS_support >=.001`, and a nonzero interior interface K change with peak `.01 K_in`. The independent outer epsilon bump, supported away from H, has `MAX |delta epsilon_ind|/(A_in/r)=.01`, outer RMS >=1e-4 and support RMS >=.001. These are final represented-field checks, not nominal coefficients.
6. Recompute the wall root/cutoff/Ew from the perturbed K and fields, then lift the sole dependent epsilon-H trace. With the retained amplitude relation, `Ew_new(H)/Ew_base(H)=(K_new/K_base)^(3/2)`: a +1% K change relative to the same base scale produces about +1.5037% trace change. Do not force it back to 1%. Report induced trace peak and require <=2% relative to local unperturbed Ew(H); report the total epsilon norm separately from its independent 1% bump. If the developed local K differs from K_in, use the exact ratio rather than the 1.5037% illustrative number.
7. Require nonzero changes in velocity gradients, k, epsilon, interior K and cutoff above ten times independently assessed representation noise. Admissibility is checked on the final represented and trace-consistent state. Adjust support/shape, preserving these fixed peak/RMS/support bounds, if the first shape fails; do not silently lower the amplitude. Add the opposite sign or a second independently specified shape as a targeted robustness test.

The support and RMS lower bounds are prospective engineering choices, not paper-derived disturbance spectra. If two signed seeds cannot satisfy both these amplitude/support bounds and the revised constitutive/scalar gates, report a finite-perturbation feasibility failure with its location/margin. Do not replace it with a quiet developed-state cost measurement.

## 7. Required integration evidence before closure

- Preserve exact inlet source/calibration identities; separately hash new boundary calibration/reference columns and fixed top/outlet loads.
- Verify corner U/k/epsilon values with scaled tolerances and dependent trace on every initial/restart path. Keep current scalar history, owner and cutoff identities; never import an H25 restart as an H27 state.
- Implement component-dependent essential masks: inlet U,W; ground U,W; top W only; outlet free. Verify actual restricted pressure constant action and the absence/presence of a gauge from those spaces.
- Assemble full symmetric stress and total boundary tractions, current scalar diffusion and conservative advective flux. Include known Q3-k mechanical-pressure forcing with sign -2/3 and outer one-sided H ownership.
- Pin reference-load interpolation error and ensure the new outlet reference tensor is common at H. Reject stale calibration/current-field mixing. Adding molecular viscosity in volume while leaving stress/flux boundary operators high-Re is not the selected model.
- Verify the full independently integrated scalar and mechanical-energy budgets, including inlet reactions, top epsilon loss, outlet advection/work, ground reaction and sampled roughness work.
- Distinguish three objects in all outputs: retained imposed inlet; pinned counterfactual boundary/reference calibration; resolved stationary/developed 2D H27 solution. No stationary claim, numerical convergence, physical-envelope PASS or Task3 readiness is established by this source audit.

## 8. Follow-up: applicable wall principal action and retained-inlet diagnostic

The historical principal gate does not establish admissibility of the H27 wall tensor. In `RANSWallEvolutionDiagnostics.cpp:752-760`, an inner observation evaluates the old inner covariance and immediately continues; `:762-797` applies the full WJ derivative/principal search only to outer observations. The helper at `:37-70` includes normal-stress derivatives in an incompressible x-z rank-one search, with a separate transverse polarization. It is not a global 3D certificate. H26's adapter contract (`wj-h26-contract.md:87-107`) checks the full new tensor, trace, covariance, stress work and collar; its `h26/checks.h` does not supply a wall full-principal admissibility certificate. In particular, positive auxiliary shear slope does not establish the sign of the H27 tensor's normal-stress derivative.

For the independent 2D diagnostic let `tau=k/Ew`, `S=tau*(G+G^T)/2`, `O=tau*(G-G^T)/2`, `s=tr(S^2)`, `w=tr(O^2)`, and `F=N^3-cp*N^2-(2w+2.7s)N+2cp*w`, with `cp=1.8`. For a trace-free directional gradient H,

\[
 \delta N=\frac{2.7N\delta s+2(N-c')\delta w}
 {3N^2-2c'N-2w-2.7s},\qquad
 R_W=k\{2I/3+\beta_1S+\beta_4(SO-OS)\},
\]
\[
 D=N^2-2w,\quad\beta_1=-6N/(5D),\quad\beta_4=-6/(5D).
\]

Differentiate both coefficients, the commutator and the normal entries. For `n=(cos(theta),sin(theta))`, `a=(-sin(theta),cos(theta))`, `H=a tensor n`, compute

\[
 \mathcal Q(\theta)=H:\{\nu(H+H^T)-\delta R_{wall}[H]\}.
\]

At chi=0, only the off-diagonal derivative is replaced by `-nu_A*(H+H^T)_xz`; the full WJ normal derivatives remain. These purely x-z directions are applicable even when the supported mean flow is 2D; a negative witness with a nonzero spanwise propagation component alone would not establish a 2D obstruction. Scalar variables are frozen in this local momentum principal test; a positive result alone would still not certify the full coupled system.

The owned `h27/inlet_principal.py` diagnostic loads the unchanged historical `NeutralReference` and exact stored calibration (no table refit), scans physical inner radii and cutoff `(exp(-12)H,H]`, includes supplied reference-budget root cutoffs, checks the full-height outer inlet and both H owners, and compares analytic derivatives with four central-difference steps at witnesses. A negative tail witness at positive r with `Ew=A/r` applies to every cutoff `c<=r`. The distinct `D_inner->infinity`, `k_reference->K` control keeps reference shear and the zero-endpoint integrated budget, to separate a reference-profile problem from a fixed-inlet stress problem. None of these diagnostic controls changes the application inlet or production model.

Ground alone cannot prove failure for every positive cutoff: at fixed retained ground k and gradient, `q=k*G*c/A->0` as `c->0`; the WJ diagonal derivative tends to zero, and the planar action tends to `nu+nu_A*cos(2theta)^2>0`. Such tiny cutoffs may fail the required reference budget; that separate condition must be assessed rather than assuming every positive cutoff is a physically admissible root. Finite radius/cutoff scans are witness searches, not an all-cutoff theorem. The verifier is prepared for an admitted single-thread Linux run and was not executed numerically on the Mac.
