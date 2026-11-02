**M5 closure for the conservative curved CLS/EP-AWF patch**

M5 is closed on 3 November 2026 on `codex/cls-epawf-revised`.
The full numerical qualification passed at
`42b8ba83897cca1202b91032a97492c7d463c503`; supplementary note checks and
25 regressions passed at `e6bb342afaafeb05cf1ef9c195e1ab0d74263055`.
The cached-probe replay fix and targeted revalidation passed at
`387e1e31d1ca38b77e59b0ea78c090aec8c2de90`. All 28 registered tests are
covered, with the original failed probe receipt retained. The final audit
closes all ten M5 equation-owner IDs and verifies 22 actual stage/output
reports. The maximum accepted scaled stage norm is 9.52e-13; the largest
solve has 2,016 unknowns. Native all-mode metric errors at order16 are below
3.56e-15. Independent physical VTU coordinates/values agree within
2.23e-15/3.56e-15. M6 and later milestones remain open.

**Material Passport**

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: implementation and independent computational verification.
- Origin Date: 2026-11-02.
- Verification Status: VERIFIED for the conservative mapped patch, geometry,
  full equations, constraints, native stages, MPI and complete-history restart.
  Physical CLS stability and general inner-k positivity remain FAILED AS SPECIFIED.
- Version Label: cls-epawf-revised-m5-v1.

The governing specification is the unchanged
[revised note](/Users/jay/Documents/nektar-v2026/docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex),
its normal-coordinate, conservative curved transport, tilted matching-surface,
epsilon-reaction and metric-identity equations, and the unchanged
[M0 contract](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-contract-v1.json).
The original CLS tensor, Yap law, rough AWF, molecular conventions and cutoff
constants are retained. No constitutive retuning, transported-variable floor,
production clip, or new inner epsilon equation is introduced.

**Physical geometry and native spaces**

`RANSWallNormalMap` supplies analytic flat/inclined-plane, outward cylinder,
Gaussian graph and piecewise cosine-hill maps in 2D or full 3D. The Gaussian
surface is the explicitly smooth initial hill verification; the cosine hill
retains its C1 feet and one-sided curvature. A native patch face is required
at each foot. H may vary with both wall coordinates.

The computational XML remains a canonical affine tensor chart `(q1,q2,rho)`.
Native Modified_A CG bases and assembly maps are pulled back through
`Phi=X(q)+r*n(q)`, with `r=rho*H(q)/H0` and matching chart `rho=H0`.
This is an analytic mapped native-CG method. It does not import arbitrary
Cartesian curved XML geometries. The outer region uses the same certified
normal tube as the inner layer. Cartesian velocities, independent regional
pressures and their common volume-mean gauge retain the M4 spaces/backend.

With `B=grad_Gamma(n)`, positive outward circular curvature is `+1/R`.
The legacy graph curvature `dt/ds=kappa_old*n` has the opposite sign.
The map uses `A=I+r*B`, `J=det_T(A)` and `B_r=B*A^{-1}`. Every gradient,
physical mass matrix, volume-mean pressure constraint, source and face load
uses these metrics. Inner ray weights are `sqrt(g0)*J*dr*dq`; chart volume
quadrature uses `sqrt(g0)*J*H/H0` once. The physical r_K and r_m are fixed
when H or the computational partition changes.

Whole-domain certificates supplement point tests. Graphs use a global C1,1
Hessian/reach bound and a positive H bound, including the complete outer
normal range. The quadratic-remainder bound makes every admitted normal
foot the unique closest foot and bounds both principal stretches away from
zero. The outward cylinder uses a nonwrapping arc and positive outward
radius. These certificates assume an isolated wall; clearance from other
walls is a separate geometric requirement. Positive Jacobians at vertices
alone are insufficient. Conservative certificates can reject a map whose
validity they cannot prove.

All interface and lateral loads use the physical cofactor area vector.
At matching, `dS_H/dq=sqrt(g0)*J_H*(n-q_H)` with
`q_H=A_H^{-T}*grad_Gamma(H)`. Consequently the per-ground mass flux is
`m_H=J_H*(w_H-v_H dot grad_Gamma(H))`. The implementation retains this
height-slope contribution and all Cartesian components of the rotating frame.

**Conservative equations and constraints**

The weak operator retains full Cartesian convection, stress, production,
normal momentum, scalar transport, physical storage and one common pressure
gauge. Its one-sided outer forcing is evaluated at constant physical r
before taking the H trace. Differentiating the tangential velocity projection
includes `grad(n)=B_r`; freezing the frame would omit terms. The roughness
force is constant along each ray and tangent to its own base normal, including
its Cartesian vertical component. It acts only below the physical roughness
height and enters the resolved energy work.

Velocity and k have common native trace unknowns. Separate regional residuals
allow direct checks of equal and opposite loads on free matching rows after
an actual coupled stage. A generic unsolved state need not have equal regional
tractions or diffusive k fluxes; those jumps are reported as diagnostics.
No additional pressure Dirichlet trace or separate pressure correction is added.

Inner epsilon remains the algebraic current-K profile. The outer projected
H trace is recomputed in every residual, including finite-difference columns,
so its current sampler/H dependencies enter the coupled Jacobian. The diagnostic
inner epsilon flux is declared as `U*Ew` with source zero. Its nonzero residual
and the outer-minus-inner epsilon flux reaction are recorded; neither is
silently imposed as another inner PDE or prescribed outer flux.

Mapped prescribed velocities receive a physical surface-L2 projection subject
to zero total boundary mass flux. Side/edge/corner values are fixed and the
minimum-norm correction is distributed over Cartesian top-face interior modes.
Its reported size is the surface RMS velocity change. This makes finite
trace projections compatible with incompressibility; it adds no volume mass
source. Actual accepted stages must satisfy both the original weak continuity
rows and a negligible pressure-gauge multiplier. Initial mapped function data
are lifted onto essential/dependent traces before native histories are created;
an initial weak-divergence projection error remains a reported warm-start
error, not an accepted stage result.

The unsolved polynomial operator controls deliberately start with incompatible
net boundary flow. Their top projections have RMS changes .238–.395; they
exercise the algebra and are not physical manufactured solutions. In the
actual manufactured CFD cases, the correction is zero or below 6.7e-18.
The zero final flux therefore records the declared boundary projection, rather
than asserting that arbitrary supplied continuous boundary data were already
compatible.

Curved and variable-H operational tests use `chi=eta=0`. Nonplanar chi and
curved anchoring are explicitly rejected in this M5 path. Isolated inclined
constant-H neutral verification may retain chi=eta=1, and repeats the flat
reference's rotation and actual native IMEX checks. The prescribed transition
and sensors belong to M6.

**Manufactured data and independent verification**

The geometry ladder progresses through an inclined plane, circular tube,
constant-H smooth graph, variable H, full 3D Gaussian surface and explicit
cosine-foot patches. Independent NumPy analytic geometry, Cartesian finite
differences, Jacobi bases and indexed dimensional CLS/AWF formulas reconstruct
the physical equations. Native maps/signs identify the FE space; exported
stress, production and volume weights are not oracle inputs.

For scalar manufactured cases, `U=0`, constant-in-space k(t), `pi=c*z` and
outer `E=A(t)/r`, with `A=k^{3/2}/2.55`. Thus
`Div(D_e*grad(E))=Cmu(0)*k^2/1.3*(1/r^2-tr(B_r)/r)`.
The retained Yap ratio is exactly `2.55/2.5=1.02`; its active source is included.
The cosine curvature is piecewise one-sided, rather than smoothed across feet.

The circular manufactured case additionally has inner tangential speed
`f=.3*(1+.1*sin(t))*(r/H)^2*(1-r/H)^2` and outer velocity zero.
It retains `-(f^2/(R+r))*n` in acceleration and the full tangential viscous
term `(D*s)' + 2*D*s/(R+r)`, where `s=f'-f/(R+r)`.
The prescribed strong sources compensate these terms, storage, production,
original dissipation and scalar diffusion independently of the weak assembly.
Pressure/velocity/epsilon errors and physical normal velocity/strong divergence
are measured under refinement. Weak pressure constraints and strong divergence
are reported separately.

All-mode weak Piola tests use actual native volume derivatives and physical
face loads at quadrature orders 8,12,16 with identically zero prescribed volume
source. Matching and lateral cancellation is checked to the isolated-interface
1e-10 tolerance. This is verified finite quadrature accuracy, not a claim of
an exact collocation curl identity for arbitrary geometry.

The complete residual oracle independently reconstructs the declared finite
outer Gauss rule. A separate higher-order diagnostic retains the first
12-versus-16-point discrepancy for a generic state crossing the CLS invariant
maximum. This prevents conflating two different finite operators with an
algebraic implementation error. Smooth manufactured quadrature refinement is
checked separately. No universal quadrature-convergence claim is made for
arbitrary nonsmooth constitutive states.

The lateral partition test uses identical global chart polynomials represented
exactly on one and two physical patch columns, comparing physical fields and
actual weak loads against the same polynomial tests. Spatial refinement of a
solution is a different comparison. MPI 1/2/4 tests reassign physical cell
ownership and compare all coupled fields, the common gauge and every native
history slot. Restart changes both rank count and ownership while retaining
the complete native IMEX history.

**Output, execution and limits**

Mapped `.chart.fld` files explicitly retain chart geometry and map metadata.
Per-field `.vtu` output contains the actual Cartesian points and scalar values,
including independent regional pressures. It is suitable for physical
visualization without mislabelling chart coordinates as Cartesian terrain.

The M4 dense backend remains bounded to 3,000 coupled solve unknowns. This
milestone qualifies conservative geometry and a bounded reference solver;
scalable distributed sparse assembly, arbitrary unstructured terrain imports,
M6 transitions and robust long-time terrain validation remain later work.
The previously demonstrated physical CLS and inner-k failures remain unchanged.

Linux validation uses only Tailscale (`jay@100.105.247.110`). The source is
synchronized into the isolated run checkout, with the unrelated HIT main
checkout and shared installation preserved. Compiler concurrency is 16;
solver/tests are sequential, with MPI 1/2/4 and OMP/OpenBLAS/MKL/BLIS threads1.
Heavy artifacts remain below
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised`.

The [integration receipt](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m5/integration.json)
binds all source pins, passing coverage, equation delivery, field audits,
retained failed attempts and the Linux raw-artifact index. The
[evidence directory](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m5/README.md)
contains the complete qualification and note/replay/regression receipts.
All original source hashes remain unchanged. The unrelated Linux HIT checkout
remains clean at f4b70baec10016dd098005bc5a5507b87e4a57bf.

**Completed main qualification**

The complete M5 qualification passed at
`42b8ba83897cca1202b91032a97492c7d463c503`, with 4,352 independent checks
and 49 recorded commands. The native geometry probe adds 5,178 checks and
26 expected rejections; its independent NumPy audit adds 4,608 checks over
13 geometry cases. The supplementary note-example probe is registered
separately at `e6bb342afaafeb05cf1ef9c195e1ab0d74263055` and passes
454 independent checks. A fingerprint confirms that all 688 pre-existing
qualification source inputs and 41 runtime executables/linked libraries are
identical across these two revisions; only the supplementary test inputs
and CMake registration were added.

The complete legacy run passed 25 tests, including full M3 and M4, and
exposed one cached-session probe replay failure. Existing `case.opt` and
`.opt` files made a clone inherit an optimizer sidecar as its first input;
the installed session reader then attempted an invalid short-name suffix
comparison. At `387e1e3`, both affected probes explicitly reopen the original
XML. The only code edits are that filename selection and comments. The
production solver, model kernels and libraries remain identical. Cached
repetition, targeted M4 spatial verification and all six curved operator
oracles pass 63,612 replay checks for these two probe executables; the original failed receipt
is retained.

| Check | Measured result |
| --- | --- |
| Complete Cartesian momentum, k, outer epsilon and pressure weak loads | Maximum absolute native/independent block discrepancy 3.11e-15 |
| Matching mass/momentum, full traction, k flux and epsilon reaction | Maximum scaled interface/integrated-budget discrepancy 5.44e-15 |
| Conservative physical surface-L2 projection | Maximum independently integrated remaining boundary flux 5.28e-16 |
| Identical fields on one/two lateral patch columns | Maximum scaled weak-moment differences 6.12e-7, 1.15e-8, 1.45e-11 at outer Gauss orders 12,24,48 |
| Smooth variable-H manufactured quadrature, orders12→16 | Maximum scaled physical state/native-history change 2.68e-12 |
| Fixed-geometry h refinement, variable H | Epsilon RMS error 6.01e-5→4.81e-6 |
| Full versus condensed curved stages | All coupled fields, pressure gauge and complete native histories agree within 1e-8 |
| MPI/ownership | Serial, ranks2/4 with shift1 agree within 1e-8; each physical cell is owned once |
| Restart/ownership | Three uninterrupted IMEX2 steps equal one saved plus two resumed steps on ranks4 with shift2 |
| Policy/geometry negative controls | Curved eta1, curved chi1, folded between-node geometry and a straddled C1 foot are rejected |

The fixed-physical-geometry p sequences give the following physical RMS errors.
The pressure column is the Euclidean norm of the two regional physical RMS
errors, with one common gauge alignment.

| Geometry / degree | Outer epsilon | Pressure | Strong div(U) |
| --- | --- | --- | --- |
| Circular moving layer, p3 | 6.16e-5 | 3.96e-5 | 4.80e-6 |
| Circular moving layer, p5 | 1.37e-6 | 8.82e-7 | 3.04e-8 |
| Circular moving layer, p6 | 2.11e-7 | 8.90e-7 | 2.55e-8 |
| Variable H, p3 | 6.01e-5 | 4.04e-6 | 4.58e-7 |
| Variable H, p5 | 1.34e-6 | 3.84e-8 | 1.06e-8 |
| Full 3D surface, p3 | 6.07e-5 | 3.97e-6 | 4.39e-7 |
| Full 3D surface, p4 | 8.94e-6 | 3.84e-7 | 6.32e-8 |
| Explicit C1 feet, p3 | 6.67e-5 | 2.14e-10 | 5.39e-11 |
| Explicit C1 feet, p5 | 1.87e-6 | 1.46e-11 | 2.75e-12 |

The circular pressure error plateaus near 8.9e-7 at the fixed timestep .005;
no asymptotic pressure order is claimed from this sequence. Strong divergence
is a sampled physical diagnostic and is nonzero at finite p, despite passing
weak pressure constraints. These distinctions are retained in the reports.

The generic max-switch control has `tr(G^2)` ranging from -0.248 to +0.312.
Changing its outer rule from12 to16 changes one momentum block by 2.98e-5.
The independently reconstructed *same finite rule* agrees with the native
operator to roundoff. The separate patch and smooth-MMS studies establish
quadrature convergence for their declared cases; the generic control does not
establish universal constitutive quadrature accuracy.

The note's exact manufactured circle has curvature .05, `a=1+.1*cos(.4*s)`
and `H=2+.2*sin(.3*s)`. The supplementary native map/face calculation and
independent physical-coordinate integration reproduce the negative matching
mass flux -0.08237640447 at s=1.3, while the omitted-slope value is positive
+0.04283039671. They also reproduce the pressure increment 0.12407192173 and
metric cutoff 0.009458614009. Independently formed top, ground and lateral
fluxes close the mass, arbitrary-constant-Cartesian-stress momentum and
variable-diffusivity scalar budgets. This standalone example is a geometric
manufactured-field check; the native-CG evolving CLS cases above are separate.


**Reproduction and final evidence**

The isolated Release build uses the production incompressible-solver source
list, the source SpatialDomains library, installed MPI Nektar dependencies,
and strict warning-as-error checks for the new code. The shared installation
is not overwritten. The selected geometry library is recorded by `ldd` and
all runtime artifacts are fingerprinted.

The exact main commands, source revisions, timeouts, logs and thread bindings
are recorded in the receipts. Their essential forms are:

```bash
cmake -S "$source/solvers/IncNavierStokesSolver/Tests/RANS" -B "$build" \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$build" --parallel 16
python3 "$source/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m5.py" \
  --build "$build" --output "$run/oracle"
ctest --test-dir "$build" --parallel 1 --output-on-failure \
  -E '^RANSCLSEPAWFM5(Oracle|NoteOracle)$'
```

Here source/build/run are the isolated directories below the Linux run root.
Set each of `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and
`BLIS_NUM_THREADS` to1,
`PYTHONDONTWRITEBYTECODE=1`, `TMPDIR=$run_root/tmp` and
`LD_LIBRARY_PATH=$build/geometry-lib:...` as in the recorded commands.
MPI checks use explicit ranks1/2/4 with core binding. The note oracle and
cached-probe replay have separate source-bound command receipts. A composed
regression XML records the 25 unchanged passing tests and the corrected test;
it does not replace the original failed XML.

The source/runtime audit permits only the exact original-XML filename
substitution in the two replayed probes. The solver executable and all other
runtime artifacts remain byte-identical to the complete numerical
qualification. Both probes are then independently rechecked, including cache
reuse, rather than treating their changed binaries as already tested.

![M5 physical refinement and conservative patch convergence](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m5/figures/curved-verification.png)

The [vector figure](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m5/figures/curved-verification.pdf)
and reproducible plotting script accompany the numeric tables.
