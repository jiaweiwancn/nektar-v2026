# M7 — Admissibility and perturbation diagnosis

## Material Passport

- Mode: experiment-agent run and reproducibility validation.
- Scope: unchanged revised-note CLS and equilibrium-preserving AWF.
- Branch: `codex/cls-epawf-revised`; baseline M6 `f41c7e46986928d380866f16df284c791b4d5dac`.
- Authority: approved implementation plan, frozen M0 contract, revised note.
- Execution: isolated Linux source/build and run directories, Tailscale SSH
  `jay@100.105.247.110`; numerical library threads one; compiler jobs 16.
- Status: implementation/qualification in progress. No stability qualification.

## Experiment contract fixed before execution

The original CLS and scalar AWF source and the frozen note/contract stay
unchanged. Diagnostics are opt-in (`RANSWallDiagnostics=1`). The accepted
state must have nonnegative sampled inner k and strictly positive sampled K,
outer k and outer epsilon. A rejected candidate is recorded as a rejection,
never a valid trajectory. Covariance and tangent signs are measured independently.
No clipping, filtering, SVV, constitutive retuning or positivity intervention
is introduced. There is no claim about polynomial positivity everywhere.

Every initial/accepted native stage is audited at the actual split ray source
quadrature, outer source quadrature, fixed prescribed-source quadrature,
physical mass quadrature, both matching traces, physical K and transition
samplers, native GLL reconstruction and a denser equispaced reconstruction.
Each element uses `2*p+3` points per coordinate for the dense set.
The independent Jacobi oracle reconstructs every accepted k polynomial on that
set. A deliberately negative between-node polynomial must be rejected even
though native GLL and source GL subsets are positive.

The native nonlinear campaign separately perturbs velocity, k and epsilon,
then their combination, a near-interface disturbance and a disturbance crossing
the prescribed upstream transition on a cosine-foot map. Amplitudes, wavelengths,
dt, p and mesh are changed independently. Isolated flat verification uses
eta=1 with the original frozen defect; unanchored paired controls and all
transition/terrain runs use eta=0. An unperturbed native benchmark spans five
flow-through times; bounded nonlinear perturbations are short diagnostic runs,
not arbitrary-time stability tests. A separate matching-height study fixes
physical roughness, nu, rK and rm and reports the resulting wall-model dependence.

The three-component solenoidal transverse test uses a clearly labelled reduced
linearized Fourier/native-CG system supported on `4<z<12`, strictly above H=2.
The compiled unchanged constitutive law supplies its tangent. The complete
transverse stress block and base advection are implicit in native IMEX1/2;
the explicit prescribed body term is zero. Independent indexed stress,
quadrature, matrix exponential and time recurrences check this experiment.
Wavenumbers 1/2/4, amplitudes 1e-5/1e-4/1e-3 and independent dt/p/h changes
test the adverse growth described in the note. This is not a nonlinear 3D
terrain trajectory.

Acceptance remains the frozen contract: scaled coupled residual <=1e-9,
isolated identities <=1e-10, MPI physical differences <=1e-8. Strong divergence
and reconstructed traction/flux jumps are measured separately from weak
constraints. Native stage energy uses the actual physical mass operator,
the unmasked velocity/constant-k test, prescribed and roughness work,
essential-boundary reactions and the separate anchor work. Trapezoidal
integrals of sampled accounting defects are labelled as such. They do not
replace the native IMEX storage identity. Negative principal diffusion,
growth, or loss of sampled admissibility leaves robust-model qualification
failed even when implementation checks pass.
