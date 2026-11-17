# Current full velocity tangent for steady preparation

The source-v10 physical steady continuation reduces its true fixed-scale norm
from 0.003037822745 to 0.000836441294, then exhausts 320 Krylov iterations at
the next base. The achieved linear norm is 0.000274354069, above its requested
forcing target. Accepted history and a separately marked unqualified candidate
are preserved. This is a preconditioning failure, not proof of nonexistence.

Source v11 uses the exact current frozen-(k,E,nu_A,r) derivative of V with
respect to all four supported velocity-gradient entries in the STEADY inverse
approximation. It includes all cross terms and the full three-dimensional WJ
deviatoric trace terms for weakly divergent planar fields. Owned native signed
rules assemble the entire matrix, including essential reaction rows. The
matrix is allowed to be nonsymmetric/indefinite; it is not clipped or replaced
by a positive stress. Current advection, pressure/AWF, scalar and dependent
Gamma blocks remain as before. Actual residuals and finite-difference Jacobian
actions still use the original complete nonlinear equations. The advancing
algorithm is unchanged. `WJSteadyVelocityPC=PositiveFrozenV1` retains the older
inverse approximation; `FullGradientV1` selects this one.

Native/independent probe qualification precedes the physical reference: exact
directional versus assembled actions, two finite-difference levels with an
independent full-WJ value implementation, independent Jacobi matrix assembly,
both constitutive backends, signed-basis transformations, weak divergence,
collar endpoints, negative-tangent retention and malformed input rejection.
Neither this approximate block inverse nor the probe certifies the coupled
physical reference.

The scalar Krylov loop and final check now form the same explicit `A*x-rhs`
vector before the fixed physical dual norm. The older fused `rhs-A*x` and
`A*x-rhs` evaluations rounded differently: 1.0032e-12 versus a 1e-12 target in
a physical control and 3.1799e-13 versus a 3e-13 target in the tighter MMS.
All tolerances remain unchanged. This consistency repair is tested with the
tighter MMS, rather than relaxing the acceptance threshold.

The preliminary physical split comparison passed full-row actions below
2e-14 but its rigorous pressure field bound was 2.888e-8, above 1e-8, at
coupled tolerance 3e-11. Repeat at 3e-12 coupled tolerance; retain that stricter
field gate. These first-step timings remain diagnostics, not production costs.

Prospective v11 admission in the continuing H27 account: build <=600 s;
native/independent tangent <=150 s; focused scalar checks <=120 s; seeded
steady MMS and analytic reference <=240 s; physical steady continuation plus
independent assessment <=1730 s; tighter physical split checks <=450 s.
The remaining MMS tranche keeps its original cumulative 1200 s cap, carrying
forward the already charged 99.5993 s from v10 and preserving every failure.
Use four build workers, one science rank, six numerical-library thread counts
one, and the existing 72-worker/RSS/storage safeguards. No long pilot or full
Task3 horizon is admitted by this diagnostic tranche.
