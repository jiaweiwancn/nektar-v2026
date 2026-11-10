# Nonlinear duration and response qualification

Starting revision: `dfd94161fd9d97f7f1a8ce2966e9cc9461b6d686`, clean implementation
worktree, no newer local or remote commits. The previous evidence and eight
unresolved epsilon temporal-order observations remain immutable. Specifications
retain SHA256 `056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64`
and `17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40`.

New authoritative Linux root:
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-nonlinear-20261110`.
Source, build, copied immutable runtime, diagnostic runtime, tmp and trajectories
are isolated. The primary Linux checkout and concurrent CLS/WJ jobs are untouched.
The initial relevant thread census is 39; numerical libraries use one thread,
initial native concurrency is two, and the aggregate limit is 72. Resource
availability is checked again at each launch. The existing 3000-unknown guard
and unsupported-route guards remain unchanged.

## Contract declared before new evolution

No governing turbulence/wall equations, coefficients, sources, molecular policies,
one-sweep stage order, pressure convention or native histories change. Eta and
corrections remain zero. No reference cancellation, resets, trace repair,
radicand clipping or new regularization. B5/B11 retain 2355/2997 unknowns, p5/p3
spaces, fixed model-specific meshes and 48/12/12/12/12 integration.

Admit A93 and A95 Natural Accepted independently. Validate both .0016/.016
short comparisons, both model-specific common-test/spatial receipts, scalar
observation coverage, actual consumed inputs and byte-level runtime equivalence.
E94E95 success cannot satisfy either admission. New four-member cohorts use
.002/.001 and milestones .016,.1,.5. If uncertainty fails while admissible,
allow one fresh .001/.0005 cohort; further refinement needs a recorded limiting
observable and a new bounded contract. No checkpoint timestep is edited.

At every milestone all four field/completion/identity/source gates and
`max(delta_t_B5,delta_t_B11) <= max(1e-7,.1*delta_h)` must pass before continuing.
Neutral field RMS limits remain U,k,epsilon,Pi <=1e-3, w<=1e-4, mean wall/full
interface traction <=1e-3, on the existing fixed physical scales. Full native
and fixed common-test momentum-x/inner-k residuals, and non-floor nonlinear-model
vertical momentum, require fine/coarse <=.85 above the existing 1e-8 residual
floor and five times evaluation/quadrature uncertainty (absolute 1e-10).
Analytical balance, reconstruction, reduction, field/time accuracy and execution
are separate statuses. Other equations and new fine modes remain visible.

Final requirements include both spaces' alpha-zero audits, all four trajectories'
scalar observations through .5, frozen quadrature on both, and independent full
R/Q/production/normal-traction checks on actual native samples. Accepted states
are counted uniquely as well as with repeated restart endpoints. The one-column
epsilon constant test remains a prescribed-row reaction identity.

After a model passes neutral duration, use StreamfunctionPolyV1 amplitude 2e-5,
and a selected half-amplitude 1e-5 control, at .016 then conditionally .1. The
streamfunction and first normal derivative vanish at the essential boundaries
and interface. Bind the actual perturbed initialized vector/history separately
from the neutral reference and boundary map. Same-space dt controls require
identical actual initialization. Compare response to a matching unperturbed
numerical trajectory, never to a regenerated neutral initializer.

Response spatial differences must be <=max(5e-7,.1*response RMS); timestep
differences <=max(1e-7,.05*response RMS) on the existing physical scales. A
component below these absolute floors is unresolved. The half-amplitude test
reports departure from linear scaling with limit max(5e-7,.1*half-response RMS),
without assuming exact linearity or monotone decay. Retain the existing hard
physical and absolute field/traction accuracy gates. A93/EP is the preselected
short wall-source/production comparison, with all four companions at .016;
extend only for a specific resolved wall question. No temporal order, global
stability or open-outflow claim follows from these experiments.

Pressure, vertical velocity, Qxz and full-production responses must be resolved
for extension to .1; individual smaller scalar/stress responses retain their
unresolved labels. Stress response uses scale u_star^2 and production response
uses epsilon(H). Initial coefficients and native zero-step history are sealed
separately from the neutral reference; their stable artifact identities are
checked on every same-space temporal comparison and exact resume.

For exactly divergence-free planar flow, A93's Q:S vanishes because the two
nonzero strain eigenvalues are opposite and Omega^2 is a multiple of the
identity on that plane. A95's commutator has zero contraction with S generally.
This special kinematic fact does not permit dropping full production: the
variable viscosity responds, Qxz can be nonzero, and the discrete raw gradient
can have divergence. Report the nonlinear and isotropic production contributions
separately and keep the general three-dimensional constitutive evidence intact.
