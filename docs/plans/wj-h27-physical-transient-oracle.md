# H27 independent accepted physical timestep audit

This source-only extension leaves the frozen steady oracle and its controls
unchanged. `physical_transient_oracle.py` imports their independent Modified_A
basis, sparse Riesz operators, full three-dimensional WJ tensor, amended wall
stress, reference, geometry ownership, current scalar sources, and physical
boundary calculations. It does not import the native constitutive or matrix
assembly implementation. The native report's coefficients and signed maps are
inputs; native residual numbers and cached advection never enter the independent
equation.

For accepted new state `q`, input history states `a,b`, and unchanged `dt`, the
time derivative is `(q-a)/dt` at startup, then
`[1.5(q-a)-0.5(a-b)]/dt`. The complete explicit transport load is `C(a)` at
startup, then `2C(a)-C(b)`. `C` includes domain velocity advection, the sampled
convective roughness force with its factor two and negative residual sign, and
conservative scalar advection with all exterior fluxes, including the outer
epsilon H face. Current stress, pressure, known-k roughness force, sources,
diffusivity, natural data, and dependent epsilon trace remain current. Pressure
has neither a time mass nor an advection load.

The storage action is integrated independently from stable coefficient
differences. Transport at each historical state uses that state's H trace and
cutoff panels. Every saved advection coefficient, including constrained rows,
is compared against independent reconstruction using full sparse Riesz norms.
Input current/advection hashes bind the output's previous hashes. All times,
steps, inventories, basis identity, and unchanged dt are checked. Startup must
contain duplicate fresh histories. The earlier BDF state is hashed and checked
scientifically but is not falsely claimed bound to a separate earlier output
that was not supplied.

Scalar budgets include BDF inventory rate, explicit conservative outflow,
current production/destruction, prescribed natural flux, and essential
reactions. Reactions make an algebraic weak-budget identity; separately reported
owner-gradient flux budgets require spatial resolution evidence. The actual
integral of production minus destruction is never constrained to zero.

Mechanical budgets test the BDF/EX equation by the **accepted current** velocity.
The explicit history kinetic-energy boundary flux is reported separately; it
cannot replace current-velocity work against old advection. An independent
cross-advection integration identity checks this distinction. For BE,
`u·D_t u = D_t(|u|²/2) + |u-a|²/(2dt)`. For BDF2 it equals the difference of
`(|u|²+|2u-a|²)/4` over dt plus `|u-2a+b|²/(4dt)`. The physical kinetic-energy
BDF rate is also reported; it is not identified with this BDF2 G-energy rate.
Current raw turbulent production cancels algebraically in the mechanical-plus-k
account, with weak divergence and sampled roughness work retained.

Default independent Gauss orders are 20 and 24 on owner-aware split panels.
Their residual difference must be at most 20% of the unchanged residual limits:
momentum/continuity `1e-9`, k/epsilon `1e-8`. Cache error must be at most `1e-10`
in the fixed full Riesz units. Accounting identities must be at most `1e-10`
in the frozen physical scales. Physical profiles, retained-inlet envelopes,
interface owner gradients, and optional same-time refinement comparison reuse
the independent steady implementations. Instantaneous steady residuals remain
explicit observations, not transient convergence gates.

The new control program independently creates all signed mass and advection
coefficients from exact spatial polynomials and quadratic time factors with
unsplit Gauss8 integration, splitting only the sampled force's support. It tests
BE and BDF2/EX2, all coefficient moments, storage inventories, mechanical time
work, missing mass/convection/AWF, erroneous BE/EX1 substitutions,
nonconservative scalar transport, malformed history rejection, and a corrupted
cache that must fail comparison while leaving the independent residual
unchanged. The fixture deliberately violates physical boundary conditions and
must remain `NOT_MET`; it is neither a physical solution nor a manufactured
forced calculation.

Linux commands, admitted and launched only by the root resource runner:

```sh
python3 h27/physical_transient_oracle_test.py --table TABLE --out NEW-controls.json --seconds 60
python3 h27/physical_transient_oracle.py --native RUN.result.json --table TABLE --out NEW-audit.json --seconds 300
```

The audit also accepts a standalone coefficient-bearing accepted scientific
state. Optional `--compare-finer OTHER.result.json` requires the same physical
time and independently reconstructs that state's time-discrete residual and
caches. It is a comparison, not a temporal order proof. Outputs are immutable,
single-threaded, with an explicit wall cap and an 8 GiB address-space default.
A measured case can request up to 20 GiB through `--address-space-gib` under
its named resource admission, subject to stricter inherited limits and a
separate live RSS watch. The two orders run sequentially. Only the coarse
actions and residuals remain during the fine solve. Exactly equal history
coefficients share independently evaluated transport. All saved caches still
get their own comparison. The outlet-B reaction rows use the complete
accepted time-discrete momentum action. Exceptions preserve a partial
`NOT_MET` report.

No numerical execution occurred on the Mac during source preparation. Overall
H27 status remains `NOT_MET` regardless of a successful accepted-step equation
subgate: continuum admissibility, physical flux-budget refinement, temporal
order, finite perturbations, and the long-time objective are separate evidence.
