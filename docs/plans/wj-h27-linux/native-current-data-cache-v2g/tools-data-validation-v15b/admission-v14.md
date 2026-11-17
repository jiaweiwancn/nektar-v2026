# Native terminal-output foundation admission

The active approved-v2 physical equations, reference, inlet, boundary laws and
current-state guards remain unchanged. The separately proposed v4 transport
model is inactive. Source13 now supplies the numerical inverse-cache/history
foundation. The next bounded experiment addresses terminal output and an
independent reduction of actual evolving physical observations. It cannot close
the unresolved model/reference, long-pilot, h/p/q/dt or production-cost gates.

Add an opt-in `WJTerminalOutput=True` with finite `WJTerminalBegin>=0` and
`WJTerminalEnd>WJTerminalBegin`. Preserve existing periodic scientific output
when disabled. Force both accepted-time brackets of each terminal endpoint,
including endpoints that are not integer multiples of dt. Retain the exact
requested interval; never move an endpoint to the output grid or fabricate an
accepted state at an interpolated time. Include final output and emit a
create-only scientific index with exact source/reference/signature, interval,
frame clocks, file identities, reasons and provenance. Separate fresh initial
observations and imported histories from steps accepted in the current process.
Require complete coefficients for this first full-state output implementation.
Count all serialization/index work in complete initialization/step/export costs.

The independent Linux reducer reads the indexed original signed Q3/Q2 owner
maps, verifies byte hashes, reconstructs physical fields and full-v2 stress
without a native constitutive call, and retains owner gradients at coincident
coordinates. Use tensor Gauss5 plus endpoint-inclusive samples on every cell;
give the additional endpoints zero quadrature weight. Report finite sampled
MAX separately from spatial RMS and make no continuum-supremum claim. Include
U/W, k, modified pressure, full Reynolds stress, velocity gradient, outer
epsilon, interface K/Gamma/cutoff and ground traction. Keep the open-traction
pressure level.

For each observed quantity interpolate the observation linearly between its
actual scientific frames. Clip each segment to the unchanged requested
interval. Integrate the resulting linear observation and its square exactly;
verify both against an independent three-point time-Gauss calculation. Compute
variance by integrating squared deviations from the already computed mean,
avoiding a subtraction of nearly equal raw moments. Reject gaps or incomplete
coverage. Serialize double-precision reduced arrays and bind their geometry,
owners, scales and complete time weights. Reload and compare the reduction;
do not silently discard the raw qualification evidence.

Compare scheduled frames against output at every accepted step on the same
physical trajectory. Require exact state/history hashes at common clocks.
Apply the unchanged H27 numerical RMS/MAX scales to both mean and standard
deviation differences: U/W1e-4, k1e-3, Pi/R/traction1e-3, gradU1e-2 and outer
epsilon/Ew1e-2. GradU uses the fixed retained-inlet local shear magnitude
tau0/(nu+nuA(Kin,r)) below H and u_tau/(kappa*r) above H for every component.
K/Gamma/cutoff use1e-3 of Kin/Ain-over-H/H. This is a temporal
output-sampling foundation, not a spatial or time-integration accuracy proof.

Use the retained700-cell physical case, dt=.005, CoupledJFNKV1,
OperatorDriftV1/.02/eight bases, nonlinear3e-12, linear1e-12, no Anderson and
all current acceptance guards. Advance to .06 with the nonaligned terminal
interval[.0173,.0537]. Run periodic science stride3, dense stride1, and a
stride3 six-step prefix plus cold six-step continuation. The segment
manifest authorizes prefix through step6 and continuation after step6. Any
duplicate bracket at the restart boundary must have identical scientific
state/history hashes. Do not combine abandoned future branches. Bind at least
one independent current-equation endpoint audit. Add a tiny MMS scheduling
case with aligned endpoints to cover exact-grid behavior, a zero-begin MMS
case to identify fresh initial observations, and one bounded failing MMS step
to verify an incomplete index does not label a fresh history as accepted.

This prospective validation revision is made before generating/executing any
source14 fixtures. Moving the restart from3 to6 puts it strictly inside the
requested test interval and exercises nonzero integration contributions from
both processes. The compiled source14 capsule is unchanged; the final checker,
fixture generator and revised admission are retained as a separate hash-bound
validation overlay. No scientific gate or physical equation changes.

Freeze a new immutable source14 capsule after committing source13. Admit
materialization60s, four-worker build600s, fixtures30s, four physical native
operations760s, tiny MMS30s, early independent MMS scheduling20s,
reduction/comparison600s, one independent physical
current-equation audit180s and evidence90s: <=2370 additional aggregate Linux
operation seconds, under the same43200s account and3600s protected evidence
reserve. Failures and overlapping operations remain charged. Use one science
rank, all six numerical-library thread counts1, at most72 actual shared-host
workers, science RSS8GiB/build16GiB, the existing8GiB source/build and8GiB
science limits and15GiB free reserve. Resolve every unknown process outcome
before retrying. No full Task3 run is admitted.

Full D still needs representative developed/perturbed states, current-model
accuracy-qualified cadence and storage allocation, a tested monitoring
escalation policy, complete restart provision and a complete coupled-cost
forecast. The original Task3 horizon505.02466521278404 and terminal interval
[447.8483832973273,505.02466521278404] remain unchanged.
