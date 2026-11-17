# Unchanged-v2 implicit Newton qualification admission

This is a numerical algorithm change only. Active v2 stress, inner diffusivity,
cutoff/reference, physical inlet/boundaries and conservative BE/BDF2-EX2
equations remain unchanged. The separate v4 scientific decision is pending.
Do not implement v4 or infer any v4 qualification from this work.

`WJImplicitSolver=PicardV1` remains the default control. The opt-in
`CoupledJFNKV1` solves the same fully current implicit state with the tested
Newton/FD-GMRES engine. Actual alpha is 1/dt at BE startup and 1.5/dt later;
complete accepted essential/Gamma history enters the stable mass differences.
EX advection remains fixed accepted history and has no current transport
Jacobian. The inverse approximation contains the complete frozen-scalar
velocity tangent, current positive scalar diffusion/reaction, exact physical
mass, transported-Q3-k pressure action and dependent Gamma trace derivative.
It is not a new time scheme, pseudo-time source or frozen-current residual.

Uncommitted Newton/FD states retain polynomial positivity and current typed
K-range/covariance/backflow guards. Only typed invalid trials are rejected;
unrelated root/range/input errors remain fatal. Every physical acceptance
recomputes the full planar principal certificate and current residual with the
gauge multiplier removed. No Newton candidate is persisted as an accepted
history. A typed linear failure alone permits one cold retry from the original
history at the same dt; the injection marker is preserved. A single cooperative
step deadline includes forcing, both possible solves and final validation;
retry passes only the remaining allowance. No history is committed after the
deadline. Complete operation watchdogs remain the hard external wall cap.

Current residual/PC callback counts include interrupted attempts. If an
exception prevents returned engine statistics, `counts_complete=false` marks
the missing Krylov/Jacobian counts while retaining actual wall time, rebuilds,
callback counts and disjoint costs. Failure evidence exports these facts in
the failed attempt, separately from the prior accepted state's diagnostics.
Newton trace records are explicitly uncommitted bases/Armijo candidates,
exclude FD probes and are not accepted updates. Inclusive trace observations
are never added to disjoint cost categories. Native coefficient/gradient/source
sampling is timed exclusively from separately measured WJ/guard work;
Newton/GMRES/Armijo orchestration is timed after subtracting disjoint callbacks.

## Initial bounded Linux matrix

Freeze source-v12 from source-v11 plus the complete current Mac overlay and
hash every file. Retain all older capsules, binaries, failures and accounts.
Build the full normal-main factory executable with four workers; admit
materialization <=60 s and build <=600 s (measured prior build about 218 s).

Create fresh paired Picard/Newton fixtures. Use identical fields, sources,
quadrature, tolerances, masks, clocks and histories; anderson_depth=0 for both.
Only the implicit solver differs. Initial tests are:

* PolynomialScalarsV1/body=1 continuous MMS, nx=nr=2, dt=.02, four steps,
  nonlinear 3e-12 and scalar-control linear 3e-13. Exercise BE plus BDF2,
  sampled body, changing Gamma/cutoff, gauge and complete scalar histories.
* Historical MMS, dt=.0025, eight steps, nonlinear 3e-11 and control linear
  1e-11; retain the previously documented micro-dt roundoff policy.
* Physical retained-inlet startup on the original 700 cells, dt=.005, three
  steps, nonlinear 3e-12, linear 1e-12. Do not call it a qualified pilot or
  production cost. Compare complete fields/current histories at the same end.
* A Newton MMS retry at step 2, a cold continuation from step 1, and a
  one-Newton-iteration failure and tiny cooperative-deadline failure. Verify
  accepted-history preservation, compatibility/algorithm rejection and actual
  source/signature/cache binding. No broad scientific exception retry.

Independent continuous/physical transient oracles assess complete current
equations and both integration rules. The separately frozen
`implicit_algorithm_check.py` enforces <=1e-8 field/history differences with
fixed reference scales, original owners and no pressure-level shifting.
This matrix is a prerequisite; full order/physical h/p/q/dt and warm cost gates
remain open even if it passes.

Admit fixtures <=30 s, paired MMS native operations <=180 s total, physical
control/Newton <=360 s total, retry/restart/negative controls <=180 s total,
independent equations/comparisons <=720 s total and evidence <=90 s. Together
with materialization/build this is <=2220 aggregate operation seconds in the
same 43200 s account, with 3600 s protected for evidence. All measured failed
and overlapping operations carry forward; no new allowance is created.
Use one science MPI rank, six library thread counts one, <=72 actual host
workers, science RSS 8 GiB/build 16 GiB and existing storage/free-space caps.
No full Task3 horizon is authorized.
