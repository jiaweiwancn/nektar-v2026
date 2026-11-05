# M8-R12 frozen Lf completion and trial-scale comparison

Work in progress; numerical conclusions will use the designated Linux system.
Reviewed baseline: 6ab80ffe1895811e958e703c2a4a51e3e0496999, no newer remote commits.
The separate Picard worktree and running jobs are preserved.

The physical base, J, h0, residual and source bindings remain R11's immutable
R9-initial / R8 iterate-10 inputs. E=2.911641510023811. No physical equation,
quadrature, certificate, reserve, radius, history or tolerance is changed.

Before execution: the frozen-QP investigation/reference/repair budget is 900 s,
excluding builds. One instrumented 96-admission corrected replay is permitted.
At most two substantive strategies are allowed: corrected dual exchange with
an extension justified by the trace, and a feasible-start method only if needed.
Any extension will be declared before that attempt. The subsequent matched
native phase has a 600 s cooperative cap, nine trials per policy, three Nf
tangents, ten Armijo fractions per tangent and 50 new Nf source values. N1 is
not rerun. Original gates are 5e-12 for original rows, 1e-7 actions, 2e-12
native root and 1e-9 for stationarity/active equality/complementarity. Before
results: normalized objective upper-bound allowance 2e-9, reference/native
objective agreement and primal-dual gap 2e-8. These optimizer accuracy tests
do not change physical or original-row gates.

Build allocation: 16 compiler workers. Numerical runs: one MPI rank, all six
numerical thread variables explicitly one, with one allocated CPU. Existing
three Picard solver jobs (three threads each, plus their MPI/runtime wrappers)
were observed and preserved, leaving the combined allocation below 72.
