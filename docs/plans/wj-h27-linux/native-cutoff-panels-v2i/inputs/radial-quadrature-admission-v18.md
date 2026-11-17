# Approved-v2 radial quadrature diagnosis and refinement

Source18 repairs the known-material H trace without changing physics or the
2e-13 check. Both continuous MMS families and the700-cell12-step physical
case have exactly matching source15/source18 coefficients in all seven
checked state/history roles, passing the independent current equations,
original-owner/dependent/lift/mass-dual comparisons and unchanged1e-8 gates.

The28x200 radial steady state converges natively (four accepted Newton steps,
77 linear iterations; true norm6.882641149697246e-12). Its independent
20x20 audit finds momentum1.6028390980667496e-12, continuity
2.657574229526232e-15, epsilon1.7324355799200504e-12, but
k1.8293736150710952e-7, above1e-8. Boundary checks pass. It is not a qualified
reference. Its inner-k physical RMS/MAX.1327946821/.440731146 still miss the
fixed.05/.10 envelope. Do not attribute this equation discrepancy to physical
development or promote a same-quadrature algebraic zero.

Diagnose and refine quadrature while keeping all physical equations, native
source18, mesh/spaces, reference table, inlet/roughness, current residual
limits and guards fixed. A moving-cutoff crossing of a newly inserted radial
face is a hypothesis, not yet an established cause. Retain the prior8x12
native state and every comparison even when an independent subgate fails.

1. Audit that8x12 state independently at32x32 to distinguish oracle-order
   sensitivity from the native quadrature discrepancy.
2. Solve the same steady equations with native16x12 integration. Use the
   converged8x12 field only as a fresh initial guess; rebuild all dependent
   quantities/explicit history, keep time0, nonlinear1e-10/linear1e-12 and
   FullGradientV1. Require native current convergence, then independent32x32
   and48x48 current-equation/boundary/accounting checks. All failed gates stay
   visible. Compare complete current signed-space fields and derived state
   against8x12 under the frozen numerical budgets.
3. Only if16x12 misses an independent equation gate, admit a32x12 native
   refinement from the16x12 guess, with the same independent32/48 audits and
   signed-space comparisons. Check measured16x12 memory before admission;
   never exceed8GiB to obtain a result. Any failure unrelated to quadrature
   stops dependent execution for diagnosis.

Prior radial8x12 native time454.015566s, max tree RSS5,286,068,224bytes;
equation audit113.149694s. Higher native x order increases sample storage;
forecast16x12 at6--7GiB and32x12 at7--8GiB, keeping the hard watchdog.
The better initial guess may reduce nonlinear work; it is not a cost promise.
Forecast native16 at400--900s and conditional32 at700--1500s; independent
order32 at250--400s and48 at500--850s. Use actual whole operation costs.

Prospective upper bounds: original-state order32 audit400s, fixtures30s,
native16 1400s, its two independent audits400+900s and same-grid assessment120s;
conditional native32 1800s, its audits400+900s and assessment120s; evidence120s.
Total<=6590 additional operation seconds, cumulative radial/trace/quadrature
phase<=11570s including the previous4980s admission. This is inside the same
43200s campaign account, carried forward from OPEN12204.860811439692s before
the radial/trace phase, with3600s protected evidence and no allowance reset.
Do not spend the conditional branch until its stated trigger and memory
review apply. All serial/thread, <=72 actual-worker,8GiB science RSS,
source/build8GiB/science8GiB/free15GiB caps remain binding. New isolated output
prefixes only. V4 remains inactive and no full Task3 horizon is admitted.
