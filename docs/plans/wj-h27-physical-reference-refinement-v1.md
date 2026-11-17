# Same-inlet v2 reference refinement, prospective admission

Source-v11's 700-cell steady calculation converges after eight additional
Newton steps and 432 Krylov iterations from the explicitly marked earlier
candidate. It advances no physical time. Independent momentum/continuity/k/
epsilon residuals are 3.214e-12 / 2.185e-14 / 3.916e-10 / 2.155e-11 and the
weak budget identities pass. Native covariance/k and principal/mu minima are
0.278765 and 0.00400050. This establishes one current discrete stationary state,
not a resolved physical reference or positive H27 closure.

The independent retained-inlet comparison fails the fixed k envelope. At
x=-7.5 the inner k RMS/MAX departures are 0.09965/0.40790 in Kin units,
versus fixed limits 0.05/0.10; x=-5 gives 0.13875/0.26489. U/W and outer
epsilon have separate reports. Do not change those limits, replace the inlet,
or interpret a steady algebraic zero as mesh convergence.

Construct a finer tensor mesh by bisecting EVERY original XML x interval and
physical-elevation interval. Derive the inventory from the retained XML rather
than guessing uniform dimensions; the resulting mesh has 2800 cells. Keep
Q3/Q2, exact geometry/roughness, physical boundaries, table, full v2 stress,
all scalar equations and 8x12 native split integration. Project the converged
coarse fields only as an initial guess, rebuild the dependent Gamma trace and
fresh histories, and solve the fine current steady equations. Use the same
source-v11 executable and full-gradient steady preconditioner, Krylov dimension
80, linear tolerance 1e-12 and nonlinear tolerance 1e-10.

Compare full current independent equations and actual physical budgets, plus
U/k/epsilon/R/Pi/gradient/traction fields at common physical owner points and
every original station. Apply the existing 1e-4 / 1e-3 / 1e-2 etc numerical
budgets without posthoc loosening. Record both meshes' whole/inner/outer
physical-development envelopes, including failures. If needed, separately
refine quadrature or the inlet region from these results; no production pilot
is admitted before resolving a scientific or numerical reference failure.

Admit <=30 s generation, <=2200 s fine native execution (1800 s internal
steady deadline), <=900 s independent comparison and <=120 s compact
assessment. This bounded matrix is <=3250 aggregate operation seconds in the
continuing 43200 s account, with 3600 s protected for evidence. All prior
charges and failures carry forward. Use one science rank, six library thread
counts one, RSS<=8GiB and shared actual workers<=72. Do not execute Task3.
