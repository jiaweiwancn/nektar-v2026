# Spatial MMS inner linear tolerance control

Source-v11 passes the complete refined temporal study, including its tighter
3e-13 linear / 3e-12 coupled run. U orders are 2.11355/2.05316, k
2.01494/2.01092, epsilon 2.01863/2.01267, and pressure 1.98973/1.98865.
Mesh/quadrature/tolerance separation and owner-gradient budget convergence
all pass on the finer 32-cell temporal base.

The separate Historical spatial MMS at dt=0.0025 reaches a scalar linear
roundoff plateau: 1.60704e-12 after 400 iterations versus its requested 1e-12.
The complete accepted prior step had residual 7.20e-12. The canonical stopping
repair correctly rejects that unmet inner target; no failed state is accepted.
Preserve the failed 8/32-cell study attempts and their accepted histories.

Repeat all four spatial variants freshly on the SAME source-v11 executable
with inner linear tolerance 1e-11 and the SAME coupled tolerance 3e-11. Keep
T=0.02, dt=0.0025, nx=nr=2/4/8 and the fine dt=0.00125 control, quadrature,
manufactured equations and independent acceptance limits unchanged. This is
an internal linear-solver setting, not a relaxation of any H27 current-equation
or physical-development gate. Every accepted endpoint still needs the full
current fixed-scale residual, interface and conservation checks. The complete
temporal tolerance study above and the spatial matched-time contamination
control remain required; unresolved rates are not declared successful.

Continue the same 1200-second MMS tranche with 444.1090821 s already charged,
not a fresh allowance. Admit each native variant <=120 s sequentially and the
independent analysis <=400 s only while the cumulative remaining reservation
permits it. Retain all failures. One rank, numerical-library threads one,
shared actual workers<=72 and existing RSS/storage limits apply.
