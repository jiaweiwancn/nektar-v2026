# M8-R2 initial representation preflight

Outer6 passes the declared **finite-set initial coverage** checks for both
inner laws and both source quadratures. Outer9 is not needed. The mesh has
2165 unknowns, retains the six R1 inner cells, p3/p1 and native GLL8, and
changes only the outer three intervals to six. All four actual initial
coefficient vectors are bitwise identical, SHA256
`f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca`.
This is an admissible initialization result, not stationarity, an implicit
stage root, positivity throughout polynomial elements, or trajectory stability.

The [compact summary](preflight-002-summary.json) binds each XML/report/audit,
the producer `9f20b3053c303e8805f5183f668a8ec468d9656b`, corrected independent
analysis `dd7d7cf2b71b532b56288aa9d525ae5955602698`, and the before/after
runtime fingerprints. All producer C/C++ hashes are identical in the analysis
checkout. Source, runtime, tools, inputs and both clean revisions remained
unchanged. Root ran CPU49 with all six numerical thread settings equal to one;
the read-only curator used CPU50 and reconstructed no fields.

| Case | Source orders tangential/normal | Independent checks | Inner covariance minimum | Saturation source measure, of 12 | Native startup s | Domain s | Residual/constraints/norms s |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Linear base | 8/16 | 823 | .0772933168271 | 0 | 10.015, retained execution | .190670 | .131344 |
| Linear high | 12/24 | 943 | .0772933168271 | 0 | 23.029 | .305670 | .241293 |
| SaturatedR2 base | 8/16 | 827 | .110802957528 | .0421042732205 | 10.015 | .191877 | .130931 |
| SaturatedR2 high | 12/24 | 947 | .110802957528 | .0431203736176 | 23.029 | .305189 | .243464 |

All 3540 independent checks pass. The full finite outer covariance minimum is
`.176807110004236`; the projected dependent trace has the same positive minimum.
The worst normalized gradient-entry departure is `.000500194924459` in all
four audits, leaving margin `.000499805075541` to the unchanged `.001` bound.
Its witness is the first outer-cell upper endpoint, r=7/3; both sides of the
old r=8/3 witness are also retained. R1's maximum was `.003412338694780`.
Local and sampled K/K0 remain one to roundoff. The stricter sampled-K
implementation guard has margin `.0100`. Actual outer epsilon/epsilon0 is
within `[.983252473144720,1.000005401761861]`; the projected matching trace
is within `[.983252473144721,1.000000000000001]`. The limiting epsilon-box
margin is `.003252473144720`, with relative trace coefficient residual
`5.32e-15` or less. The independent inventory includes every exported actual
scalar observation, freshly constructed volume/baseline-trace/current-trace
rays and the extra constitutive diagnostic rays.

The global finite inner a/c maximum is `1.09303567444`, above rho=.8.
Thus the projected field activates the new law despite the analytic neutral
reference's strict inactive bound. Its selected source production and mean
stress work both change by `-.0123957428841` at base quadrature and
`-.0123918162236` at high quadrature; the exchange defect is zero.
The capped source remains on its plateau: lambda=1, limiter/cap active measures
zero, and E=D=`8.43616705170709` (base) or `8.43625113409570` (high).
The redistribution integral is below `4.7e-17`. No cap-induced epsilon
attenuation is observed at initialization.

Quadrature sensitivity is recorded at identical coefficients, without extra
native solves. For Linear, increasing normal16 to24 changes complete weak
loads only by roundoff; changing tangential8 to12 gives a k mass-dual load
difference `.00834036`. For SaturatedR2, normal16 to24 gives streamwise/k
difference norms `.000270838`/`.0000558726`; changing both orders gives a k
difference `.00833786`. The active-set indicator measure varies with the
quadrature and is reported at each rule. These measurements are not continuum
quadrature-error bounds. The largest native/independent implementation error
is below `1.5e-14` in the reported comparison scales.

The initial raw unanchored norms remain substantial. At base quadrature the
streamwise/k norms are `9.86324`/`62.59848` for Linear and `8.22728`/`62.60960`
for SaturatedR2; the outer-epsilon norm is `.215669` for both. No anchor or
field replacement makes these vanish. All common physical transition anchors
retain exactly the R1 chi values. The epsilon-sensor departure drops from
maximum `.00413633327` to zero because the refined first outer cell no longer
supports the H-trace basis at the fixed r=2.5 sensor; the state gate remains
one, so this representation effect does not change the observed initial chi.
Sensor definitions, transition, roughness, H and model coefficients are fixed.

The first independent audit stopped on an outer-H element-side lookup defect;
the [stop receipt](preflight-001-stop.json) and all original outputs are kept.
The corrected verifier explicitly selects outer side+1. The completed
Linear/base native initialization was reused after byte-identical XML,
binary/source and before/after runtime checks; it was not rerun. Completion
made exactly three new native initializations and independently audited all
four reports. No implicit stage, Jacobian probe or evolution was executed by
this preflight. Fresh-J cost is measured only inside the bounded stage family;
the current 10–20 minute planning estimate is not a runtime conclusion.
M8 remains open.
