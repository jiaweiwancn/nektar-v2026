# Frozen steady Historical MMS spatial study

This prospective test-only study uses the unchanged source-v11 registered
solver and approved v2 model at exactly t=0, alpha=0, with current advection.
The native steady MMS route and independent `core.Oracle.forcing(...,
steady=True)` remove the analytic time derivatives from the continuous
manufactured loads. No timestep is advanced and no fictitious BE/BDF2 or
EX1/EX2 history is constructed. The advancing spatial attempts remain
evidence: their reported matched-time control differences contaminate k by
58% and U by 33%, so they cannot establish spatial rates. This study removes
that source of error analytically; it does not replace temporal qualification.

Use the following fixed matrix, with Q3 velocity/k/epsilon and Q2 pressure.
The x domain is [0,1]; each of [0,H] and [H,.5] has nr equal intervals.

| Role | Label | nx=nr | Cells | Native qx x qr | Coupled nonlinear tolerance |
| --- | --- | ---: | ---: | --- | ---: |
| h | h2 | 2 | 8 | 12 x 18 | 3e-12 |
| h | h4 | 4 | 32 | 12 x 18 | 3e-12 |
| h | h8 | 8 | 128 | 12 x 18 | 3e-12 |
| quadrature | caller label | 8 | 128 | 16 x 24 | 3e-12 |
| tolerance | caller label | 8 | 128 | 12 x 18 | 3e-13 |

Keep the actual linear tolerance identical across h and quadrature runs;
the tolerance control may tighten it. Match the constitutive backend, stress
split, steady preconditioner, geometry, body amplitude, historical MMS
variant, source and reference identities. A stored dt is only an identity
field here and must match; it never enters the independently evaluated
steady equations. This document admits no numerical operations by itself;
the parent forecasts and admits native runs and analysis in the existing
H27 ledger, preserving all previous charges and failures.

`h27/steady_mms_spatial.py` reads five completed normal-main steady artifacts,
requires native convergence, zero physical time and genuine step-zero state,
and checks current/previous state and advection hashes. The manifest is:

```json
{
  "schema": "WJH27SteadyMMSSpatialManifestV1",
  "source_identity": "<exact source-v11 SHA256>",
  "runs": [
    {"label": "h2", "native": "h2.result.json", "role": "h"},
    {"label": "h4", "native": "h4.result.json", "role": "h"},
    {"label": "h8", "native": "h8.result.json", "role": "h"},
    {"label": "h8q", "native": "h8q.result.json", "role": "quadrature"},
    {"label": "h8tol", "native": "h8tol.result.json", "role": "tolerance"}
  ]
}
```

Paths are relative to the manifest directory or absolute. Comparisons are
predefined as h2/h4, h4/h8, h8/quadrature and h8/tolerance. An optional
`comparisons` array must name exactly those four pairs in that order, with
`kind`, `coarse`, `fine` keys; it cannot omit a control or change a gate.

Analyze each endpoint independently at 20x24 and 24x32 split Gauss rules.
Use the existing signed Jacobi reconstruction, continuous SymPy fields,
independent full WJ root/tensor and physical owner gradients. Require both
current equations to meet the original momentum/continuity/gauge 1e-9 and
k/epsilon 1e-8 limits. Their quadrature difference must be <=20% of those
limits. Require actual essential/dependent-Gamma trace and cached current
advection agreement <=1e-10. Verify scalar constant-test and mechanical
velocity-test weak budgets against independently integrated current terms,
retaining essential reactions; identity errors must be <=1e-10 in fixed
reference units. Record both owners' physical interface-flux jumps separately.
Require positive native scalar/principal margins, nonnegative covariance,
and independently sampled positive scalars/nonnegative covariance. Native
principal coverage is identified as such; no independent continuum principal
certificate is claimed by this study.

Report scaled physical L2/RMS and gradient H1-seminorm/RMS errors, sampled MAX
and locations, both independent quadratures, and common owner-aware field
differences. Smooth Q3 scalar expectations are L2 order4 and H1 order3.
The prospective fine-pair h4/h8 requirements are U/k/epsilon L2 order>=3.5
and gradU/gradk/gradepsilon order>=2.5 when both errors exceed 1e-12. If the
fine error is within that fixed band, mark the rate unresolvable and make no
rate claim; a coarse-in-band/fine-out-of-band reversal fails. The h2/h4 rate
is diagnostic. Pressure has expected Q2 orders3/2 and remains diagnostic;
do not shift its gauge to improve the error.

Each quadrature and tolerance control's common-point RMS difference must be
<=max(1e-12,10% of the h8 RMS error), for U/k/epsilon and their gradients.
The independent error and difference integrals must agree under quadrature
refinement within max(1e-12,0.1% of the finer integral norm). MAX sampling
and physical owner-flux jumps are diagnostics, not substituted for weak
identities or claimed pointwise continuity.

Run only on Linux under the admitted resource runner:

```
python3 steady_mms_spatial.py --manifest MANIFEST.json --table TABLE.json \
  --out NEW-REPORT.json --seconds ADMITTED_ANALYSIS_CAP
```

The checker creates immutable output, preserves partial failures, uses one
numerical thread and caps address space at8GiB. PASS means only this frozen
steady MMS spatial study passed. It does not qualify advancing histories,
time accuracy, the same-inlet physical reference, Task3 readiness or all H27.
