# Fixed-state residual resolution and scalar contact

The native probe and Python audit are diagnostics on saved coefficients. They
do not call `Solve`, advance a time history, alter a source law, or change the
original `2e-12` coupled root criterion. The parent orchestration owns Linux
execution, source/binary fingerprints, worker allocation and fresh run paths.
All six numerical thread settings are one; use one allocated worker.

## Native and independent commands

For each of the exact saved Linear and SaturatedR2 auxiliary packets, retain
the original h0/history, `time=.0025`, `alpha=3200`, physical session and source
quadrature. The native probe accepts:

```text
RANSWallResolutionProbe SESSION.xml INPUT_PACKET.json FRESH_NATIVE.json
```

`require_history_h0` defaults to true and checks rebuilt native h0 bitwise.
The report exports separate physical terms, the baseline-order reconstruction,
compensated sums, long-double products/sums as decimal strings, actual native
residuals, constraints and full scalar evaluation inventories.

The independent audit is read-only and refuses Mac numerical execution,
tracked changes in its analysis source, and overwriting an existing result:

```text
taskset -c 50 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r3_resolution.py --report FRESH_NATIVE.json --session SESSION.xml --packet INPUT_PACKET.json --source SOURCE --producing-revision NATIVE_SOURCE_SHA --native-receipt NATIVE_RUNNER.json --output FRESH_ANALYSIS.json
```

The producing native revision and clean analysis revision are separate epochs.
The oracle checks saved state/history identity, target time/alpha, indexed
stress and production, independent source-panel/Legendre/Jacobi reconstruction,
each physical term, complete weak sum, actual projected constraint targets,
and actual native residual. Earlier failed candidates remain failed candidates.

## Compact extraction

After the two auxiliary and one declared contact analyses, the following
read-only extractor binds native-report, oracle-report and coefficient hashes.
It imports no closure, quadrature, oracle or solver. Use distinct fresh output
paths; add further retained positive boundary-path diagnostics only when the
parent has already authorized them.

```text
taskset -c 50 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 SOURCE/solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r3_resolution_summary.py --case Linear3200 LINEAR_NATIVE.json LINEAR_ANALYSIS.json --case SaturatedR2_3200 SAT_NATIVE.json SAT_ANALYSIS.json --case ConstructedContact CONTACT_NATIVE.json CONTACT_ANALYSIS.json --output FRESH_SUMMARY.json --markdown FRESH_TABLES.md
```

The JSON retains complete term-by-block norms and compact source/hash pointers.
The Markdown contains the main residual table, per-block and per-term tables,
trace mismatch, and conditional physical contact-rate tables. Original interval
endpoints remain in JSON; displayed intervals round outward to ten significant
digits. Failed independent analyses remain explicitly failed and supply no
promoted independent gate.

## Interpretation of the arithmetic tables

The implementation already subtracts state minus history before applying the
mass matrix. Its mass dot products already use long-double accumulation. The
audit measures those existing paths rather than attributing an absent defect.
The compensated comparison uses the actual exported compensated component
totals, including the native final binary64 additions. Extended totals retain
native rounded physical primitives and stored mass entries; they do not
recompute all constitutive inputs at arbitrary precision.

The summation enclosure is
`|fl(sum x_i)-sum x_i| <= gamma_n sum |x_i|`, with `u=2^-53` and
`gamma_n=n*u/(1-n*u)`, for frozen binary64 increments. It excludes basis,
geometry, profiles, projections, constitutive products and mass-entry errors.
The rational mass certificate instead fixes the declared rounded residual
coefficients and exact polynomial physical mass, and uses exact rational
inverse and solve residuals to bound the quadratic form. It does not enclose
errors in those residual inputs or native stored-mass entries. Neither bound
alone establishes a nonlinear stagnation floor or its cause.

Free physical term norms exclude constrained rows. Strong-row physical term
loads are boundary reactions, retained separately. Only the fully constrained
residual has a constraint-residual norm. Native generic `term_*` norm arrays in
the raw report must not be interpreted as physical constraint defects.

Outer isotropic cancellation is an alternative algebraic evaluation of the
same reduced-pressure stress and production. The full production is `-R:G`,
including `-(2k/3) div U`; no pointwise nonnegative-production assertion follows.
Net inner destruction appears exactly once as `+D` in the k residual. The
separate `E` and `-C` entries are auxiliary checks, not extra physical sources.

## Contact and trace scope

`compute_scalar_rate=true` requests a consistent free/boundary mass solve on
the supplied saved state. The single predeclared construction instead sets
`construct_contact=true`, defines its polynomial from original h0, and lifts
the dependent EPS trace without replacing any k coefficient:

```text
k=K0*((1-f)+f*g)
g=(r/.01)^2 below the existing first inner interface .01; g=1 above
f=1 on [-3,1], with the declared quadratic ramps to0 at x=-4 and2
```

The exact `0/+-K0` coefficient identity proves k nonnegative throughout these
polynomial elements. The probe and oracle separately check native scalar
inventories, dependent-trace consistency and finite covariance/outer-box/lambda
conditions. A successful reconstruction audit does not imply that every
restricted physical gate passed. There is no shape or model search if this
one diagnostic fails.

Always report the original fixed GLL location `(1.6046496089512394,0,0)` and the
constructed ground zero `(1,0,0)` separately. The original location lies outside
the constructed zero plateau and is positive there. A rate at that positive
point is not a zero-contact inward-invariance test. The historical terminal
near `8e-15`, with unresolved sign and nonzero trace mismatch, is not a certified
positivity counterexample.

The physical k rate is obtained from
`M_FF kdot_F=-L_k,F-M_FB kdot_B`, followed by the actual polynomial basis-row
evaluation. No mass lumping or load-coefficient sign substitutes for that rate.
U and stored reduced Pi are fixed in this scalar-subsystem diagnostic, while
mechanical `p/rho=Pi-2k/3` varies with k in forcing and sensors.

Full scalar-subsystem solvability is assessed only when all actual baseline
and current trace projection nodes verify unit state/admissibility gates,
strict `b>delta`, `lambda=1` and H above both profile cutoffs. On that verified
branch the separate baseline/current projected derivative is
`1.5*Epre(H)/K*sampledKdot`; designated downstream correction modes remain zero,
`g_E=0`, and the induced EPS boundary mass is included in the second SPD mass
solve. The oracle independently reconstructs this derivative and projection.
Otherwise the result remains an explicitly conditional k-only rate; it does
not assert a complete constrained k/EPS evolution. Lambda approaching zero
remains outside the frozen positive-EPS restricted contract.

A negative rate at a verified zero could support a scoped discrete
inward-invariance diagnosis only with its arithmetic and admissibility
qualifications. It cannot prove absence of a backward implicit root,
continuum instability, or a trajectory positivity failure. A nonnegative rate
in this one construction is not a positivity theorem either.
