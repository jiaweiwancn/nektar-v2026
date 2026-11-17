# Native terminal-output and reduction checkpoint

This checkpoint follows `f740934c4bbe14b4e1c51b7011a0822a15129dc3`.
The approved `FiniteNuIsotropicWallFullWJCollarV2` equations, reference,
retained inlet, boundaries, BDF2/EX2, actual production and current acceptance
guards are unchanged. V4 remains inactive. `TASK3_READY=NO`; H27-C/D and the
full coupled production-cost gate remain open. No full Task3 run was launched.

## Implementation and pinned execution

`WJTerminalOutput=True`, `WJTerminalBegin` and `WJTerminalEnd` opt into native
scientific output at both accepted-time brackets of a nonaligned endpoint,
or the exact accepted state of an aligned endpoint. Periodic output remains
available. A frame is written once per owned step, with merged reason bits.
The solver never moves the requested interval or invents an accepted state
at an interpolated time. A create-only `.science-index.json` binds source,
reference, full signature, clocks, complete/failed run status, file bytes/SHA,
reasons, write times and observation provenance. Fresh initial observations,
imported histories and steps accepted in this process are separate origins.

Linux namespace: `/home/jay/build/wj-h27-coupling-r1-20261117`.
Source14 identity:
`d89a803bf15a60bddfe6a0d3fe60f44e07e35d7a4c209b261c735388b90b537f`.
Runtime ELF SHA256:
`c58e4899f0e7880256f806944ea7612f1e3630cca0b93e0f38df875caafa259e`.
The reference-table SHA256 is unchanged:
`b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d`.
All1717 locally present C++/header/CMake files in the executed capsule match
this Mac checkpoint. The standard-main Release build took205.200580s with
four compiler jobs and at most15 actual shared-host workers. Science used
one serial process/rank and all six numerical-library thread counts1.

The compiled capsule precedes the final Python checker/fixture overlay and
the prospective cadence refinement. Both immutable validation overlays and
admissions are retained separately. Final checker SHA256:
`64471e899287c1659f893d0de23800e28be404e1764044790c4fb8d2d5f89e1e`.
Exact commands, XML/option files, runtime libraries and receipts accompany
the [evidence inventory](wj-h27-linux/native-terminal-output-v2f/inventory.json).

## Actual evolving output experiment

Use the retained14x50=700-cell Q3/Q2 physical case to T=.06, dt=.005,
CoupledJFNKV1, OperatorDriftV1/.02/eight bases, nonlinear3e-12 and linear1e-12.
The test terminal interval is exactly[.0173,.0537], measure.0364. This short
interval tests output handling; it does not replace the physical-duration,
integration-accuracy or original Task3 window requirements.

The stride3 first attempt retains steps3/4/6/9/10/11/12. It passes native
bracketing and every cold-segment statistic, but misses the inner-stress mean
sampled MAX: .0013810776285185788 against the fixed .001 stress-over-tau limit.
Its inner-R RMS is .00031979532367094226. All raw inputs, failed reports,
reduced arrays and the separate numeric failure detail remain on Linux.

A prospectively admitted stride2 refinement retains2/3/4/6/8/10/11/12.
The stride1 trajectory supplies all12 accepted states. Common-clock current,
previous and both explicit-history hashes are exactly equal. A separate prefix
advances6 steps and a cold process advances6 more. Its checkpoint links all
four complete histories, and the first cold accepted step binds that import.
The window integrates .0127 time units from the prefix and .0237 from the
continuation. Prefix owns through6 and continuation owns new accepted states
after6; there is no conflicting or silently selected future branch.

## Independent reduction

The Linux checker reconstructs original signed Q3/Q2 owners and full 3D WJ
stress without a native basis/stress call. It preserves each owner gradient
at coincident coordinates. Its Gauss5 tensor samples plus zero-weight exact
endpoints give spatial RMS and finite sampled MAX; they do not prove a
continuum supremum or spatial quadrature qualification. Include U/W, k, Pi,
full R, full gradU, outer epsilon, inner Ew, interface K/Gamma/cutoff and
ground traction. Ground traction is `(V_xz,V_zz-Pi)` on the upward normal.
Pressure retains its physical open-traction level. All scales are fixed in
the [admission](wj-h27-terminal-output-admission-v1.md).

Interpolate each observation linearly between its actual frames, clip each
segment to the unchanged requested interval, and integrate the observation
and its square exactly. An independent three-point time-Gauss calculation
checks both moments with a declared floating-arithmetic reserve. Centered
variance is integrated directly as a nonnegative quadratic; no negative
variance floor or cancellation subtraction is used. Incomplete coverage is
rejected. Exported double arrays round-trip exactly through compressed NPZ.
They preserve sampled statistics and their owners/scales, not full state fields.

Both means and standard deviations pass the unchanged numerical budgets.
The table reports the most demanding RMS/MAX across those statistics and
the applicable inner/outer/interface/ground groups; full rows remain in the
independent report.

| Observable | Scaled RMS | Sampled MAX | Limit |
| --- | ---: | ---: | ---: |
| U/W | 8.79208e-8 | 6.75394e-7 | 1e-4 |
| k | 8.83343e-5 | 3.93138e-4 | 1e-3 |
| Pi | 2.93350e-6 | 1.38213e-5 | 1e-3 |
| Full R | 1.99455e-4 | 8.87647e-4 | 1e-3 |
| gradU | 2.14013e-5 | 2.04377e-4 | 1e-2 |
| Outer epsilon | 3.18966e-9 | 1.11112e-7 | 1e-2 |
| Inner Ew | 9.80696e-8 | 1.27606e-7 | 1e-2 |
| Interface K | 6.23718e-8 | 6.99025e-8 | 1e-3 |
| Interface Gamma | 9.35902e-8 | 1.04894e-7 | 1e-3 |
| Interface cutoff | 3.98080e-9 | 4.46089e-9 | 1e-3 |
| Ground traction | 5.51116e-7 | 4.99922e-6 | 1e-3 |

All cold-segment/uninterrupted means and standard deviations meet1e-8.
The largest Pi difference is7.68845e-11 and traction3.66921e-11; the other
observations are closer. The independent new cold endpoint residuals are
momentum9.25497e-15, continuity6.60855e-16, k1.25910e-13 and
epsilon1.16821e-14. Its current-equation/history-boundary/quadrature subgates
are `MET`; its overall physical report remains `NOT_MET`.

Aligned [.02,.06] and zero-begin [0,.06] PolynomialScalarsV1/body1 MMS runs
also pass their independent scheduling checks and have identical endpoint
histories. A real one-Newton-step failure leaves step0 intact, publishes an
incomplete index, and labels its step0 frame as a fresh observation. The
failed native receipt and complete preserved-history journal remain.

## Timings and retained data

These costs include actual native output. Every physical final step still
uses four Newton updates and seven Krylov iterations. Per-step variation,
rebuild/reuse events and disjoint categories remain in the assessment.

| Physical run | Whole operation, s | Complete advancement, s | Scientific frames, s |
| --- | ---: | ---: | ---: |
| Stride2 scheduled12 | 98.680523 | 92.995783 | 1.046206 |
| Stride1 dense12 | 98.088882 | 93.371302 | 1.571532 |
| Stride2 prefix6 | 54.081512 | 48.343401 | .522910 |
| Stride2 cold6 | 52.254160 | 45.687749 | .525912 |

Frame times are subsets of the reported advancing times. Index writes take
about.00006–.00010s within final export. Whole-operation variation is retained;
this table is not a production-speed or qualified-dt forecast.

The scheduled sampled-statistics NPZ is8,485,993bytes. Full native frames and
all three reduced artifacts remain in Linux run storage. The larger array
payload stays outside the source checkout; `mac-source-binding.json` records
its Linux identity and the Mac review copy at
`/private/tmp/wj-h27-terminal-scheduled-statistics-v14.npz`. The full immutable
archive remains at `/private/tmp/wj-h27-native-terminal-output-v2f.tar.gz` on
Mac and `evidence/native-terminal-output-v2f.tar.gz` in the Linux namespace.
It is25,792,665bytes, SHA256
`c975b9742374beb37872d75e4f63913ee69b445bba12e1506a86398d2841b4c3`.

Publication took8.841094s. The same campaign account is OPEN at
10630.33897288586 aggregate operation seconds, including both failed controls
and the failed cadence3 assessment. No allowance was reset;43200s total and
3600s protected evidence remain the governing limits.

## Remaining work

This qualifies a short terminal-output/reduction foundation. Raw scientific
frames have deliberately been retained; a bounded production reduction/retention
policy and its storage forecast are still required. Zero-`NumSteps` terminal
output also needs its frame serialization moved under the export timer before
that timing label is qualified; the measured advancing cases above are unaffected,
and whole-operation accounting already includes all work. The next numerical
checkpoint will repair and test that path.

Source15 now repairs and verifies that zero-step timing path; see the
[current-data and cold-restart checkpoint](wj-h27-current-data-cache-checkpoint-v2.md).

The unresolved v2 same-inlet reference, evolved physical pilots, h/p/q/dt
accuracy, monitoring escalation and complete per-trajectory/campaign cost
remain open. V4 awaits its separate scientific decision. The original Task3
horizon505.02466521278404 and terminal interval
[447.8483832973273,505.02466521278404] remain unchanged; this checkpoint supplies
no qualified full-Task3 launch command.
