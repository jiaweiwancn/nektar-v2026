# CT16: separate bounded production-path engineering

Status: **IN PROGRESS; resource admission remains UNRESOLVED**. The user
explicitly authorized this separate phase after CT15 closed UNRESOLVED at
`5168ed5fba02a4c92a3847827be3a2dbb6ca071b`. CT15's report, ledger, archives and
failed attempts remain historical evidence. No CT16 history inherits their
component identity.

The [prospective contract](cao-tamura-accepted-16-evidence/contract.json) fixes
CT16's start at 2026-11-16 13:30:25 UTC (21:30:25 Asia/Shanghai). Its allocation
is 24,000 process CPU seconds, 28,800 elapsed wall seconds and 16 GiB of new
live data, with 16 GiB per-process address space. One numerical family/rank,
six numerical-library thread settings of one, two compiler jobs, at most 16
project workers and 72 aggregate Linux workers apply. The clock has not reset.
The external supervisor counts closed and live work even while the shared
worker guard suspends the owned family. Every failed attempt remains charged.
Heavy data and private runtimes are under
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-16-20261116`; shared installations,
other jobs and historical archives are untouched.

All CT8 physical and CT14 numerical criteria remain fixed. In particular the
manifest SHA256 is
`004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296`;
Original/Natural A95, 62,533 slots/961 elements, P4/P2, HA/128, actual current K,
transported inner k, dependent epsilon and the full signed stress/covariance/
radicand inventories remain. IMEX1 is one flow/epsilon/k sweep, with dt
0.00025/0.000125 and its measured post-k epsilon lag. VectorCOLAMDV1 and the
consistent free-mass pressure action remain. The Off oracle and ExactStepV1
open-route guard stay in place. No physical-preservation execution beyond the
already qualified short window has been launched.

| Gate | Current status and limit |
|---|---|
| G1 | V5–V10 small controls and actual selected first-two-step comparisons PASS; affected selected quadrature/alpha controls remain incomplete |
| G2 | CT15 V4 baseline pair PASS; fresh changed-component .002 pair UNATTEMPTED |
| G3 | Timing-off 0/1/2 native measurements available; complete workflow and frozen coarse6-to8 restart/retention heldout UNATTEMPTED |
| G4 | Two selected native retirement/recovery cycles UNATTEMPTED |
| G5 | UNADMITTED; all measured native-only entry extrapolations already exceed the complete CPU ceiling |
| G6 | Source/contracts and sealed progress published; final closure/accounting still pending |

Physical T=.407 remains **UNATTEMPTED**. Full-transit production and experimental
replication remain **UNQUALIFIED**. The passing short baseline, prospective
resource admission and longer-time physical qualification are separate claims.

The [sealed progress snapshot](cao-tamura-accepted-16-evidence/progress-v2/results.json)
binds completed work through V10 and the exact trace-index experiments. It
records 9,556.785661073 closed-family CPU seconds before its own collector;
subsequent failed/active attempts and publication costs still count. Its
[artifact inventory](cao-tamura-accepted-16-evidence/progress-v2/artifact-inventory.json)
hashes 3,707 payloads kept outside Git. The
[scientific receipts](cao-tamura-accepted-16-evidence/progress-v2/scientific-evidence.json)
bind actual runtime bytes, configurations, numerical histories and commands.
Only the deliberate paused-supervisor negative has an incomplete worker audit;
it is retained explicitly and never treated as clean numerical admission.

Timing-off native commands use the actual selected case, dt=.00025 and matching
scientific work. Every V5–V10 accepted vector in these experiments is binary64
identical to the corresponding CT14 vector; maximum coefficient difference is
zero. The nonnegative envelope is
`startup=max(0,C0)` and
`step=max(0,C1-C0,(C2-C0)/2,C2-C1)`. Detailed timers below only diagnose costs.

| Policy/component | C0 / C1 / C2 native CPU seconds | Repeated-step envelope |
|---|---:|---:|
| CT15 V4 | 56.598205 / 189.308833 / 350.148692 | 160.839859 |
| V5 StreamScalarV1 | 49.903598 / 157.136224 / 273.466898 | 116.330674 |
| V6 ScopedScatterV1 | 53.181892 / 154.906970 / 253.475760 | 101.725078 |
| V7 ScopedAuditPointsV1 | 46.089977 / 137.980507 / 234.472438 | 96.491931 |
| V8 ExactFlatMetricV1 | 44.271072 / 127.618765 / 212.064898 | 84.446133 |
| V9 ImmutableRowLRUV1 | 41.024357 / 116.387451 / 196.404393 | 80.016942 |
| V10 StreamGradientV1 | 34.413750 / 92.753797 / 148.738315 | 58.340047 |

V5 streams scalar FE values with the original double products and centred
long-double summation, and bounds the per-call exact constitutive point map at
2,097,152 entries. V6 removes repeated hash lookups while scattering the same
chronological matrix additions through a bounded, assembly-local pointer plan.
V7 gives the independent native audit its own current-state point calculations;
it never borrows domain validity. V8 specializes only the validated canonical
identity Flat map and retains all general-map fallbacks. V9 uses individual
least-recent eviction at the unchanged 1,048,576-row limit. V10 streams both
actual velocity gradients through bounded 3x128 scratch with the exact original
arithmetic. No policy caches a current coefficient, factor or constitutive
validity across states.

All four small operators, current flow/K/epsilon/time/data/alpha/mask/space
controls, Off updates, genuine continuation, full rejected-step rollback,
foreign histories and continuation with both row/point quotas reduced to 1024
pass. V8–V10 include 44 ordered dependency/geometry queries. Cold-row scalar,
flat-metric, LRU and gradient verifiers compare binary64 values, metadata and
signed zero where applicable. Each component has a new execution/history
identity; previous histories are rejected rather than relabelled.

The tested V10 native source is
`5b4e3cf40520cfb361eb0f8eb13b5bd6a561b88b`. Its private governing executable
SHA256 is `b1fdae6a08c245b9af894dcf5fadce4463d27c81420409670f81e721b2ee803f`;
the probe is `130718db7aa38627e616d64bc1d502bbee89d7887183ebe67c092167b35f15ec`.
The runtime is `build-gradient-v1/runtime` under the external root. Later
controller/publication commits do not imply rebuilt native executables.

The [future ledger](cao-tamura-accepted-16-evidence/progress-v2/forecasts.json)
reports fresh-zero and conditional own-.002 starts for .407 and 0.1/1/3 transits,
using decimal member-wise ceiling and explicit endpoints. V10's fresh entry
native-only extrapolation is 362,833.6509975 CPU seconds with the 25% headroom
applied once: 15.118 times the entire 24,000-second production ceiling. It
includes 153 length-32 segments and two separate zero initializations, but
omits independent audits, milestones, retention/replay, controllers and
publication. It is an unvalidated short extrapolation, not a mathematical lower
bound for all future states or implementations. Every horizon is UNADMITTED;
complete campaign CPU/wall and qualified retention peak remain unknown. The
zero-other-cost 3.931203931 CPU/step allowance falls to 2.839039466 after only
this measured startup model, before any other required workflow costs.

V10's broad one-step diagnostic assigns 18.336474 exclusive CPU seconds to
21 residual evaluations, 15.402599 to current constitutive inventories, 8.416919
to physical audits, 6.559479 to flow assembly and 5.501287 to repeated-point
native audits. It makes 21,870,516 current scalar evaluations and 7,538,384
paired gradient evaluations. These scopes identify further work; they are not
acceptance speedups. Basis construction itself is only .663 CPU seconds, so a
basis-quota sweep is not selected.

The exact `CT16-ExactTraceIndexV1` independent-audit experiment preserves every
complete observation and decision for the sealed actual initial/coarse/fine
states. Full bundle CPU changes from 16.325126/16.964765/16.221196 to
13.737221/13.807352/13.712643. Its 187-entry index uses 19,132 bytes, holds only
immutable first-match H-plus addresses and preserves signed-zero membership.
Missing/nonfinite queries, foreign packets, stale geometry and failed-scope
cleanup pass their controls. Actual Python and mapped-library bytes are bound
in the receipts. This still exceeds the entire per-step CPU allowance even if
native evolution were free.

The [first bounded batch experiment](cao-tamura-accepted-16-evidence/independent-batch-contract-v1.json)
failed the unchanged gradient-extremum geometry criterion. Its cost and source
remain recorded; it is not adopted. The
[narrowed experiment](cao-tamura-accepted-16-evidence/independent-batch-contract-v2.json)
keeps gradient witnesses on original PreparedRows arithmetic and batches only
physical-value norms. The
[within-audit residual candidate](cao-tamura-accepted-16-evidence/attempt-7-contract.json)
reuses one caller-owned raw residual only for identical immediate audit
consumers, with fresh original-vector bitwise controls. Its first build failed
a namespace qualification before solver execution (96.976697 CPU seconds);
the correction is pushed and the failure retained. The
[coverage candidate](cao-tamura-accepted-16-evidence/attempt-8-contract.json)
changes only exact membership representation, with a bounded hash and complete
ordered overflow fallback. Both require their own controls and measurements.

Earlier retained failures include the legacy probe-name guard and the original
cache-quota history mismatch; the fixes preserve component identity while
keeping disposable capacities as resource metadata. The
[packed contraction prototype](cao-tamura-accepted-16-evidence/primitive-contraction-v1/results.json)
is bitwise exact for 37,914 rows but shows no consistent benefit and is not
adopted. The
[flat-metric primitive](cao-tamura-accepted-16-evidence/primitive-flat_metric-v1/results.json)
is bitwise exact at 28,288 actual coordinates and supported the subsequently
measured V8 specialization. Failed and insufficient attempts remain in CT16's
single ledger. Work continues within the original finite allocation; none of
these local results closes G1–G6 or authorizes Task 3 physical execution.
