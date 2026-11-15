# H19 reproduction and scope

The completed root is `/home/jay/build/wj-h19-target-timestep-20261115`.
Its 5400 s ledger is CLOSED. Commands below describe retained execution, not
permission to reopen the allocation or repeat a native reservation.

Native source capsule b2365e52 and scientific Python 4108e3da are separately
pinned; publication sources/late helper are in the separate source supplement.
Use the exact private GNU 11.4/Release runtime and loaded-library closure from
runtime-build/runtime-qualified/closure records. Do not substitute host libc.
Six numerical-library thread controls were 1, compile parallelism 8, one serial
native plus bounded observer; total scientific allocation stayed below 72.

## Bootstrap versus ordinary restart

`RANSWallResolutionProbe SESSION.xml INPUT.json OUTPUT.json` accepts the explicit
`step_zero_bootstrap` helper operation. The new XML has dt=.005 and otherwise
the exact scientific semantics of the manifest-selected fine original session.
The packet contains `source_dt=.0025`, the fine source XML pin, its original
report/checkpoint pins and an owned `checkpoint_base`. The source report is
strictly imported under the source's own signature before any new history exists.

`StepZeroIMEX1` requires IMEX1, absolute step/time zero and a distinct dt-bearing
signature. It creates a new checkpoint from preserved physical arrays and native
`Initialise`/`Capture`; it never transfers fine GLM slots. Its implicit callback
throws if asked to solve. `step_zero_bootstrap_checks` uses a small sealed input,
tests time/step/order and ordinary wrong-dt refusals, and writes no new fixture CP0.
The ordinary strict import/restart APIs are unchanged.

The helper explicitly audits its new CP0; it does not run the governor loop.
The restarted governor uses:

```xml
<I PROPERTY="RANSWallRetention" VALUE="BoundedV1" />
<I PROPERTY="RANSWallMonitoringPolicy" VALUE="CheckpointAuditV1" />
<I PROPERTY="RANSWallEvaluationReuse" VALUE="ExactStepV1" />
<I PROPERTY="RANSWallImmutableRowReuse" VALUE="Off" />
```

Its immutable intent uses recovery/scientific stride1, maximum absolute step2,
dt bits `3f747ae147ae147b`, the new coarse CP0 as `original_qh0_checkpoint`, and
the fine CP0 separately as `spatial_lineage_anchor`. Time stays0 in XML and
NumSteps2 is additional work from its own new CP0. TimeStep is never excluded
from the signature.

## Recorded command forms

Let R denote the closed H19 root and S its source Tests/RANS directory.
The bootstrap reserves exactly one slot:

```sh
/usr/bin/python3 -B "$S/run_wj_h19.py" --root "$R" \
  --label coarse-native-CP0-v1 --seconds 240 --phase qualification --kind bootstrap -- \
  /usr/bin/python3 -B "$S/study_wj_h19.py" --root "$R" --mode bootstrap
```

After identity, focused/native component gates and measured cost admission:

```sh
/usr/bin/python3 -B "$S/run_wj_h19.py" --root "$R" \
  --label coarse-target-v1 --seconds 750 --phase target --kind target \
  --attempts 2 --workers 2 -- \
  /usr/bin/python3 -B "$S/study_wj_h19.py" --root "$R" --mode target
```

The only protected final read used phase target_readback/seconds 120 and
`study_wj_h19.py --mode strict`, after native exit. No fine trajectory or new
fine strict reader was run. Exact argv, threads, worker allocations, input pins,
maps, peak memory and failed/unknown charges are retained in operation receipts.

The final focused Python invocation selected test_wj_h19_timestep,
test_wj_h18_checkpoint_audits, test_wj_h18_readback and the single existing
RepairedMonitor child-exit test:18 distinct methods. The prior misspelled selector
receipt remains failed. Native RANSWallMonitorChecks and the small history helper
component ran separately. CTest registration alone is not executed evidence.

## Comparison and failure semantics

H18 exact replay is unchanged. H19 compares coarse1/fine2 at T=.005 and
coarse2/fine4 at .01, with separate producers, signature/dt, callbacks/alpha and
serialized time bits. Its static cache compatibility explicitly proves spatial/
profile independence of dt; historical pins remain unchanged. State/closure/
profile/history-dependent results are recomputed.

Physical PASS, provisional screen status and continuation are separate. A
REVIEW_STOP ACK may commit the current physically qualified checkpoint, then
the governor throws a distinct TEMPORAL_REVIEW_STOP before another advance.
It does not relabel a numerical-review exceedance as a constitutive-domain exit.
Missing/invalid ACKs and native hard exits retain their original protocol/physics
failures. No reservation may be retried or substituted.

The realized common inventory includes both face owners, current native witnesses
and both states' moving-cutoff neighborhoods. All row identities and signed arrays
remain Linux in the hash-bound NPZ/inventory files. The compact publication retains
regional RMS/max screens, signed means/components, principal/direction witnesses,
all 780 Q1 vectors, interface/weak/residual evidence and explicit unavailable data.

Final native strict packets preserve projected tests and genuine completed-stage
history; alpha0 spatial norms are separate. No new24/32 quadrature campaign,
governing factorization, pressure regauge or model/mesh/boundary change belongs to
this comparison. H11 temporal adequacy and duration/application admission remain
unresolved. This provisional screen is not a production timestep selection.
