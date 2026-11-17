# H27 persisted scientific statistics and actual owned retention

**Bounded foundation passed; positive H27 closure remains open.
TASK3_READY=NO.** This checkpoint uses the approved v2 equations and the
source15 executable from `9c843a8e0b6543b46832cd9ee7963879e650eefd`.
It adds Python analysis and retention tools. No C++ change or new build was
made. The proposed v4 scalar/reference amendment remains inactive.

## Executed recipe and identities

Authoritative Linux root:
`/home/jay/build/wj-h27-coupling-r1-20261117`.
New experiment: `science/scientific-retention-v16`.
Immutable helper overlays: `tools-scientific-retention-v16` and
`tools-scientific-retention-v16b`.

| Item | Exact identity / configuration |
| --- | --- |
| Active model | `FiniteNuIsotropicWallFullWJCollarV2` |
| Executed native source | `47fb3597c4e56348361cfab6eadf9a34ba8d89027e57ed01153ca4d9c5ab10d1` |
| Executed ELF SHA256 | `be3dd8fa87236e3c25a5e6b890077090587dd6f85fa83cc309a31a46415b1a05` |
| Reference table SHA256 | `b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d` |
| Mesh / spaces | Retained700-cell14x50; Modified_A Q3 velocity/k/outer-epsilon, Q2 pressure |
| Assembly quadrature | x8/r12 with existing model/cutoff splitting |
| Time / solver | dt=.005, BE startup then BDF2/EX2, CoupledJFNKV1, no Anderson |
| Current tolerances | nonlinear3e-12, linear1e-12 |
| Inverse cache | CurrentDataDriftV1, relative threshold.02, maximum8 bases, full events |
| Scientific observations | Every accepted step, full coefficients, terminal test window[.0173,.0537] |
| Recovery | Global accepted steps3/6/9/12, across all three native epochs |
| Execution | One serial science process/rank, all six numerical-library thread limits1 |

The active tensor is `Rwall=(1-chi)Riso+chi*RW`, where
`Riso=2kI/3+bA*n^T+n*bA^T` and `bA=-nuA*(I-n*n^T)*M*n`.
The geometry-only C2 collar reaches full WJ below H. Physical viscosity,
full3D DevStrainRawWorkV1, raw production, independent Cd, inner
`Dk=nu+Cd*k^2/(sigmaK*Ew)`, actual `P-Ew`, the v2 counterfactual cutoff,
dependent outer trace A(K)/H, original inlet/roughness and all current guards
are unchanged. No floor, clipping, EP or equilibrium forcing was added.

Mac verification binds all1717 locally present C++/header/CMake inputs to the
executed source15 capsule, and the three new helpers byte-for-byte to v16b.
Mac ran no numerical validation. The binding is
[mac-source-and-transport-binding.json](wj-h27-linux/native-scientific-retention-v2h/mac-source-and-transport-binding.json).

## Persisted moments and independent comparison

The [admission](wj-h27-scientific-retention-admission-v1.md) freezes the
linear-segment, positive-weight centered-moment formula, source/owner/clock
checks and2e-12 streaming versus independent RMS/MAX limits. Full owner and
scale inventories remain distinct at coincident coordinates. Observations
are the independent source14 signed Gauss5-plus-endpoint inventory:16 groups
cover inner/outer U,k,Pi,R,G, outer epsilon/inner Ew, interface K/Gamma/cutoff
and ground total traction. Sampled MAX is not a continuum bound.

An uninterrupted12-step physical reference is compared with three separate
four-step native processes. Epoch2 imports the full step4 binary; epoch3
imports step8. Each V1 statistical generation was committed before launching
the next native epoch. V2 removes cumulative-list copying: each JSON stores
only that batch's frames/segments, cumulative counts and the exact previous
generation identity. The same native data were reprocessed in three cold V2
statistics processes without rerunning native trajectories. All three V1/V2
NPZ payloads are byte-identical, so their persisted moments and last
observations are exactly equal. Both the ordered V1 experiment and immutable
V2 replay are retained in the evidence.

| Generation | Last step / time | Covered terminal measure | Complete window |
| --- | --- | ---: | --- |
| 1 | 4 / .02 | .002700000000000001 | No |
| 2 | 8 / .04 | .0227 | No |
| 3 | 12 / .06 | .0364 | Yes |

Each NPZ is written create-only, reloaded and compared exactly. A separately
computed two-pass centered variance also checks analytic first/second moments
against three-point time Gauss integration. All mean/standard-deviation
RMS/MAX comparisons pass: streaming versus the same split trajectory has
largest MAX1.7763568394002505e-15 (limit2e-12); split versus uninterrupted has
largest MAX3.1809803183089075e-11 (limit1e-8, outer pressure standard deviation).
The stale step4 generation supplied to the step8-to12 batch is rejected as
not the immediate predecessor, without publishing any output generation.

All six endpoint history roles pass fixed1e-8 coefficient, two-resolution
original-owner field and dependent K/Gamma/cutoff/Ew comparisons. The largest
coefficient difference is previous pressure1.5674735406497862e-9 in tau0
units. No pressure alignment is used. The separate full current-equation
oracle passes momentum1.0014471150431699e-14, weak continuity
6.459298532511859e-16, k1.2170955770430932e-13 and epsilon
1.165117609507495e-14, with current history, boundary, accounting and refined
quadrature subgates MET. Its overall physical status remains NOT_MET.

## Actual archive, deletion and cold recovery

After the independent checks pass, a protected qualification-only tar archive
is written and every member hash is checked. A durable verified-reduction
record and explicit owned deletion plan precede actual unlink operations.
No path outside the new experiment is owned. Original native scientific
indices remain unchanged and continue to describe pre-retention writes;
the separate receipt records which payloads are now absent.

| Payload | Written | Retained | Actually removed |
| --- | --- | --- | --- |
| Full split scientific frames | 1 through12 | 11,12 | 1 through10 |
| Scheduled recovery binaries | 3,6,9,12 | 9,12 | 3,6 |
| Native final binaries | 4,8,12 | 8,12 | 4 |
| V2 statistics NPZ | generations1,2,3 | 2,3 | 1 |

The14 removed payloads total56,672,766bytes. The8 retained selected payloads
total29,898,925bytes. Qualification-only evidence in the protected-retention
directory totals32,083,570bytes at assessment time, counted separately from
working data. It is not credited as a production storage saving. Earlier
experiments and the V1 prototype payloads remain intact.

All four scheduled binary writes independently match the four complete
current/previous/explicit history hashes of their accepted scientific frames.
The removed3/6 bytes are audited directly from the protected archive;9/12
are audited from retained files. After rotation,9 and12 each import through
a separate normal native zero-step process. Each exported final binary is
byte-identical to its complete input binary; all four exported histories
bind independently. Both logs contain zero accepted steps, there are no
inverse refreshes/reuses, and no advancing residual claim is made. The final
scientific generation still reloads after the older NPZ is removed.

See [retention assessment](wj-h27-linux/native-scientific-retention-v2h/science/scientific-retention-v16/retention-assessment-v2.json),
[actual receipt](wj-h27-linux/native-scientific-retention-v2h/science/scientific-retention-v16/protected-retention-v2/retention-receipt.json),
[four-write audit](wj-h27-linux/native-scientific-retention-v2h/four-recovery-write-audit.json)
and [post-rotation cold imports](wj-h27-linux/native-scientific-retention-v2h/science/scientific-retention-v16/post-retention-check-v2.json).

## Measured costs, storage and resources

| Native case | Whole operation s | Complete advancing s | Accepted steps |
| --- | ---: | ---: | ---: |
| Uninterrupted | 81.955984 | 76.543287 | 12 |
| Epoch1, fresh | 32.386005 | 26.876817 | 4 |
| Epoch2, cold from4 | 30.653572 | 24.771446 | 4 |
| Epoch3, cold from8 | 30.697549 | 24.766697 | 4 |
| Retained recovery9, cold zero-step | 6.519620 | 0 | 0 |
| Retained recovery12, cold zero-step | 6.419298 | 0 | 0 |

The uninterrupted initial block1–3 mean is6.711461s/step; warm blocks4–6,
7–9 and10–12 mean6.443475,5.948488 and6.411005s/step. The raw assessment
retains every Newton/Krylov count, full step timing and serialization timing.
Every accepted step takes4 Newton iterations. Krylov counts are
9,7,8,7,7,7,7,7,7,7,7,7 uninterrupted; the split processes respectively use
9,7,8,7 / 7,7,7,7 / 7,7,7,7. Zero-step imports perform neither iteration.
The three V1 statistics operations cost6.616354/6.586299/6.534136s; V2 costs
6.535232/6.430429/6.488751s. Independent moments/history plus actual retention
cost41.681207s. These external scientific observations/reductions/audits are
not relabeled equation work. Native frame/index timers are subsets of native
complete step/export totals, not additions to them. Each cold9/12 process
has approximately1.10s construction,3.62s initialization and.042s final
export; whole-process costs above also include normal process overhead.

At assessment, the entire new experiment uses216,101,101bytes; native epoch
outputs plus V2 generation working data use44,467,742bytes. Inputs, FLD,
logs, JSON and qualification evidence are included in the run-directory
measurement. Maximum measured SCI/recovery/statistics payloads are
4,709,595/859,058/8,542,328bytes. An illustrative500-step every-step batch
at this700-cell recipe projects2,537,601,650bytes for raw science, recovery,
final and statistics payloads, plus2,875,875bytes of current-batch log/index
data. This is a measured fixture extrapolation, not a rigorous future-byte
bound or selected production storage qualification. Metadata, FLD, fixed
inputs and cumulative JSON history remain separately charged.

Likewise, dt=.005 would require101005 steps for the unchanged T. Applying
the largest measured warm block to that count gives650,823.196s of
advancement, or813,528.995s with1.25 provision, before other required costs.
Four such trajectories total2,603,292.785s before the provision and fixed
costs. These are **unqualified short-fixture extrapolations**, not complete
production forecasts or a substituted timestep selection. They do not pass
the per-trajectory86,400s complete-cost gate. A larger dt requires actual
physical accuracy evidence; quiet-state/kernel timing is insufficient.

This phase charges331.783368565841s including6.365610752953s publication,
within the1050s admission. The same aggregate account remains OPEN at
12,204.860811439692s, leaving30,995.139188560308s of the43200s ceiling;
3600s remains protected for evidence. There are no RESERVED operations.
It is a progress checkpoint, not a closed H27 resource account.

Maximum actual shared-host workers12/72 and process-tree RSS770,551,808bytes
are within limits. Source/build/runtime uses2,623,931,230bytes and all H27
science1,666,090,661bytes, each below8GiB; free ext4 space70,955,089,920bytes
exceeds15GiB. The general Linux checkout's exact head and original dirty
developer-guide entry remain preserved. The main Mac worktree and unrelated
untracked pilot remain untouched.

Evidence archive:1,183,442bytes, SHA256
`96c064d25a23e75e2076c102ba34b1140d95ddf9e20955383362e0d6d3e04317`.
All265 archived inventory files were verified after transport. Large raw
coefficient/owner/profile/NPZ arrays and the protected qualification archive
remain under the Linux experiment; they are not copied into Git. Duplicated
older helper/source exports are also omitted from the compact Git selection;
their exact archive paths and hashes are recorded in
[repository-evidence-selection.json](wj-h27-linux/native-scientific-retention-v2h/repository-evidence-selection.json).
The full transport archive remains available locally and on Linux. See
[inventory](wj-h27-linux/native-scientific-retention-v2h/inventory.json) and
[publication receipt/account](wj-h27-linux/native-scientific-retention-v2h/publication-receipt-and-account.json).

## Remaining closure and executable commands

| Required status | Current result |
| --- | --- |
| H27_MODEL_AND_REFERENCE | NOT_MET: latest5600-cell inner-k physical RMS/MAX.1326434/.4236560 versus.05/.10; latest x-only numerical refinement.00802872/.13088236 versus.001/.001 |
| H27_COUPLED_NUMERICS | Full qualification NOT_MET; independent kernel, continuous-MMS, temporal/spatial, current-equation and short-history foundations documented in prior checkpoints |
| H27_PHYSICAL_PILOTS_AND_ACCURACY | NOT_MET: no selected h/p/q/dt and no required>=26.13350407707765 baseline/startup/both-sign evolving matrix |
| H27_RESTART_MONITORING_STORAGE | NOT_MET: this short cold/statistics/rotation foundation passes; selected-case monitoring escalation and complete production storage remain open |
| H27_COMPLETE_COST | TARGET_NOT_MET: no complete accuracy-qualified per-case/campaign forecast |
| TASK3_READY | NO |

The finite-resolution reference results do not establish a resolved continuum
barrier. The concrete next scientific decision is the separately documented
[v4 proposal](wj-h27-historical-inner-transport-proposal-v4.md), not a silent
change to scalar diffusion/cutoff. The original physical horizon
505.02466521278404 and terminal interval
[447.8483832973273,505.02466521278404] remain binding. No qualified Task3
launch command exists and the full horizon was not launched.

The executed fixture/controller and every session/options file are in the
evidence. For a fresh reproduction use a separate experiment directory and
new operation labels/output prefixes; all existing outputs are immutable.
The normal solver invocation was:

```bash
R=/home/jay/build/wj-h27-coupling-r1-20261117
"$R/runtime-v15/IncNavierStokesSolver" --use-opt-file \
  "$R/science/scientific-retention-v16/control.xml.opt" \
  "$R/science/scientific-retention-v16/control.xml"
```

That exact case has already executed and must not be overwritten. A bounded
read-only recheck of the retained cold imports can use a new immutable output
name, inside the existing aggregate account:

```bash
R=/home/jay/build/wj-h27-coupling-r1-20261117
python3 -B "$R/tools/run.py" --root "$R" \
  --label scientific-post-retention-review-v1 --phase analysis --seconds 30 -- \
  python3 -B "$R/tools-scientific-retention-v16b/h27/scientific_recovery_check.py" \
    --manifest "$R/science/scientific-retention-v16/post-retention-manifest-v2.json" \
    --out "$R/science/scientific-retention-v16/post-retention-review-v1.json" \
    --seconds 25
```

This recheck does not advance physical time or qualify Task3. Mutating
retention must not be rerun in the already rotated namespace.
