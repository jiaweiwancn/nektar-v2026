# H20 resumption: actual commands and immutable roles

The original root and new root are both CLOSED. These commands document executed
receipts; they are not instructions to relaunch their consumed reservations.

Mac worktree: `/Users/jay/Documents/nektar-v2026/.worktrees/wj-earsm-awf`.
New Linux root: `/home/jay/build/wj-h20-three-level-resume-20261115`.
The original root `/home/jay/build/wj-h20-three-level-20261115` was imported
read-only. Native binaries/private closure remained in the manifest-selected
H19 runtime. No packages, mounts or unrelated campaigns were modified.

All Linux imports/launches used:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

`resume_wj_h20.py` initialized the linked ledger, verified the old sealed prefix,
performed bounded fsync/readback and copied the source capsule. Corrected source
was synchronized from the Mac. `run_wj_h20.execute` atomically reserved and
metered every following operation. Exact argv, caps, inputs, workers and elapsed
receipts are in the publication/operations tables.

In `source/.../Tests/RANS`, actual operations were:

```sh
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode setup
python3 -B -m unittest -v test_wj_h20_three_level test_wj_h19_timestep test_wj_h18_checkpoint_audits
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode preflight
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode admission --level M
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode target --level M
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode strict --level M
```

M's governor succeeded. The passive serialized-step assertion was then resolved
using `repair_wj_h20_postprocess.py` on saved native/commit/strict bytes, without
another M launch. Original source/config/receipt files were retained. A versioned
`source-v2` and active configuration/cache-compatibility record were frozen; the
same test modules ran there, including halt-preservation and type regressions.

In that versioned source, actual remaining native operations were:

```sh
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode admission --level C
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode target --level C
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode strict --level C
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode bootstrap
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode admission --level F
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode target --level F
python3 -B study_wj_h20.py --root /home/jay/build/wj-h20-three-level-resume-20261115 --mode strict --level F
```

The existing `step_zero_bootstrap` reader operation created only F's CP0, with
source XML/checkpoint under dt=.0025 and destination under.00125. It generated
native history, not a copied/rescaled history or relabelled header. All governing
restarts used own-dt-compatible histories and original qh0 anchors. NumSteps was
4/2/16 additional work; absolute ceilings8/4/16; XML Time stayed0.

After native exit, runtime maps were checked using
`study_wj_h16.runtime_qualify`. Publication ran from versioned publication-code
directories. The first passive export hit its time cap; the exact-byte-tested
once-per-key NPZ reader produced publication-v2 in a new charged operation.
`close_wj_h20_resumed.py` sealed linked accounting and the essential archive.
No native work, reader or temporal trajectory was repeated for this export fix.
