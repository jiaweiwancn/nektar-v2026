# H20 executed commands and source roles

Source/Git worktree:
`/Users/jay/Documents/nektar-v2026/.worktrees/wj-earsm-awf`.
Linux root: `/home/jay/build/wj-h20-three-level-20261115`.
The root/ledger is closed; the failed native reservation must not be relaunched.
Commands below describe the executed receipts, not permission to reopen it.

All Linux numerical imports/launches used:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

After selected-input bootstrap and source transfer, the existing atomic execution
helper admitted these commands (absolute paths are recorded in each request):

```sh
python3 -B /home/jay/build/wj-h20-three-level-20261115/source/solvers/IncNavierStokesSolver/Tests/RANS/study_wj_h20.py --root /home/jay/build/wj-h20-three-level-20261115 --mode setup
python3 -B -m unittest -v test_wj_h20_three_level test_wj_h19_timestep test_wj_h18_checkpoint_audits
python3 -B /home/jay/build/wj-h20-three-level-20261115/source/solvers/IncNavierStokesSolver/Tests/RANS/study_wj_h20.py --root /home/jay/build/wj-h20-three-level-20261115 --mode admission --level M
python3 -B /home/jay/build/wj-h20-three-level-20261115/source/solvers/IncNavierStokesSolver/Tests/RANS/study_wj_h20.py --root /home/jay/build/wj-h20-three-level-20261115 --mode target --level M
python3 -B /home/jay/build/wj-h20-three-level-20261115/source/solvers/IncNavierStokesSolver/Tests/RANS/study_wj_h20.py --root /home/jay/build/wj-h20-three-level-20261115 --mode strict --level M
```

The target command reserved four attempts but its one governor process exited
with `invalid bounded monitoring limits` before restoration. The strict command
selected the original H18 M anchor, not a newly completed H20 state. No C/F,
bootstrap, extra governing trajectory or build command ran.

The same three unittest modules were then run through a separate
`source-correction-v1/.../Tests/RANS` import path, using Python `-c` to place that
path first. Its exact command is in the post-stop receipt. This executed the two
new preflight/halt methods, for 26 methods; the first run had 24. Frozen execution
files at `source/...` were not replaced. Runtime maps were qualified with
`study_wj_h16.runtime_qualify` on the retained H20 receipts.

Native governor/reader wrappers were the manifest-selected H19 wrappers under
`/home/jay/build/wj-h19-target-timestep-20261115/runtime`, preserving the private
loader/library closure. No host-package changes or alternative runtime occurred.

The original publisher at `source/.../publish_wj_h20.py` failed on an absent local
build directory. The corrected publisher/closer ran from `publication-code-v2`;
both prefixes and receipts are retained. Final transport verified the archive
SHA256 and every selected member before copying compact evidence into Git.

The current source fixes ACK timeout/preflight and first-failure reporting, but
those changes have no new governing execution evidence. Any subsequent numerical
qualification requires explicit new bounded admission and must preserve this
closed campaign and its failure.
