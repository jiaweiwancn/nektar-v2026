# WJ-H26 fresh reproduction

The executed source is `8aeb00a2f37437cd17e043a9ba6f86a6f208e1d6`. The original
H26 journal and capsules are closed. This recipe uses a new isolated parent,
a fresh independent operation journal and version1 outputs; it reads the
original immutable source/reference geometry and qualified stock dependencies.
Choose an unused parent and an explicitly allocated new execution allowance.
The runner retains the H26 limits:21600 operation seconds,1800 protected for
evidence,72 total host workers, four compiler jobs, one MPI rank, six numerical
thread controls1,8/16GiB science/build RSS,8GiB source/build and science storage,
and15GiB free-space reserve.

Run on the authoritative Linux server through `ssh jay@100.105.247.110`:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export H26_BASE=/home/jay/build/wj-h26-native-flow-20261117
export H26_REPRO=/home/jay/build/wj-h26-reproduction-new/wj-h26-native-flow-20261117
python3 -B - <<'PY'
import hashlib,json,os,shutil
from pathlib import Path
base=Path(os.environ['H26_BASE']);root=Path(os.environ['H26_REPRO'])
assert not root.exists(), 'choose a new parent; never overwrite a capsule'
root.mkdir(parents=True)
for name in ('tmp','science','operations','inputs'):(root/name).mkdir()
source=root/'source-v1'
shutil.copytree(base/'source-v8',source,ignore=shutil.ignore_patterns('__pycache__'))
def identity(p):
    return dict(path=str(p),bytes=p.stat().st_size,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest())
(root/'source-v1.json').write_text(json.dumps(dict(
    revision='8aeb00a2f37437cd17e043a9ba6f86a6f208e1d6',
    copied_from=identity(base/'source-v8.json'),
    files=[identity(p)for p in sorted(source.rglob('*'))if p.is_file()]),indent=2)+'\n')
(root/'accounting.json').write_text(json.dumps(dict(
    schema='WJH26FiniteCampaignV1',closed=False,ceiling=21600,
    evidence_protected=1800,entries={}),indent=2)+'\n')
template=root/'science/profile-v3';template.mkdir()
shutil.copy2(base/'science/profile-v3/case.xml',template/'case.xml')
PY
export TMPDIR="$H26_REPRO/tmp"
T="$H26_REPRO/source-v1/solvers/IncNavierStokesSolver/Tests/RANS/h26"
U="$H26_BASE/tools-v10/solvers/IncNavierStokesSolver/Tests/RANS/h26"
python3 -B "$T/run.py" --root "$H26_REPRO" --label fixtures-reproduce --phase analysis --seconds 60 -- python3 -B "$T/fixtures.py" --out "$H26_REPRO/science/fixtures-v2"
python3 -B "$T/run.py" --root "$H26_REPRO" --label build-reproduce --phase build --seconds 300 -- python3 -B "$T/build.py" --root "$H26_REPRO" --version 1
python3 -B "$T/run.py" --root "$H26_REPRO" --label prepare-reproduce --phase analysis --seconds 60 -- python3 -B "$T/campaign.py" --root "$H26_REPRO" --version 1 --mode prepare
python3 -B "$T/campaign.py" --root "$H26_REPRO" --version 1 --mode verify
python3 -B "$T/campaign.py" --root "$H26_REPRO" --version 1 --mode constant_limit
python3 -B "$T/campaign.py" --root "$H26_REPRO" --version 1 --mode benchmark
python3 -B "$T/run.py" --root "$H26_REPRO" --label analyze-reproduce --phase evidence --seconds 120 -- python3 -B "$U/analyze.py" --root "$H26_REPRO" --version 1 --stage C --output "$H26_REPRO/publication-reproduce"
```

The source/build/runtime configuration and all original loaded-library hashes
are in [the identity record](wj-h26-linux/checkpoint-c/identities.json) and
[runtime audit](wj-h26-linux/checkpoint-c/runtime-resource-audit.json). The build
helper reads the pinned old stock objects and verifies their dependencies
against the new copied source; it does not overwrite those objects or the shared
Nektar++ installation. It produces new executables and loader wrappers.

Each scientific child has its own immutable request/receipt/log. Existing passed
labels are checked/reused; a failed attempt requires a new explicit label after
diagnosis. Preserve failures and check all unknown remote outcomes before
launching replacements. Do not weaken current-system residual targets or change
the source/geometry/parameters to reproduce a timing decision.

The original final audit additionally verifies old-account/source preservation,
all original mapped libraries, resource bounds and publication transport hashes.
A new reproduction needs its own end-of-campaign audit/close, with its new
namespace and original-reference paths recorded explicitly. The original H26
ledger must remain closed.
