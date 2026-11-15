# H24 authoritative compact evidence

Status: **CLOSED PASS**,28/28 advances from H20 .02 to .04, three terminal strict
reads. The scientific execution and closed ledger are final. No campaign should
be relaunched from these reproduction commands without a separate authorization.

The authoritative healthy-ext4 Linux root is
`/home/jay/build/wj-h24-flat-continuation-20261115`.
[authoritative-index-v1.json](authoritative-index-v1.json) pins publication-v1,
closure-v1, checkpoints A/B and the unchanged historical H23 CLOSED ledger.
Its SHA256 is
`b638402b15e8242bbdf2ab825f3753050fb241a03022a40f7aa25dc30fd356bc`.

Read [the evidence report](../wj-h24-flat-continuation-evidence.md),
[frozen contract](../wj-h24-flat-continuation-contract.md) and
[exact executed commands](../wj-h24-reproduction.md). Scientific Python capsule
`75457717ebb89a42bcc1cb2155c8cd0fff9ec94f`, native capsule
`057e6b7465660cac5ec8617abf5c876fa8bff71f` and passive publication source
`3133cfd3` have separate per-file/runtime pins. The documentation commit does not
replace compiled-source provenance.

Useful compact entry points:

| Path | Meaning |
|---|---|
| `closure-v1/verification.json` | Six scientific closure criteria, operation counts and capacity checks |
| `closure-v1/accounting-final.json` | Authoritative CLOSED final spend; includes setup/closure charges |
| `publication-v1/component-gates.json` | Scientific gates and the then-active publication reservation snapshot |
| `publication-v1/startup-restored-bindings.json` | Actual H20 .02 producer/callback/commit and own-CP0 chains |
| `publication-v1/state-history-bindings.json` | All28 new recovery/completion/callback/qualification bindings |
| `publication-v1/audit-role-status.json` | M-only, C-pair and F-triplet availability and physical assessments |
| `publication-v1/all-signed-Q1-moments.csv` | Complete255+525 moments for every participating assessment |
| `publication-v1/all-pair-regional-sensitivity.csv` | All-pair, all-region RMS/MAX absolute screens |
| `publication-v1/contraction.csv` |200 declared resolved/unresolved records; raw differences/floors retained |
| `publication-v1/fixed-probe-triplets.json` | Signed fixed-owner values and common actual directions |
| `publication-v1/representation.csv`, `reference-drift.csv` | q_eq/original-qh0/.02-change separation |
| `publication-v1/terminal-strict-residuals.json` | Genuine completed-stage versus spatial consistent-mass dual norms |
| `publication-v1/all-timing-scopes.json.gz` | Inclusive/exclusive native timings, lossless copy |
| `publication-v1/actual-storage-categories.json`, `conditional-cost.json` | Measured bytes and unadmitted future cadence scenarios |
| `publication-v1/heavy-manifest.json` | Original Linux heavy arrays/checkpoints remain in place |
| `local-storage-aliases.json` | Local representations and original hash/size bindings |
| `transport-verification-v1.json` | Exact archives,192 sealed files,136 lossless transformations and figure QA |

Both compact archives were transported and verified; these are lossless transport
containers, not additional scientific payloads:

| Archive | Bytes | SHA256 |
|---|---:|---|
| `h24-essential-v1.tar.gz` |2,435,098|`45157c23cda852b360149d43fb1dd7ab524b0d1694e72c1589256c2d875d370b`|
| `h24-compact-v1.tar.gz` |4,156,535|`1e30f9cbf040d59c715abc9793312237598478f41d4058b672fbaf3e0e7c2c7c`|

`publication-v1/copied-evidence.json` records source pins and gzip transforms;
original remote publication-manifest bytes are retained. The local manifest
covers the local representations, including A/B, and excludes itself. Verify it
from this directory without importing numerical modules:

```bash
python3 - <<'PY'
import hashlib, json
from pathlib import Path
p = Path('.')
m = json.loads((p / 'local-compact-manifest.json').read_text())
for item in m['files']:
    b = (p / item['path']).read_bytes()
    assert len(b) == item['bytes'], item['path']
    assert hashlib.sha256(b).hexdigest() == item['sha256'], item['path']
print('PASS', len(m['files']), 'local files')
PY
```

The publication snapshot conservatively charges an in-flight180-second operation;
closure replaces that reservation with its actual24.692938-second receipt and
separately seals final accounting at3052.381070 seconds. No mutable journal is
used as an immutable accounting seal. H23's failure/charges remain separate.

Numerical finite-support PASS and short-window contraction do not admit production
dt, the unchanged505.02466521278404 horizon, terrain or experimental use. The
physical7 mm observations remain inside the9.5 mm matching plane.
