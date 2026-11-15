# CT14 evidence and reproduction

Authoritative Linux archive:

`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-14-20261115`

The mount was verified as `/dev/sda2`, `fuseblk`. Source/Git development uses
the dedicated Mac worktree under `/Users/jay/Documents/nektar-v2026/.worktrees/`.
The Linux main checkout, its unrelated dirty work, shared installation, older
archives and historical physical failures are preserved.

## Bindings and roles

| Archive | Role |
|---|---|
| build-v1 / capacity-v1 | Scalar Scotch ND; original factor40,862,187 rejects |
| build-v2 / capacity-v2 | Native-vector Scotch ND; original factor40,119,892 rejects |
| build-v3 / capacity-v3 | Preserve external ND order; same40,119,892 rejects |
| build-v4 / capacity-v4 | Directed native-vector COLAMD; all four frozen controls pass |
| build-runtime-v1/runtime | Immutable qualified governing/diagnostic executables and dependencies |
| component-v1 | Four complete small operators, independent actions, original full-domain update, restart/rollback and incompatible identities |
| collect-component-v1 | Independent primitive/history/resource recomputation |
| static-v1 | Complete62533-slot selected static admission, two canonical exports and four initialized frozen scalar controls |
| collect-static-v1 | Independent static verification plus regional initial absolute pressure |
| pair-v1 | Fresh .001/.0005 histories to the commonT=.002; physical decision is separate from stepping |
| resources-* | Closed owned-family CPU/wall/worker receipts, including failed attempts and controllers |

Native source for the admitted component is `6477d9d9`. Governing executable
SHA256 is `9ae7d3c63492c82a4c9f979616f673d8e4909a5c41ee8069fc1592f2fb23da06`;
diagnostic probe SHA256 is
`95d6a0258926ab250b69ea0940976f43a993bf4bdccc3001a5819f4f36835824`.
Component controls use that source; the selected static controller uses
`087fa720`; static collection uses `c0224508`; first physical-pair dispatch is
`eca8f3ab`. Subsequent documentation/publication changes do not rebuild or
relabel those runtimes. Request manifests bind individual controller SHA256s,
contracts, inputs, runtime/library bytes and predecessor seals.

`capacity-component-v1/` contains complete compact numerical results and
provenance. `static-v1/` contains publication views. Their
`compact_publication` metadata binds the complete original external bodies;
large numerical vectors are replaced only by count/hash/archive references.
Those views are not checkpoint or replay inputs. Full CSR/state/history/stage
files remain at their external hash-bound paths. The initial static publication
copied a full initialized vector; the normal subsequent correction removes it
from the current view without rewriting the earlier checkpoint or native seal.

The original CT8 physical manifest remains SHA256
`004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296`.
The frozen experimental dataset and completed digitization audit are unchanged.
The original CT11 pressure failure and every rejected ordering remain failures.

## Commands

Use fresh destinations. Never rerun into a sealed directory. Preserve complete
requests and their `completed-artifacts.json` before collecting results. Each
command below is run through the owned worker launcher:

```text
python3 -B TESTS/ct_accepted3_command.py --control NEW/resources-LABEL -- COMMAND
```

Set all six thread variables to one:
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
`BLIS_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`.
Use an isolated `TMPDIR` on the verified large mount. Retain16-project/72-aggregate
worker limits, two compiler jobs, one numerical family, and16GiB per-process
address space. Every failed family remains charged to the same allocation.

Here TESTS is the pushed source's
`solvers/IncNavierStokesSolver/Tests/RANS`; NEW/source is its clean Linux
checkout. CT9--CT13 are the preserved sibling archives named above. RT is the
new private runtime. The actual requests record full paths and executable/
library hashes; these templates omit only repeated path prefixes.

```text
ct_accepted14_build.py --root NEW --ct9 CT9 --output NEW/build-capacity --phase capacity
ct_accepted14_capacity.py --ct13 CT13 --runtime NEW/build-capacity/runtime \
  --contract SOURCE/docs/plans/cao-tamura-accepted-14-evidence/contract.json \
  --policy VectorCOLAMDV1 --output NEW/capacity
ct_accepted14_build.py --root NEW --ct9 CT9 --output NEW/build-runtime --phase runtime
ct_accepted14_component.py --ct9 CT9 --ct10 CT10 --ct11 CT11 --runtime RT \
  --capacity NEW/capacity --contract CONTRACT --output NEW/component
ct_accepted14_static.py --ct12 CT12 --ct13 CT13 --runtime RT \
  --component NEW/component --capacity NEW/capacity --contract CONTRACT --output NEW/static
ct_accepted14_pair.py --runtime RT --static NEW/static --component NEW/component \
  --capacity NEW/capacity --contract NEW/admitted-pair-contract.json --output NEW/pair
ct_accepted14_collect.py --root NEW --runtime RT --capacity NEW/capacity \
  --component NEW/component --static NEW/static --pair NEW/pair \
  --self-control NEW/resources-collection --output NEW/collection
```

For replication, the prospective admitted-pair contract must bind the newly
qualified static/component/capacity result SHA256s and actual runtime/initial
identities. Keep every numerical/physical limit, horizon, scale and classifier
unchanged. The original contract is intentionally rejected for different new
receipts; an old history must never be relabelled to match a new component.
Keep actual first-pair .001/.0005 request identities. Publication compaction
does not change these execution inputs or their seals.

The original failed ordering controls can be reproduced individually with their
preserved utility/runtime, exact canonical A/RHS and recorded work allocation.
Their commands and hashes are in capacity-v1/v2/v3 request/result bodies. No
timeout is classified as measured fill. Do not repeat known failures as a blind
ordering sweep.

The receipt classifier's discriminating controls run with:

```text
python3 -B -m unittest discover -s TESTS -p test_ct_accepted14_receipts.py -v
```

Actual evolution retains optional factor/solve receipts with
`NEKTAR_RANS_SPARSE_RECEIPTS`; detailed timing is off. Every original-system
residual, refinement count, factor fill and canonical binding remains required.
Scalar/preconditioner COLAMD, consistent-mass LDLT and the full
`D[(alpha*M_free)^-1 G]-C` pressure action are unchanged.

## Interpretation

Static/operator/observer success does not establish a physical-pair pass.
Only the unchanged complete physical classifiers on the fresh common-time
histories decide short qualification. Initial representation error, own-initial
drift, absolute reference error, common-time differences, physical current-K
tractions, scalar budgets and regional pressure remain separate quantities.
Imposed epsilon trace and refreshed post-k lag are separate; no repair occurs.

A passingT=.002 pilot establishes this short numerical qualification only.
Long flat preservation still needs a passing longer-time contract and measured
CPU/wall/output affordability with qualified retention/checks. The first hill
also needs mapped region-aware numbering, Cartesian/cofactor gradients,
conservative transmission, outstanding curved spatial/temporal gates and
physical roughness-sublayer/thin-layer/approach-flow admission. No terrain,
other-model displaced route, EP, MPI/3D or variable-height production is launched.
