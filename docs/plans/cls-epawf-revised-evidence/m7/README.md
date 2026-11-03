# M7 evidence

These receipts distinguish verification of the faithful implementation from
qualification of the original physical model. Passing tests reproduce the
adverse transverse growth; robust-model stability remains failed. Stopped
unforced trajectories remain unqualified.

Start with `integration-audit.json` and the [M7 report](../../cls-epawf-revised-m7.md).
The integration receipt maps all 33 registered gates to their actual runs and
revisions. The 31 earlier gates use the complete CTest/JUnit receipt. The two
M7 gates use their actual native and strict-transverse driver receipts. This
is combined source-bound coverage, not a claim that all gates ran in one
CTest invocation.

- `native-validation.json`: final native campaign, including independent
  physical mass norms/energy/covariance, MPI and perturbed history restart,
  H dependence and three required initial-state rejection controls.
- `independent-replay.json`: independent requalification of the first complete
  native campaign, including actual trace-transfer GLL sampler locations.
- `transverse-validation.json`: independent indexed tensor, complete weak
  matrices, native IMEX1/2, covariance and independent dt/p/mesh comparisons.
- `dense-guard.json`: all other audited admissibility sets remain positive;
  dense reconstruction alone rejects the negative initial k polynomial.
- `active-energy.json`: independent initial mixed-chi and final unforced
  Upstream energy, physical sensor/gate and matching-traction accounting.
- `unforced-flow-diagnosis.json`: complete control/pulse and dt-half observed
  prefixes, preserved solver stops, independent covariance and energy replay.
  Collection success does not imply passage qualification.
- `regression-receipt.json` and `regression-junit.xml`: all 31 earlier gates.
- `qualified-source-sync.json` and `qualified-runtime-fingerprint.json`:
  matching compiled physics, runtime images, loader bindings, numerical
  environment, installed numerical packages and tools after source sync.
- `artifact-index.json`: hashes/locations of full raw attempts on Linux,
  including failed setup and unqualified physical-flow attempts.

The actual compiled numerical source is `bf202be`; the final audit drivers
are qualified at `dea75da`. There are no C/C++ changes between those revisions.
The frozen revised note, contract, original CLS and original scalar AWF hashes
are checked independently. The documentation-only closure is separately
compared against the full qualified numerical/runtime fingerprint. Its final
receipt is stored on Linux under `m7-final-sync/` to avoid a self-referential
Git commit hash in the committed receipt.

PNG/PDF scientific figures and their exact plotting/data receipt are included.
Runner/postprocessing scripts are retained byte-for-byte where their hashes
are cited. Generated meshes, full states, field files, MPI partitions and
checkpoints remain outside the source checkout under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised`.
