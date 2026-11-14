# CT-Accepted-9: pressure consistency and serial scalable one-sweep

The reviewed and actual starting revision is
`0c28610e685aef95f7e6f21b9a4fdcd305dbf319`, clean in the implementation
worktree. The separate main/CLS/WJ worktrees and their untracked files are
untouched. Mac HTTPS fetch timed out; the Linux SSH-over443 remote check
confirmed the same head. The active storage is `/media/jay/新加卷1`.

The [contract](cao-tamura-accepted-9-evidence/contract.json) precedes new
trajectories. Task A uses the sealed CT8 dense runtime and exactly one fresh
.0005 history. Task B independently qualifies `AcceptedSerialSparseV1`; neither
task can supply a missing gate for the other. All CT4--CT8 results remain sealed.

The physical input is the CT8 central three-point displaced A95 reference at
7.5mm, physical floor0, displacement3.5mm. The original coefficients, Natural
q/kappa, surface roughness, viscosity conventions, full-height inlet, RN roof,
absolute outlet traction and one-sweep ordering remain fixed. No reference
restoration, pressure reset/regauging, clipping, source cancellation or extra
governing correction is allowed.

Initial dependency inspection found Eigen sparse headers and the installed
Nektar GMRES/matrix-action interfaces. PETSc remains OFF. The new adapter will
assemble sparse rows directly, retain consistent mass, and apply the pressure
Schur operator without storing M-inverse-G columns. Small bounded dense trace
and radial projection systems remain explicit. Results, failures and exact
admission scopes will be appended after authoritative Linux qualification.
