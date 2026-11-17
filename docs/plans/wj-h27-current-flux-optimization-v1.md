# Equivalent current-flux solver changes, source v10

This is a numerical optimization of the user-approved v2 equations. It changes
neither stress, raw production, dependent dissipation, scalar source/diffusion,
physical boundary data, nor any acceptance threshold.

Let `M=G+G^T`, `s=P M n`, `a=(1-chi)nu_A`, and
`b=nu+chi nu_W`. Use the positive frozen operator
`B=b M+a(s tensor n+n tensor s)` and the complete current remainder `V-B`.
Its work is `B:G=b|M|^2/2+a|s|^2>=0`. Below the collar `B=V` including weakly
divergent planar fields. The old isotropic-mu split remains selectable as
`WJStressSplit=IsotropicMuV1`; the new default on this opt-in solver is
`PositiveShearV1`. Existing constitutive `mu` remains the admissibility scale.
Both b and a enter measured preconditioner drift. Molecular viscosity is used
as the scale for zero shear entries in this refresh-only diagnostic.

Complete current residuals directly integrate `grad(test):V`, the unchanged
pressure/AWF and transported-Q3-k actions, and each scalar's own signed
diffusion/reaction/source actions. This avoids rebuilding sparse matrices just
to multiply them by their current coefficient vectors. Transient residuals form
full BDF coefficient differences before mass action, including essential/Gamma
history. The implicit linear equations continue to assemble the current matrix.
`WJCheckCurrentActions=True` additionally compares all full rows, including
essential reactions, to the assembled route with a fixed physical-scale 1e-10
bound. The steady route checks the augmented continuity/gauge rows as well;
an explicit uncommitted weak-divergence/nonzero-lambda probe exercises them.

Steady Newton finite-difference samples retain scalar polynomial positivity,
cutoff interval, root and covariance checks. The full continuum principal audit
is applied to the initial/base states and every nonlinear line-search candidate
before recording a best state or accepting it. It is also repeated explicitly
at final physical acceptance. These uncommitted derivative probes are not
accepted physical states. Only exactly identical unknown vectors within a fixed
steady solve may reuse a fully validated closure. Domain and full nonlinear
validation counts/times are reported separately.

Prospective Linux qualification: generic guard/exception controls; independent
full/outer scalar direct actions; mixed-boundary regressions; seeded nonlinear
polynomial steady MMS with all-row and augmented checks; evolving polynomial
MMS; and same-inlet physical first-step comparison of both splits, including
independent transient equations and <=1e-8 scaled field differences. No speed
claim is made before these measurements. A bounded physical steady continuation
may then reuse the explicitly pinned v9 candidate only as an initial guess and
must solve the v10 current equations afresh.

Admission against the continuing 43200 s account (3998.5351 s used at publication,
3600 s protected for evidence): source materialization/fixtures <=90 s;
build <=600 s; generic controls <=60 s; focused mixed/scalar checks <=300 s;
tiny coupled checks/oracles <=780 s; two physical diagnostic prefixes <=420 s
and independent split/current checks <=420 s; steady continuation and its
independent assessment <=1730 s; finer MMS temporal/spatial studies <=1200 s.
The conservative sum is <=5600 s, before any longer physical pilot admission.
Most tiny operations are expected to take substantially less than their caps.
The amended forecast is frozen while the first v10 build is running and before
its dependent science is admitted. Build uses four workers; science uses one
rank and numerical-library threads one, subject to the shared 72-actual-worker
ceiling and existing RSS/storage watchdogs. Failed attempts remain charged.
No full Task3 run is admitted here.
