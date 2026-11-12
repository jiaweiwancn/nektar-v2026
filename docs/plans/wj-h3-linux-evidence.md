# WJ-H3: corrected initialization and bounded principal diagnosis

H3 corrects the physical reduced-pressure initializer and qualifies order12
among the tested outer integration orders. It does **not** remove the refined
dt=.02 principal obstruction: a fresh corrected run rejects T=.06 at an
interior point. Its dt=.01 companion and all four predetermined B controls
accept T=.08 without correction. These are short diagnostic transients, not a
developed or numerically converged rough-flat precursor. H2's T=8 denial and
WJ-7's RESOURCE_LIMIT/T3 NOT_RUN remain unchanged.

The clean fetched baseline was `b83ec8979cb029ab39136ef89ff14f88a5031a75`.
Mac source/Git work used the existing `codex/wj-earsm-awf` worktree; each commit
was pushed normally. Other worktrees, dirty main checkouts and concurrent
campaigns were preserved. The [contract](wj-h3-contract.md) was committed before
numerical execution. The additional matched base-mesh seed comparison is
zero-step only and consumes no governing attempt.

| Role | Exact source identity |
|---|---|
| Frozen contract | `d3e91c5355d82bb7c057577a330cfcbf0e60c38b` |
| Pressure correction/default seed and label correction | `612cf0d16da6573a1d679296fd4bfa183a54302d` |
| Producing H3 generator/control/passive-integration source | `c1f9652495913ab67ae5114907f2a0374fc4cb19` |
| Versioned all-state assessment | `12eb50ed27722f398d1a61226e62f8a6221a44c2` |
| Matched zero-step comparison and publication | `d35aee73b97d4c3c9dd229c61edd790070ec2673` |
| Reused governing C++ producer | `c344dc2209f030e426e2ee5501d5e5b72a6d2c46` |
| Reused native observer producer | `e62b5f4d265f068dd5af3f090bfb0e385c3ae53e` |

Linux is authoritative. USB SSH was unavailable; the documented Tailscale
fallback `ssh jay@100.105.247.110` was used. After the server reboot, the same
NTFS volume UUID `5AB8A669B8A642FD` was mounted at `/media/jay/新加卷1`.
The old mount directory was not remounted or populated. Historical input paths
were translated read-only to the new mount; their file contents were unchanged.
The new isolated root is
`/media/jay/新加卷1/nektar-v2026-runs/wj-h3-20261112`.
It contains separate producing `source`, `observer-source-v2`, and
`publication-source-v3` checkouts and immutable versioned outputs.

No native C++/header ABI input changed. All340 governing compilation-source
dependencies and all133 files in the sealed H2 runtime closure were verified
before reuse. Governing ELF SHA-256 is
`9ae37a1aba62755b10ad99c9484c498cef9d927a75fb950f00344f4ffc5d06fe`;
the passive observer is
`6f35b837f0dd13eb155e0fe34f0dc670c609ee6cbe739609a5de288dd49ed06d`.
The selected source-built SpatialDomains library is
`0c23e0de180a233bd88f70b265ae9865aeac456e6482f8ff0fa1b4be84d69dab`.
[Actual process mappings](wj-h3-linux/actual-loaded-runtime.json) contain129
files:123 immutable runtime/loader files match the sealed H2 closure, and six
are process-specific OpenMPI shared-memory files, not replacement libraries.
The [metadata comparison](wj-h3-linux/mac-runtime-closure-verification.json)
distinguishes those roles.

The reused native build is GNU11.4/C++17 Release, `-O3 -DNDEBUG
-ffp-contract=off -fno-fast-math`, `RANS_BUILD_SOLVER=ON`. H3 configured a
separate external CMake test directory with the same arithmetic flags and
source-built SpatialDomains selection; it did not rebuild or install native
code. Every run used one serial rank and all six OMP/OpenBLAS/MKL/BLIS/NumExpr/
vecLib thread controls set to1. Admission snapshots and five-second resource
guards accounted for other campaigns and MPI helpers under the72-worker cap.

Original Mac, worktree and Linux specification bytes agree:

```
main     1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124
appendix 0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c
```

The literal `quad`/`\quad` typography erratum and correct `(d+a*h)/a^2`
coefficient were retained. No specification was changed.

The physical sign correction follows directly from the complete stress tensors:
`T_w,nn=-pi_w` and `T_o,nn=2K/3-pi_o-R_o,nn` for the initial parallel flow.
Thus `pi_w-pi_o=R_o,nn-2K/3`. At the unchanged continuous H2 completion,
K=.012618181994557477 and Rzz=.0056162152115941805, so the corrected inner
offset is **-.002795906118110805**. The old positive offset gives the predicted
normal-traction error .00559181223622161. Before editing, retained native base
and refined CP0 maxima reproduced .005591598549377207 and
.005591162480535834. [Pre-edit verification](wj-h3-linux/preedit-verification.json)
binds those reports. The independent regression derives both physical tensors;
the old expression fails, while the corrected continuous balance is within
2e-18. The analytic neutral initializer's already-correct negative sign was
not reversed.

`SignCorrectedV1` changes only this physical pressure seed. Every nonpressure
Function value string in A is byte-identical to the retained H2 refined input.
`VerticalBalanceV1` separately sets C=Rzz(H)-2K/3, pi_w=C, and
pi_o(z)=2k(z)/3-Rzz(z)+C, using the same completion and physical gradient.
Continuous Tzz=-C throughout both regions. Its pressure and frozen-traction
expressions are approximated on declared physical panels; splitting these
panels fixed an initially detected2.80e-11 approximation error without changing
U/k/epsilon. Final pressure/traction expression tests retain the2e-11 limit.
The physical Q1 pressure projection is a separate error source.

The represented CP0 is not the continuous profile. In A, the corrected
interface normal maximum on native interface reporting points is about3.00e-5,
with independently reconstructed interface stresses agreeing within4.3e-17.
In the base-mesh balanced CP0 the maximum is about4.98e-5. No strong H pressure
row was added. At identical base23/q8 inputs the two corrected initializers have
exactly identical differential coefficients; their uncancelled vertical
consistent-mass dual residuals are7.19215390e-4 (constant corrected offset)
and5.75677575e-4 (balanced seed). The latter is not machine-zero equilibrium.
Startup has no completed-stage packet. [Matched result](wj-h3-linux/matched-startup.json)
and [figure](wj-h3-linux/startup-balance.png) keep these distinctions.

Only completion labels changed: L=f*delta remains the input Blackadar
parameter, while ell(infinity)=L+kappa*H. At f=.09,H=.5,delta=6.25 the actual
asymptote is .12055342908373864*delta. `WJFiniteDepthCompletionV2` reports both
quantities. The epsilon function, U/k completion and historical H2 files were
not changed to make the old label true. Historical H2 wording identifying
.09*delta as the actual asymptote is superseded by this correction.

The retained-state diagnosis reconstructs actual native coefficients, including
both normal-face owners and points x=7.5,7.9,7.99,7.999999,8. At the old failing
face z=.6720111978440669, owner17's original represented CP0 at x=8 has
Gxz=.2384288826349788, tau=11.926411670922962 and principal/nu_eff=.06166322.
The continuous prescribed values .24047947338930892/11.92700416 are therefore
not substituted for native CP0. By rejected T=.06, owner17 has
Gxz=.30425583808613166, Gzz=.0651082852991303 and principal/nu_eff=-.01239024.
It is already negative at x=7.5 and7.9; this is not just an isolated endpoint.
The adjacent owner18 at x=8 has Gxz=.23598320910954942 and positive .06531218.
The polynomial response extends through the element, with differing traces.
[Owned-face data](wj-h3-linux/failing-face.csv) and the
[figure](wj-h3-linux/retained-face.png) retain both owners.

Independent local algebra and the compiled shared WJ action agree:

| Unapplied diagnostic state | Principal / nu_eff |
|---|---:|
| Original rejected state | -.01239023667349 |
| Rejected G, prescribed-profile k/epsilon | -.01220802347033 |
| Prescribed-profile G, rejected k/epsilon | +.05855628243464 |
| Continuous prescribed profile | +.05876076287334 |
| Remove tr(G)I/3 from rejected G | -.01239023667349 |
| Keep only rejected Gxz | -.01799627515994 |
| Planar tracefree diagnostic | -.02188232748098 |

These are diagnostic substitutions, never state corrections. WJ already uses
deviatoric strain; removing the three-dimensional trace leaves stress/principal
form unchanged. Covariance and raw production remain positive at the original
failure. [Retained summary](wj-h3-linux/retained-summary.json),
[supplement](wj-h3-linux/closure-v1/retained-supplement.json), and
[adverse packets](wj-h3-linux/adverse-witnesses.json) retain directions, physical
positions, covariance/root diagnostics, states and genuine callback hashes.
The tests also preserve the q=3.5 adverse principal witness. No universal
q>3.5 admission rule was introduced.

Passive outer integration used8/12/16/24/32 on identical represented fields and
test functions, split by actual elements/regions. No order48 was needed.
The24/32 reference differences are below the original component/regional floors.
The old native8 versus independent24 outer-stress discrepancies are reproduced
as3.81343321e-9 (base T=.08) and4.35353287e-9 (refined half-dt T=.08), within
2e-13 of the published3.81343314e-9/4.35353290e-9. Identical-native-node algebra
still passes independently. The old refined regional-sum failure stays recorded.

The largest8/32 discrepancy is mode1, the x/5 test, block1 vertical momentum:
global3.81335485e-9/4.35345723e-9. Its leading element contributions are
4.36170094e-9 in base geometry21 and4.36066025e-9 in refined geometry41.
Both cover physical z=[7.507811773998694,12.5], the finite-depth outer tail,
not the near-H principal witness. Full signed per-element arrays remain in
the indexed heavy views. The additional base T=.08 confirmation is a labelled
post-control passive supplement; it did not change q_star or rewrite the
original selection/publication.

Order12 is the smallest demonstrated adequate member of8/12/16 on retained
and corrected startup fields. All B endpoints also pass against the converged
reference: their outer-stress12/32 errors are1.26e-10–2.82e-10 and24/32 regional
changes below8e-14. Native full inner/outer/exterior three-mode checks pass the
unchanged2e-9+2e-8*abs(component) and4e-9+2e-8*abs(regional sum) limits.
TangentialOrder changes outer x/normal integration and inner tangential ray
consumers; NormalOrder remains16 and exterior order max(16,q)=16.
H2 checkpoints keep their original signatures; these are explicit passive
integration views. B uses fresh native q12 CP0s. See
[selection](wj-h3-linux/quadrature-selection.json),
[table](wj-h3-linux/quadrature.csv), [weak loads](wj-h3-linux/weak-budgets.json),
and [figure](wj-h3-linux/quadrature.png). A retains its original q8 limitation.
An absolute weak-load error was not compared to a dimensionless principal
margin as evidence of causation. No quadrature-only governing contrast was run.

The original WJ constants/tensors, DevStrainRawWorkV1, direct epsilon, raw
P=-R:G and independent Cd diffusion remain unchanged. So do positive AWF nu,
the current-K counterfactual cutoff, actual uncapped P_A-E_c, dependent A(K)/H,
outer-only epsilon, shared U/k, regional pressures, Q3/Q1 spaces and open
boundary roles. Cases use nu=6.423173803526447e-5,H=.5,h_s=.03721697831176141,
d0=.0875,top20 and physical domain[0,8]x[d0,20]. They do not inherit the old neutral-ABL
dimensional parameters. No force, restoration, clipping, cap, stabilization,
extra H condition, pressure gauge or constitutive regularization was added.

| Fresh case | Normal cells / Nx / states | q / dt | Accepted T | Attempts / accepted | Endpoint or rejected principal / nu_eff | Governing invocation seconds |
|---|---|---|---:|---:|---:|---:|
| A sign only |43 /1 /1966|8 /.02|.04|3 /2|-.0003461673 at rejected .06|16.8332|
| A half dt |43 /1 /1966|8 /.01|.08|8 /8|+.0003279840|33.1484|
| B zero traction |23 /1 /1062|12 /.01|.08|8 /8|+.0047438908|26.5340|
| B frozen seed traction |23 /1 /1062|12 /.01|.08|8 /8|+.0046755680|40.2470|
| B zero traction |23 /3 /2630|12 /.01|.08|8 /8|+.0045578544|69.8859|
| B frozen seed traction |23 /3 /2630|12 /.01|.08|8 /8|+.0047832443|84.1750|

The base Nx2 zero-step native count is1846. Refined Nx2 would require3418
and was not evolved; the3000 guard remains. Nx3 uses exact knots
0,8/3,16/3,8. A's half-dt companion was needed because the sign-corrected
coarse run still exits at the same time as historical H2, whose half-dt run
also reached .08. There was no additional halving or unchanged failed retry.
A and B differ in pressure profile, normal mesh, timestep and quadrature;
their difference is not attributed to one factor.

The new A failure is **interior**, x=3,z=1.1947077409162705, geometry24:
Gxz=.17325767788115604, divU=.0009831072015387939,
k=.010339055157814254, epsilon=.0005216044058824632,
N=4.085694338106877, F_N=31.224348232636256, root residual6.94e-18,
nu_eff=.017640160477040822 and nu_d=.017866786359557348.
Its covariance minimum .0028894757711102987 and Praw=.0005231069629597746
are positive, yet the principal form is negative. The first hard exit was
retained before commit with predecessor, history and failed state. The sign
fix removes a real initial normal-traction error and changes the limiting
state substantially; it does not establish constitutive ellipticity.

B isolates outlet and horizontal factors at common times .00,.01,...,.08.
The frozen outlet is the complete continuous seed traction, separately in
each region, with no second gauge or pressure Dirichlet condition. Its actual
native expression values agree with independent seed tensors within2e-11.
The same initial differential data, scalar data, dt, normal mesh and q_star
are held fixed within each comparison. No regional regauging was applied.
Pointwise outlet defects at T=.08 fall from .0049270 to .00075264 for Nx1
and .0045412 to .00063356 for Nx3 when frozen seed traction is applied.
These are represented natural-boundary defects, not strong constraints.

| Paired T=.08 comparison, all common owned samples | max abs delta Gxz | max abs delta(principal/nu_eff) | max abs delta Ux |
|---|---:|---:|---:|
| Outlet, Nx1 |.0146148|.0276098|.000532046|
| Outlet, Nx3 |.0400937|.0619923|.00123685|
| Nx1/Nx3, zero traction |.0391762|.0608397|.000723747|
| Nx1/Nx3, frozen traction |.0136973|.0264572|.000445354|

These extrema can occur at different points. They do not bound the error of
the minimum by a bulk field percentage. Full
[paired local records](wj-h3-linux/paired-local.csv) include k, epsilon, pi,
normal strain/divergence, both normal owners, horizontal-face owners and
intersection ownership. The positive B minima show short-run admissibility
on the declared inventory; substantial local gradient sensitivity remains.

Weak continuity and represented strong divergence are separate. At T=.08 the
outer continuous-Q1 projected divergence L2 norms are6.80e-17–2.45e-16, while
the complements have norms .00358486–.00377506. Sampled strong divergence is
.0059620–.0078620. Broken-element-Q1 complements are .00258881–.00276273;
they are not the same space as the continuous pressure-test complement.
This offline diagnostic mass/projection algebra never changes the governing
velocity or solve counts. It provides no evidence by itself of an unstable
pressure pair or of a unique cause of the principal exit.

Native readbacks retain actual completed callback histories/time/alpha.
For B Nx3 zero/frozen at .08, uncancelled spatial (u,v,k,full outer epsilon)
mass-dual norms are respectively
(.000831794,.000815486,.048312632,.005931698) and
(.000801984,.000211514,.048304467,.005924930).
Completed-stage norms are instead
(.000067586,.000255217,.000067683,.000013396) and
(.000067585,.000125028,.000067655,.000013402).
Both regional continuity blocks are near2e-17; the gauge is absent.
Neither stage defects nor spatial residuals are described as stationarity.
Readbacks have zero initialization Lift work and zero governing diagnostic solves.

The signed mean wall stress at B endpoints is .00269106–.00269247.
Actual integrated inner production is .0101254–.0101281; current-K reference
production/destruction is .0199591–.0199620. Their reference budget difference
is below2.3e-16, while actual production is deliberately not equated to it.
External mass flux is within4.8e-14. Complete signed mode loads, scalar sources,
diffusion, external reactions and interface readbacks are retained in
[weak budgets](wj-h3-linux/weak-budgets.json),
[interface history](wj-h3-linux/interface-history.csv), and
[native observations](wj-h3-linux/native-observations.json). Pointwise traction
and k-flux defects remain visible; no strong condition was added to remove them.

The configured inventory has75 registrations; exactly these three CTests
were selected, and all passed (12 unittest cases in total):

```
ctest --test-dir "$R/test-config" -j 1 --output-on-failure \
  -R '^(WJH3PressureSeeds|WJH3PassiveAlgebra|WJH2MeasuredProfiles)$'
```

The tests cover the failing old sign, both corrected seed conventions, the
unchanged analytic seed, expression errors, unchanged epsilon formula/labels,
x-knot validation, full3D batch algebra and the adverse rank-one action.
Native fresh runs and strict readbacks exercise the production session path.
A supported changed pressure-seed Function restart fails specifically at
checkpoint-signature validation before governing integration; original CP bytes
remain unchanged. Ordinary accepted CP0 restarts are exercised by every short
trajectory. No historical EP suite or large SVD campaign was repeated.

Executed driver commands are represented below; their exact argv/XML/code/log
hashes, exits and times are in
[final command receipts](wj-h3-linux/closure-v1/final-command-receipts.json).
They were wrapped with `run_wj_h3.py` reservations. Existing immutable output
labels must not be rerun in place.

```
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
R=/media/jay/新加卷1/nektar-v2026-runs/wj-h3-20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export LD_LIBRARY_PATH=/media/jay/新加卷1/nektar-v2026-runs/wj-h1-20261111/build-v1/geometry-lib
python3 -B study_wj_h3.py --root "$R" --mode retained
python3 -B study_wj_h3.py --root "$R" --mode setup
python3 -B study_wj_h3.py --root "$R" --mode select-q
python3 -B study_wj_h3.py --root "$R" --mode controls
python3 -B assess_wj_h3.py --root "$R" --mode analyze
python3 -B assess_wj_h3.py --root "$R" --mode observe
python3 -B assess_wj_h3.py --root "$R" --mode count-and-signature
python3 -B assess_wj_h3.py --root "$R" --mode startup-pair
python3 -B publish_wj_h3.py --root "$R"
```

There were43 scheduled governing attempts:42 accepted and one rejected.
Each executed exactly one U predictor, pressure, k and epsilon solve and one
bound call: counts(43,43,43,43,43). No material correction occurred; all43 raw/
corrected k snapshots are byte-identical. The ledger conservatively charges44
units; the unexecuted fourth coarse-A attempt is not reclaimed. No four-attempt
native-code regression reserve was needed. The closed conservative charge is
1945.004793s, including setup120s and closure300s, against3600s. Nested parent/
child times are conservatively both charged. Peak recorded RSS is3090665472
bytes (2.878GiB); all consuming source/output roots total about2.157GB, below
16GiB. No resource ceiling was crossed. Build time is0; configuration and
observer/plotting/failed setup costs are included. The failed first expression
test, incorrect test cwd, CMake quoting error and intended native A rejection
remain in the command receipts; none was overwritten as a PASS.

[The original publication manifest](wj-h3-linux/manifest.json) binds21 compact
files and its heavy index. Linux closure reverified all356 indexed heavy files
and57 new initial/accepted packet bindings. The
[closure bridge](wj-h3-linux/closure-v1/seal.json) retains the original manifest,
supplement, final ledger and command identities without rewriting it.
[Mac transport verification](wj-h3-linux/mac-transport-verification.json) checks
compact bytes independently; heavy field/checkpoint bytes were verified on
Linux, not copied or claimed independently read on Mac. CSV CRLF bytes are
preserved with local `.gitattributes`. Heavy fields/checkpoints remain outside Git.

The first limiting H3 gate is A's in-plane principal failure, not an unresolved
integration gate. Initializer correction and quadrature are PASS; A dt=.02 is
OBSTRUCTED; A dt=.01 and B are bounded admissible transients; spatial/temporal
accuracy and developed calibration remain UNRESOLVED/NOT_QUALIFIED. The observed
outlet and mesh effects do not identify a unique constitutive/numerical cause.
T=.08 means T*U_infinity/H_h=.08 and only .01 flow-throughs of this strip.

A later package should first test the **same balanced/q12 refined43-cell Nx1**
problem at .01/.005 to the same .08 horizon, then a fixed base23/Nx3 .01/.005
pair, retaining current hard predicates and checking limiting-point gradients/
principal changes, not only2% field norms. This isolates the still-uncombined
normal/seed/timestep effects before increasing duration. Zero complete traction
is the physical candidate; frozen seed traction remains a startup diagnostic.
An eventual corrected finite-horizon T=.4 pair would require40/80 steps;
at the measured Nx3 cost up to10.52s/step, a conservative25% headroom plus
600s observation/closure is about2180s. This is a proposal, not authorization
or a launched extension; it requires resolved local-error/admission gates first.

Developed rough-flat calibration subsequently needs an upstream finite-depth
inlet whose flat evolution reaches Fig.3 at the true virtual crest x=0,
for example a separately qualified [-4,12] fetch/outlet model with the physical
ceiling retained. It must not silently move Fig.3 to x=-4. The downstream
buffer and scalar-backflow model need explicit checks. Current important
refinements exceed capacity: refined43/Nx2=3418 and base23/Nx4=3414 exceed3000.
The present two-level gradient sensitivity therefore motivates separately
qualified mapped-route capacity/operator work before a resolved developed run.
Even T=16 at .005 would require3200 sweeps, roughly9.4h at the measured small
case rate before refinement/observation costs; a short dense prototype is not
production feasibility. No new pressure pair or stabilization is justified by
these data. No longer proposal, curved physical run, hill fit, EP, sparse/MPI
work, exact equilibrium, global stability/realizability or stationarity claim
is included in H3.
