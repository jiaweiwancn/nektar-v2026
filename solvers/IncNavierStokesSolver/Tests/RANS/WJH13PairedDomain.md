# WJ-H13 paired domain and full-height neutral reference

H13 prepares the actual `[-10,15] × [effective wall,20]` target. It does not
run the old `[0,8]` strip again. The user's selected inlet is the full-height
neutral WJ reference, with the calibrated finite-molecular-viscosity Budget
wall column below `H=.15`; it supersedes H13 finite-depth/ambient designs.

`wj_h13_cases.py` reuses the frozen H5 reference and H6 bisection, and accepts
no new wall height, intercept, roughness or mesh search. The H5 generator's
default remains the old single-column case; its new optional `xknots` argument
is used only by this paired generator. The complete target has 14 horizontal
columns, 16 inner and 34 outer cells, Q3 differential fields and regional Q1
pressures. `wj_h13_domain.py` counts entities and matrix payload sizes using
integers/sets, without assembling global matrices.

For reproduction on Linux, set the six numerical-library thread controls to
one before importing Python. With the RANS test directory on `PYTHONPATH`:

```bash
python3 -B wj_h13_cases.py --source /path/to/matching/source --output /new/offline/case
python3 -B wj_h13_cases.py --source /path/to/matching/source --preflight /new/offline/case/paired-domain.json
```

The second command is expected to return **2 / OBSTRUCTED**. This external
preflight validates the exact flat XML, shared topology and source capacity
predicates before any native execution. It never calls a native executable.
`NumSteps=0` alone provides no safety against dense constructor/projection
allocations. Do not run the generated 24,688-state session on the current
mapped dense route (capacity 3000; the 4331 exception excludes mapped cases).
Moving a session requires an explicit path alias while preserving its bytes;
the checked preflight function accepts a supplied alias path without rewriting
the original manifest. A changed native capacity/source requires fresh review.

Artifacts include:

- `flat/case.xml`: complete physical Functions, `FlatInletTraceV1`,
  `ReferenceTractionV1`, zero steps, and an explicit NOT_RUN warning.
- `paired-domain.json`: shared chart/connectivity/roles, actual binary64 knot
  vectors, physical reference, separate map identities, counts and capacity.
- `flat-geometry.json`, `terrain-geometry.json`: physical vertex specializations
  of the same chart, with the normal displaced wall and independent outer map.
- `terrain-chart-topology.xml`: **geometry only**, without Conditions or a
  fabricated runnable curved OpenTunnel boundary problem.
- `observations.json`: Fig.3 remains at its actual hill-absent virtual crest;
  original-floor 7 mm observations remain inner AWF predictions.
- `preservation-protocol.json`: distinct representation/drift/reference/fetch
  errors, physical norm/scalar scopes, stations, transit and diffusion clocks,
  fresh CP0 requirements and later admission/duration/sensitivity phases.
- `mapped-backend-prerequisite.json`: source-located dense call graph and the
  smallest complete execution-route conversion needed for target admission.

The top uses fixed reference outward epsilon flux `+Cd*K²/(sigma_e*rT)` and
tangential traction `+u_tau²`. It is not a prescribed top epsilon value,
fixed gradient, dynamic Robin law or SlipV1. Regional complete outlet tractions
fix pressure; no second gauge or interface pressure equation is added. Physical
gradients use `r=physical_z-d0`, not the outer chart coordinate. Inner molecular
traction is included; individual Reynolds normal stresses need not match.

The next implementation is a separately qualified mapped block/operator path
covering **startup and observation as well as governing solves**. Preserve
consistent mass, full nonsymmetric `D*(alpha*M_free)^-1*G-C`, actual pressure
actions and restrictions, k-then-epsilon, dependent epsilon identity rows and
columns, current-K cutoff, all bounds/predicates, rollback and history. Qualify
local actions/solves against existing small dense mapped controls and adverse
witnesses before a fresh target CP0. A session selector or guard increase alone
cannot remove the dense allocations. See the H13 evidence and machine-readable
plan for exact source locations and block dimensions.

The full-height logarithm reaches about `1.09003 U_infinity` at 800 mm. It is
an explicitly selected idealized inflow, not a finite-depth tunnel fit. No
matching height is physically optimal here. Putting every physical 7 mm stress
target outside the AWF remains incompatible with the retained roughness inverse.
No native preservation, duration, temporal, hill or experimental PASS follows
from these artifacts.
