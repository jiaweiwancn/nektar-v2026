# Frozen B3 true linear residual versus correction size

This is the fixed-base follow-up to the [B3 direction accuracy and curvature study](wj-h27-c1-b3-direction-accuracy-curvature.md).
The sweep used the 5,238-cell accepted B3 field at hill height fraction .0024 with coefficient digest
`ced15ecef5037f0e638664999a6ec82f46e02585a95895f907aa14ed3704c600`.
The native physical residual was 6.438421472031332e-6 and the row-scaled L2
residual was 1.4277267361975784e-4. Both original 0.3 forcing thresholds were
active, 1.9315264416093993e-6 for the physical block maximum and
4.283180208592735e-5 for the row-scaled L2 norm.

| Scaled correction infinity norm | Fresh physical error | Fresh scaled L2 error | Measured result |
| ---: | ---: | ---: | --- |
| 0.00020533 | 6.36946702e-6 | 1.41249871e-4 | Both exceed target |
| 0.025 | 4.69907260e-6 | 1.04403727e-4 | Both exceed target |
| 0.1 | 3.88823685e-6 | 8.64088808e-5 | Both exceed target |
| 1,000 | 3.66203163e-6 | 8.12114307e-5 | Both exceed target |
| 5,000 | 2.75624002e-6 | 6.11435824e-5 | Both exceed target |
| 8,500 | 1.97033484e-6 | 4.37444991e-5 | Both exceed target |
| 8,750 | 1.91489082e-6 | 4.25113269e-5 | Smallest sampled pass at the primary spacing with the existing action check |
| 9,250 | 1.80344586e-6 | 4.00446526e-5 | Both pass at the primary and half spacings |
| 10,000 | 1.63739992e-6 | 3.64220977e-5 | Smallest sampled pass at all three spacings with the existing action check |
| 13,000 | 9.92633308e-7 | 2.21737770e-5 | Smallest sampled pass with the empirical allowance added to each error |

Each row was evaluated at an unchanged base using the unshifted current
action and the original block norms. Sizes use the stored unknown scaling
`||D_x^-1 d||_infinity`. The candidates minimize the larger forcing-ratio
prediction in one fixed rank-16 Krylov subspace under the indicated size
cap. A fresh native application on the resulting vector supplies each
reported error. The same factor built the whole subspace. Its colored FD
spacing was 4e-6 and its inverse-only mass was 0.009539982703780335.
Correction mass, correction diffusion, and residual regularization were
zero. The current hybrid action used split quadrature 10 by 16.

All 105 native actions used central stencils at the three prescribed
relative spacings 512e-6, 256e-6 and 1024e-6 with zero feasibility reductions.
The action check compares twice the physical and scaled primary/half
spacing differences to half of the corresponding forcing target. At 8,750,
the two allowances were 5.53601968e-7 and 8.56985879e-6. The half-spacing
physical error was 1.93504748e-6, slightly above target. At 13,000 the
error plus allowance was 1.90774034e-6 physical and 3.60504018e-5 scaled.
These are empirical spacing allowances.

The lowest passing checked size among the sampled scalings of the archived
R9 vector was 12,278.482. Its repeated raw actions were bitwise equal.
On the 8,750 candidate the pressure and k block sizes were 8,750 and
7,974.32. This finite subspace and archived-vector experiment advanced zero
physical time. All final packets retain the exact coefficient digest.

## Identities and cost

The isolated tested source revision was
`3349353ad247ef982e47fe9dcc8871bc8fb08057`. The Linux native executable
SHA-256 was `2d4f238e7ea00c1304e78551fb3adaa3e3f538aea57edc21ba0d8a3ec4eb6f0c`.
The authoritative directory is

```text
/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-frozen-size.pnvR6Uif
```

The cold basis run cost 2,033.91 wall seconds, including 1,488.91 s of
colored assembly and 125.44 s of factorization, with peak RSS 7.89 GiB.
Fresh action replay cost 474.81 s for the initial sweep, 218.68 s for
the crossing check, 135.90 s for the allowance check, and 325.55 s for
the archived direction. These process durations overlap. The campaign
elapsed time through evidence collection was 3,511 s.

[Numerical results](wj-h27-linux/b3-frozen-size/combined-results.csv) ·
[Residual curves](wj-h27-linux/b3-frozen-size/frozen-size.png) ·
[Machine-readable action summaries](wj-h27-linux/b3-frozen-size/combined-results.json) ·
[Input and run identities](wj-h27-linux/b3-frozen-size/run-evidence.json)

## Reproduction

The focused Linux `WJH27SteadyKrylovProbe` verifies the basis and two-direction
batch on a known scaled linear action and checks that accepted coefficients
stay at the seed. The optimizer control checks a constrained problem with a
known solution. `WJDiagnosticFrozenSizeColumns` selects the opt-in basis
capture. `WJDiagnosticActionProbeCount` selects the digest-bound batch
evaluation. The experiment retains the original equations and diagnostic
state custody. B3 original-root qualification remains open.

From a Linux Git checkout containing this report, copy the tested source
bundle and the [replay scripts](wj-h27-linux/b3-frozen-size/) into a fresh run
directory. The scripts use the archived R9 input and existing Linux
installation. The build uses eight workers and each solver uses one rank
with single-threaded numerical libraries.

```bash
campaign_runs='/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs'
retained="$campaign_runs/wj-c1-b3-frozen-size.pnvR6Uif"
helpers="$(git rev-parse --show-toplevel)/docs/plans/wj-h27-linux/b3-frozen-size"
replay="$(mktemp -d "$campaign_runs/wj-c1-b3-frozen-size-replay.XXXXXXXX")"
mkdir "$replay/tmp"
cp "$helpers"/*.sh "$helpers"/*.py "$helpers/source-bundle.sha256" "$replay/"
cp "$retained/source.bundle" "$replay/"
git clone --shared --no-checkout "$retained/source" "$replay/source"
cd "$replay"
sha256sum -c source-bundle.sha256
bash workers.sh build-workers.ps 62
bash build.sh
bash basis.sh
bash verify.sh
bash refine.sh
bash margin.sh
```

Run the archived-direction control and regenerate the combined evidence in
the same replay directory.

```bash
r="$PWD"
source threads.sh
export PYTHONPATH="$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27"
python3 frozen_size.py self-test
archive='/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/runs/wj-c1-b3-dir-curvature.SfYckSXM/hybrid-fine'
mkdir ray
bash workers.sh ray-workers.ps 68
python3 ray.py "$r/ray" "$archive/state.action-merit-unresolved-base-1-proposal-0-sequence-3.json"
bash verify-ray.sh "$r/ray"
python3 combine.py "$r"
python3 plot.py combined-results.json frozen-size.png --ray ray/frozen-size-results.json
```

The retained Linux directory stores the raw vectors, coefficient packets,
and full logs. The input manifest and each action receipt record their
identities. The Linux source bundle serves the exact numerical revision
listed above. The integration and evidence commits add the published
reports and replay material to that tested numerical source.
