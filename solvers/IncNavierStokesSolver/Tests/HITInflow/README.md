# HIT inflow contracts and experimental fixtures

This directory contains the Milestone 1 configuration examples and the principal
5.08 cm grid dataset from Comte-Bellot and Corrsin (1971),
[Journal of Fluid Mechanics 48, 273-337](https://doi.org/10.1017/S0022112071001599).
The inlet implementation and solver regressions belong to later milestones.
Neither XML fragment includes a mesh or a complete boundary definition, and the
current solver does not implement `HITInflow`.

The normative interface and integration contracts are in
[the Milestone 1 specification](../../../../docs/plans/hit-inflow-milestone-1.md).
The original project sequence is in
[the implementation plan](../../../../docs/plans/hit-inflow-implementation-plan.md).

## Source and units

`cbc1971-source.json` preserves the checked source numbers as decimal strings
in centimetre units, including their printed precision and missing entries.
The complete Table 2(a), Table 3, and principal Table 4 rows were checked
visually against the local original PDF on 2026-10-28. Their printed pages are
296, 298, and 299 (PDF pages 24, 26, and 27). The source record contains that
PDF's SHA-256 digest. The source PDFs are supplied separately under
`docs/literature/` and are not required to regenerate the numeric fixtures.

The measurement-region speed is 12.7 m/s, from the footnote on printed page
289 (PDF page 17). The approach speed is 10 m/s. Stations use
`theta = U0*t/M`, with travel time through the upstream contraction, rather
than simply `x/M`. Taking station 42 as the computational reference gives
station offsets of 3.612896 m and 8.322564 m at stations 98 and 171.

Wavenumbers are angular wavenumbers. Multiply source `cm^-1` by 100 to obtain
`rad/m`. Multiply either spectrum in `cm^3/s^2` by `1e-6` to obtain `m^3/s^2`.
No extra factor of `2*pi` belongs in this conversion.

## Numeric files

| Files | Meaning | Role |
| --- | --- | --- |
| `cbc1971_theta42_spectrum_3d.dat` | Table 3 shell spectrum, 19 points | Experimental inlet input |
| `cbc1971_theta98_spectrum_3d.dat` | Table 3 shell spectrum, 19 points | Downstream validation |
| `cbc1971_theta171_spectrum_3d.dat` | Table 3 shell spectrum, 18 points | Downstream validation |
| `cbc1971_theta*_longitudinal_1d.dat` | Table 2(a) one-sided longitudinal spectra, 22 points each | Virtual-probe validation |
| `cbc1971_stations.dat` | Table 4 converted to SI, with relative times and distances | Gross statistical targets |

The first comment in each spectrum declares its kind and units. Both kinds have
two numeric columns, but their meanings differ. The future inlet loader must
accept `KIND=Shell3D` and reject `KIND=Longitudinal1D`. Shell energy integrates to
one half of total velocity variance; the one-sided longitudinal spectrum
integrates to streamwise velocity variance.

Missing source entries remain missing: individual station files omit them.
There is no smoothing, extrapolation, taper, or normalization in these files.
The small nonmonotonic sequence in Table 3 at station 98 is printed in the
original and is preserved. The source spectrum has finite coverage, and matching
its measured part does not establish complete energy or dissipation integrals.
The runtime contract makes extensions explicit.

The Table 4 consistent initial viscosity is about `1.5e-5 m^2/s`. Figure 9's
caption instead uses `2*nu=0.28 cm^2/s`; this source discrepancy remains recorded.
`lambda` retains the paper's Taylor transverse microscale convention. The last
station-file column is the dimensionless product `(lambda/Lg)*Re_lambda`.

## Reproduction and checks

From the repository root, using only Python's standard library:

```bash
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/prepare_cbc1971.py
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/prepare_cbc1971.py --check
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/check_contract.py
```

The preparation command rewrites only the seven generated `.dat` files and
`cbc1971-manifest.json`. It uses exact decimal scaling and writes a deterministic
manifest with source, generator, and output SHA-256 digests. `--check` performs
no writes and fails on changed data or digests. `--output-dir` permits comparison
against a separate directory. To verify a particular supplied source copy too:

```bash
python3 solvers/IncNavierStokesSolver/Tests/HITInflow/prepare_cbc1971.py --check \
  --source-pdf 'docs/literature/Simple Eulerian time correlation of full-and narrow-band velocity signals in grid-generated isotropic turbulence.pdf'
```

`check_contract.py` checks these examples and their fixture-level constraints.
It is not the future C++ parser and does not certify a solver run. The numerical
generator, pressure convergence, MPI sampling, and physical comparisons have
their own later milestone checks.
