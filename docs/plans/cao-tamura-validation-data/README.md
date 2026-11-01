# Cao-Tamura 2006: frozen experimental data v1.1

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent; PDF visual audit
- Origin Mode: CV1 source audit, extraction and data freeze
- Origin Date: 2026-11-01
- Verification Status: source/data checks passed; no CFD or experimental-accuracy claim
- Version Label: cao-tamura-2006-v1.1

This is the CV1 dataset for the
[calibrated validation plan](../cao-tamura-calibrated-validation.md).
It contains **98 plotted series, 1,151 digitized points and 39 reported scalar
facts** from Cao and Tamura (2006),
[doi:10.1016/j.jweia.2005.10.001](https://doi.org/10.1016/j.jweia.2005.10.001).
The [publisher record](https://www.sciencedirect.com/science/article/abs/pii/S0167610505001108)
confirms the article identity. These are transcribed numerical facts and
derived plot coordinates, not author-native raw observations or time series.
Repeated plots are not independent experiments.

The source is the user-supplied 497,768-byte, 19-page PDF in `docs/literature/`,
identified by SHA256
`adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df`.
PDF pages and printed article pages coincide. The PDF, original page renders,
overlays and the adjacent error-prone OCR Markdown are not redistributed here.
No data from the separate *Wind-Tunnel Experiments* paper, later low-hill
experiments, or other authors' comparison marks are merged into this dataset.

## Files and roles

| File | Contents |
|---|---|
| [source_manifest.json](source_manifest.json) | Citation, fingerprint, facility/geometry, coordinate definitions, dimensional convention, unresolved facts and source discrepancies |
| [figure_inventory.csv](figure_inventory.csv) | All fourteen figures, pages, axes, use/exclusion and qualifications |
| [reported_scalars.csv](reported_scalars.csv) | Geometry, inflow and held-out scalar targets, with source locator and reported/assumed/inferred status |
| [series.json](series.json) | Each plotted series' quantity, case, station, coordinate scale, role and point count |
| [profiles.csv](profiles.csv) | Published axis values, physical-coordinate conversions, extraction envelopes, source symbol identity and fixed height weights |
| [extraction.json](extraction.json) | Audited tick/bar anchors, log-axis flags, per-station offsets and selected drawing ranges |
| [scores.json](scores.json) | Frozen calibration/holdout split, primary/secondary endpoints, normalization, coverage policy and engineering targets |
| [uncertainty.json](uncertainty.json) | Extraction bounds and separately identified missing experimental/facility uncertainties |
| [snapshot.json](snapshot.json) | SHA256 of the nine versioned data/definition files, including this README |

The executable [data auditor](../../../solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_data.py)
generates the three extraction products and verifies the freeze. Compact
[CV1 evidence](../cao-tamura-validation-evidence/cv1.json) records the executed
checks, source consistency, tools, commands and source hashes.

The role split is six calibration series (Fig. 3a/b/c, S and R), two reserved
flat-flow cross-check series (Fig. 3d), 50 validation series, 22 duplicate
cross-check series, fourteen wake-mean cross-check series and four dimensional
holdout series. The reported flat-floor scalar constraints accompany the
calibration. Fig. 3d is reserved as the required unused flat-profile check;
Table 1's vertical-rms ratios accompany it. Hill data never enter calibration.

## Extraction and checks

The plots contain vector symbols. The primary extraction uses pdfplumber:
circle centres come from their bounding boxes; triangle centres are the
three-vertex centroid. Taking the triangle bounding-box midpoint would bias
the plotted height. Only the filled symbol path is selected; its outline is
not another point. Source drawing indices and anchors are bound to the PDF
hash. Exact repeated drawing layers in Figs. 8 and 10 are not counted twice.
Connecting curves, legends, coloured dividing/peak guides and Raupach's Fig.
3c comparison crosses are excluded from observations.

A separate pypdf content-stream interpreter tracks transformations, fill
colour and paint operators without using pdfplumber's geometry. It agrees on
all 1,151 selected centres, shapes and fills (maximum centre difference zero
at stored source precision).

The initial unpublished v1 draft misclassified four Fig. 6a upwind/downwind
series as duplicate crest curves. Final CV1 review corrected their roles in
v1.1, retaining the draft snapshot locally; all numerical point values and
acceptance bounds are unchanged, and no fitting or CFD had begun.

A Poppler 300 dpi raster pass estimates centres from opposite outer-ink edges
of 472 isolated circle symbols: maximum
distance 0.179 PDF points, below the unchanged 0.5-point extraction envelope.
Crowded symbols and triangles are covered by the complete vector comparison
and visual page overlays, rather than claiming independent raster centres
for obscured marks. All seven plotted pages (7, 9, 11-15) were visually checked.
The method uses independent software backends; it is not two independent
human digitizers or two independent experiments.

Eighteen symbols encoded underneath white legend boxes are retained with
`quality=legend_occluded` and **weight zero**. They are not scored or fitted.
The remaining 1,133 points are marked visible; proximity/overlap with another
symbol still limits a blind raster reading, as the receipt records.

The fixed extraction envelope is **0.5 PDF points per plotted axis**, transformed
through the recorded axis map. Logarithmic abscissae use the larger of the two
transformed deviations. This is a conservative digitization envelope, not a
measurement standard deviation or 95% interval. Axis semantics, physical scale
and source inconsistencies require separate treatment even when a vector
centre is known precisely.

Repeated-figure checks use physical global height and linear interpolation
only over the overlapping measured range. Fig. 4 versus Fig. 8 maximum
differences are 0.0104 (SS) and 0.00761 (RR) in U/Uinf; Fig. 5 versus Fig. 10
differences are 0.00185 (SS) and 0.00106 (RR) in sigma_u/Uinf. All fit within
the recorded combined extraction comparison bounds. The four Fig. 6b markers
at yprime/Hh=0.25 agree with the text to within 0.0011 in S. These checks test
transcription and axis offsets, not CFD accuracy.

## Coordinate and column contract

Cases S/R denote undisturbed smooth/rough references; SS/RR denote corresponding
inflow and hill surfaces. RS is rough inflow with a smooth hill; SR is smooth
inflow with a rough hill. Paper x is streamwise and y is global vertical;
`yprime=y-yh(x)` is **vertical above the physical base surface**, not the
closest-normal coordinate used in an AWF strip. Hh=0.04 m; the crest is x=0.
Solver 3D vertical is z, so the paper's vertical variance/covariance involves
R_zz/R_xz there, not the spanwise R_yy/R_xy.

| Figure | Published height | Quantity and interpretation |
|---|---|---|
| 3a | y in mm | U/Uinf |
| 3b | y/delta | -uv/Uinf^2; use labelled 0.5/1 ordinate ticks, not frame bottom as zero |
| 3c/d | (y-d)/delta | Longitudinal/vertical sigma divided by u_tau |
| 3e | log10(y-d in mm) on horizontal axis | U/u_tau on vertical axis; duplicate consistency check |
| 4, 5, 8, 10 | **Global y/Hh** | Offset U/Uinf or sigma_u/Uinf profiles; use the per-station zero ticks and explicit scale bar |
| 6 | yprime/Hh | Fractional speed-up S |
| 7 | yprime/delta | Same-vertical-height ratios to the matching flat reference |
| 9, 12 | Fixed physical probe height | U or sigma_u in m/s against x/Hh; smooth 2 mm, rough 7 mm above floor |
| 11 | y/delta75 | -uv/Uinf^2 at x/Hh=7.5; numerical delta75 is unspecified |

The nine Fig. 4/5/10 stations are
`[-2.5,-1.25,0,1.25,2.5,3.75,5,6.25,7.5] Hh`.
Fig. 8 covers `[0,1.25,2.5,3.75,5,6.25,7.5] Hh`, with the source-label
qualifications below. Every retained height, rather than an invented uniform
probe array, is enumerated in the CSV. Multi-profile origins use actual station
ticks: the nearby printed `0` text and broken reference profiles are not
sufficiently precise to define each ordinate's velocity origin.

`axis_x/axis_y` retain the published axes; `value` is the physical plotted
quantity named in `series.json` (for Fig. 3e it is axis_y, otherwise the
profile abscissa; traverses use axis_y). `px/py` are PDF points with origin at
the page top-left. `source_curve` is a zero-based pdfplumber curve index.
`radius_x/radius_y` support the rendering audit, not physical error estimates.
`axis_*_bound` and `value_bound` are extraction-only envelopes.

`x_over_H`, `y_over_H` and `yprime_over_H` identify physical probe coordinates.
Conversions involving delta use **nominal** delta=0.25 m and rough d=0.0035 m;
their status is marked in `series.json`. They do not turn approximate/assumed
source parameters into exact measurements. Fig. 11's physical height columns
are deliberately empty because delta75 is not numerically specified. Its
published normalized coordinate remains usable with the declared conditional
shape-comparison convention in `scores.json`.

`weight` is a per-series normalized Voronoi-height weight on the visible
published coordinate (streamwise position for traverses). Equal-height
overplots share a cell; they do not gain independent weight. Independent
coordinates are grouped only at six decimal places to avoid merging distinct
low-height observations. An unavailable CFD observation cannot silently be
removed and the rest reweighted to pass: report the omitted weight/range and
label any subset result conditional. Below-interface observations require
the CV2 reconstruction/availability decision.

## Decisions and unresolved source issues

- **Dimensional scale:** Fig. 3a gives U_H/Uinf about 0.79351 (S), 0.68670 (R).
  With Re_H=12000/10500 these imply Uinf/nu about 378067/382261 per metre.
  Rough Re_s=3.97 instead implies 389216 per metre; the Re_H-derived rough
  Re_s is 3.899, about 1.8% lower. Treat the paper's Reynolds numbers as
  approximate; do not manufacture exact simultaneous agreement. The frozen
  analysis convention uses Hh and the case's U_H as scales, nu*=1/Re_H.
  For illustration only, assuming air nu=1.5e-5 m^2/s gives Uinf about
  5.671/5.734 m/s; neither speed nor that viscosity is reported by the paper.
  Do not calibrate a dimensional scale using held-out hill amplitudes.
- **Fig. 3e normalization:** applying reported u_tau/Uinf to the repeated
  log-axis plot differs from Fig. 3a by up to 0.05014 (S), 0.04351 (R) U/Uinf.
  It is a source consistency discrepancy, beyond subpixel extraction error.
  Fig. 3a/b/c take calibration priority; Fig. 3e is excluded from fitting.
  Do not refit u_tau just to remove the discrepancy.
- **Reattachment:** the quoted crest-referenced x_R/Hh=5.4 (SS), 6.5 (RR) is
  obtained by ground extrapolation of dividing points. Finite-height Fig. 9
  zero crossings and wall-shear zeros are distinct observations. The source
  does not specify its extrapolation formula. The frozen comparison requires
  a declared tail method and sensitivity, not a fitted extraction method that
  forces the quoted values. A 10 mm traverse spacing is a resolution indicator,
  not a complete uncertainty interval on x_R.
- **Source labels/ranges:** Fig. 8 literally labels stations `2.53` and `3.0`;
  2.5 and 5.0 are inferred from the regular station pattern, hill geometry and
  Figs. 4/10. Fig. 8 stays a qualified cross-check. Text says the near-floor
  traverse ends at 320 mm (8 Hh), while Fig. 9 markers extend to about 8.75 Hh
  and Fig. 12 to about 8.5 Hh. Both source statements are preserved. Table 1
  prints sigma_tau in its denominator headers; interpreting this as u_tau is
  supported by the surrounding text and Fig. 3 labels. Fig. 10's crowded
  denominator typography is interpreted using its text and repeated Fig. 5.
- **Facility:** the reported 5% blockage equals Hh/roof height; frontal area
  over the full 1 m tunnel width with a 0.8 m hill gives 4%. End plates may
  define an effective test width, but their missing dimensions cannot be
  filled in as facts. A central-plane/periodic-span proxy needs facility
  sensitivities and cannot be called a fully specified tunnel reproduction.
- **Mixed cases:** RS/SR identify upstream and hill conditions, not exact
  downstream coverage or material junctions. CV7 must bracket a hill-only
  surface patch with retained surroundings and an upstream-foot change
  persisting downstream, fixed before examining mixed-case scores.
- **Missing turbulence/error data:** spanwise variance, epsilon, integral
  scales, raw time series and a complete uncertainty budget are absent.
  Table 1's near-surface ratios have no unique measurement height. Do not
  convert them into pointwise boundary conditions or infer k uniquely from
  two variances. The ±1.5% probe pitch-voltage calibration statement is not
  a 95% velocity/stress uncertainty. No author-native data were supplied or
  found in the inspected paper/publisher preview; this is a bounded search
  statement, not proof that such data do not exist elsewhere.

The source measurement uncertainty remains unquantified. Source ambiguities
are recorded for later sensitivity/coverage decisions, not declared solved
by the data freeze. This satisfies CV1's explicit requirement to record
unresolved facility, units and mixed-surface details. It does not close CV2's
model/reference/observation work or validate the turbulence models.

## Reproduce and preserve the freeze

The standard-library check can run without the PDF:

```bash
python3 solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_data.py \
  --data docs/plans/cao-tamura-validation-data --report /tmp/ct2006-check.json
```

For source reconstruction, use Python with `pdfplumber`, `pypdf`, `Pillow`
and the `pdftoppm` executable (Poppler). CV1 used the bundled Mac runtime;
the exact versions and commands are in the evidence. No Nektar build, MPI
job or CFD result is involved in this source audit.

```bash
CV1_PY=/Users/jay/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
CV1_PDF="docs/literature/Experimental study on roughness effects on turbulent boundary layer flow over a two-dimensional steep hill.pdf"
"$CV1_PY" solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_data.py \
  --data docs/plans/cao-tamura-validation-data --pdf "$CV1_PDF" \
  --raster /tmp/ct2006-raster-check --report /tmp/ct2006-source-check.json
"$CV1_PY" solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_data.py \
  --data /tmp/ct2006-fresh-extraction --pdf "$CV1_PDF" --extract \
  --report /tmp/ct2006-fresh-check.json
```

Compare the fresh `profiles.csv`, `series.json` and `extraction.json` byte for
byte against this directory. Fresh extraction never overwrites a sealed
snapshot. The nine definition/data files are hashed in `snapshot.json`;
alteration is rejected by ordinary verification. Any later correction to
axes, source facts or scores requires a new data version, an explicit reason
and reassessment of affected calibration/validation. Add future evidence to
the sibling evidence directory, not by changing frozen v1.1 tables.
