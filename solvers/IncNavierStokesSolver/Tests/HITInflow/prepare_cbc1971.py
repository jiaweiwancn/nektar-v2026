#!/usr/bin/env python3
"""Reproduce or check the small CBC 1971 fixtures using only the standard library.

This prepares measurements, not extrapolated inflow spectra or solver results.
Run without --check to regenerate the seven .dat files and their manifest.
"""

import argparse
from decimal import Decimal, InvalidOperation, getcontext
import hashlib
import json
from pathlib import Path
import sys


getcontext().prec = 40
HERE = Path(__file__).resolve().parent
SOURCE = HERE / "cbc1971-source.json"
STATIONS = ("42", "98", "171")
D = Decimal


def require(condition, message):
    if not condition:
        raise ValueError(message)


def positive(value):
    number = D(value)
    require(number.is_finite() and number > 0, "Expected a finite positive value")
    return number


def digest(data):
    return hashlib.sha256(data).hexdigest()


def numeric_line(values):
    return " ".join(format(value, ".12E") if value else "0.000000000000E+0"
                    for value in values)


def validate_source(source):
    require(source["schema"] == "nektar.hit.cbc1971.source.v1", "Unknown source schema")
    require(tuple(source["stations"]) == STATIONS, "Unexpected principal stations")
    for table_name, expected_rows in (("table2a", 23), ("table3", 20)):
        rows = source[table_name]["rows"]
        require(len(rows) == expected_rows, "Unexpected row count in " + table_name)
        previous = D(0)
        for row in rows:
            require(len(row) == 4, "Expected four source columns in " + table_name)
            k = positive(row[0])
            require(k > previous, "Source wavenumbers must increase")
            previous = k
            for value in row[1:]:
                if value is not None:
                    positive(value)
    require(len(source["table4"]["rows"]) == 3, "Unexpected Table 4 row count")
    for theta, row in zip(STATIONS, source["table4"]["rows"]):
        require(len(row) == 10 and row[1] == theta, "Incorrect Table 4 layout")
        for value in row:
            positive(value)
    experiment = source["experiment"]
    for key in ("grid_mesh_m", "approach_speed_m_s", "measurement_speed_m_s",
                "reference_theta", "initial_viscosity_m2_s"):
        positive(experiment[key])
    require(experiment["reference_theta"] == "42", "Unexpected reference station")
    require(len(source["source"]["pdf_sha256"]) == 64, "Invalid source PDF digest")


def products(source, source_bytes):
    """Return deterministic output bytes and metadata without writing files."""
    validate_source(source)
    output = {}
    records = {}
    doi = source["source"]["doi"]
    expected_counts = {"table2a": (22, 22, 22), "table3": (19, 19, 18)}
    for table, kind, suffix in (
            ("table2a", "Longitudinal1D", "longitudinal_1d"),
            ("table3", "Shell3D", "spectrum_3d")):
        page = source["source"][table + "_printed_page"]
        for column, theta in enumerate(STATIONS, 1):
            rows = [(positive(row[0]) * D(100), positive(row[column]) * D("1e-6"))
                    for row in source[table]["rows"] if row[column] is not None]
            require(len(rows) == expected_counts[table][column - 1],
                    "Missing-value layout changed in " + table)
            name = "cbc1971_theta{}_{}.dat".format(theta, suffix)
            lines = [
                "# NEKTAR_HIT_SPECTRUM_V1 KIND={} K_UNITS=rad/m E_UNITS=m3/s2".format(kind),
                "# DOI={} TABLE={} PRINTED_PAGE={} THETA={}".format(doi, table, page, theta),
                "# Measured table entries only; no extension, smoothing, or normalization.",
                "# columns: k_rad_per_m spectrum_m3_per_s2",
            ] + [numeric_line(row) for row in rows]
            output[name] = ("\n".join(lines) + "\n").encode("utf-8")
            records[name] = {
                "table": table, "printed_page": page, "theta": int(theta),
                "kind": kind, "rows": len(rows),
                "k_min_rad_per_m": str(rows[0][0]),
                "k_max_rad_per_m": str(rows[-1][0]),
                "role": "inlet" if table == "table3" and theta == "42" else "validation",
            }

    experiment = source["experiment"]
    mesh = positive(experiment["grid_mesh_m"])
    approach = positive(experiment["approach_speed_m_s"])
    convection = positive(experiment["measurement_speed_m_s"])
    reference = positive(experiment["reference_theta"])
    rows = []
    for row in source["table4"]["rows"]:
        values = [positive(value) for value in row]
        require(values[0] * D("0.01") == mesh, "Grid mesh inconsistent with Table 4")
        elapsed = (values[1] - reference) * mesh / approach
        rows.append((values[1], elapsed, convection * elapsed,
                     values[2] * D("0.01"), values[3] * D("0.0001"),
                     values[4] * D("0.01"), values[5] * D("0.01"),
                     values[6] * D("0.01"), values[7] * D("0.01"),
                     values[8], values[9]))
    require(rows[1][2] == D("3.612896") and rows[2][2] == D("8.322564"),
            "Station mapping changed; retain post-contraction convection speed")
    name = "cbc1971_stations.dat"
    lines = [
        "# NEKTAR_HIT_STATIONS_V1 UNITS=SI",
        "# DOI={} TABLE=table4 PRINTED_PAGE=299".format(doi),
        "# theta=U0*t/M; U0=10 m/s; Ubar=12.7 m/s; M=0.0508 m; reference_theta=42",
        "# columns: theta dt_s dx_m urms_m_s epsilon_m2_s3 eta_m lambda_m Lg_m Lf_m Re_lambda lambda_over_Lg_times_Re",
    ] + [numeric_line(row) for row in rows]
    output[name] = ("\n".join(lines) + "\n").encode("utf-8")
    records[name] = {"table": "table4", "printed_page": 299, "rows": 3,
                     "role": "validation", "units": "SI"}
    for name, contents in output.items():
        records[name]["sha256"] = digest(contents)
    manifest = {
        "schema": "nektar.hit.cbc1971.manifest.v1",
        "source_sha256": digest(source_bytes),
        "source_pdf_sha256": source["source"]["pdf_sha256"],
        "generator_sha256": digest(Path(__file__).read_bytes()),
        "source_doi": doi,
        "conversion": {"k_cm_inv_to_rad_per_m": "100", "spectrum_cm3_s2_to_m3_s2": "1e-6"},
        "files": records,
    }
    output["cbc1971-manifest.json"] = (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Compare fixtures without writing")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--source-pdf", type=Path, help="Also verify the supplied local PDF digest")
    args = parser.parse_args()
    try:
        source_bytes = SOURCE.read_bytes()
        source = json.loads(source_bytes)
        if args.source_pdf:
            require(digest(args.source_pdf.read_bytes()) == source["source"]["pdf_sha256"],
                    "Supplied PDF differs from the visually checked source")
        expected = products(source, source_bytes)
        if args.check:
            for name, contents in expected.items():
                require((args.output_dir / name).read_bytes() == contents,
                        "Generated fixture differs: " + name)
            print("PASS: six spectra, station mapping, missing values, SI conversions, and SHA-256 manifest")
        else:
            args.output_dir.mkdir(parents=True, exist_ok=True)
            for name, contents in expected.items():
                (args.output_dir / name).write_bytes(contents)
            print("Wrote seven numeric fixtures and their SHA-256 manifest to " + str(args.output_dir))
    except (ValueError, KeyError, InvalidOperation, OSError) as error:
        parser.exit(1, "error: {}\n".format(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
