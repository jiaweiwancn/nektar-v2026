#!/usr/bin/env python3
"""Collect the immutable frozen-size study summaries."""
import csv
import json
from pathlib import Path
import sys

from frozen_size import pin, write_json

root = Path(sys.argv[1])
sources = [root/name/"frozen-size-results.json" for name in (".", "refine", "margin", "ray")]
study, refinements, margin, archived = [json.loads(path.read_text()) for path in sources]
keys = ("targets", "physical_time_advanced")
for other in (refinements, margin, archived):
    assert all(study[k] == other[k] for k in keys)
    for k in ("base_coefficient_sha256", "source_identity", "equation_identity", "free_id_sha256", "current_fd"):
        assert study["numerical_identity"][k] == other["numerical_identity"][k]
joint = sorted(study["candidates"] + refinements["candidates"] + margin["candidates"], key=lambda row: row["native_scaled_infinity"])
report = {
    "schema":"WJH27FrozenCorrectionSizeCombinedV1", "sources": [pin(path) for path in sources],
    "targets": study["targets"], "basis_rank":study["basis_rank"], "numerical_identity":study["numerical_identity"],
    "physical_time_advanced":0., "repeat_archived_raw_equal":archived["repeat_archived_raw_equal"],
    "smallest_nominal_pass":min((row for row in joint if row["both_targets_pass"]),
                               key=lambda row:row["native_scaled_infinity"], default=None),
    "smallest_verified_size":min((row for row in joint if row["pass_and_resolved"]),
                                key=lambda row:row["native_scaled_infinity"], default=None),
    "smallest_pass_all_three_spacings":min((row for row in joint if row["action_resolved"] and
                                             all(e is not None and all(a <= b for a, b in zip(e, study["targets"]))
                                                 for e in row["residuals"])),
                                      key=lambda row:row["native_scaled_infinity"], default=None),
    "smallest_pass_with_allowance":min((row for row in joint if row["both_targets_with_allowance"]),
                                     key=lambda row:row["native_scaled_infinity"], default=None),
    "archived_smallest_verified_size":archived["smallest_verified_size"],
    "candidates":joint + archived["candidates"],
}
write_json(root/"combined-results.json",report)
with (root/"combined-results.csv").open("x") as file:
    writer = csv.writer(file)
    for idx, source in enumerate(sources):
        with source.with_suffix(".csv").open() as data:
            reader = csv.reader(data)
            header = next(reader)
            if idx == 0:
                writer.writerow(header)
            writer.writerows(reader)
