#!/usr/bin/env python3
"""Build a factor-free size sweep along the archived R9 direction."""
from pathlib import Path
import hashlib
import re
import sys

import numpy as np

from c1_action_gap import read_audit
from frozen_size import BINDINGS, array_spec, pin, write_json

root = Path(sys.argv[1])
archive = Path(sys.argv[2])
meta, arrays = read_audit(archive)
checks = meta["measurements"]
assert checks["forcing"] == .3 and meta["fixed_k_tangent"] and meta["correction_mass"] == 0
assert meta["base_coefficient_sha256"] == "ced15ecef5037f0e638664999a6ec82f46e02585a95895f907aa14ed3704c600"
source, = set(re.findall(r"[a-f0-9]{64}", (root.parent/"source-fingerprint.txt").read_text()))
d = arrays["direction"][:, 0]
sx = arrays["unknown_scale"][:, 0]
size = float(np.max(np.abs(d/sx)))
radii = [.00020533, .025, 1000.] + [a*size for a in (.72, .75, .77, .79, .82, .85, .87, .9)]
pieces = [("scaled_archived_R9_direction", radius, d*(radius/size)) for radius in radii]
pieces.extend((label, None, d) for label in ("archived_R9_direction", "repeat_archived_R9_direction"))
candidates = []
output = root/"candidates.bin"
with output.open("xb") as stream:
    for index, (label, radius, vector) in enumerate(pieces):
        payload = vector.astype("<f8", copy=False).tobytes()
        stream.write(payload)
        candidates.append(dict(index=index, label=label, radius=radius, size=float(np.max(np.abs(vector/sx))),
                               proposed_scaled_l2=float(np.linalg.norm(vector/sx)),
                               direction_sha256=hashlib.sha256(payload).hexdigest()))
common = {key:meta[key] for key in BINDINGS}
common.update(source_identity=source,current_fd=.000512,factor=None)
write_json(root/"candidates.json", {
    "basis": None, "archived_direction":pin(archive), "input":pin(output),
    "forcing":.3, "targets":[checks["physical_target"],checks["scaled_target"]],
    "basis_receipt_identity":common, "basis_rank":0,
    "frozen_base":array_spec(meta,"raw_base_residual"), "candidates":candidates})
