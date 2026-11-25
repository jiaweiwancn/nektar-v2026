# SPDX-License-Identifier: MIT
"""Required physical chart locations for an original normal-hill release."""

from math import isfinite
from sys import float_info


def locations(geometry):
    """Keep support controls on the exact reference h, independent of the mesh."""
    H, h = geometry["H"], geometry["h"]
    L = geometry["hill_half_length"]
    left, top = geometry["left"], geometry["top"]
    return {
        "upstream_collar": (left + .228*(-L-left), H*(1111./1500.)),
        "downstream_outer": (L, H+min(.957899, .5*(top-H))),
        "support_original": (.75*L, h),
        "support_refined": (.5625*L, h),
        "support_downstream_foot": (L, h),
    }


def tolerance(x, y, H):
    # Chart face bounds come from decimal mesh vertices and a displacement
    # subtraction. The exact reference h is left untouched at either owner.
    return 64*float_info.epsilon*max(1., abs(x), abs(y), H)


def coincides(x, y, H):
    return isfinite(x) and isfinite(y) and abs(x-y) <= tolerance(x, y, H)


def contains(element, a, r, H):
    xl, xr, rl, rr = element["bounds"]
    return (xl-tolerance(a, xl, H) <= a <= xr+tolerance(a, xr, H)
            and rl-tolerance(r, rl, H) <= r <= rr+tolerance(r, rr, H))


def valid_owners(data, margins, points):
    """Check the full fixed ownership set and the exact reported support witness."""
    geometry = data["geometry"]
    H = geometry["H"]
    cells = {e["geometry"]: e for e in data["elements"]}
    required = locations(geometry)
    expected = {(label, i) for label, (a, r) in required.items()
                for i, e in cells.items() if contains(e, a, r, H)}

    def valid(point):
        cell = cells.get(point["geometry"])
        return cell is not None and contains(cell, point["a"], point["rho"], H)

    seen = set()
    for point in points.get("fixed_controls", []):
        label = point.get("control")
        if not valid(point):
            return False
        if label is None:
            continue
        key = (label, point["geometry"])
        if key not in expected or key in seen:
            return False
        a, r = required[label]
        if not (coincides(point["a"], a, H) and coincides(point["rho"], r, H)):
            return False
        seen.add(key)
    if expected != seen or {label for label, _ in seen} != required.keys():
        return False

    records = points["local_minima"]
    if not all(valid(value) for record in records for value in record["values"]):
        return False
    support = [record for record in records if record["region"] == "forcing_support"]
    witness = margins["regions"]["forcing_support"]["witness"]
    return len(support) == 1 and any(
        point["side"] == "reported-minimum" and point["geometry"] == witness["geometry"]
        and coincides(point["a"], witness["a"], H) and coincides(point["rho"], witness["rho"], H)
        for point in support[0]["values"])


def positive_point(point):
    spread = point["differentiation_spread"]
    minima = point["minima"]
    return (isfinite(spread) and spread >= 0. and len(minima) >= 3
            and all(isfinite(point[key]) and point[key] > bound for key, bound in
                    (("k", 0.), ("epsilon", 0.), ("covariance_over_k", 1e-10), ("mu", 0.)))
            and all(isfinite(m["over_mu"]) and m["over_mu"]-spread > 1e-8 for m in minima))


def positive_region(region):
    uncertainty = region["max_directional_difference"]
    return (region["sample_pass"] and region["samples"] > 0
            and isfinite(uncertainty) and uncertainty >= 0.
            and all(isfinite(region[key]) and region[key] > limit for key, limit in
                    (("minimum_k", 0.), ("minimum_epsilon", 0.), ("covariance_over_k", 1e-10),
                     ("principal_over_mu", 1e-8+uncertainty))))
