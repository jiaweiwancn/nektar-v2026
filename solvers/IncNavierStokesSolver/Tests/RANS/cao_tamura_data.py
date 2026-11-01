#!/usr/bin/env python3
"""Freeze/audit published Cao--Tamura data; does not run or calibrate CFD.

Extraction requires pdfplumber and pypdf. Verification of the frozen tables
requires only Python's standard library. The source PDF is user supplied and
is never copied into the repository. See the data README for uncertainty,
source ambiguities, exclusion decisions and reproduction commands.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
PDF_SHA = "adffb5c1d740be791211d82ce35778a37a3315a6dc72c56c80b70cc6642a43df"
DATA_VERSION = "cao-tamura-2006-v1.1"
STATIONS = [-2.5, -1.25, 0., 1.25, 2.5, 3.75, 5., 6.25, 7.5]
HILL_HEIGHT = .04
DELTA_NOMINAL = .25
COORDINATE_BOUND = .5  # PDF point, extraction envelope, not experimental CI.
FROZEN_FILES = ["README.md", "extraction.json", "series.json", "profiles.csv",
                "reported_scalars.csv", "figure_inventory.csv", "source_manifest.json",
                "uncertainty.json", "scores.json"]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def axis(a, b, logarithmic=False):
    return dict(anchors=[a, b], logarithmic=logarithmic)


def transform(spec, position):
    (p0, v0), (p1, v1) = spec["anchors"]
    if spec["logarithmic"]:
        v0, v1 = math.log10(v0), math.log10(v1)
    value = v0 + (position-p0)*(v1-v0)/(p1-p0)
    return 10**value if spec["logarithmic"] else value


def definitions():
    """Manually audited tick/bar anchors and drawing selections, PDF points.

    Curve indices are zero based pdfplumber page.curves indices. Select filled
    symbol paths only; outlines, connecting lines, legends and coloured
    dividing/peak guides are not additional observations. Multi-profile zero
    origins use station ticks, not the manually positioned text '0' labels.
    """
    panels = []

    def add(name, page, selection, ax, ay, quantity, height, cases,
            role="validation", stations=None, origins=None, unit=None,
            orientation="profile", duplicate_of=None):
        panels.append(dict(id=name, page=page, selection=selection,
                           axis_x=ax, axis_y=ay, quantity=quantity,
                           height_kind=height, cases=cases, role=role,
                           stations=stations, origins=origins, unit=unit,
                           orientation=orientation, duplicate_of=duplicate_of))

    # Each selection = [first curve, last curve, case, shape].
    add("f3a", 7, [[2, 30, "R", "circle"], [33, 67, "S", "circle"]],
        axis([106.554, 0], [222.474, 1]),
        axis([222.797, 0], [93.638, 300]), "U_over_Uinf", "y_mm", ["S", "R"],
        role="calibration")
    # Fig. 3(b): frame bottom is NOT the labelled y/delta=0 position.
    add("f3b", 7, [[294, 328, "S", "circle"], [331, 359, "R", "circle"]],
        axis([273.171, 0], [378.452, .004]),
        axis([160.403, .5], [95.809, 1]), "minus_uv_over_Uinf2", "y_over_delta",
        ["S", "R"], role="calibration")
    add("f3c", 7, [[148, 182, "S", "circle"], [185, 213, "R", "circle"]],
        axis([106.019, 0], [233.999, 3]),
        axis([406.592, 0], [290.298, 1]), "sigma_u_over_utau",
        "y_minus_d_over_delta", ["S", "R"], role="calibration")
    add("f3d", 7, [[106, 140, "S", "circle"], [75, 103, "R", "circle"]],
        axis([272.324, 0], [399.528, 2]),
        axis([405.6, 0], [289.95, 1]), "sigma_vertical_over_utau",
        "y_minus_d_over_delta", ["S", "R"], role="flat_crosscheck")
    add("f3e", 7, [[220, 254, "S", "circle"], [256, 284, "R", "circle"]],
        axis([185.228, 1], [342.853, 1000], True),
        axis([569.938, 0], [448.793, 30]), "U_over_utau", "y_minus_d_mm",
        ["S", "R"], role="duplicate_crosscheck", orientation="log_profile",
        duplicate_of="f3a")

    multiprofiles = [
        ("f4a", 9, 232, 458, "SS", "U_over_Uinf", 143.541, 77.148, 16.582,
         [150.327, 171.101, 191.793, 212.568, 233.237, 254.012, 274.789, 295.477, 316.254]),
        ("f4b", 9, 3, 199, "RR", "U_over_Uinf", 271.484, 204.497, 16.773,
         [150.581, 171.486, 192.495, 213.505, 234.408, 255.416, 276.429, 297.329, 318.341]),
        ("f5a", 9, 526, 750, "SS", "sigma_u_over_Uinf", 452.485, 385.986, 83.125,
         [151.806, 172.647, 193.378, 214.219, 234.973, 255.814, 276.653, 297.384, 318.227]),
        ("f5b", 9, 796, 992, "RR", "sigma_u_over_Uinf", 573.730, 507.232, 83.120,
         [151.588, 172.407, 193.138, 213.978, 234.710, 255.552, 276.370, 297.123, 317.942]),
        ("f10", 14, 25, 448, "both", "sigma_u_over_Uinf", 144.375, 78.293, 84.935,
         [137.476, 158.741, 179.966, 201.252, 222.537, 243.801, 265.028, 286.312, 307.597]),
    ]
    for name, page, first, last, case, quantity, bottom, top, scale, origins in multiprofiles:
        sel = ([[25, 221, "RR", "circle"], [224, 448, "SS", "circle"]]
               if case == "both" else [[first, last, case, "circle"]])
        add(name, page, sel, axis([0, 0], [scale, 1]),
            axis([bottom, 0], [top, 4]), quantity, "y_over_H",
            ["SS", "RR"] if case == "both" else [case],
            role="duplicate_crosscheck" if name.startswith("f5") else "validation",
            stations=STATIONS, origins=origins, unit=scale,
            duplicate_of="f10" if name.startswith("f5") else None)

    add("f6a", 11, [[1, 25, "SS_upwind", "circle"], [27, 49, "RR_upwind", "circle"],
                    [51, 69, "SS", "triangle"], [71, 85, "RR", "triangle"],
                    [87, 111, "SS_downwind", "triangle"], [113, 135, "RR_downwind", "triangle"]],
        axis([68.521, -1], [218.146, 1]), axis([242.121, 0], [71.182, 8]),
        "speedup", "yprime_over_H", ["SS", "RR"], role="duplicate_crosscheck",
        duplicate_of="f6b_crest")
    add("f6b", 11, [[150, 168, "SS", "circle"], [170, 184, "RR", "circle"],
                    [186, 202, "SR", "triangle"], [204, 218, "RS", "triangle"]],
        axis([271.841, 0], [421.688, 1]), axis([242.204, 0], [71.182, 8]),
        "speedup", "yprime_over_H", ["SS", "RR", "RS", "SR"])
    add("f7a", 12, [[1, 19, "SS", "circle"], [21, 35, "RR", "circle"]],
        axis([88.376, .2], [208.835, 1]), axis([245.38, 0], [97.047, 1]),
        "sigma_u_over_flat_same_height", "yprime_over_delta", ["SS", "RR"])
    add("f7b", 12, [[42, 60, "SS", "circle"], [62, 76, "RR", "circle"]],
        axis([263.798, 0], [413.251, 2]), axis([244.863, 0], [97.794, 1]),
        "uv_over_flat_same_height", "yprime_over_delta", ["SS", "RR"])
    add("f8", 13, [[2, 138, "SS", "circle"], [140, 250, "RR", "circle"]],
        axis([0, 0], [40.822, 1]), axis([161.785, 0], [80.138, 2]),
        "U_over_Uinf", "y_over_H", ["SS", "RR"], role="wake_mean_crosscheck",
        stations=[0., 1.25, 2.5, 3.75, 5., 6.25, 7.5],
        origins=[68.353, 119.382, 170.410, 221.438, 272.467, 323.524, 374.552],
        unit=40.822, duplicate_of="f4")
    add("f9", 13, [[863, 903, "SS", "circle"], [777, 817, "RR", "circle"]],
        axis([175.832, 4], [272.206, 8]), axis([321.238, 0], [251.042, 2]),
        "U_m_per_s", "probe_fixed_height", ["SS", "RR"],
        role="dimensional_holdout", orientation="traverse")
    add("f11", 15, [[2, 30, "SS", "circle"], [33, 55, "RR", "circle"]],
        axis([156.705, 0], [335.465, .01]), axis([282.264, 0], [106.064, 1]),
        "minus_uv_over_Uinf2", "y_over_delta75", ["SS", "RR"])
    add("f12", 15, [[107, 145, "SS", "circle"], [66, 104, "RR", "circle"]],
        axis([156.693, 3], [342.079, 9]), axis([497.172, .5], [366.348, 1]),
        "sigma_u_m_per_s", "probe_fixed_height", ["SS", "RR"],
        role="dimensional_holdout", orientation="traverse")
    return panels


def marker(curve):
    command = "".join(q[0] for q in curve["path"])
    if not (curve["fill"] and 1 < curve["width"] < 8 and 1 < curve["height"] < 8):
        return None
    if command == "mcccc":
        shape = "circle"
        x = (curve["x0"]+curve["x1"])/2
        y = (curve["top"]+curve["bottom"])/2
    elif command == "mlll":
        shape = "triangle"
        vertices = [q[1] for q in curve["path"][:3]]
        x, y = (sum(v[k] for v in vertices)/3 for k in range(2))
    else:
        return None
    return dict(px=x, py=y, shape=shape,
                radius_x=curve["width"]/2, radius_y=curve["height"]/2,
                white=min(curve["non_stroking_color"]) > .9)


def groups_by_station(points, stations):
    groups = [[]]
    for point in points:
        if groups[-1] and point["py"] > groups[-1][-1]["py"] + .5:
            groups.append([])
        groups[-1].append(point)
    assert len(groups) == len(stations), (len(groups), stations)
    # All these panels include a highest-height point in each station.
    groups.sort(key=lambda group: min(group, key=lambda p: p["py"])["px"])
    return groups


def physical_height(panel, case, x, y, station):
    d = .0035 if case.startswith("R") else 0.
    wall_height = 0. if case in ["S", "R"] else hill(station)
    kind = panel["height_kind"]
    if kind == "y_mm":
        yy = y*.001/HILL_HEIGHT
    elif kind == "y_over_delta":
        yy = y*DELTA_NOMINAL/HILL_HEIGHT
    elif kind == "y_minus_d_over_delta":
        yy = (y*DELTA_NOMINAL+d)/HILL_HEIGHT
    elif kind == "y_minus_d_mm":
        yy = (x*.001+d)/HILL_HEIGHT
    elif kind == "y_over_H":
        yy = y
    elif kind == "yprime_over_H":
        yy = y + wall_height
    elif kind == "yprime_over_delta":
        yy = y*DELTA_NOMINAL/HILL_HEIGHT + wall_height
    elif kind == "probe_fixed_height":
        yy = .175 if case.startswith("R") else .05
    else:
        # delta at x/H=7.5 is unspecified: preserve the published coordinate.
        return None, None
    return yy, yy-wall_height


def hill(x_over_h):
    return math.cos(math.pi*x_over_h/5)**2 if abs(x_over_h) < 2.5 else 0.


def height_weights(rows, independent="axis_y"):
    """Voronoi height cells; duplicate-height markers share one cell's weight."""
    levels = collections.defaultdict(list)
    # Quantisation groups nearly coincident overplots, not independent repeats.
    for row in rows:
        levels[round(row[independent], 6)].append(row)
    ordered = sorted(levels)
    if len(ordered) == 1:
        for row in rows:
            row["weight"] = 1/len(rows)
        return
    total = ordered[-1]-ordered[0]
    for i, level in enumerate(ordered):
        left = ordered[0] if i == 0 else (ordered[i-1]+level)/2
        right = ordered[-1] if i == len(ordered)-1 else (level+ordered[i+1])/2
        for row in levels[level]:
            row["weight"] = (right-left)/total/len(levels[level])


def extract(pdf, output):
    import pdfplumber
    assert sha(pdf) == PDF_SHA, "Source PDF fingerprint changed"
    assert not (output/"snapshot.json").exists(), "Extract to a fresh directory; data are frozen"
    output.mkdir(parents=True, exist_ok=True)
    panels = definitions()
    rows, series = [], []
    with pdfplumber.open(pdf) as doc:
        for panel in panels:
            page = doc.pages[panel["page"]-1]
            legends = [r for r in page.rects if r["fill"] and
                       min(r["non_stroking_color"]) > .9 and
                       r["width"] > 35 and 4 < r["height"] < 60]
            for first, last, case, shape in panel["selection"]:
                pts = []
                for index in range(first, last+1):
                    point = marker(page.curves[index])
                    if point and point["shape"] == shape:
                        assert point["white"] == case.startswith("R"), (panel["id"], case, index)
                        pts.append(dict(point, source_curve=index))
                groups = (groups_by_station(pts, panel["stations"])
                          if panel["stations"] else [pts])
                for group_index, points in enumerate(groups):
                    station = panel["stations"][group_index] if panel["stations"] else 0.
                    if panel["id"] == "f6a":
                        station = -2.5 if "upwind" in case else 2.5 if "downwind" in case else 0.
                    if panel["id"] == "f11":
                        station = 7.5
                    offset = panel["origins"][group_index] if panel["origins"] else 0.
                    tag = (f"_x{station:g}" if panel["stations"] else "")
                    sid = f"{panel['id']}_{case}{tag}"
                    role = panel["role"]
                    duplicate = panel["duplicate_of"]
                    if panel["id"] == "f6a" and ("upwind" in case or "downwind" in case):
                        role, duplicate = "validation", None
                    meta = dict(id=sid, panel=panel["id"], case=case,
                                station_x_over_H=(None if panel["orientation"] == "traverse" else station),
                                quantity=panel["quantity"],
                                height_kind=panel["height_kind"], role=role,
                                page=panel["page"], offset_pdf_pt=offset,
                                shape=shape, white=case.startswith("R"),
                                coordinate_conversion=("delta75_unspecified" if panel["height_kind"] == "y_over_delta75"
                                    else "nominal_delta_0.25m" if "delta" in panel["height_kind"]
                                    else "paper_probe_height" if panel["height_kind"] == "probe_fixed_height"
                                    else "physical_Hh_0.04m"),
                                duplicate_of=duplicate, count=len(points))
                    series.append(meta)
                    current = []
                    for point in points:
                        x = transform(panel["axis_x"], point["px"]-offset)
                        y = transform(panel["axis_y"], point["py"])
                        xbound = max(abs(transform(panel["axis_x"], point["px"]-offset+sgn*COORDINATE_BOUND)-x)
                                     for sgn in [-1, 1])
                        ybound = max(abs(transform(panel["axis_y"], point["py"]+sgn*COORDINATE_BOUND)-y)
                                     for sgn in [-1, 1])
                        xx = x if panel["orientation"] == "traverse" else station
                        yy, yp = physical_height(panel, case, x, y, xx)
                        value = y if panel["orientation"] in ["log_profile", "traverse"] else x
                        bound = ybound if panel["orientation"] in ["log_profile", "traverse"] else xbound
                        quality = ("legend_occluded" if any(r["x0"] < point["px"] < r["x1"] and
                                   r["top"] < point["py"] < r["bottom"] for r in legends)
                                   else "visible")
                        row = dict(series_id=sid, point_id=len(current),
                                   source_curve=point["source_curve"], px=point["px"], py=point["py"],
                                   radius_x=point["radius_x"], radius_y=point["radius_y"],
                                   axis_x=x, axis_y=y, axis_x_bound=xbound, axis_y_bound=ybound,
                                   x_over_H=xx, y_over_H=yy, yprime_over_H=yp,
                                   value=value, value_bound=bound, weight=0., role=role,
                                   quality=quality)
                        current.append(row)
                    height_weights([r for r in current if r["quality"] == "visible"],
                                   "axis_x" if panel["orientation"] in ["log_profile", "traverse"] else "axis_y")
                    rows.extend(current)
    write_json(output/"extraction.json", dict(data_version=DATA_VERSION, source_sha256=PDF_SHA,
               coordinate_extraction_bound_pdf_pt=COORDINATE_BOUND, panels=panels))
    write_json(output/"series.json", series)
    with (output/"profiles.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: format(v, ".10g") if isinstance(v, float) else v
                             for k, v in row.items()})
    return rows, series


def independent_paths(pdf):
    """Second parser: interpret pypdf paint operators/CTM, no pdfplumber API."""
    from pypdf import PdfReader
    from pypdf.generic import ContentStream
    result = collections.defaultdict(list)
    reader = PdfReader(pdf)
    for number in [7, 9, 11, 12, 13, 14, 15]:
        page = reader.pages[number-1]
        ctm = (1., 0., 0., 1., 0., 0.)
        colour = (0., 0., 0.)
        stack, path = [], []

        def mapped(x, y):
            a, b, c, d, e, f = ctm
            return a*x+c*y+e, float(page.mediabox.height)-(b*x+d*y+f)

        for operands, operator in ContentStream(page.get_contents(), reader).operations:
            if operator == b"q":
                stack.append((ctm, colour))
            elif operator == b"Q":
                ctm, colour = stack.pop()
            elif operator == b"cm":
                a,b,c,d,e,f = map(float, operands)
                A,B,C,D,E,F = ctm
                ctm = (a*A+b*C, a*B+b*D, c*A+d*C, c*B+d*D, e*A+f*C+E, e*B+f*D+F)
            elif operator in [b"scn", b"rg", b"g"]:
                values = tuple(map(float, operands))
                colour = values if len(values) == 3 else values*3
            elif operator in [b"m", b"l", b"c"]:
                values = list(map(float, operands))
                path.append((operator.decode(), [mapped(*values[i:i+2]) for i in range(0, len(values), 2)]))
            elif operator == b"h":
                path.append(("h", []))
            elif operator == b"re":
                path.append(("rect", []))
            elif operator in [b"f", b"F", b"f*", b"B", b"B*", b"b", b"b*"]:
                commands = "".join(q[0] for q in path)
                points = [point for _, group in path for point in group]
                if commands == "mcccc":
                    xx, yy = zip(*points)
                    if 1 < max(xx)-min(xx) < 8 and 1 < max(yy)-min(yy) < 8:
                        result[number].append(("circle", min(colour) > .9,
                                               (max(xx)+min(xx))/2, (max(yy)+min(yy))/2))
                elif commands == "mlll":
                    vertices = [path[i][1][0] for i in range(3)]
                    result[number].append(("triangle", min(colour) > .9,
                                           *(sum(p[k] for p in vertices)/3 for k in range(2))))
                path = []
            elif operator in [b"S", b"s", b"n"]:
                path = []
            elif operator == b"Do":
                raise AssertionError("Unexpected form/image in audited plot pages")
    return result


def parser_check(pdf, rows, series):
    paths = independent_paths(pdf)
    meta = {s["id"]: s for s in series}
    errors, counts = [], collections.Counter()
    for row in rows:
        panel = meta[row["series_id"]]
        candidates = paths[panel["page"]]
        error = min(math.hypot(px-row["px"], py-row["py"])
                    for shape, white, px, py in candidates
                    if shape == panel["shape"] and white == panel["white"])
        assert error < 1e-6, (row["series_id"], row["point_id"], error)
        errors.append(error)
        counts[panel["page"]] += 1
    return dict(method="pypdf paint-operator/CTM interpreter vs pdfplumber symbols",
                points=len(errors), maximum_centre_difference_pdf_pt=max(errors),
                points_by_page=dict(counts), status="PASS")


def source_consistency(rows, series):
    """Independent repeated-figure and reported-number checks, not CFD scores."""
    grouped = collections.defaultdict(list)
    for row in rows:
        if row["quality"] == "visible":
            grouped[row["series_id"]].append(row)

    def interpolate(data, key, value, ordinate="value"):
        points = sorted((r[key], r[ordinate]) for r in data)
        for (x0, y0), (x1, y1) in zip(points, points[1:]):
            if x0 <= value <= x1 and x1 > x0:
                return y0+(value-x0)/(x1-x0)*(y1-y0)
        return None

    result = {"scope": "data/source consistency only; no CFD comparison"}
    repeats = []
    for first, second in [("f4a_SS", "f8_SS"), ("f4b_RR", "f8_RR"),
                          ("f5a_SS", "f10_SS"), ("f5b_RR", "f10_RR")]:
        differences = []
        for name in grouped:
            if name.startswith(second+"_x"):
                other = first+"_x"+name.split("_x")[-1]
                for row in grouped[name]:
                    expected = interpolate(grouped[other], "y_over_H", row["y_over_H"])
                    if expected is not None:
                        differences.append(abs(row["value"]-expected))
        maximum = max(differences)
        # These use different print packing; the combined 0.5-pt extraction
        # envelopes are >0.04 for U and >0.01 for sigma/Uinf.
        limit = .04 if first.startswith("f4") else .013
        assert maximum < limit, (first, second, maximum)
        repeats.append(dict(first=first, second=second, points=len(differences),
                            maximum_absolute_difference=maximum, comparison_bound=limit))
    result["duplicate_figure_checks"] = repeats
    inflow = []
    for case, utau, reynolds in [("S", .033, 12000), ("R", .051, 10500)]:
        uh = interpolate(grouped["f3a_"+case], "axis_y", 40.)
        ratio = reynolds/(HILL_HEIGHT*uh)
        differences = []
        for row in grouped["f3e_"+case]:
            expected = interpolate(grouped["f3a_"+case], "y_over_H", row["y_over_H"])
            if expected is not None:
                differences.append(abs(utau*row["value"]-expected))
        inflow.append(dict(case=case, U_H_over_Uinf=uh,
                           Uinf_over_nu_per_m_from_Re_H=ratio,
                           Uinf_m_per_s_IF_nu_1p5em5=ratio*1.5e-5,
                           f3e_vs_f3a_max_U_over_Uinf_difference=max(differences),
                           f3e_role="duplicate cross-check, excluded from calibration"))
    rough_ratio = 3.97/(.051*.0002)
    result["Reynolds_consistency"] = dict(inflows=inflow,
        rough_Uinf_over_nu_per_m_from_Re_s=rough_ratio,
        rough_Re_s_from_Re_H=inflow[1]["Uinf_over_nu_per_m_from_Re_H"]*.051*.0002,
        interpretation="About 1.8% rough Re_s discrepancy; approximate source values, no unique absolute pair (Uinf, nu)")
    crest = []
    for case, reported in [("SS", .457), ("RR", .875), ("RS", .918), ("SR", .245)]:
        row = min(grouped["f6b_"+case], key=lambda r: abs(r["axis_y"]-.25))
        assert abs(row["axis_y"]-.25) < row["axis_y_bound"]
        assert abs(row["value"]-reported) < .01
        crest.append(dict(case=case, reported=reported, digitized=row["value"],
                          digitized_yprime_over_H=row["axis_y"]))
    result["text_vs_crest_markers"] = crest
    result["traverse_ranges"] = {name: [min(r["axis_x"] for r in grouped[name]),
                                          max(r["axis_x"] for r in grouped[name])]
                                   for name in ["f9_SS", "f9_RR", "f12_SS", "f12_RR"]}
    return result


def raster_check(pdf, rows, series, directory):
    """Poppler raster check on isolated, visible circle symbols.

    The location search is local; it checks rendering/selection, not a blind
    second assignment of axis semantics. All points additionally receive the
    independent vector-parser check. Crowded/triangle points are reviewed in
    page overlays rather than pretending that an obscured symbol has an
    independent raster centre.
    """
    import subprocess
    import pdfplumber
    from PIL import Image, ImageDraw
    directory.mkdir(parents=True, exist_ok=True)
    meta = {s["id"]: s for s in series}
    pages = collections.defaultdict(list)
    for row in rows:
        pages[meta[row["series_id"]]["page"]].append(row)
    dpi = 300
    scale = dpi/72
    checked, by_page, skipped = [], collections.Counter(), collections.Counter()
    with pdfplumber.open(pdf) as source:
        clutter = {n: [m for c in source.pages[n-1].curves if (m := marker(c))]
                   for n in pages}
    for page, points in sorted(pages.items()):
        prefix = directory/f"p{page:02d}"
        target = prefix.with_suffix(".png")
        subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page),
                        "-singlefile", "-r", str(dpi), "-gray", "-png", str(pdf), str(prefix)],
                       check=True, timeout=60, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        raster = Image.open(target).convert("L")
        overlay = raster.convert("RGB")
        draw = ImageDraw.Draw(overlay)
        for row in points:
            info = meta[row["series_id"]]
            px, py, rx, ry = (row[k] for k in ["px", "py", "radius_x", "radius_y"])
            colour = "#cc3300" if not info["white"] else "#0055cc"
            if row["quality"] != "visible":
                colour = "#dd00dd"
            draw.ellipse(((px-rx-.4)*scale, (py-ry-.4)*scale,
                          (px+rx+.4)*scale, (py+ry+.4)*scale), outline=colour, width=1)
            if info["shape"] != "circle" or row["quality"] != "visible":
                skipped["triangle_or_legend"] += 1
                continue
            if any((math.hypot(px-other["px"], py-other["py"]) > 1e-5 or
                    other["white"] != info["white"]) and
                   math.hypot(px-other["px"], py-other["py"]) <
                   max(rx,ry)+max(other["radius_x"],other["radius_y"])+.6
                   for other in clutter[page]):
                skipped["crowded"] += 1
                continue
            if py-ry < 0:
                skipped["edge"] += 1
                continue
            def grey(x, y):
                x0, y0 = math.floor(x), math.floor(y)
                fx, fy = x-x0, y-y0
                return sum(raster.getpixel((x0+ix, y0+iy))*wx*wy
                           for ix, wx in [(0, 1-fx), (1, fx)]
                           for iy, wy in [(0, 1-fy), (1, fy)])

            def edge(angle):
                ex, ey = math.cos(angle), math.sin(angle)
                radius = 1/math.sqrt((ex/rx)**2+(ey/ry)**2)
                previous = None
                crossings = []
                for i in range(51):
                    distance = (.35+1.4*i/50)*radius
                    value = grey((px+ex*distance)*scale-.5, (py+ey*distance)*scale-.5)
                    if previous and previous[1] < 128 <= value:
                        crossings.append(previous[0]+(distance-previous[0])*(128-previous[1])/(value-previous[1]))
                    previous = distance, value
                # The first outer crossing belongs to the symbol. A later
                # crossing may be a different series' connecting stroke.
                outer = [distance for distance in crossings if distance > .75*radius]
                return outer[0] if outer else None

            pairs = []
            for k in range(16):
                angle = math.pi*k/16
                a, b = edge(angle), edge(angle+math.pi)
                if a is not None and b is not None:
                    pairs.append((math.cos(angle), math.sin(angle), (a-b)/2))
            if len(pairs) < 8:
                skipped["insufficient_visible_edges"] += 1
                continue
            # Opposite outer-ink edge midpoints estimate the centre. Reject
            # angular outliers from connecting strokes, not coordinate errors.
            for iteration in range(3):
                aa = sum(x*x for x,y,z in pairs)
                ab = sum(x*y for x,y,z in pairs)
                bb = sum(y*y for x,y,z in pairs)
                az = sum(x*z for x,y,z in pairs)
                bz = sum(y*z for x,y,z in pairs)
                determinant = aa*bb-ab*ab
                dx, dy = (az*bb-bz*ab)/determinant, (bz*aa-az*ab)/determinant
                if iteration < 2:
                    pairs = sorted(pairs, key=lambda v: abs(v[0]*dx+v[1]*dy-v[2]))[:max(8, len(pairs)-2)]
            error = math.hypot(dx, dy)
            checked.append(dict(series=row["series_id"], point=int(row["point_id"]),
                                centre_difference_pdf_pt=error, usable_edge_pairs=len(pairs)))
            by_page[page] += 1
        overlay.save(directory/f"overlay-{page:02d}.png")
    assert len(checked) >= 200
    # A 0.5 pt total extraction envelope; one raster pixel is 0.24 pt.
    maximum = max(r["centre_difference_pdf_pt"] for r in checked)
    failed = [r for r in checked if r["centre_difference_pdf_pt"] > COORDINATE_BOUND]
    write_json(directory/"raster-points.json", checked)
    assert not failed, failed
    return dict(method="Poppler 300 dpi / opposite outer-ink edge midpoints on isolated circles",
                status="PASS", checked_points=len(checked), checked_by_page=dict(by_page),
                skipped=dict(skipped), maximum_centre_difference_pdf_pt=maximum,
                minimum_edge_pairs=min(r["usable_edge_pairs"] for r in checked),
                point_receipt_sha256=sha(directory/"raster-points.json"),
                image_sha256={f"p{page:02d}.png": sha(directory/f"p{page:02d}.png") for page in pages},
                overlays="All selected symbols overlaid; legends magenta; S solid red / R open blue")


def load_rows(data):
    rows = list(csv.DictReader((data/"profiles.csv").open()))
    for row in rows:
        for key in row:
            if key not in ["series_id", "role", "quality"]:
                row[key] = None if row[key] == "" else float(row[key])
    return rows


def verify(data):
    series = json.loads((data/"series.json").read_text())
    rows = load_rows(data)
    metadata = {s["id"]: s for s in series}
    assert len(metadata) == len(series)
    counts, weights = collections.Counter(), collections.defaultdict(float)
    checks = ["source binding", "series/count integrity", "finite values",
              "positive extraction envelopes", "vertical/physical-height transform",
              "normalised per-series height weights", "canonical extraction recipe"]
    for row in rows:
        assert row["series_id"] in metadata
        assert row["role"] == metadata[row["series_id"]]["role"]
        assert row["quality"] in ["visible", "legend_occluded"]
        counts[row["series_id"]] += 1
        weights[row["series_id"]] += row["weight"]
        assert row["axis_x_bound"] > 0 and row["axis_y_bound"] > 0
        assert (row["weight"] == 0) == (row["quality"] == "legend_occluded")
        assert all(math.isfinite(value) for value in row.values() if isinstance(value, (float, int)))
        if row["y_over_H"] is not None:
            surface = metadata[row["series_id"]]["case"]
            wall_height = 0. if surface in ["S", "R"] else hill(row["x_over_H"])
            assert math.isclose(row["y_over_H"]-row["yprime_over_H"], wall_height, abs_tol=1e-8)
            assert row["yprime_over_H"] >= -.02, row
    for s in series:
        assert counts[s["id"]] == s["count"]
        assert math.isclose(weights[s["id"]], 1, abs_tol=2e-8)
    recipe = json.loads((data/"extraction.json").read_text())
    assert recipe["source_sha256"] == PDF_SHA
    assert recipe["data_version"] == DATA_VERSION
    assert recipe["panels"] == definitions(), "Extraction recipe differs from audited data version"
    if (data/"snapshot.json").exists():
        frozen = json.loads((data/"snapshot.json").read_text())
        assert frozen["source_sha256"] == PDF_SHA
        assert frozen["version"] == DATA_VERSION
        assert sorted(frozen["files"]) == sorted(FROZEN_FILES)
        for name, digest in frozen["files"].items():
            assert Path(name).name == name
            assert sha(data/name) == digest, f"Frozen data changed: {name}"
        checks.append("frozen snapshot file hashes")
    if (data/"reported_scalars.csv").exists():
        scalars = list(csv.DictReader((data/"reported_scalars.csv").open()))
        assert len({s["id"] for s in scalars}) == len(scalars)
        assert all(math.isfinite(float(s["value"])) for s in scalars)
        assert json.loads((data/"source_manifest.json").read_text())["pdf"]["sha256"] == PDF_SHA
        score = json.loads((data/"scores.json").read_text())
        assert set(s["id"] for s in series if s["role"] == "calibration") == set(
            name for name in score["data_roles"]["calibration"] if name in metadata)
        assert not any(s["case"] in ["SS", "RR", "RS", "SR"] for s in series if s["role"] == "calibration")
        assert set(s["id"] for s in series if s["role"] == "flat_crosscheck") == {"f3d_S", "f3d_R"}
        inventory = list(csv.DictReader((data/"figure_inventory.csv").open()))
        assert [int(row["figure"]) for row in inventory] == list(range(1,15))
        checks.extend(["calibration/holdout separation", "figure inventory and scalar identifiers"])
    return dict(status="PASS", data_version=DATA_VERSION, series=len(series), points=len(rows),
                quality=dict(collections.Counter(r["quality"] for r in rows)),
                roles=dict(collections.Counter(s["role"] for s in series)),
                checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--pdf", type=Path)
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--raster", type=Path,
                        help="Local output for Poppler cross-check/overlays (requires --pdf)")
    parser.add_argument("--seal", action="store_true", help="Write immutable versioned data snapshot hashes")
    args = parser.parse_args()
    if args.extract:
        assert args.pdf, "--extract requires --pdf"
        rows, series = extract(args.pdf, args.data)
    else:
        rows = load_rows(args.data)
        series = json.loads((args.data/"series.json").read_text())
    if args.seal:
        assert not (args.data/"snapshot.json").exists(), "Snapshot already exists"
        verify(args.data)
        write_json(args.data/"snapshot.json", dict(version=DATA_VERSION,
                   source_sha256=PDF_SHA, files={name: sha(args.data/name) for name in FROZEN_FILES}))
    report = verify(args.data)
    report["source_consistency"] = source_consistency(rows, series)
    if args.pdf:
        assert sha(args.pdf) == PDF_SHA
        report["independent_parser"] = parser_check(args.pdf, rows, series)
    if args.raster:
        assert args.pdf, "--raster requires --pdf"
        report["raster"] = raster_check(args.pdf, rows, series, args.raster)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        write_json(args.report, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
