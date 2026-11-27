#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect compact R10 evidence, retaining paths and hashes of the raw Linux run."""
import argparse
import json
from pathlib import Path
import re
import time
import xml.etree.ElementTree as ET

from resources import identity
from c1_support_controls import positive_point, positive_region, valid_owners

BUILDS = ("finite", "bounded", "recycled")
CASES = ("frozen-trials", "stale-base", "stale-scale", "stale-vector", "bounded-main",
         "cold-resume", "cold-reanchor", "changed-resume", "bounded-recycled", "bounded-enriched",
         "enriched-resume", "model-refinement", "bounded-active", "endpoint-rich")
LINEAGES = {"bounded-main":"main-lineage.json", "bounded-recycled":"recycled-lineage.json",
            "bounded-enriched":"enriched-lineage.json", "enriched-resume":"enriched-resume-lineage.json",
            "bounded-active":"active-lineage.json"}


def read(path):
    return json.loads(path.read_text())


def present(directory, names):
    return {name:identity(directory/name) for name in names if (directory/name).is_file()}


def timing(path):
    if not path.is_file():
        return None
    lines = path.read_text().splitlines()
    elapsed = next((line.split(": ",1)[1] for line in lines if "Elapsed (wall clock)" in line),None)
    peak = next((line.split(": ",1)[1] for line in lines if "Maximum resident set size" in line),None)
    if elapsed is None or peak is None:
        return None
    value = 0.
    for number in elapsed.split(":"):
        value = value*60.+float(number)
    return {"wall_seconds":value,"max_rss_kib":int(peak),"receipt":identity(path)}


def endpoint_audit(directory):
    result = {}
    packet = directory/"state.diagnostic-steady.json"
    if not packet.is_file():
        return result
    native = read(packet)
    result["native_domain"] = {
        key: native[key] for key in ("minimum_k", "minimum_epsilon", "minimum_covariance_over_k",
                                    "original_principal_over_mu", "physical_whole_cell_scalar_lower_bounds")
    }
    result["native_domain"]["interface_K_table_margin_over_Kin"] = (
        native["state"]["interface_K_table_margin_over_Kin"])
    result["native_quadrature"] = [native["state"]["quadratureX"], native["state"]["quadratureR"]]
    session = ET.parse(directory/"run.xml")
    option = session.find(".//SOLVERINFO/I[@PROPERTY='WJCheckCurrentActions']")
    result["native_action_crosschecks_requested"] = option is not None and option.get("VALUE") == "True"
    result["dependent_trace_direction_relative_error"] = native["dependent_trace_direction_relative_error"]
    equations, margins, points = (directory/(f"{name}.json") for name in ("equations", "margins", "points"))
    if equations.is_file():
        equations = read(equations)
        result["equations"] = {key: equations[key] for key in (
            "residuals", "residual_pass", "limits", "order", "radial_order", "epsilon_trace",
            "integration_minimum_jacobian", "raw_boundary", "power", "mass_integration_identity_error")}
        result["equations"]["interface"] = {
            key: equations["interface"][key] for key in (
                "shared_trace_max", "physical_face_measure_max_difference", "weak_free_interface_rows_max",
                "raw_traction_jump_rms", "raw_k_diffusion_jump_rms", "K_over_Kin")}
    if margins.is_file():
        margins = read(margins)
        result["one_sided_faces"] = margins["one_sided_faces"]
        result["regions"] = margins["regions"]
        result["positive_regions"] = all(positive_region(region) for region in margins["regions"].values())
        if points.is_file():
            points = read(points)
            controls = points["fixed_controls"]
            minima = [value for region in points["local_minima"] for value in region["values"]]
            result["owners"] = {
                "fixed_samples": len(controls), "local_samples": len(minima),
                "valid": valid_owners(native["state"], margins, points),
                "positive": all(positive_point(point) for point in controls + minima),
            }
            smallest = min(controls + minima, key=lambda point: min(
                value["over_mu"] for value in point["minima"])-point["differentiation_spread"])
            result["owners"]["smallest_adjusted_sample"] = {
                key: smallest[key] for key in ("a", "rho", "geometry", "side", "differentiation_spread", "minima")}
    admission = directory/"admission.txt"
    result["admission"] = admission.read_text().strip() if admission.is_file() else None
    return result


def main(root, output):
    inputs=read(root/"frozen-trials/frozen-trial-inputs.json")
    original=read(root/"frozen-trials/frozen-trial-results.json")
    builds={}
    for stage in BUILDS:
        at=root/"bin"/stage
        fingerprint, = re.findall(r"[0-9a-f]{64}",(at/"source-fingerprint.txt").read_text())
        builds[stage]={"revision":(at/"source.txt").read_text().strip(),
                       "solver":identity(at/"IncNavierStokesSolver"),
                       "probe":identity(at/"WJH27SteadyKrylovProbe"),
                       "source_fingerprint":fingerprint,
                       "configure_cost":timing(at/"configure.cost"),
                       "build_cost":timing(at/"build.cost"),"probe_cost":timing(at/"probe.cost"),
                       "files":present(at,("controls.json","executable.sha256","source-fingerprint.txt",
                                           "bundle.sha256","workers.ps"))}
    cases={}
    for name in CASES:
        at=root/name
        if not at.is_dir():
            continue
        record={"wall":timing(at/"cost.txt"),"equation_cost":timing(at/"equations.cost"),
                "margin_cost":timing(at/"margins.cost"),"point_cost":timing(at/"points.cost"),
                "files":present(at,("run.xml","run.log","exit-status","inputs.json","identities.sha256",
                                     "state.diagnostic-steady.json","state.steady.json","state.diagnostic-uncommitted.json",
                                     "equations.json","margins.json","points.json","admission.txt"))}
        if (at/"exit-status").is_file():
            record["exit_status"]=int((at/"exit-status").read_text())
        uncommitted=at/"state.diagnostic-uncommitted.json"
        if uncommitted.is_file():
            data=read(uncommitted)
            record["uncommitted"]={"custody":data["iterate_custody"],"original_residuals":data["original_residuals"],
                                   "coefficient_sha256":data["state"]["state_sha256"]}
        if (at/"state.steady.json").is_file():
            ledger=read(at/"state.steady.json")
            record["solver"]={name:ledger[name] for name in ("status","reason","last_accepted_norm",
                         "physical_time_advanced","total_seconds","accepted_newton_steps","jacobian_actions")}
        if (at/"state.diagnostic-steady.json").is_file():
            packet=read(at/"state.diagnostic-steady.json")
            record["equation_identity"]=packet["equation_identity"]
            record["coefficient_sha256"]=packet["state"]["state_sha256"]
            record["original_residuals"]=packet["original_residuals"]
            record["sparse_inverse_work"]=packet["sparse_inverse_work"]
        cases[name]=record
    lineages={}
    for case,name in LINEAGES.items():
        path=root/name
        if path.exists():
            run=read(path)
            for row in run["trials"]:
                stencils=row.pop("stencils")
                row["stencil_summary"]=[{"purpose":s["purpose"],"method":s["method"],"relative":s["relative"],
                                         "reduction":s["reduction"]} for s in stencils]
            lineages[case]={"receipt":identity(path),"run":run}
    candidates=[(max(row["residuals"]),name,row)
                for name,value in lineages.items() for row in value["run"]["accepted"]]
    best_norm, best_case, best=min(candidates,key=lambda item:item[0])
    qualifier=root/"endpoint-rich"
    independent=endpoint_audit(qualifier)
    qualified=independent.get("admission")=="QUALIFIED"
    profiles=root/"paired-level3"/"paired-comparison.json"
    final_refinement=[]
    for path in sorted((root/"model-refinement").glob("state.action-frozen-size-finite-trial-*.json")):
        d=read(path)
        final_refinement.append({"receipt":identity(path),"payload":d["raw_payload"],"measurements":d["measurements"],
                                 "notes":d["notes"],"proposal":d["proposal"]})
    frozen=original["trials"]
    frozen_useful=any(t["metrics"]["fraction"]==1 and t["notes"]["decision"]=="DECREASE_OBSERVED" for t in frozen)
    started=int((root/"start.epoch").read_text())
    elapsed=int(time.time())-started
    audit_end=qualifier/"finish.epoch"
    gates={
        "BEST_B3_STATE_REPRODUCED":{"status":"PASS","scope":"exact archived coefficients, scales, bindings and six native blocks"},
        "FROZEN_BOUNDED_NONLINEAR_TRIALS_ASSESSED":{"status":"PASS" if frozen_useful else "OPEN","scope":"four canonical archived directions on one base"},
        "BOUNDED_DIRECTION_POLICY_VERIFIED":{"status":"PASS" if (root/"resume-controls.json").is_file() else "OPEN",
                                            "scope":"focused native controls and actual original-equation trial custody"},
        "SUSTAINED_ORIGINAL_PROGRESS_ASSESSED":{"status":"PASS" if len(candidates)>=6 else "OPEN",
                                               "scope":"accepted nonlinear progress and measured tail"},
        "B3_ORIGINAL_ROOT_QUALIFIED":{"status":"PASS" if qualified else "OPEN",
                                      "scope":"independent exact-coefficient original admission"},
        "PAIRED_AXIAL_LEVEL_THREE_COMPLETE":{"status":"PASS" if profiles.is_file() else "OPEN",
                                             "scope":"retained A2/A3 and B2/B3 plus A3/B3 on frozen and owner coverage"},
    }
    for gate in ("NEXT_SPATIAL_DISCRIMINATOR_COMPLETED","POSITIVE_ORIGINAL_ANCHOR_RESOLUTION_SUPPORTED",
                 "ADVERSE_ORIGINAL_BRANCH_RESOLUTION_SUPPORTED","HEIGHT_CONTINUATION_RELEASED",
                 "FULL_HEIGHT_ORIGINAL_PHYSICAL_PILOT_PASSED","C1_CLOSED","C2_READY",
                 "LONG_FLAT_PRODUCTION_COST_CLOSED"):
        gates[gate]={"status":"OPEN","scope":"separate qualification remains pending"}
    gates["REMOTE_PUSH_VERIFIED"]={"status":"OPEN","scope":"await verification of the evidence publication"}
    snapshots=dict(initial=(root/"initial-workers.log").read_text(),
                   sampled_counts=[int(m) for path in root.glob("*.log") if path.is_file()
                   for m in re.findall(r"counted workers (\d+)",path.read_text(errors="replace"))])
    result={
        "schema":"WJH27B3BoundedNonlinearProgressV1",
        "reviewed_remote_head":"441f854236f10d96e514d48118978ddf5a757c33","campaign_root":str(root),
        "campaign":{"start":(root/"start.iso").read_text().strip(),"budget_seconds":10800,
                    "elapsed_evidence_seconds":elapsed,"extension_seconds":0,
                    "elapsed_to_audit_finish_seconds":int(audit_end.read_text())-started if audit_end.is_file() else None,
                    "workers":snapshots,"build_workers":16,"rank_threads":1},
        "fixed_original_inputs":inputs,"frozen_trials":original,
        "physical_time_advanced":0,"native_root_limit":1.e-10,"builds":builds,"cases":cases,
        "resume_control":read(root/"resume-controls.json") if (root/"resume-controls.json").is_file() else None,
        "lineages":lineages,"best_accepted":{"case":best_case,"native_max":best_norm,**best},
        "selected_endpoint":present(root,("selected-endpoint.json","endpoint-rich.setup.cost","final-audit-workers.ps")),
        "recycled_comparison":read(root/"recycled-comparison-inputs.json"),
        "optimizer_refinement":{"receipts":final_refinement,"files":present(root/"model-refinement",
                                 ("refinement-inputs.json","cost.txt","state.steady.json"))},
        "active_model":{"files":present(root,("active-fixed-model.json","active-helper-source.txt",
                                             "active-helper-bundle.sha256","active-model.cost")),
                        "proposal":read(root/"active-fixed-model.json") if (root/"active-fixed-model.json").is_file() else None,
                        "helper_revision":(root/"active-helper-source.txt").read_text().strip()},
        "independent_endpoint":independent,
        "spatial":{"A2_A3_retained_gap":.02551944209,"threshold_ceiling":.0007594305446050912,
                   "fixed_core":"73c264cf3d237311049fec71f07c32fa44a5b38c2cd2cd308c85103c9f62ab98",
                   "paired_result":identity(profiles) if profiles.is_file() else None},
        "replay_files":present(root,("build-stage.sh","run-finite.sh","run-finite-controls.sh",
                            "run-bounded.sh","run-resume-controls.sh","run-recycled.sh","run-recycled-pair.sh",
                            "run-model-refinement.sh","start-active.sh","run-active.sh","run-endpoint.sh","audit-best.sh",
                            "threads.sh","workers.sh","archive.sha256","paired-archive.sha256",
                            "collect_bounded.py","refine_model_candidate.py","package_evidence.py","plot_progress.py",
                            "pending-admission.txt")),
        "gates":gates,
    }
    with output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    main(args.root.resolve(),args.out)
