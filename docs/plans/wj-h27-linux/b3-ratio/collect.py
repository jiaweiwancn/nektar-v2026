#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect compact WJ-H27 ratio-policy evidence and verify accepted custody."""
import argparse
from pathlib import Path
import re
import sys
import time

import numpy as np

from c1_action_gap import read_audit
from c1_bounded_trials import SETTINGS, digest, load, packed, save
from resources import identity
from replay_ratio import SEED


def receipts(case):
    native = []
    for path in sorted(case.glob("state.action-*.json")):
        info = load(path)
        if info.get("schema") == "WJH27SteadyActionAuditV1":
            native.append((path,info))
    return sorted(native,key=lambda row:row[1]["capture_sequence"])


def compact(path, info):
    return {
        "receipt":identity(path), "payload":info["raw_payload"], "base":info["newton"],
        "proposal":info["proposal"],
        "base_coefficients":info["base_coefficient_sha256"], "factor_serial":info["factor_serial"],
        "factor":info["factor"], "policy":info["bounded_policy_sha256"],
        "span":info["notes"].get("matched_residual_span"),
        "measurements":info["measurements"], "notes":info["notes"],
        "blocks":{a["name"]:a["physical_blocks"] for a in info["arrays"] if "physical_blocks" in a},
        "stencils":[{k:s[k] for k in ("purpose","method","relative","reduction")} for s in info["stencils"]],
    }


def assemble(case, seed):
    origin = load(seed)
    previous = origin
    initial = max(origin["original_residuals"])
    start = int((case/"start.epoch").read_text())
    bases, trials, axes = [], [], []
    active_basis = None
    previous_receipt = None
    preceding_products = {}
    for path, info in receipts(case):
        assert info["correction_mass"]==0 and info["fixed_k_tangent"]
        assert info.get("bounded_acceptance","WJH27BoundedDirectionV2")=="WJH27BoundedDirectionV2"
        stage = info["stage"]
        row = compact(path, info)
        published = path.stat().st_mtime
        if stage=="bounded-current-basis":
            bases.append(row)
            active_basis = (path,info)
        elif stage=="bounded-direction":
            basis_path, basis_info = active_basis
            assert all(info[name]==basis_info[name] for name in (
                "base_coefficient_sha256","bounded_policy_sha256","equation_identity","free_id_sha256"))
            assert info["factor_serial"]==basis_info["factor_serial"]
            result = Path(str(basis_path)+".proposal-"+str(info["proposal"])+".json")
            model_settings = load(result)
            assert model_settings["basis_sha256"]==identity(basis_path)["sha256"]
            assert model_settings["radius"]==row["measurements"]["radius"]
            assert model_settings["declared_iterations"]==1000
            row["optimizer"] = {"receipt":identity(result),"result":model_settings}
            row["wall_from_preceding_receipt"] = published-previous_receipt
            measurements = row["measurements"]
            cumulative = measurements["current_products_at_base"]
            previous_products = preceding_products.get(info["newton"],0)
            row["incremental_current_products"] = cumulative-previous_products
            if "physical_actual" in measurements:
                decrease = measurements["physical_base"]-measurements["physical_actual"]
                row["reduction_per_incremental_current_product"] = decrease/max(1,cumulative-previous_products)
                row["reduction_per_trial_wall_second"] = decrease/max(.001,row["wall_from_preceding_receipt"])
            preceding_products[info["newton"]] = cumulative
            trials.append((path,info,row))
        elif stage=="bounded-residual-axes":
            assert info["factor_serial"]==1
            _, raw = read_audit(path)
            residual = raw["raw_base_residual"][:,0]
            offsets = np.r_[0,np.cumsum(info["free_counts"])]
            selections = [(slice(offsets[4],offsets[5]),"k"),
                          (slice(offsets[5],offsets[6]),"outer_epsilon")]
            if "scaled_momentum_rhs" in raw:
                selections.append((slice(0,offsets[2]),"momentum"))
            for sl, name in selections:
                expected = np.zeros_like(residual)
                expected[sl] = -residual[sl]
                assert np.array_equal(expected,raw["raw_negative_"+name+"_residual"][:,0])
                assert np.array_equal(expected*raw["row_scale"][:,0],raw["scaled_"+name+"_rhs"][:,0])
                coupled = raw["scaled_"+name+"_coupled_correction"][:,0]
                row[name+"_scaled_correction_by_unknowns"] = [
                    float(np.linalg.norm(v)) for v in np.split(coupled,offsets[1:-1])]
                row[name+"_scaled_factor_backsolve_error"] = float(np.linalg.norm(
                    raw["raw_"+name+"_factor_backsolve_error"][:,0]*raw["row_scale"][:,0]))
            axes.append(row)
        previous_receipt = published
    checkpoints, products = [], 0
    for path in sorted(case.glob("state.accepted-*.json")):
        packet = load(path)
        state, controller = packet["state"], packet["merit_controller"]
        radius = controller["state"]["radius"]
        free, ordering = packed(packet,state["state_sha256"])
        free_prev, _ = packed(previous,previous["state"]["state_sha256"])
        assert all(state[key]==origin["state"][key] for key in SETTINGS)
        assert packet["attempted_time"]==0 and packet["iterate_custody"]=="LAST_ACCEPTED"
        free_sha = digest(free.tobytes())
        matched = [(p,i,r) for p,i,r in trials if r["base_coefficients"]==previous["state"]["state_sha256"]
            and r["notes"]["decision"] in ("MATCHED_ELIGIBLE","BOUNDED_ACCEPTED")
            and any(a["name"]=="trial_free" and a["sha256"]==free_sha for a in i["arrays"])
            and r["measurements"]["radius"]==radius["used"]]
        assert matched
        attempted = max(r["measurements"]["current_products_at_base"]
                        for _,i,r in trials if i["newton"]==matched[0][1]["newton"])
        trial_path, meta, row = matched[0]
        assert ordering==meta["free_id_sha256"] and row["policy"]==controller["policy_sha256"]
        _, vectors = read_audit(trial_path)
        d = np.asarray(controller["accepted_displacement"],dtype="<f8")
        assert digest(d.tobytes())==controller["accepted_displacement_sha256"]
        assert np.array_equal(free-free_prev,d)
        assert np.array_equal(vectors["trial_free"][:,0],free)
        assert np.array_equal(vectors["raw_remainder"],vectors["raw_actual_residual"]-vectors["raw_model_residual"])
        norm = max(packet["original_residuals"])
        actual = row["measurements"]["physical_actual"]
        assert np.isclose(norm,actual,rtol=1e-12,atol=1e-16)
        assert row["measurements"]["model_ratio_accepted"]==1
        products += int(attempted)
        elapsed = path.stat().st_mtime-start
        checkpoints.append({"packet":identity(path),"coefficients":state["state_sha256"],
            "residuals":packet["original_residuals"],"selected_trial":identity(trial_path),"radius":radius,
            "actual_to_predicted_decrease":controller["state"]["actual_to_predicted_decrease"],
            "displacement_sha256":controller["accepted_displacement_sha256"],
            "policy_sha256":controller["policy_sha256"],"stagnation_updates":radius["window_updates"],
            "controller_accepted_total":controller["state"]["accepted"],
            "margins":{k:packet[k] for k in ("original_principal_over_mu","minimum_covariance_over_k",
                                           "physical_whole_cell_scalar_lower_bounds")},
            "K_over_Kin":state["interface_K_over_Kin"],"cumulative_current_products":products,
            "elapsed_to_checkpoint":elapsed,"cumulative_reduction_per_product":(initial-norm)/max(1,products),
            "cumulative_reduction_per_second":(initial-norm)/max(1.,elapsed)})
        previous = packet
    final = load(case/"state.diagnostic-steady.json")
    outcome = load(case/"state.steady.json")
    assert outcome["physical_time_advanced"]==0 and outcome["accepted_newton_steps"]==len(checkpoints)
    assert final["state"]["state_sha256"]==previous["state"]["state_sha256"]
    return {"seed":identity(seed),"initial_max":initial,"best_max":max(previous["original_residuals"]),
        "current_bases":bases,"trials":[r for _,_,r in trials],"residual_axes":axes,
        "accepted":checkpoints,"terminal":{k:outcome[k] for k in (
            "status","reason","physical_time_advanced","total_seconds","jacobian_actions")},
        "inverse_preparation":final["sparse_inverse_work"],"last_original_blocks":final["original_residuals"]}


def collect(root, output):
    sys.path.insert(0,str(root/"source/docs/plans/wj-h27-linux/b3-bounded"))
    from package_evidence import endpoint_audit, timing

    def present(path, names):
        return {n:identity(path/n) for n in names if (path/n).is_file()}

    old = root.parent/"wj-c1-b3-bounded.xjpRylk7"
    paired = root.parent/"wj-c1-paired-axial.UXeQ6Noe"
    start_epoch = int((root/"start.epoch").read_text())
    audit_finish = root/"endpoint-rich/finish.epoch"
    frozen = {name:load(root/name/"frozen-trial-results.json") for name in ("terminal","refined")
              if (root/name/"frozen-trial-results.json").is_file()}
    builds = {}
    for stage in ("ratio","coupled","coupled-probe-compile","coupled-probe-incomplete","momentum"):
        at = root/"bin"/stage
        if not (at/"finish.epoch").is_file():
            if (at/"build.cost").is_file():
                builds[stage]={"revision":(at/"source.txt").read_text().strip(),
                    "build_cost":timing(at/"build.cost"),
                    "log":identity(root/("build-"+stage+".log"))}
            continue
        fingerprint, = re.findall(r"[0-9a-f]{64}",(at/"source-fingerprint.txt").read_text())
        builds[stage] = {"revision":(at/"source.txt").read_text().strip(),"source_fingerprint":fingerprint,
            "solver":identity(at/"IncNavierStokesSolver"),"probe":identity(at/"WJH27SteadyKrylovProbe"),
            "build_cost":timing(at/"build.cost"),"receipts":present(at,(
                "bundle.sha256","controls.json","workers.ps","executable.sha256","source-fingerprint.txt"))}
    cases = {}
    for name in ("terminal","refined","retained","resume-control","changed-policy","scalars","momentum","endpoint-rich"):
        at = root/name
        if (at/"exit-status").is_file():
            cases[name] = {"wall":timing(at/"cost.txt"),"exit_status":int((at/"exit-status").read_text()),
                "equations_cost":timing(at/"equations.cost"),"margins_cost":timing(at/"margins.cost"),
                "files":present(at,("run.xml","run.log","inputs.json","identities.sha256","state.steady.json",
                    "state.diagnostic-steady.json","state.diagnostic-uncommitted.json","equations.json","margins.json",
                    "points.json","admission.txt"))}
    lineages = {}
    if "retained" in cases:
        lineages["retained"] = assemble(root/"retained",old/SEED)
    if "scalars" in cases:
        seed = Path(load(root/"matched-seed.json")["packet"]["path"])
        lineages["scalars"] = assemble(root/"scalars",seed)
    if "momentum" in cases:
        seed = Path(load(root/"momentum-decision.json")["packet"]["path"])
        lineages["momentum"] = assemble(root/"momentum",seed)
    best = min((item for run in lineages.values() for item in run["accepted"]),
               key=lambda a:max(a["residuals"]),default=None)
    endpoint = endpoint_audit(root/"endpoint-rich")
    matched = {(t["span"],t["measurements"]["radius"]) for t in
               lineages.get("scalars",{}).get("trials",[]) if t["base"]==0}
    complete_panel = all((span,radius) in matched for span in ("retained","k","outer_epsilon","joint")
                         for radius in (.00625,.025))
    reduced = {}
    for path in sorted((root/"reduced").glob("*.json")):
        reduced[path.stem] = {"data":load(path),"receipt":identity(path),
                             "cost":timing(path.with_suffix(".cost")) if path.name=="legacy-025.json"
                             else timing(Path(str(path)[:-5]+".cost"))}
    gates = {
        "BEST_B3_STATE_REPRODUCED":{"status":"PASS" if "terminal" in frozen else "OPEN",
                                    "scope":"signed base, reference and original coefficients and six-block replay"},
        "BOUNDED_ACCEPTANCE_SCALE_CONSISTENT":{"status":"PASS" if "terminal" in frozen and
            (root/"resume-controls.json").is_file() else "OPEN",
            "scope":"model-ratio finite replay, native controls and unchanged/changed cold resume"},
        "REDUCED_MODEL_ACCURACY_ASSESSED":{"status":"PASS" if "refined" in frozen else "OPEN",
            "scope":"same terminal current span and matched fresh nonlinear trials"},
        "CURRENT_RESIDUAL_ENRICHMENT_ASSESSED":{"status":"PASS" if complete_panel else "OPEN",
            "scope":"one shifted inverse, raw RHSs and eight held-base matched trials"},
        "SUSTAINED_ORIGINAL_PROGRESS_ASSESSED":{"status":"PASS" if "scalars" in lineages else "OPEN",
            "scope":"accepted decrease and measured stop, without root qualification"},
        "B3_ORIGINAL_ROOT_QUALIFIED":{"status":"PASS" if endpoint.get("admission")=="QUALIFIED" else "OPEN",
            "scope":"native and independent original tolerances and domain"},
    }
    for gate in ("PAIRED_AXIAL_LEVEL_THREE_COMPLETE","NEXT_SPATIAL_DISCRIMINATOR_COMPLETED",
                 "POSITIVE_ORIGINAL_ANCHOR_RESOLUTION_SUPPORTED","ADVERSE_ORIGINAL_BRANCH_RESOLUTION_SUPPORTED",
                 "HEIGHT_CONTINUATION_RELEASED","FULL_HEIGHT_ORIGINAL_PHYSICAL_PILOT_PASSED","C1_CLOSED","C2_READY",
                 "LONG_FLAT_PRODUCTION_COST_CLOSED","REMOTE_PUSH_VERIFIED"):
        gates[gate]={"status":"OPEN","scope":"separate prerequisite remains pending"}
    save(output,{"schema":"WJH27B3ModelRatioRecoveryV1",
        "reviewed_remote_head":"a80496e6a9fd57053e825c46d995511ad268b18b","campaign_root":str(root),
        "campaign":{"start":(root/"start.iso").read_text().strip(),"budget_seconds":10800,
                    "elapsed_evidence_seconds":int(time.time())-start_epoch,
                    "elapsed_through_audit_seconds":int(audit_finish.read_text())-start_epoch
                        if audit_finish.is_file() else None,
                    "build_workers":{"ratio":2,"coupled":2,"momentum":4},
                    "threads_per_rank":1,"extension_seconds":0},
        "physical_time_advanced":0,"native_original_limit":1e-10,"builds":builds,
        "base_archive":present(root/"terminal",("frozen-trial-inputs.json",)),
        "reused_r10_receipts":present(old/"endpoint-rich",(
            "state.diagnostic-steady.json","equations.json","margins.json","points.json")),
        "legacy_optimizer":identity(old/"optimizer-source/solvers/IncNavierStokesSolver/Tests/RANS/h27/"
                                    "c1_bounded_active_optimizer.py"),
        "r10_replay_executable":identity(old/"bin/recycled/IncNavierStokesSolver"),
        "helper_identities":present(root/"source/solvers/IncNavierStokesSolver/Tests/RANS/h27",
            ("c1_bounded_active_optimizer.py","c1_bounded_optimizer.py","c1_recycled_basis.py")),
        "joint_model":identity(root/"source/docs/plans/wj-h27-linux/b3-frozen-size/frozen_size.py"),
        "execution_inputs":present(root,("threads.sh","capacity.sh","build.sh","phase-one.sh","phase-two.sh",
                                         "audit.sh","solve.sh","matched-seed.json","endpoint-selection.json",
                                         "momentum-decision.json","phase-momentum.sh")),
        "worker_receipts":{p.name:identity(p) for p in root.glob("*workers.ps")},
        "frozen":frozen,"reduced":reduced,"cases":cases,"lineages":lineages,
        "best_accepted":best,"endpoint":endpoint,"restart":present(root,("resume-controls.json",)),
        "spatial":{"fixed_core_sha256":"73c264cf3d237311049fec71f07c32fa44a5b38c2cd2cd308c85103c9f62ab98",
                   "retained_A2_A3_gap":.02551944209,"tolerance_ceiling":.0007594305446050912,
                   "A3":present(paired/"A3-rich",("state.diagnostic-steady.json","margins.json","equations.json"))},
        "gates":gates})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=Path)
    p.add_argument("output",type=Path)
    args=p.parse_args()
    collect(args.root.resolve(),args.output)
