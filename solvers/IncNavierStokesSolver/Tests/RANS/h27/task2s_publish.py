#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal compact 2R-S evidence and a complete prospective radial-run admission."""

import argparse
import json
from pathlib import Path
import tarfile

from resources import atomic, identity
from task2s_contract import CAPS


def read(path):
    return json.loads(path.read_text())


def first_radial_forecast(root, out, book):
    floor = 300
    caps = CAPS["S2_first"]
    total_reservation = sum(seconds+7 for seconds in caps.values())
    settled = sum(v["seconds"] for v in book["entries"].values() if v["status"] != "RESERVED")
    account = dict(charged=book["opening_seconds"]+settled,
                   staged=book["tranche_seconds"]-settled,
                   nonprotected=book["ceiling"]-book["protected_seconds"]-book["opening_seconds"]-settled)
    records = root/"science/task2rs-s1"
    replay = read(records/"observation-replay.json")
    wedge = read(records/"wedge-fit.json")
    assert replay["archive_replay"]["status"] == "PASS" and len(replay["rows"]) == 690 and len(replay["failures"]) == 28
    assert wedge["fitted_shape_pencil"]["status"] == "PASS"
    witnesses = [
        root/"operations/r3-terminal-steady/receipt.json",
        root/"operations/r3-terminal-current/receipt.json",
        root/"operations/r3-outlet-refine-steady/receipt.json",
        root/"operations/r3-outlet-lean-steady/receipt.json",
        root/"operations/s1-spatial-replay-corrected/receipt.json",
        root/"operations/s1-boundary-and-admissibility/receipt.json",
        root/"operations/s1-wedge-and-fit-serialization-retry/receipt.json",
    ]
    want = total_reservation+floor
    atomic(out,dict(schema="WJH27Task2RSRadialAdmissionV1",
        live_before_publication_charge=account, capped_operations_seconds=caps,
        total_reservation_seconds=total_reservation, evidence_hold_seconds=floor,
        staged_requirement_including_evidence_hold=want,
        minimum_new_stage_before_publication_charge=max(0.,want-account["staged"]),
        memory=dict(measured_terminal_peak=8005373952,science_rss_cap=8*1024**3,
            estimated_peak_bytes=round(8005373952*62/60+.25*1024**3),
            basis_dimension_previous=40,basis_dimension_next=30,
            largest_unknown="Sparse factor fill and peak coexistence on the changed wall mesh."),
        case=dict(cells=[114,62],polynomial_degrees=[3,2],native_quadrature=[6,8],
                  terminal_x_width=.00244140625,first_radial_width=.000125,
                  internal_solver_cap=980,original_exterior_contract=True,time_code=0.),
        evidence=[identity(p) for p in [*witnesses,records/"observation-replay.json",
                  records/"boundary-audit.json",records/"wedge-fit.json"]],
        interpretation=dict(contracting="Measure all rows against both neighbors and test a bounded sharper x/r level.",
           rss_stop="Use measured quiet strips to fit the local splits, prechecking inlet representation.",
           persistent_strong_divergence="Check near-wall and collar gradients, and require common resolved amplitude before claiming a singular obstruction.",
           worker_preflight="Defer and retry on the same ledger when worker headroom returns.")))


def joint_forecast(root):
    base = root/"science/task2rs-s2-radial"
    design = root/"science/task2rs-s2-joint-design/forecast.json"
    scope = read(design)
    last = root/"operations/s2-radial-steady/receipt.json"
    peak = read(last)["result"]["peak_tree_rss_bytes"]
    original = read(base/"inputs.json")["cells"]
    nx,nr = scope["cells"]
    predicted = round(peak*nx*nr/(original[0]*original[1])+.25*1024**3)
    if predicted >= 8*1024**3 or not scope["collar_h_split"]:
        raise ValueError("joint candidate exceeds the estimated science footprint or loses the collar control")
    caps = dict(inputs=45,steady=1100,current=255,radial_pair=155,terminal_pair=155,
                local_diagnostics=100,wedge=33,radial_projection=35,joint_projection=35,publication=100)
    required = sum(t+7 for t in caps.values())+300
    if read(base/"terminal-comparison.json")["comparison_status"] != "NOT_MET":
        raise ValueError("joint diagnostic admission depends on the observed radial-only shortfall")
    atomic(root/"evidence/task2rs-s2-joint-forecast.json",dict(
        schema="WJH27Task2RSJointAdmissionV1", capped_operations_seconds=caps,
        staged_requirement_including_evidence_hold=required,evidence_hold_seconds=300,
        selected_geometry=identity(design),previous_memory_receipt=identity(last),
        memory=dict(peak_bytes_measured=peak,predicted_with_contingency=predicted,
                    contingency_bytes=round(.25*1024**3),limit_bytes=8*1024**3,
                    uncertainty="Small-cell aspect ratio can change pivoting and sparse factor fill. Runtime guard is mandatory."),
        dependencies=[identity(base/p) for p in ("terminal-comparison.json","parent-comparison.json",
                    "independent.json","boundary-audit.json","wedge-fit.json")],
        planned_interpretation="Compare both previous terminal owners and preserve every quiet-strip control. "
            "Use three distinct x/r levels for the predeclared amplitude criterion. "
            "If regular-like contraction is absent, distinguish bandwidth sensitivity, collar strong divergence "
            "and the conditional lower-wall covariance limit before advancing a diagnosis."))


def tip_forecast(root, lean=False):
    base = root/"science/task2rs-s2-radial"
    joint = root/"science/task2rs-s2-joint"
    name = "tip-lean" if lean else "tip"
    design = root/f"science/task2rs-s2-{name}-design/forecast.json"
    scope = read(design)
    check = read(joint/"radial-comparison.json")
    wedge = read(joint/"wedge-fit.json")
    if wedge["status"]!="INCONCLUSIVE" or check["comparison_status"]!="NOT_MET":
        raise ValueError("isolated tip admission requires unresolved joint results")
    original = read(base/"inputs.json")["cells"]
    if scope["collar_h_split"] or not scope["relocated_tip"] or scope["cells"][1] != original[1]:
        raise ValueError("tip test must retain the parent radial grid outside the first two intervals")
    last = root/"operations"/("s2-tip-steady" if lean else "s2-radial-steady")/"receipt.json"
    previous = read(last)["result"]
    if lean:
        if previous["reason"]!="RSS_LIMIT" or not scope["removed"] or scope["cells"][0]>=original[0]:
            raise ValueError("lean tip retry requires the measured RSS stop and fewer x cells")
    elif scope["cells"] != original or scope["removed"]:
        raise ValueError("the isolated tip control retains the parent grid population")
    predicted = previous["peak_tree_rss_bytes"]*scope["cells"][0]/original[0]+.125*1024**3
    if predicted >= 8*1024**3:
        raise ValueError("insufficient estimated memory headroom")
    caps = dict(inputs=45,steady=1000,current=255,radial_pair=155,joint_pair=155,
                local_diagnostics=110,wedge=40,projection=35,publication=95)
    if lean:
        caps.update(steady=1090,radial_pair=160,joint_pair=160,near_wedge=55)
    decision = ("Retest both tip widths with their original aspect ratio and unchanged collar. "
                "Reduce the x cell count through measured quiet-strip merges after the RSS stop. "
                "Include each strip in the comparison and repeat the narrower angular analysis."
                if lean else
                "Halve the last x and first r widths at the same aspect ratio by moving just their adjacent grid lines. "
                "This preserves the collar and all upstream x edges of the 114 by 62 parent. "
                "Measure whether the new widespread collar differences from the joint test recur, "
                "and whether the retained problem has a stable angular mode across three resolved states.")
    atomic(root/f"evidence/task2rs-s2-{name}-forecast.json",dict(
        schema="WJH27Task2RSIsolatedTipAdmissionV1", selected_geometry=identity(design),
        capped_operations_seconds=caps,staged_requirement_including_evidence_hold=sum(t+7 for t in caps.values())+300,
        evidence_hold_seconds=300,previous_memory_receipt=identity(last),
        memory=dict(estimated_peak_bytes=round(predicted),limit_bytes=8*1024**3,
                    uncertainty="Scale the measured prior peak by the x population. Small-cell pivot fill may vary."),
        measured_dependencies=[identity(joint/name) for name in ("radial-comparison.json",
                               "boundary-audit.json","wedge-fit.json")],
        decision=decision))


def resume_forecast(root):
    base = root/"science/task2rs-s2-tip-lean"
    run = root/"operations/s2-tip-lean-steady/receipt.json"
    peak = read(run)["result"]["peak_tree_rss_bytes"]
    steady = read(base/"state.steady.json")
    seed = base/"state.steady-candidate.json"
    packet = read(seed)
    if ("DEADLINE" not in steady["reason"] or steady["accepted_newton_steps"]<1
            or steady["true_norm"]>1e-7 or packet["qualified_reference"]
            or packet["state"]["step"]!=0 or packet["state"]["time"]!=0):
        raise ValueError("continuation needs a contracting, typed zero-time candidate")
    upper_unknowns = sum(packet["state"]["sizes"])+packet["state"]["sizes"][0]
    basis_increment = 2*(64-30)*upper_unknowns*8
    predicted = peak+basis_increment+.125*1024**3
    if predicted >= 8*1024**3:
        raise ValueError("Krylov augmentation exceeds estimated RSS headroom")
    caps = dict(inputs=45,steady=890,current=255,radial_pair=160,joint_pair=160,
                local_diagnostics=110,wedge=40,near_wedge=55,projection=35,publication=95)
    atomic(root/"evidence/task2rs-s2-tip-resume-forecast.json",dict(
        schema="WJH27Task2RSTypedSteadyContinuationV1",previous=identity(run),steady=identity(base/"state.steady.json"),
        seed=identity(seed),fixed_mesh=read(base/"inputs.json")["geometry"],same_current_equations=True,
        capped_operations_seconds=caps,staged_requirement_including_evidence_hold=sum(t+7 for t in caps.values())+300,
        evidence_hold_seconds=300,memory=dict(prior_peak_bytes=peak,maximum_unknowns=upper_unknowns,
        krylov_dimensions=[30,64],two_dense_basis_increment_bytes=basis_increment,
        contingency_bytes=round(.125*1024**3),estimated_peak_bytes=round(predicted),limit_bytes=8*1024**3),
        reason="Restart Newton from the typed, unqualified same-model iterate. "
               "The deadline left a 236-fold residual contraction and a long restarted inner Krylov solve. "
               "Increase the basis space with measured RSS headroom. Rebuild lifts, current data and histories."))


def nested_forecast(root,resume=None):
    base = root/"science/task2rs-s2-tip-resume"
    new_inputs = read(base/"inputs.json")
    old_nx,old_nr = new_inputs["cells"]
    next_cells = [old_nx+1,old_nr+1]
    run = root/"operations/s2-tip-resume-steady/receipt.json"
    peak = read(run)["result"]["peak_tree_rss_bytes"]
    bound = peak*next_cells[0]*next_cells[1]/(old_nx*old_nr)+.2*1024**3
    pair = read(base/"radial-comparison.json")
    wedge,near=read(base/"wedge-fit.json"),read(base/"wedge-near.json")
    if (bound>=8*1024**3 or pair["comparison_status"]!="NOT_MET"
            or pair["redistributed_strip_controls"]["failures"]
            or wedge["status"]!="INCONCLUSIVE" or near["status"]!="INCONCLUSIVE"):
        raise ValueError("no basis for an additional admitted nested control")
    caps = dict(inputs=45,steady=1030,current=250,pair=160,local_diagnostics=105,
                wedge=40,near_wedge=55,projection=35,publication=95)
    evidence = [identity(base/f) for f in ("radial-comparison.json","wedge-fit.json","wedge-near.json")]
    name = "nested"
    if resume:
        name = "nested-finish" if resume=="finish" else "nested-resume"
        previous = "nested-resume" if resume=="finish" else "nested"
        prev = root/f"science/task2rs-s2-{previous}"
        report = read(prev/"state.steady.json")
        seed = prev/"state.steady-candidate.json"
        state = read(seed)
        if ("DEADLINE" not in report["reason"] or report["accepted_newton_steps"]<1
                or report["true_norm"]>=(1e-8 if resume=="finish" else 1e-7) or state["qualified_reference"]
                or state["state"]["step"]!=0 or state["state"]["time"]!=0):
            raise ValueError("nested continuation lacks meaningful contraction on a typed zero-time iterate")
        run = root/f"operations/s2-{previous}-steady/receipt.json"
        peak = read(run)["result"]["peak_tree_rss_bytes"]
        bound = peak+.125*1024**3
        if bound >= 8*1024**3:
            raise ValueError("measured nested peak leaves insufficient estimated headroom")
        caps["steady"] = 670 if resume=="finish" else 790
        new_inputs = read(prev/"inputs.json")
        evidence += [identity(seed),identity(prev/"state.steady.json")]
    atomic(root/f"evidence/task2rs-s2-{name}-forecast.json",dict(
        schema="WJH27Task2RSNestedTipAdmissionV1",current_mesh=new_inputs["geometry"],
        next_cells=next_cells,capped_operations_seconds=caps,
        staged_requirement_including_evidence_hold=sum(t+7 for t in caps.values())+300,
        evidence_hold_seconds=300,prior_receipt=identity(run),
        memory=dict(prior_peak_bytes=peak,predicted_peak_bytes=round(bound),limit_bytes=8*1024**3),
        dependencies=evidence,typed_continuation=resume is not None,
        interpretation="Add one bisection in each terminal direction while keeping every existing edge and the upstream inventory. "
                       "Compare the scale response directly with the preceding corner-only pair. "
                       "Repeat both broad and narrow angular checks on a sequence whose tip contracts in both directions. "
                       "Use the actual nested comparison and resolved fits to distinguish continued under-resolution from a supported mode."))


def comparison_summary(path, original_failures):
    value = read(path)
    rows = value["rows"]
    worst = {}
    for field in sorted({name.split("/")[-1] for name in rows}):
        group = {k:v for k,v in rows.items() if k.endswith("/"+field)}
        peak = max(group,key=lambda name:group[name]["MAX"])
        rms = max(group,key=lambda name:group[name]["RMS"])
        worst[field] = dict(MAX=group[peak]["MAX"], MAX_row=peak, RMS=group[rms]["RMS"],
                            RMS_row=rms, limits=group[peak]["limits"], peak=group[peak].get("peak"))
    extra = value.get("redistributed_strip_controls",{})
    return dict(witness=identity(path),first=value["first"],second=value["second"],
        count=len(rows),status=value["comparison_status"],failed=value["failures"],extrema=worst,
        original_28={name:{key:rows[name][key] for key in ("RMS","MAX","limits","status")}
                     for name in original_failures},
        quiet_strip_inventory=len(extra.get("rows",{})),quiet_strip_failures=extra.get("failures",[]))


def s2_summary(root, out, book):
    prior = root/"science/task2rs-s1/observation-replay.json"
    old = read(prior)["failures"]
    cases, comparisons = {}, {}
    for folder in sorted((root/"science").glob("task2rs-s2-*")):
        xml = folder/"inputs.json"
        if not xml.exists():
            continue
        inputs = read(xml)
        result = dict(inputs=identity(xml), cells=inputs["cells"],quadrature=inputs["quadrature"],
                      mode_counts=inputs["modes"],tolerances=inputs["tolerances"],time_code=inputs["time_code"])
        for name,key in (("state.steady.json","native_steady"),("state.result.json","native_state"),
                         ("state.steady-candidate.json","typed_unqualified_seed"),
                         ("independent.json","current"),("projection.json","projection"),
                         ("boundary-audit.json","local"),("wedge-fit.json","wedge"),("wedge-near.json","near_wedge")):
            path = folder/name
            if path.exists():
                result[key] = identity(path)
        if "native_steady" in result:
            s = read(folder/"state.steady.json")
            result["steady"] = {key:s[key] for key in ("status","reason","true_norm","total_seconds",
                                            "accepted_newton_steps","jacobian_actions")}
        if "current" in result:
            independent = read(folder/"independent.json")
            core = independent["core"]
            result.update(physical=independent["physical"],equations=independent["equations"],
                current_residuals=core["residuals"],physical_envelopes=core["physical_envelope"],
                ground_against_Tref=core["ground_shear_diagnostic"],
                weak_identity_error=core["accounting_identity_max_error"],stations=independent["observations"],
                owner_flux={key:item for key,item in independent["interfaces"]["norms"].items() if key.startswith("H/")})
        if "local" in result:
            src = read(folder/"boundary-audit.json")
            state = src["cases"][folder.name]
            result.update(continuity=state["continuity"],vertex=state["weak_natural_boundary"]["raw_vertex"],
                          admissibility_samples=state["samples"],collar=src["collar"])
        for witness,path in (("wedge","wedge-fit.json"),("near_wedge","wedge-near.json")):
            if witness not in result:
                continue
            w = read(folder/path)
            result[witness+"_outcome"] = dict(status=w["status"],reason=w["reason"],
                resolved_fit_point_counts={name:s["common_points"] for name,s in w["studies"].items()},
                fits=[dict(study=name,window=window, amplitude=fit["amplitude"],
                           uncertainty=fit["sampled_mesh_amplitude_bound"],
                           heldout={g:(v["singular_and_regular"]["budget_RMS"],v["regular_only"]["budget_RMS"])
                                    for g,v in fit["heldout"].items() if g!="training"})
                      for name,s in w["studies"].items() for window,info in s["fits"].items()
                      for fit in info["levels"] if fit["status"]=="FIT"])
        for path in sorted(folder.glob("*comparison.json")):
            comparisons[folder.name+"/"+path.name] = comparison_summary(path,old)
        cases[folder.name] = result
    base = read(root/"evidence/s0-baseline.json")
    relevant = {label:dict(phase=v["phase"],status=v["status"],seconds=v["seconds"],
                                 peak_tree_rss_bytes=v["peak_tree_rss_bytes"],
                                 peak_host_workers=v["peak_host_workers"],reason=v["reason"],
                                 receipt=identity(root/"operations"/label/"receipt.json"))
                for label,v in book["entries"].items() if label.startswith("s2-") and v["status"]!="RESERVED"}
    finest = cases["task2rs-s2-nested-finish"]
    pair = comparisons["task2rs-s2-nested-finish/nested-comparison.json"]
    if (finest["equations"]!="MET" or finest["physical"]!="MET" or pair["status"]!="NOT_MET"
            or finest["wedge_outcome"]["status"]!="INCONCLUSIVE"
            or finest["near_wedge_outcome"]["status"]!="INCONCLUSIVE"):
        raise ValueError("the final measured disposition requires review")
    nested_solver_labels = ("s2-nested-steady","s2-nested-resume-steady","s2-nested-finish-steady")
    final_analysis_labels = ("s2-nested-finish-current-verified-path","s2-nested-finish-pair",
        "s2-nested-finish-local","s2-nested-finish-wedge","s2-nested-finish-near","s2-nested-finish-projection")
    elapsed = sum(relevant[k]["seconds"] for k in (*nested_solver_labels,*final_analysis_labels))
    previous_peak = max(relevant[k]["peak_tree_rss_bytes"] for k in nested_solver_labels)
    nx,nr = finest["cells"]
    spent = sum(entry["seconds"] for entry in book["entries"].values() if entry["status"]!="RESERVED")
    live = book["ceiling"]-book["protected_seconds"]-book["opening_seconds"]-spent
    next_bound = dict(next_nested_cells=[nx+1,nr+1],measured_solver_and_final_analysis_seconds=elapsed,
                      initial_setup_and_future_publication_are_additional=True,evidence_hold_seconds=300,
                      live_nonprotected_seconds_before_publication=live,
                      worst_measured_nested_peak_bytes=previous_peak,
                      cell_scaled_peak_before_margin_bytes=round(previous_peak*(nx+1)*(nr+1)/(nx*nr)),
                      science_cap_bytes=8*1024**3)
    atomic(out,dict(schema="WJH27Task2RSResolutionEvidenceV1",source=base["compiled_source"],
        application=base["artifacts"]["application"],reference=base["artifacts"]["table"],geometry=base["geometry"],
        original_28=old, original_690_witness=identity(prior),cases=cases,comparisons=comparisons,operations=relevant,
        next_like_for_like_admission=next_bound,
        checkpoint_status=dict(S0="PASS",S1="PASS",S2="INCONCLUSIVE",S3="NOT_RUN",S4="NOT_RUN"),
        gates=dict(DIAGNOSIS="INCONCLUSIVE",ORIGINAL_PROBLEM_SPATIAL_QUALIFICATION="NOT_MET",
                   RECALIBRATED_NATIVE_REFERENCE="NOT_MET",SHORT_SCREEN="NOT_RUN",COMPLETE_COST="NOT_RUN",
                   TASK2_COMPLETE="NO",TASK3_READY="NO")))


def package(root, stage):
    book = read(root/"accounting.json")
    active = {k for k,v in book["entries"].items() if v["status"] == "RESERVED"}
    if any(not name.startswith(stage+"-publish") for name in active):
        raise ValueError(f"still reserved at publication: {sorted(active)}")
    evidence = root/"evidence"
    if stage == "s1":
        first_radial_forecast(root,evidence/"task2rs-s2-first-forecast.json",book)
    if stage == "s2":
        s2_summary(root,evidence/"task2rs-s2-summary.json",book)
    paths = {Path(__file__).resolve(),evidence/"s0-baseline.json"}
    if stage == "s1":
        paths.add(evidence/"task2rs-s2-first-forecast.json")
    if stage == "s2":
        paths.add(evidence/"task2rs-s2-summary.json")
    paths.update((root/"tools").glob("s-support*/h27/task2s*.py"))
    paths.update((root/"tools").glob("s-support*/h27/physical_oracle.py"))
    paths.update((root/"tools").glob("s-support*/h27/task2r_physical-audit.py"))
    paths.update((root/"tools").glob("s-support*/h27/task2r_verify.py"))
    paths.update((root/"tools").glob("s-support*/cao-tamura-fig3/rough-fig3.csv"))
    paths.update(p for p in evidence.glob("task2rs-*-forecast.json"))
    for label,event in book["entries"].items():
        if label.startswith((stage+"-","s0-")) and event["status"]!="RESERVED":
            directory = root/"operations"/label
            paths.update(directory/name for name in ("request.json","receipt.json","command.log"))
    for directory in (root/"science").glob("task2rs-"+stage+"*"):
        for p in directory.iterdir():
            if p.is_file() and p.suffix in (".json",".xml",".opt") and p.name not in (
                    "state.result.json","state.steady-candidate.json","state.failure.json"):
                paths.add(p)
    paths = sorted(paths)
    manifest = evidence/f"task2rs-{stage}-manifest.json"
    atomic(manifest,dict(schema="WJH27Task2RSArchiveV1",stage=stage,
          files={str(p.relative_to(root)):identity(p) for p in paths},
          field_storage="The native states and raw histories stay in the named ext4 case directories."))
    archive = evidence/f"task2rs-{stage}.tar.gz"
    with tarfile.open(archive,"x:gz",compresslevel=5) as stream:
        for p in [*paths,manifest]:
            stream.add(p,arcname=str(p.relative_to(root)),recursive=False)
    atomic(evidence/f"task2rs-{stage}-bundle.json",dict(archive=identity(archive),manifest=identity(manifest)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=Path)
    p.add_argument("stage",choices=("s1","s2","s3","joint-forecast","tip-forecast","tip-lean-forecast",
                                    "tip-resume-forecast","nested-forecast","nested-resume-forecast",
                                    "nested-finish-forecast"))
    args=p.parse_args()
    if args.stage=="joint-forecast":
        joint_forecast(args.root)
    elif args.stage=="tip-forecast":
        tip_forecast(args.root)
    elif args.stage=="tip-lean-forecast":
        tip_forecast(args.root,lean=True)
    elif args.stage=="tip-resume-forecast":
        resume_forecast(args.root)
    elif args.stage=="nested-forecast":
        nested_forecast(args.root)
    elif args.stage=="nested-resume-forecast":
        nested_forecast(args.root,resume="resume")
    elif args.stage=="nested-finish-forecast":
        nested_forecast(args.root,resume="finish")
    else:
        package(args.root,args.stage)


if __name__=="__main__":
    main()
