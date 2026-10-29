#!/usr/bin/env python3
"""Linux acceptance of full CLS transport and the conservative original AWF."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {k:"1" for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS")}
os.environ.update(THREADS)
if sys.platform.startswith("linux"):
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[-16:])
from generate_m4 import awf_reference, generate_cls, info, parameter


def require(condition, message):
    if not condition: raise RuntimeError(message)


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    require(sys.platform.startswith("linux"), "Authoritative validation requires Linux")
    root = Path(__file__).resolve().parents[4]
    build = args.build.resolve();output = args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    solver = build/"IncNavierStokesSolver";audit = build/"RANSFieldAudit"
    binary_hashes = {p.name:digest(p) for p in (solver,audit)}
    records = {}; commands = []

    def execute(command,cwd,label,reject=None,timeout=900):
        command = list(map(str,command));print("Running "+label,flush=True)
        commands.append({"command":command,"cwd":str(cwd),"timeout":timeout})
        process = subprocess.run(command,cwd=cwd,capture_output=True,text=True,timeout=timeout,env=os.environ.copy())
        log = process.stdout+process.stderr;(cwd/(label+".log")).write_text(log)
        require((process.returncode != 0 and reject in log) if reject else process.returncode == 0,
                f"{label}: inspect {cwd/(label+'.log')}")
        return log

    def run(name,ranks=1,modify=None,reject=None,**options):
        directory = output/name;directory.mkdir(exist_ok=True)
        session = directory/(name+".xml");generate_cls(session,**options)
        tree = ET.parse(session)
        if ranks > 1:info(tree.getroot(),"GlobalSysSoln","IterativeStaticCond")
        if modify:modify(tree.getroot())
        tree.write(session,encoding="utf-8",xml_declaration=True)
        fingerprint = {"fixture_sha256":digest(session),"binaries":binary_hashes,"ranks":ranks,"rejection":reject}
        saved = directory/"case.json"
        if saved.exists() and json.loads(saved.read_text()).get("fingerprint") == fingerprint:
            record = json.loads(saved.read_text());print("Reusing "+name,flush=True)
        else:
            command = (["mpirun","--bind-to","none","--map-by","slot","-np",str(ranks)] if ranks>1 else [])+[solver,session.name]
            log = execute(command,directory,name,reject=reject)
            record = {"fingerprint":fingerprint,"options":options,"ranks":ranks,"solver_command":list(map(str,command))}
            if reject:record["expected_rejection"] = reject
            else:
                require("Total Computation Time" in log,"Solver did not finish")
                raw = execute([audit,session.name,name+".fld",name+".fld"],directory,"audit")
                record["field_audit"] = json.loads(raw[raw.index('{"time"'):])
                rows = [{k:float(v) for k,v in row.items()} for row in csv.DictReader((directory/(name+".rans.csv")).open())]
                require(len(rows) == options.get("steps",4),"Incomplete stage diagnostics")
                require(all(math.isfinite(v) for row in rows for v in row.values()),"Nonfinite stage diagnostics")
                require(all(row["min_k"]>0 and row["min_eps"]>0 and row["guard_activations"]==0 for row in rows),"Unexpected state guard")
                record["stages"] = rows
                require(max(row["live_global_systems"] for row in rows)<=rows[0]["live_global_systems"]+2,"Growing matrix cache")
                wall = directory/(name+".awf.csv")
                if wall.exists():
                    wall_rows = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(wall.open())]
                    require(len(wall_rows)==len(rows),"Incomplete wall diagnostics")
                    require(all(math.isfinite(v) for row in wall_rows for v in row.values()),"Nonfinite wall diagnostic")
                    require(all(row["min_kP"]>0 and row["min_epsP"]>0 and row["chi"]==0 and row["guard_activations"]==0 for row in wall_rows),"Original wall state or branch failure")
                    require(max(row["epsilonH_error"] for row in wall_rows)<1e-9,"Dissipation algebraic boundary failed")
                    require(max(row["mass_error"] for row in wall_rows)<1e-7,"Integrated strip continuity failed")
                    require(all(row["budget_available"]==1 for row in wall_rows),"Missing constant-test budget")
                    for row in wall_rows:
                        momentum_scale = max(1.,abs(row["traction_x"]),abs(row["traction_y"]))
                        k_scale = max(1.,2*row["area"]*(row["mean_production"]+row["mean_dissipation"]))
                        require(max(abs(row["momentum_budget_x"]),abs(row["momentum_budget_y"]))/momentum_scale<2e-7,"Momentum inventory/flux budget failed")
                        require(abs(row["k_budget"])/k_scale<2e-7,"Integrated k inventory/source budget failed")
                    record["wall_stages"] = wall_rows
                record["artifacts_sha256"] = {p.name:digest(p) for p in directory.iterdir() if p.is_file() and p.suffix in (".xml",".csv",".his",".log")}
            saved.write_text(json.dumps(record,indent=2)+"\n")
        records[name] = record
        (output/"progress.json").write_text(json.dumps(records,indent=2)+"\n")
        if "field_audit" in record:
            print(name+": "+str({k:v["l2_error"] for k,v in record["field_audit"]["fields"].items()}),flush=True)
        return record

    temporal = []
    for dt in (.025,.0125,.00625):
        temporal.append(run(f"cls_dt{dt:g}",dt=dt,steps=round(.1/dt),modes=5,points=12))
    rates = {}
    for field in ("u","v","k","eps_tilde","p"):
        errors = [r["field_audit"]["fields"][field]["l2_error"] for r in temporal]
        rates[field] = [math.log(a/b,2) for a,b in zip(errors,errors[1:])]
        require(min(rates[field])>1.65,f"CLS temporal convergence failed: {field} {rates[field]}")
    for name in ("rans_Yap","rans_quadratic_norm","rans_cubic_norm"):
        require(temporal[-1]["field_audit"]["fields"][name]["max"]>1e-5,"Manufactured flow does not exercise "+name)
    run("cls_model_output_p8",dt=.00625,steps=16,modes=9,points=16)
    require(records["cls_model_output_p8"]["field_audit"]["fields"]["rans_Yap"]["l2_error"]<
            .1*temporal[-1]["field_audit"]["fields"]["rans_Yap"]["l2_error"],"Derived-field polynomial projection failed to refine")
    run("cls_molecular",dt=.0125,steps=8,modes=5,points=12,diffusion="Molecular")
    run("cls_quadrature",dt=.0125,steps=8,modes=5,points=16)
    run("cls_mpi",ranks=4,dt=.0125,steps=8,modes=5,points=12)
    run("cls_3d",ranks=4,dimension=3,dt=.01,steps=2,modes=5,points=10)
    require(max(v["l2_error"] for k,v in records["cls_3d"]["field_audit"]["fields"].items() if not k.startswith("rans_"))<1e-4,"Full 3D manufactured CLS failed")
    run("cls_missing_cubic",dt=.00625,steps=16,modes=5,points=12,omit_cubic=True)
    require(records["cls_missing_cubic"]["field_audit"]["fields"]["p"]["l2_error"]>5*temporal[-1]["field_audit"]["fields"]["p"]["l2_error"],"Missing cubic stress forcing not detected")
    run("cls_missing_yap",dt=.0125,steps=8,modes=5,points=12,modify=lambda r:parameter(r,"RANSYapCoefficient",0))
    require(records["cls_missing_yap"]["field_audit"]["fields"]["eps_tilde"]["l2_error"]>10*temporal[1]["field_audit"]["fields"]["eps_tilde"]["l2_error"],"Missing Yap not detected")
    run("cls_shifted_distance",dt=.0125,steps=8,modes=5,points=12,modify=lambda r:parameter(r,"RANSWallOrigin",.5))
    require(records["cls_shifted_distance"]["field_audit"]["fields"]["eps_tilde"]["l2_error"]>10*temporal[1]["field_audit"]["fields"]["eps_tilde"]["l2_error"],"Geometric Yap distance not exercised")
    for degree in (4,6,8):
        run(f"cls_abl_p{degree}",kind="bulk_abl",dt=.025,steps=4,modes=degree+1,points=degree+6)
    for f in ("u","v","k","eps_tilde","p"):
        require(records["cls_abl_p8"]["field_audit"]["fields"][f]["l2_error"]<records["cls_abl_p4"]["field_audit"]["fields"][f]["l2_error"],"CLS column refinement failed for "+f)
    run("awf_initial",kind="awf",dt=.025,steps=1,order=1,modes=7,points=12)
    run("awf",kind="awf",dt=.025,steps=8,modes=7,points=12)
    run("awf_mpi",kind="awf",ranks=4,dt=.025,steps=8,modes=7,points=12)
    run("awf_mpi_n8",kind="awf",ranks=4,dt=.025,steps=4,modes=7,points=12,
        modify=lambda r:parameter(r,"RANSStripPointN",8))
    require(records["awf_mpi_n8"]["wall_stages"][-1]["wall_ranks"]<4 and
            records["awf_mpi_n8"]["wall_stages"][-1]["remote_N_without_wall_samples"]>0,
            "MPI test needs N samples owned by ranks without wall faces")
    for name in ("awf_initial","awf"):
        for row in records[name]["wall_stages"]:
            ref=awf_reference(row["mean_kP"])
            tau=row["mean_uH"]/ref["j0"]-ref["j1"]/ref["j0"]*row["mean_forcing_x"]
            require(abs(tau-row["traction_x"]/row["area"])<2e-8,"Coupled original traction disagrees with independent quadrature")
            forcing=row["mean_forcing_x"]
            production=(tau*tau*ref["production_weight"]+2*tau*forcing*ref["production_weight1"]+
                        forcing*forcing*ref["production_weight2"])/2
            require(abs(production-row["mean_production"])<2e-8,"Original strip production disagrees with independent quadrature")
            require(abs(ref["epsilon_mean"]-row["mean_dissipation"])<2e-8,"Original strip mean dissipation not recovered")
            require(abs(ref["epsilon_p"]-row["mean_epsP"])<2e-8 and abs(ref["epsilon_h"]-row["mean_epsH"])<2e-8,"Original EP/EH profile mismatch")
            require(abs(2*row["mean_uH"]-row["mean_uN"]-row["mean_uP"])<2e-8,"Log matching elimination failed")
    first=records["awf_initial"]["wall_stages"][0]
    reference=awf_reference(3.190675081624288*.4**2)
    reference["mean_production"]=.4**4*reference["production_weight"]/2
    require(abs(reference["mean_production"]-.461092609)<1e-8 and abs(reference["epsilon_mean"]-.774943967)<1e-8,
            "Independent zero-forcing reference disagrees with the note")
    require(first["mean_dissipation"]>first["mean_production"] and records["awf"]["wall_stages"][-1]["mean_kP"]<first["mean_kP"],"Original AWF baseline was incorrectly balanced/reset")
    run("awf_3d",kind="awf",ranks=4,dimension=3,dt=.001,steps=2,modes=5,points=9,span_angle=.4,
        modify=lambda r:parameter(r,"RANSPressureTolerance",1e-7))
    for row in records["awf_3d"]["wall_stages"]:
        require(abs(row["traction_y"]/row["traction_x"]-math.tan(.4))<1e-5,
                "The two tangential AWF tractions lost their prescribed wind direction")
    run("awf_nonuniform",kind="awf",ranks=4,dt=.01,steps=4,modes=7,points=12,perturbation=.002)
    require(max(r["normal_velocity_max"] for r in records["awf_nonuniform"]["wall_stages"])>1e-6,"Nonuniform wall test did not exercise strip normal velocity")
    for name,h,nu in (("awf_smooth",0,1.5e-5),("awf_transitional",7e-5,1.5e-5),("awf_molecular",0,5.)):
        run(name,kind="awf",dt=.001,steps=2,modes=5,points=10,roughness=h,physical_nu=nu)
    for row in records["awf_molecular"]["wall_stages"]:
        require(abs(row["mean_production"])<1e-14 and
                abs(row["mean_epsP"]-row["mean_epsH"])<1e-12 and
                abs(row["mean_epsP"]-row["mean_dissipation"])<1e-12,
                "The entirely molecular auxiliary wall branch was not exercised")

    def body_force(r,force):
        names=[v.text for v in r.findall("CONDITIONS/VARIABLES/V")][:-1]
        f=ET.SubElement(r.find("CONDITIONS"),"FUNCTION",NAME="WallBody")
        for i,name in enumerate(names):ET.SubElement(f,"E",VAR=name,VALUE=str(force if i==0 else 0))
        forcing=ET.SubElement(ET.SubElement(r,"FORCING"),"FORCE",TYPE="Body")
        ET.SubElement(forcing,"BODYFORCE").text="WallBody"
    run("awf_reverse_traction",kind="awf",dt=.001,steps=2,modes=7,points=12,modify=lambda r:body_force(r,-1.))
    require(records["awf_reverse_traction"]["wall_stages"][0]["traction_x"]<0,"Signed reversed original traction was lost")
    def invalid_k(r):
        r.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']/E[@VAR='k']").attrib["VALUE"]="1+2*(y-2)"
    run("awf_invalid_kP",kind="awf",ranks=4,dt=.01,steps=1,modes=5,points=10,modify=invalid_k,reject="nonpositive eliminated KP")
    run("awf_invalid_N",kind="awf",ranks=4,dt=.01,steps=1,modes=5,points=10,modify=lambda r:parameter(r,"RANSStripPointN",40),reject="every physical N point")
    run("cls_invalid_distance",ranks=4,dt=.025,steps=1,modes=5,points=12,modify=lambda r:parameter(r,"RANSWallOrigin",1),reject="positive geometric distance")
    comparisons={}
    for a,b in (("cls_dt0.0125","cls_mpi"),("cls_dt0.0125","cls_quadrature"),("awf","awf_mpi")):
        raw=execute([audit,output/a/(a+".xml"),output/a/(a+".fld"),output/b/(b+".fld")],output/b,"comparison")
        comparisons[a+"__"+b]=json.loads(raw[raw.index('{"time"'):])
        require(max(v["difference_linf"] for v in comparisons[a+"__"+b]["fields"].values())<2e-7,"Full-field serial/MPI or quadrature mismatch")
    source_files=["solvers/IncNavierStokesSolver/EquationSystems/"+n for n in
        ("RANSVelocityCorrectionScheme.cpp","RANSVelocityCorrectionScheme.h","RANSAnalyticalWallStrip.cpp","RANSAnalyticalWallStrip.h")]
    source_files += ["solvers/IncNavierStokesSolver/Tests/RANS/"+n for n in ("generate_m4.py","validate_m4.py","CMakeLists.txt")]
    source_files += ["library/UnitTests/SolverUtils/RANSFieldAudit.cpp","solvers/IncNavierStokesSolver/CMakeLists.txt"]
    source_files += ["solvers/IncNavierStokesSolver/Tests/RANS/generate_m3.py","solvers/IncNavierStokesSolver/Tests/TGV3d_2_2.xml"]
    source_files += ["solvers/IncNavierStokesSolver/Tests/RANS/"+n for n in ("CLS_MMS.xml","AWFOriginal.xml","plot_m4.py")]
    source_files += ["library/SolverUtils/Turbulence/"+n for n in ("RANSModels.cpp","RANSModels.h","AnalyticalWallFunction.cpp","AnalyticalWallFunction.h","RANSWallStrip.h")]
    report={"schema":1,"milestone":"neutral_abl_rans_awf_m4","status":"passed","utc":datetime.now(timezone.utc).isoformat(),
        "source_parent":subprocess.check_output(["git","-C",root,"rev-parse","HEAD"],text=True).strip(),
        "source_sha256":{p:digest(root/p) for p in source_files},"binary_sha256":binary_hashes,"python":platform.python_version(),
        "maximum_project_workers":16,"cpu_affinity":sorted(os.sched_getaffinity(0)),"maximum_mpi_ranks":4,"concurrent_cases":1,"threads":THREADS,
        "temporal_rates":rates,"original_reference_oracle":reference,"full_field_comparisons":comparisons,"cases":records,"commands":commands}
    (output/"validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":"passed","temporal_rates":rates,"case_count":len(records)},indent=2),flush=True)


if __name__=="__main__":main()
