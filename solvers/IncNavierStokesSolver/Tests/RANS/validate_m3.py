#!/usr/bin/env python3
"""Linux standard-RANS transport acceptance: accuracy, pressure, walls and MPI."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from threading import Lock
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS={k:"1" for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS")}
os.environ.update(THREADS)
if sys.platform.startswith("linux"):
    os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[-16:])
from generate_m3 import generate, generate_3d, manufactured


def require(condition,message):
    if not condition:raise RuntimeError(message)


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def parameter(tree,key,value):
    node=tree.find("CONDITIONS/PARAMETERS")
    for item in node:
        if item.text.split("=")[0].strip()==key:
            item.text=f"{key} = {value}";return
    ET.SubElement(node,"P").text=f"{key} = {value}"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--skip-long-abl",action="store_true",help="Pilot only; cannot close Milestone 3")
    args=parser.parse_args()
    require(sys.platform.startswith("linux"),"Authoritative tests must run on Linux")
    root=Path(__file__).resolve().parents[4]
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    solver=args.build.resolve()/"IncNavierStokesSolver";audit=args.build.resolve()/"RANSFieldAudit"
    require(solver.is_file() and audit.is_file(),"Build the isolated solver and field audit first")
    binary_hashes={p.name:digest(p) for p in (solver,audit)}
    commands=[];records={};record_lock=Lock()
    def execute(command,cwd,label,timeout=240,expected_failure=None):
        print("Running "+label,flush=True)
        commands.append({"command":list(map(str,command)),"cwd":str(cwd),"timeout":timeout})
        result=subprocess.run(list(map(str,command)),cwd=cwd,env=os.environ.copy(),
                              capture_output=True,text=True,timeout=timeout)
        text=result.stdout+result.stderr
        (cwd/(label+".log")).write_text(text)
        if expected_failure:
            require(result.returncode!=0 and expected_failure in text,f"{label}: missing expected rejection")
        else:require(result.returncode==0,f"{label} failed: {cwd/(label+'.log')}")
        return text
    def run_case(name,kind="mms",ranks=1,dimension=2,modify=None,reject=None,**options):
        directory=output/name;directory.mkdir(exist_ok=True)
        session=directory/(name+".xml")
        if dimension==3:generate_3d(session,kind=kind,**options)
        else:generate(session,kind=kind,**options)
        tree=ET.parse(session)
        if ranks>1:
            for info in tree.findall("CONDITIONS/SOLVERINFO/I"):
                if info.attrib["PROPERTY"]=="GlobalSysSoln":info.attrib["VALUE"]="IterativeStaticCond"
        if modify:modify(tree.getroot())
        tree.write(session,encoding="utf-8",xml_declaration=True)
        fingerprint={"fixture_sha256":digest(session),"binaries":binary_hashes,"ranks":ranks,"rejection":reject}
        saved=directory/"case.json"
        if saved.exists() and json.loads(saved.read_text()).get("fingerprint")==fingerprint:
            record=json.loads(saved.read_text());print("Reusing exact fixture/binary result "+name,flush=True)
        else:
            command=(["mpirun","--bind-to","none","--map-by","slot","-np",str(ranks)] if ranks>1 else [])+[str(solver),session.name]
            log=execute(command,directory,name,timeout=900 if kind.startswith("abl") else 240,expected_failure=reject)
            record={"fingerprint":fingerprint,"options":options,"dimension":dimension,"kind":kind,"ranks":ranks,
                    "solver_command":command}
            if reject:record["expected_rejection"]=reject
            else:
                require("Total Computation Time" in log and "L 2 error" in log,"Solver did not finish")
                raw=execute([audit,session.name,name+".fld",name+".fld"],directory,"audit")
                record["field_audit"]=json.loads(raw[raw.index('{"time"'):])
                rows=list(csv.DictReader((directory/(name+".rans.csv")).open()))
                require(len(rows)==options.get("steps",4),"Missing implicit-stage diagnostics")
                numbers=[{k:float(v) for k,v in row.items()} for row in rows]
                require(all(math.isfinite(v) for row in numbers for v in row.values()),"Nonfinite stage diagnostic")
                require(all(row["min_k"]>0 and row["min_eps"]>0 and row["guard_activations"]==0 for row in numbers),"Unexpected state guard or nonpositive state")
                record["diagnostics"]={"stages":len(rows),"maximum_iterations":max(r["iterations"] for r in numbers),
                    "maximum_change":max(r["change"] for r in numbers),"maximum_pressure_correction":max(r["pressure_correction"] for r in numbers),
                    "maximum_strong_divergence_l2":max(r["divergence_l2"] for r in numbers),
                    "maximum_weak_divergence_l2":max(r["weak_divergence_l2"] for r in numbers),
                    "maximum_live_global_systems":max(r["live_global_systems"] for r in numbers),
                    "first_live_global_systems":numbers[0]["live_global_systems"],"last_live_global_systems":numbers[-1]["live_global_systems"],
                    "first_min_nut":numbers[0]["min_nut"],"last_min_nut":numbers[-1]["min_nut"]}
                require(record["diagnostics"]["last_live_global_systems"]<=record["diagnostics"]["first_live_global_systems"]+2,"Growing matrix-system history")
                record["artifacts_sha256"]={p.name:digest(p) for p in directory.iterdir() if p.is_file() and p.suffix in (".log",".csv",".xml",".his")}
            saved.write_text(json.dumps(record,indent=2)+"\n")
        with record_lock:
            records[name]=record
            (output/"progress.json").write_text(json.dumps(records,indent=2)+"\n")
        if "field_audit" in record:
            print(name+" errors "+str({k:v["l2_error"] for k,v in record["field_audit"]["fields"].items()}),flush=True)
        return record

    temporal={}
    # Three independent one-rank jobs, each with one numerical thread; this
    # pool closes before any MPI case starts. Prewarm symbolic forcing once.
    manufactured()
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs=[pool.submit(run_case,f"mms_o{order}_dt{dt:g}",dt=dt,
                          steps=round(.2/dt),order=order,modes=7)
              for order in (1,2) for dt in (.05,.025,.0125)]
        for job in jobs:job.result()
    for order in (1,2):
        group=[records[f"mms_o{order}_dt{dt:g}"] for dt in (.05,.025,.0125)]
        rates={}
        for field in ("u","v","k","eps","p"):
            errors=[r["field_audit"]["fields"][field]["l2_error"] for r in group]
            rates[field]=[math.log(a/b,2) for a,b in zip(errors,errors[1:])]
            require(min(rates[field])>(.75 if order==1 else 1.65),f"Incorrect order {order} for {field}: {rates[field]}")
        temporal[str(order)]=rates
    run_case("mms_highre",dt=.025,steps=8,order=2,modes=7,diffusion="HighRe")
    run_case("mms_quadrature",dt=.025,steps=8,order=2,modes=7,points=14)
    def offset_gauge(tree):
        p=tree.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']/E[@VAR='p']")
        p.attrib["VALUE"]="("+p.attrib["VALUE"]+")+1000"
    run_case("mms_gauge_offset",dt=.025,steps=8,order=2,modes=7,modify=offset_gauge)
    for field in ("u","v","k","eps","p"):
        a=records["mms_gauge_offset"]["field_audit"]["fields"][field]["l2_error"]
        b=records["mms_o2_dt0.025"]["field_audit"]["fields"][field]["l2_error"]
        require(abs(a-b)<1e-8,"Initial pressure gauge changed the solution")
    inconsistent=run_case("mms_missing_transpose",dt=.0125,steps=16,order=2,modes=7,missing_transpose=True)
    # Much of a missing transpose contribution can be absorbed into pressure.
    # Pressure consistency is therefore the decisive observable in this case.
    require(inconsistent["field_audit"]["fields"]["p"]["l2_error"]>10*records["mms_o2_dt0.0125"]["field_audit"]["fields"]["p"]["l2_error"],"Missing transpose stress was not detected in pressure")
    stiff=run_case("stiff_euler",kind="stiff",dt=.1,steps=2,order=1,modes=5)
    k,epsilon=1.,50.
    rows=list(csv.DictReader((output/"stiff_euler/stiff_euler.rans.csv").open()))
    for row in rows:
        old_ratio=epsilon/k;dt=.1
        ratio=(dt*old_ratio-1+math.sqrt((1-dt*old_ratio)**2+4*dt*1.92*old_ratio))/(2*dt*1.92)
        k/=1+dt*ratio;epsilon/=1+dt*1.92*ratio
        require(abs(float(row["min_k"])-k)<1e-8 and abs(float(row["min_eps"])-epsilon)<1e-7,"Implicit homogeneous loss differs from independent backward-Euler root")
    run_case("stiff_bdf2_guard",kind="stiff",dt=.1,steps=2,order=2,modes=5,reject="positivity guard")
    decay=[]
    for dt in (.025,.0125,.00625):
        decay.append(run_case(f"decay_dt{dt:g}",kind="decay",dt=dt,steps=round(.2/dt),order=2,modes=4))
    decay_rates={f:[math.log(a["field_audit"]["fields"][f]["l2_error"]/b["field_audit"]["fields"][f]["l2_error"],2) for a,b in zip(decay,decay[1:])] for f in ("k","eps")}
    require(all(min(r)>1.7 for r in decay_rates.values()),"Decay does not converge at second order")
    require(decay[-1]["diagnostics"]["first_min_nut"]!=decay[-1]["diagnostics"]["last_min_nut"],"Viscosity was frozen")
    run_case("decay_3d",kind="decay",dimension=3,dt=.025,steps=8,order=2,modes=4)
    run_case("decay_3d_mpi",kind="decay",dimension=3,ranks=4,dt=.025,steps=8,order=2,modes=4)
    three=run_case("mms_3d",dimension=3,ranks=4,dt=.01,steps=2,order=2,modes=5,points=11)
    require(all(f["l2_error"]<3e-5 for f in three["field_audit"]["fields"].values()),"Three-component manufactured solution failed")
    run_case("mms_mpi",ranks=4,dt=.025,steps=8,order=2,modes=7)
    comparisons={}
    for a,b in [("mms_o2_dt0.025","mms_mpi"),("mms_o2_dt0.025","mms_quadrature"),("decay_3d","decay_3d_mpi")]:
        directory=output/b
        raw=execute([audit,output/a/(a+".xml"),output/a/(a+".fld"),directory/(b+".fld")],directory,"serial-mpi-audit")
        data=json.loads(raw[raw.index('{"time"'):]);comparisons[a+"__"+b]=data
        require(all(f["difference_linf"]<2e-8 for f in data["fields"].values()),"Serial/MPI full-field disagreement")
    run_case("mpi_invalid_state",ranks=4,dt=.025,steps=2,order=2,modes=5,
             modify=lambda tree:parameter(tree,"RANSMinimumK",1.025),reject="positivity guard")
    for modes in (5,7,9):
        run_case(f"abl_p{modes-1}",kind="abl",dt=.05,steps=20,order=2,modes=modes,diffusion="HighRe",tolerance=1e-9)
    for f in ("u","k","eps"):
        require(records["abl_p8"]["field_audit"]["fields"][f]["l2_error"]<records["abl_p4"]["field_audit"]["fields"][f]["l2_error"],"ABL spatial refinement failed")
    run_case("abl_wrong_top",kind="abl",dt=.05,steps=20,order=2,modes=9,diffusion="HighRe",tolerance=1e-9,top_flux=False)
    require(records["abl_wrong_top"]["field_audit"]["fields"]["eps"]["l2_error"]>5*records["abl_p8"]["field_audit"]["fields"]["eps"]["l2_error"],"Incompatible top flux was not detected")
    run_case("abl_shifted",kind="abl_shift",dt=.05,steps=20,order=2,modes=9,diffusion="HighRe",tolerance=1e-9)
    run_case("abl_molecular",kind="abl",dt=.05,steps=20,order=2,modes=9,diffusion="Molecular",tolerance=1e-9,
             modify=lambda tree:parameter(tree,"Kinvis",.01))
    if not args.skip_long_abl:
        run_case("abl_flowthrough",kind="abl",ranks=4,dt=.05,steps=600,order=2,modes=9,diffusion="HighRe",tolerance=1e-9)
        fields=records["abl_flowthrough"]["field_audit"]["fields"]
        kappa=math.sqrt(.48*1.3*.3)
        friction_velocity=kappa*10/math.log(500)
        epsilon_top=friction_velocity**3/(kappa*30)
        require(fields["u"]["linf_error"]/10<.001 and
                fields["k"]["linf_error"]/fields["k"]["max"]<.001 and
                fields["eps"]["linf_error"]/epsilon_top<.001,"Standard ABL drift exceeds the milestone target")
    legacy=output/"legacy";legacy.mkdir(exist_ok=True)
    fixture=root/"solvers/IncNavierStokesSolver/Tests/KovaFlow_m8_short_HOBC.xml"
    shutil.copy2(fixture,legacy/fixture.name)
    legacy_log=execute([solver,fixture.name],legacy,"legacy")
    for name,target in {"u":2.30665e-8,"v":5.67848e-9,"p":1.07686e-8}.items():
        value=float(re.search(r"L 2 error \(variable "+name+r"\)\s*:\s*(\S+)",legacy_log)[1])
        require(abs(value-target)<1e-11,"Legacy VCS regression changed")
    source_files=[str(p.relative_to(root)) for p in
        [root/"solvers/IncNavierStokesSolver/EquationSystems"/n for n in ("RANSVelocityCorrectionScheme.h","RANSVelocityCorrectionScheme.cpp","IncNavierStokes.h")]+
        [root/"library/UnitTests/SolverUtils/RANSFieldAudit.cpp",Path(__file__).resolve(),Path(__file__).with_name("generate_m3.py")]]
    source_files += ["library/SolverUtils/Turbulence/RANSModels.cpp","library/SolverUtils/Turbulence/RANSModels.h","solvers/IncNavierStokesSolver/CMakeLists.txt","solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt","solvers/IncNavierStokesSolver/Tests/RANS/StandardRANS_MMS.xml","solvers/IncNavierStokesSolver/Tests/TGV3d_2_2.xml"]
    report={"schema":1,"milestone":"neutral_abl_rans_awf_m3","status":"pilot_passed" if args.skip_long_abl else "passed",
        "utc":datetime.now(timezone.utc).isoformat(),"source_parent":subprocess.check_output(["git","-C",root,"rev-parse","HEAD"],text=True).strip(),
        "source_sha256":{p:digest(root/p) for p in source_files},"binary_sha256":binary_hashes,
        "cpu_affinity":sorted(os.sched_getaffinity(0)),"maximum_project_workers":16,"maximum_mpi_ranks":4,
        "maximum_concurrent_serial_cases":3,"mpi_overlaps_serial_pool":False,"threads":THREADS,
        "time_refinement_rates":temporal,"decay_rates":decay_rates,"full_field_comparisons":comparisons,
        "cases":records,"commands":commands,"legacy_vcs_regression":"passed"}
    (output/"validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"temporal_rates":temporal,"decay_rates":decay_rates},indent=2),flush=True)


if __name__=="__main__":main()
