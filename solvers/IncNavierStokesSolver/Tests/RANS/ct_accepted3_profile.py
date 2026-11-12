#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate timing-only runtime, equivalence, and accepted-segment profile.

Primary physical cohorts stay pinned to CT2 runtime. Detailed row/law scopes
aggregate by accepted step, avoiding per-point output. No new physics cache.
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
import numpy as np
from rans_request_manifest import RequestManifest, runtime_identity, sha
from ct_accepted_open import run, require
from ct_accepted3_run import make_case
from compare_cls_epawf_pimple import same_history


def build(a):
    os.environ['CMAKE_BUILD_PARALLEL_LEVEL']='2'
    runtime=a.output/'runtime'
    m=RequestManifest(a.output,dict(schema='CT3DetailedTimingBuildV1',source_sha256=sha(Path(__file__)),
        source_revision=subprocess.check_output(['git','-c','index.threads=1','-c','core.preloadIndex=false','rev-parse','HEAD'],cwd=a.root/'source',text=True).strip(),
        compiler_jobs=2,governing_runtime=runtime_identity(a.ct2/'runtime/IncNavierStokesSolver')))
    source=a.root/'source/solvers/IncNavierStokesSolver/Tests/RANS'
    commands=[['cmake','-S',str(source),'-B',str(a.output/'build'),'-DCMAKE_BUILD_TYPE=Release','-DRANS_BUILD_SOLVER=ON',
        '-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(a.ct2/'dependency/libSpatialDomains.so.5.9.0'),
        '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'],
        ['cmake','--build',str(a.output/'build'),'--parallel','2','--target','IncNavierStokesSolver','RANSWallTimingChecks','RANSNativeDiagnosticChecks']]
    result=dict(status='fail',executions=[])
    for i,cmd in enumerate(commands):
        e=run(cmd,a.output,'build-'+str(i)+'.log');result['executions'].append(e)
        write_json(a.output/'progress.json',result);require(e['exit_code']==0,'detailed_profile_build_failed')
    runtime.mkdir();shutil.copytree(a.ct2/'runtime/lib',runtime/'lib',symlinks=True)
    for name in ('IncNavierStokesSolver','RANSWallTimingChecks','RANSNativeDiagnosticChecks'):shutil.copy2(a.output/'build'/name,runtime/name)
    # The observation adapter is intentionally the exact unchanged CT2 binary.
    for name in ('RANSPreservationProbe','RANSNewModelsProbe'):shutil.copy2(a.ct2/'runtime'/name,runtime/name)
    os.environ['LD_LIBRARY_PATH']=str(runtime/'lib')
    os.environ.update(NEKTAR_RANS_TIMING_FILE=str(a.output/'unit-timing.jsonl'),NEKTAR_RANS_DETAILED_TIMING='1')
    e=run([runtime/'RANSWallTimingChecks'],a.output,'timing-check.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE');os.environ.pop('NEKTAR_RANS_DETAILED_TIMING')
    require(e['exit_code']==0,'buffered_timing_native_test_failed')
    rows=[json.loads(x) for x in (a.output/'unit-timing.jsonl').read_text().splitlines()]
    agg={q['name']:q for q in rows if q.get('aggregated')}
    require(agg['buffered']['calls']==agg['buffered_child']['calls']==5,'buffered_timing_call_count')
    require(all(0<=q['exclusive_process_cpu_seconds']<=q['process_cpu_seconds']+1e-12 and 0<=q['exclusive_wall_seconds']<=q['wall_seconds']+1e-12 for q in rows),'nested_timing_intervals_invalid')
    result.update(status='pass',runtime={name:runtime_identity(runtime/name) for name in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNativeDiagnosticChecks')},
        unit_test=e,unit_aggregate=agg,libraries={p.name:sha(p) for p in (runtime/'lib').iterdir() if p.is_file()},
        scope='timing-only binary; detailed timings optional; no governing source change or caching; observation adapter remains CT2 byte identity')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return result


def profile(a):
    c=json.loads((a.parent/'construction.json').read_text());old=json.loads((a.reference/'case.evolution.json').read_text())
    parent=json.loads((a.parent/'case.evolution.json').read_text());dt=old['dt']
    require(abs(parent['time']-.008)<1e-13 and abs(old['time']-.02)<1e-13,'profile_requires_admitted_accepted_segment')
    m=RequestManifest(a.output,dict(schema='CT3AcceptedDetailedProfileV1',parent_sha256=sha(a.parent/'case.evolution.json'),reference_sha256=sha(a.reference/'case.evolution.json'),
        runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),start=.008,target=.02,dt=dt,source_sha256=sha(Path(__file__))))
    folder=a.output/'run';make_case(folder,a.reference/'case.xml',c,c['degree'],dt,round((.02-.008)/dt),Path(parent['checkpoint']))
    os.environ.update(LD_LIBRARY_PATH=str(a.runtime/'lib'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),NEKTAR_RANS_DETAILED_TIMING='1')
    e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE');os.environ.pop('NEKTAR_RANS_DETAILED_TIMING')
    native=json.loads((folder/'case.evolution.json').read_text())
    require(e['exit_code']==0 and not native['failure'] and same_history(old,native) and old['accepted_states']==native['accepted_states'],'detailed_timing_changed_accepted_states_or_history')
    phases={};startup={};per_step={}
    for line in (folder/'timing.jsonl').read_text().splitlines():
        q=json.loads(line);target=phases if q['step']>0 and q['time']>.008 else startup
        p=target.setdefault(q['name'],dict(exclusive_process_cpu_seconds=0.,exclusive_wall_seconds=0.,calls=0))
        p['exclusive_process_cpu_seconds']+=q['exclusive_process_cpu_seconds'];p['exclusive_wall_seconds']+=q['exclusive_wall_seconds'];p['calls']+=q.get('calls',1)
        if target is phases:
            d=per_step.setdefault(str(q['step']),dict(time=q['time'],exclusive_cpu=0.));d['exclusive_cpu']+=q['exclusive_process_cpu_seconds']
    result=dict(status='pass',bitwise_accepted_state_history_equivalence=True,execution=e,stored_slots=len(native['state']),accepted_steps_in_segment=round((.02-.008)/dt),
        start=.008,target=.02,dt=dt,accepted_exclusive_scopes=phases,startup_and_output_scopes=startup,per_step=per_step,
        native_report_sha256=sha(folder/'case.evolution.json'),timing_sha256=sha(folder/'timing.jsonl'),
        scope='sum EXCLUSIVE intervals only; nested inclusive totals not added; detailed timing overhead and output remain visible; profile does not change primary runtime attribution')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=['build','profile'],required=True)
    for n in ('root','ct2','output'):p.add_argument('--'+n,type=Path,required=True)
    for n in ('parent','reference','runtime'):p.add_argument('--'+n,type=Path)
    a=p.parse_args();r=build(a) if a.mode=='build' else profile(a);print(r['status']);return 0
if __name__=='__main__':raise SystemExit(main())
