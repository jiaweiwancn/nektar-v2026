#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Diagnose and qualify exact immutable trace addresses in independent audits."""
import argparse,cProfile,hashlib,json,os,pstats,resource,time
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted15_python import PreparedRows
from ct_accepted15_audits import bundle
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from rans_trace_index import TraceIndex,positive_h_trace_index
from rans_request_manifest import RequestManifest,sha,canonical


def controls():
    H=.1875
    packet=dict(observation=dict(grid=[[-0.,0.,H],[0.,0.,H],[1.,0.,H]],grid_sides=[1,1,1]),test_basis=[{}, {}, {}])
    oracle=[positive_h_trace_index(packet,x,H) for x in (-0.,0.,1.)]
    with TraceIndex(packet,H) as index:
        require([positive_h_trace_index(packet,x,H) for x in (-0.,0.,1.)]==oracle==[0,0,2],'trace first-match/signed-zero membership')
        foreign=dict(packet)
        require(positive_h_trace_index(foreign,1.,H)==2 and index.queries==3,'foreign packet borrowed index')
        for x,height in [(2.,H),(float('nan'),H),(1.,np.nextafter(H,1.))]:
            try:positive_h_trace_index(packet,x,height)
            except StopIteration:pass
            else:raise ValueError('missing trace accepted')
        packet['observation']['grid'][2][0]=2.
        try:positive_h_trace_index(packet,1.,H)
        except ValueError:pass
        else:raise ValueError('stale trace metadata accepted')
        packet['observation']['grid'][2][0]=1.
    try:
        with TraceIndex(packet,H):packet['observation']['grid_sides'][2]=-1
    except ValueError:pass
    else:raise ValueError('changed trace scope geometry accepted')
    packet['observation']['grid_sides'][2]=1
    require(positive_h_trace_index(packet,1.,H)==2,'trace scope was not restored after failure')
    return dict(status='pass',first_match_and_signed_zero=True,foreign_packet_fallback=True,
                missing_nonfinite_height_mismatch_rejected=True,stale_metadata_rejected=True,
                changed_scope_geometry_rejected=True,scope_restored_after_failure=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct14','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--phase',choices=['diagnostic','comparison'],required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));n.reference=lru_cache(None)(n.reference)
    scripts=Path(__file__).parent;initial=a.ct14/'static-v1/target/frozen.json';paths=[a.ct14/'pair-v2'/dt/'case.evolution.json' for dt in ('0.00025','0.000125')]
    names=['rans_trace_index.py','ct_accepted4_run.py','ct_accepted5_traction.py','ct_accepted10_compare.py','ct_accepted16_trace.py']
    m=RequestManifest(a.output,dict(schema='CT16ExactTraceIndexQualificationV1',phase=a.phase,policy=TraceIndex.policy,
        source_sha256={name:sha(scripts/name) for name in names},initial_sha256=sha(initial),actual_history_sha256=[sha(f) for f in paths],
        current_state_cache=False,immutable_address_index_only=True,physical_advance=False,
        oracle='same PreparedRows arithmetic,original first-match linear trace search,unchanged complete per-state independent bundle',
        bounds=dict(maximum_grid_rows=TraceIndex.maximum_grid_rows,maximum_entries=TraceIndex.maximum_entries,maximum_index_bytes=TraceIndex.maximum_index_bytes)))
    ini=n.read(initial);c=n.read(a.ct14/'static-v1/results.json')['case'];weights=n.inventory(c)[2];count=len(weights)
    states=[ini['observation']['state']]+[n.read(f)['state'] for f in paths]
    # This immutable analytic-reference preparation is common to both paths,
    # charged once and separately measured. Actual state physics remains fresh.
    cpu=time.process_time();wall=time.monotonic()
    for q,region in zip(ini['observation']['grid'][:count],ini['observation']['grid_regions'][:count],strict=True):n.reference(n.col(c),q[2],bool(region))
    out=dict(schema='CT16ExactTraceIndexQualificationV1',status='incomplete',phase=a.phase,policy=TraceIndex.policy,
        common_immutable_reference_startup=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall),physical_advance=False)
    def evaluate(state):return bundle(ini,state,c,weights,count)
    if a.phase=='diagnostic':
        profile=cProfile.Profile();cpu=time.process_time();wall=time.monotonic()
        with PreparedRows():q=profile.runcall(evaluate,states[-1])
        require(q['metrics']['status']=='pass','original independent endpoint scientific gate')
        profile.dump_stats(str(a.output/'bundle.pstats'));stats=pstats.Stats(profile)
        ranked=sorted(stats.stats.items(),key=lambda v:-v[1][3])[:40]
        out['diagnostic']=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall,bundle_sha256=hashlib.sha256(canonical(q)).hexdigest(),
            cumulative_functions=[dict(file=k[0],line=k[1],function=k[2],primitive_calls=v[0],calls=v[1],self_CPU=v[2],cumulative_CPU=v[3]) for k,v in ranked],
            interpretation='cProfile diagnosis only;instrumented cumulative intervals overlap and cannot be added or supply acceptance speedups')
        write_json(a.output/'original-bundle.json',q)
    else:
        out['controls']=controls();out['states']=[]
        for i,state in enumerate(states):
            cpu=time.process_time();wall=time.monotonic()
            with PreparedRows():before=evaluate(state)
            old=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall)
            cpu=time.process_time();wall=time.monotonic()
            with PreparedRows(),TraceIndex(ini,c['H']) as index:after=evaluate(state)
            new=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall,index_entries=len(index.indices),index_bytes=index.bytes,index_queries=index.queries,geometry_sha256=index.geometry_sha256)
            require(before==after,'immutable trace addresses changed an independent physical observation or decision')
            require(before['metrics']['status']=='pass' and index.queries>0,'unchanged physical gates/index coverage')
            write_json(a.output/('bundle-'+str(i)+'.json'),after)
            out['states'].append(dict(index=i,status='pass',all_observations_and_decisions_exact_equal=True,
                bundle_sha256=sha(a.output/('bundle-'+str(i)+'.json')),original=old,trace_index=new))
            write_json(a.output/'progress.json',out);print('CT16 exact trace index',i,old['CPU'],new['CPU'],flush=True)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json'])
    return 0


if __name__=='__main__':raise SystemExit(main())
