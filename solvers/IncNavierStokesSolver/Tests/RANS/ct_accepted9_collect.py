#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT9 evidence; recompute physical primitives before aggregation."""
import argparse,copy,os,json,hashlib,math
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_gate import pair,backend_controls,backend_trajectory,capacity as capacity_gate
from ct_accepted9_pressure_run import sealed
from ct_accepted_feasibility import write_json,require,parameters
from rans_request_manifest import RequestManifest,sha,canonical

LIBRARIES={}
def compact(value):
    if isinstance(value,list):return [compact(x) for x in value]
    if not isinstance(value,dict):return value
    result={}
    for key,x in value.items():
        if key=='actual_mapped_files':
            digest=hashlib.sha256(canonical(x)).hexdigest();LIBRARIES[digest]=x;result['actual_mapped_files_catalog_sha256']=digest
        elif key=='interface_rows':
            result['interface_row_count']=len(x);result['interface_rows_sha256']=hashlib.sha256(canonical(x)).hexdigest()
            result['largest_interface_row_witness']=max(x,key=lambda r:sum(v*v for v in r['traction_jump'])) if x else None
        else:result[key]=compact(x)
    return result

def timing(path):
    totals={};steps={}
    with path.open() as stream:
        for line in stream:
            q=json.loads(line);entry=totals.setdefault(q['name'],dict(calls=0,exclusive_cpu=0.,exclusive_wall=0.));entry['calls']+=q.get('calls',1)
            entry['exclusive_cpu']+=q['exclusive_process_cpu_seconds'];entry['exclusive_wall']+=q['exclusive_wall_seconds']
            s=steps.setdefault(str(q['step']),dict(exclusive_cpu=0.,exclusive_wall=0.));s['exclusive_cpu']+=q['exclusive_process_cpu_seconds'];s['exclusive_wall']+=q['exclusive_wall_seconds']
    return dict(path=str(path),sha256=sha(path),names=totals,steps=steps,scope='exclusive native process intervals only; nested inclusive totals are never added; detailed aggregation includes clock overhead when explicitly enabled')

def pressure_vectors(root,ct8):
    pre=sealed(ct8/'preflight-v6');old=sealed(ct8/'strip-T008-v1');fine=sealed(root/'pressure-finer-v1');count=pre['global_count'];weights=np.asarray(pre['weights']);ini=n.read(Path(pre['initial']['0.001']['directory'])/'initial.json');regions=np.asarray(ini['observation']['grid_regions'][:count]);scale=pre['case']['ustar']**2
    reports=[n.read(Path(old['members'][k]['directory'])/'case.evolution.json') for k in ('0.002','0.001')]+[n.read(Path(fine['members']['0.0005']['directory'])/'case.evolution.json')]
    histories=[{round(t,12):s for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)} for r in reports];out=[]
    for t in sorted(set.intersection(*(set(h) for h in histories))):
        pi=[n.observations(ini,h[t],count)[:,4] for h in histories];a,b=pi[0]-pi[1],pi[1]-pi[2];rows={}
        for label,mask in [('whole',np.ones(count,bool)),('inner',regions==0),('outer',regions==1)]:
            w=weights[mask];x,y=a[mask],b[mask];aa=float(w@(x*x));bb=float(w@(y*y));ab=float(w@(x*y));rr=float(w@((x-2*y)**2));rows[label]=dict(a_true_L2=aa**.5,b_true_L2=bb**.5,alignment=ab/(aa*bb)**.5 if aa*bb else None,a_minus_2b_true_L2=rr**.5,a_signed_mean=float(w@x/w.sum()),b_signed_mean=float(w@y/w.sum()),a_scaled_RMS=float(np.sqrt(aa/w.sum())/scale),b_scaled_RMS=float(np.sqrt(bb/w.sum())/scale),scope='absolute pressure difference vectors at equal physical times; no mean subtraction or observed-order admission')
        out.append(dict(time=t,states=[exact_vector_hash(h[t]) for h in histories],pressure=rows))
    return dict(schema='CT9PressureVectorTripletV1',dt=[.002,.001,.0005],rows=out,scale=scale,observation_rule='sealed piecewise6-point physical Gauss rule; affine P2 pressure squared is exactly polynomial-integrated; governing nonlinear quadrature studied separately')

def compact_pair(q):
    out={k:copy.deepcopy(q[k]) for k in ('schema','required_dt','target','same_actual_initial','temporal','traction_temporal')}
    out['qualification']=pair(q);out['members']={}
    for name,r in q['members'].items():
        out['members'][name]={k:copy.deepcopy(r[k]) for k in ('dt','time','directory','failure','signature','native_report_sha256','history_sha256','actual_initial_sha256','initialization_receipt','domain','constitutive','corrections','extra_solves','native_audited_extrema','final_native_audit') if k in r}
        out['members'][name]['initialization_receipt']={k:v for k,v in r['initialization_receipt'].items() if k not in ('effective_initial','live_fixed_boundary')}
        out['members'][name]['accepted_metrics']=r['accepted_metrics'];out['members'][name]['budgets']=r['budgets']
        file=Path(r['directory'])/'case.evolution.json';require(sha(file)==r['native_report_sha256'],'changed native trajectory during final export');native=n.read(file)
        keys=('time','alpha','scheduled_complete','linear_backend','predictor_linear_residual','pressure_linear_residual','native_full_G_action_defect','momentum_splitting_identity_defect','corrected_continuity_gauge','frozen_momentum_defect','pressure_krylov_iterations','mass_inverse_applications','imposed_trace_defect','refreshed_trace_defect','bound_calls','additional_governing_equation_solves')
        out['members'][name]['actual_stage_checks']=[{k:q.get(k) for k in keys} for q in native['pimple_attempts']]
    return out

def failure_witnesses(root):
    """Preserve first attempts without copying large native state packets."""
    out={}
    for label in ('build-v1','build-v2','build-v6','controls-v1','pressure-finer-v1','capacity-v1','linear-recheck-intermediate-v1','linear-recheck-p4-v1'):
        folder=root/label
        if not folder.exists():continue
        paths=[p for p in (folder/'request-manifest.json',folder/'completed-artifacts.json',folder/'results.json',folder/'progress.json',folder/'0.log',folder/'1.log') if p.exists()]
        row=dict(path=str(folder),artifacts={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in paths})
        if (folder/'results.json').exists():
            q=n.read(folder/'results.json');row.update(status=q.get('status'),failure=q.get('failure'),qualification=q.get('qualification'))
            if label=='capacity-v1':
                row['cases']={name:{k:compact(r[k]) for k in ('status','failure','native_count','zero_state_sha256','native_report_sha256','history','startup_execution','frozen_execution') if k in r} for name,r in q['cases'].items()}
            if label=='controls-v1':row['checks']={name:{k:r[k] for k in ('status','failure') if k in r} for name,r in q['checks'].items()}
        out[label]=row
    return out

def residual_refinements(folders):
    out={};prefix='RANS sparse refinement '
    for folder in folders:
        for path in folder.rglob('*.log'):
            records=[json.loads(line.split(prefix,1)[1]) for line in path.read_text(errors='replace').splitlines() if prefix in line]
            if not records:continue
            require(all(1<=r['steps']<=3 and 0<=r['final_residual']<=r['bound'] for r in records),'linear refinement log violates unchanged residual contract')
            out[str(path)]=dict(sha256=sha(path),factor_reused_solve_count=len(records),linear_residual_corrections=sum(r['steps'] for r in records),maximum_steps=max(r['steps'] for r in records),worst_final_residual_over_bound=max(r['final_residual']/r['bound'] for r in records),first_record=records[0])
    return dict(scope='linear algebra corrections within the same frozen-system solve, not extra governing equation stages or PIMPLE sweeps',logs=out)

def target_transit(ct8,capacity):
    from scipy.integrate import quad
    c=sealed(ct8/'preflight-v6')['case'];column=n.col(c);H=c['H'];d=c['displacement'];target=capacity['cases']['target-p4']['case'];Z=target['Z'];length=target['x'][-1]-target['x'][0];us=c['ustar'];kappa=parameters(c['model'])['kappa']
    require(c['model']=='A95' and Z==20. and length==55.,'foreign target/reference transit request')
    cuts=sorted({0.,H,*[z for z in (column.zv,c['roughness_input'],column.cutoff) if 0<z<H]})
    inner=sum(quad(lambda z:us**2*(H-z)/column.D(z),a,b,epsabs=1e-14,epsrel=1e-12)[0] for a,b in zip(cuts[:-1],cuts[1:],strict=True))
    primitive=lambda y:c['B']*y+us/kappa*(y*math.log(y)-y)
    flow=inner+primitive(Z-d)-primitive(H-d);bulk=flow/Z;transit=length/bulk
    return dict(schema='CT9TargetReferenceTransitV1',definition='length divided by full-height reference bulk velocity, including finite-nu inner resistance; independent quadrature plus analytic outer-log integral',length=length,roof=Z,inner_integral=inner,volume_flow_per_width=flow,bulk_velocity=bulk,transit_time=transit,stage_times={str(f):f*transit for f in (.01,.1,1.)},physical_trajectories='not_admitted: independent pressure gate failed; no target timestep or campaign is silently assigned')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct8','output'):p.add_argument('--'+k,type=Path,required=True)
    for k in ('controls','equivalence','capacity','build'):p.add_argument('--'+k,required=True,help='immutable experiment label under --root')
    a=p.parse_args();labels=[a.controls,a.equivalence,a.capacity,a.build]
    require(all(Path(x).name==x and x not in ('.','..') for x in labels),'invalid evidence label')
    required=['pressure-comparison-v2','pressure-diagnosis-v1',*labels];sources={};bodies={}
    for label in required:
        path=a.root/label
        if (path/'completed-artifacts.json').exists():sources[label]=dict(path=str(path),sha256=sha(path/'results.json'),request=sha(path/'request-manifest.json'),seal=sha(path/'completed-artifacts.json'));bodies[label]=sealed(path)
        else:sources[label]=dict(status='missing_or_incomplete',path=str(path))
    m=RequestManifest(a.output,dict(schema='CT9FinalEvidenceV1',source=sha(Path(__file__)),sources=sources,contract=sha(Path(__file__).resolve().parents[4]/'docs/plans/cao-tamura-accepted-9-evidence/contract.json')))
    pressure=bodies.get('pressure-comparison-v2');diagnosis=bodies.get('pressure-diagnosis-v1');controls=bodies.get(a.controls);equivalence=bodies.get(a.equivalence);capacity=bodies.get(a.capacity);build=bodies.get(a.build)
    binding={};expected=build['runtime']['IncNavierStokesSolver'] if build else None
    for label in (a.controls,a.equivalence,a.capacity):
        request=a.root/label/'request-manifest.json'
        actual=n.read(request)['request'].get('runtime') if request.exists() else None
        binding[label]=dict(status='pass' if expected and actual==expected else 'fail',runtime=actual)
    write_json(a.output/'runtime-binding.json',compact(binding))
    if pressure:
        write_json(a.output/'pressure-pair.json',compact(compact_pair(pressure)));write_json(a.output/'pressure-vectors.json',pressure_vectors(a.root,a.ct8))
    if diagnosis:write_json(a.output/'pressure-diagnosis.json',compact(diagnosis))
    if controls:write_json(a.output/'backend-controls.json',compact(controls))
    if equivalence:
        write_json(a.output/'backend-pair.json',compact(dict(physical_pair=compact_pair(equivalence),backend_status=backend_trajectory(equivalence),equivalence=equivalence['equivalence'],physical_continuation_admitted=False)))
    if capacity:
        cap=copy.deepcopy(capacity)
        for row in cap['cases'].values():
            if 'sparse_backend' in row:
                for k in ('predictor','pressure_increment','consistent_mass_correction'):row['sparse_backend'].pop(k,None)
                for scalar in row['sparse_backend']['scalars']:scalar.pop('diagnostic_increment',None)
            for key in ('initialization',):
                if key in row:row[key]={k:v for k,v in row[key].items() if k not in ('effective_initial','live_fixed_boundary')}
        write_json(a.output/'capacity.json',compact(cap))
    if build:write_json(a.output/'runtime.json',compact(build))
    write_json(a.output/'retained-attempts.json',failure_witnesses(a.root))
    if capacity:write_json(a.output/'target-transit.json',target_transit(a.ct8,capacity))
    profiles={str(p.relative_to(a.root)):timing(p) for p in list((a.root/a.equivalence).glob('*/timing.jsonl'))+list((a.root/a.capacity).glob('*/*-timing.jsonl'))}
    write_json(a.output/'cost.json',profiles)
    write_json(a.output/'linear-refinement.json',residual_refinements([a.root/x for x in (a.controls,a.equivalence,a.capacity)]))
    resources={}
    for file in sorted(a.root.glob('launch-*/execution.json')):
        q=n.read(file);q['high_water_receipts']={p.name:dict(sha256=sha(p),records=[json.loads(x) for x in p.read_text().splitlines()]) for p in file.parent.glob('*highwater.jsonl')};events=file.parent/'worker-guard-events.jsonl';q['suspensions']=[json.loads(x) for x in events.read_text().splitlines()] if events.exists() else [];resources[str(file.parent.name)]=q
    resource_pass=bool(resources) and all(q['resource']['aggregate_status']==q['resource']['project_status']=='pass' for q in resources.values())
    write_json(a.output/'resources.json',compact(resources));write_json(a.output/'loaded-libraries.json',LIBRARIES)
    decision=dict(schema='CT9RequiredScopeV1',physical_reference='CT8 central mathematical feasibility retained; eight failed extraction corners unchanged',pressure_reconstruction='pass' if diagnosis and diagnosis['status']=='pass' else 'incomplete',pressure_uncertainty=pair(pressure) if pressure else dict(status='incomplete'),backend_operator='pass' if controls and backend_controls(controls)['status']=='pass' else 'failed_or_incomplete',backend_trajectory='pass' if equivalence and backend_trajectory(equivalence)['status']=='pass' else 'failed_or_incomplete',capacity=capacity_gate(capacity) if capacity else dict(status='incomplete'),resources='pass' if resource_pass else 'failed_or_incomplete',target_preservation='not_admitted',terrain_readiness='unsupported mapped field numbering/operators and historical curved temporal/degree gates remain open',sources=sources)
    decision['physical_continuation_admitted']=False
    decision['runtime_binding']='pass' if all(q['status']=='pass' for q in binding.values()) else 'fail'
    decision['resource_closure']='listed numerical experiment controls are closed; final collector receipt must be bound after this process exits in the publication scope'
    # A completed investigation may have an explicit failed physical gate.
    # Its completion cannot become permission to evolve the target domain.
    decision['investigation_completion']='complete' if decision['pressure_reconstruction']==decision['backend_operator']==decision['backend_trajectory']==decision['resources']==decision['runtime_binding']=='pass' and decision['capacity']['status']=='pass' and pressure else 'incomplete'
    write_json(a.output/'decision.json',decision);files=[f.name for f in a.output.iterdir() if f.suffix=='.json' and f.name not in ('request-manifest.json','completed-artifacts.json')];m.seal(a.output,files);print(decision,flush=True);return 0 if decision['investigation_completion']=='complete' else 1
if __name__=='__main__':raise SystemExit(main())
