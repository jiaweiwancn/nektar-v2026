#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT11 primitives, exact source/runtime roles and closed scope decisions."""
import argparse,copy,json,os,subprocess
from pathlib import Path
from collections import defaultdict
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted9_pressure_run import sealed
from ct_accepted11_gate import observer,capacity
from ct_accepted11_target import target_pair
from ct_accepted_feasibility import write_json,require
from rans_request_manifest import RequestManifest,sha,runtime_identity


def read(p):return json.loads(Path(p).read_text())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct9','ct10','runtime','equivalence','capacity','controls','self-control','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--target',type=Path);a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    labels={'equivalence':a.equivalence,'capacity':a.capacity,'controls':a.controls};data={k:sealed(v) for k,v in labels.items()}
    if a.target:labels['target']=a.target;data['target']=sealed(a.target)
    requests={k:dict(directory=str(v),result_sha256=sha(v/'results.json'),request_sha256=sha(v/'request-manifest.json'),seal_sha256=sha(v/'completed-artifacts.json')) for k,v in labels.items()}
    resources={};peak_project=peak_aggregate=0
    for f in sorted(a.root.glob('resources-*')):
        if f==a.self_control:continue
        require((f/'execution.json').exists(),'unclosed resource family:'+str(f));q=read(f/'execution.json');r=q['resource'];require(r['aggregate_status']==r['project_status']=='pass','resource violation:'+str(f))
        peak_project=max(peak_project,r['maximum_project']);peak_aggregate=max(peak_aggregate,r['maximum_aggregate']);resources[f.name]=dict(execution=q,sha256=sha(f/'execution.json'),highwater={str(p.name):sha(p) for p in f.iterdir() if p.suffix=='.jsonl'})
    revision=subprocess.check_output(['git','-c','core.preloadIndex=false','-c','index.threads=1','rev-parse','HEAD'],cwd=a.root/'source',text=True).strip()
    m=RequestManifest(a.output,dict(schema='CT11CompactEvidenceV1',controller_revision=revision,inputs=requests,resources={k:v['sha256'] for k,v in resources.items()},source=sha(Path(__file__)),scope='collection completion is separate from every numerical/resource/physical status'))
    catalog={};traction_catalog={}
    def compact(v):
        if isinstance(v,dict) and {'state_sha256','ground','interface','interface_jump','mean_wall_shear'}<=set(v):
            import hashlib
            digest=hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest();traction_catalog[digest]=v
            return dict(schema='CT11TractionProfileReferenceV1',profile_sha256=digest,artifact='traction-profiles.json',state_sha256=v['state_sha256'],ground_points=len(v['ground']),interface_points=len(v['interface']),mean_wall_shear=v['mean_wall_shear'],normalized_wall_traction_error=v['normalized_wall_traction_error'],normalized_interface_traction_RMS=v['normalized_interface_traction_RMS'],wall_spatial_standard_deviation=v['wall_spatial_standard_deviation'])
        if isinstance(v,dict):
            out={}
            for k,x in v.items():
                if k=='effective_initial' and isinstance(x,list):
                    out['effective_initial_artifact']=dict(entries=len(x),native_sha256=v.get('effective_initial_sha256'),scope='full coefficient vector remains in the bound native receipt outside Git')
                elif k=='actual_mapped_files':
                    # Map contents are cryptographically bound once and shared.
                    import hashlib
                    h=hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest();catalog[h]=x;out['actual_mapped_files_catalog_sha256']=h
                else:out[k]=compact(x)
            return out
        if isinstance(v,list):return [compact(x) for x in v]
        return v
    derived_observer=observer(data['equivalence']);derived_capacity=capacity(data['capacity']);target=data.get('target');derived_target=target_pair(target) if target else dict(status='unattempted',reasons=['target prerequisites not yet admitted'])
    decision=dict(schema='CT11RequiredScopeV1',investigation_completion='complete',observer_equivalence=derived_observer,observer_capacity=derived_capacity,negative_controls=data['controls']['status'],target_pair=derived_target,actual_paired_horizon=min((r['time'] for r in target['members'].values()),default=0.) if target else 0.,per_member_horizons={k:r['time'] for k,r in target['members'].items()} if target else {},resources=dict(status='pass',maximum_project=peak_project,maximum_aggregate=peak_aggregate),first_01_transit='not_started',terrain='unattempted; mapped numbering,joined gradients/cofactors,transmission and curved temporal/spatial qualification remain separate',physical_applicability='central reference retained; numerical flat preservation is not experimental calibration,developed precursor or hill validation',historical='all CT4–CT10 seals/failures unchanged')
    if target and not target.get('first_failure'):
        plan=target.get('first_stage_resource_plan',{});families=[v['execution'] for v in resources.values() if any('ct_accepted11_target.py' in str(x) for x in v['execution']['command'])];require(len(families)==1,'target closed-family resource identity ambiguous');family=families[0];driver=plan['completed_driver_process_CPU'];require(driver>0,'missing measured driver CPU');multiplier=max(1.,family['process_cpu_seconds']/driver);remaining=plan['CPU_with_25percent_headroom']*multiplier;totalCPU=family['process_cpu_seconds']+remaining;totalWall=family['wall_seconds']+remaining*max(1.,family['wall_seconds']/family['process_cpu_seconds']);bytes_now=sum(p.stat().st_size for p in a.target.rglob('*') if p.is_file());totalStorage=bytes_now+plan['storage_bytes_with_25percent_headroom'];closed_admission=dict(complete_family_process_CPU=family['process_cpu_seconds'],measured_controller_and_children_CPU=driver,measured_family_overhead_multiplier=multiplier,total_CPU_spent_plus_remaining_with_headroom=totalCPU,CPU_limit=24000,total_wall_upper_estimate=totalWall,wall_limit=28800,measured_current_storage_bytes=bytes_now,total_storage_spent_plus_remaining_with_headroom=totalStorage,storage_limit_bytes=24*1024**3,resource_admitted=totalCPU<=24000 and totalWall<=28800 and totalStorage<=24*1024**3,scope='closed owned-family accounting includes monitor/helpers and spent work; future estimates retain25percent headroom and finite ceilings,not an accomplished horizon');decision['first_stage_resource_plan']=plan;decision['closed_campaign_resource_admission']=closed_admission;decision['continuation_admitted']=derived_target['status']=='pass' and closed_admission['resource_admitted']
        decision['next_gate']='conditional target continuation' if decision['continuation_admitted'] else ('target field/traction/budget/timestep primitives' if derived_target['status']!='pass' else 'finite target campaign resources')
    elif target:decision['continuation_admitted']=False;decision['first_failure']=target['first_failure'];decision['next_gate']='first native target rejection; no further stepping'
    else:decision['continuation_admitted']=False;decision['next_gate']='observer equivalence/capacity or required controls'
    for k,v in data.items():write_json(a.output/(k+'.json'),compact(v))
    write_json(a.output/'resources.json',compact(resources));write_json(a.output/'decision.json',decision)
    costs={}
    if target:
        for k,r in target['members'].items():
            f=Path(r['directory']);names=defaultdict(lambda:dict(calls=0,exclusive_process_cpu_seconds=0.,exclusive_wall_seconds=0.))
            for line in (f/'timing.jsonl').read_text().splitlines():
                q=json.loads(line);name=q['name'];z=names[name];z['calls']+=q.get('calls',1);z['exclusive_process_cpu_seconds']+=q['exclusive_process_cpu_seconds'];z['exclusive_wall_seconds']+=q['exclusive_wall_seconds']
            costs[k]=dict(native_execution=compact(r['execution']),zero_step_execution=compact(r['zero_step_execution']),timing_path=str(f/'timing.jsonl'),timing_sha256=sha(f/'timing.jsonl'),exclusive_categories=dict(names),accepted_steps=len(r['accepted_metrics'])-1,native_failure=r['failure'],scope='native attempted/accepted timings distinguished by report; exclusive categories only, no sum of nested inclusive timings; rejected attempts are not accepted-step cost evidence')
    write_json(a.output/'cost.json',costs)
    runtime=runtime_identity(a.runtime/'RANSPreservationProbe');old=runtime_identity(a.ct9/'build-v8/runtime/IncNavierStokesSolver');require(old['sha256']=='5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7','governing runtime changed')
    write_json(a.output/'runtime.json',dict(governing=old,diagnostic=runtime,compiled_probe=read(a.root/'build-v2/request-manifest.json'),build=compact(sealed(a.root/'build-v2')),observer_identity='NativeOwnerVolumeStreamV1',controller_revision=revision))
    historical=read(a.ct10/'qualified-v1/decision.json');require(historical['baseline_physical_uncertainty']['status']=='fail' and historical['candidate_T020']['status']=='pass' and historical['target_required_audit_capacity']['status']=='fail','historical CT10 status changed')
    write_json(a.output/'history-and-provenance.json',dict(starting_revision='8b2d56610ee12e8a2144a6becad2aa5150f263fb',controller_revision=revision,artifacts=requests,CT10_decision=dict(path=str(a.ct10/'qualified-v1/decision.json'),sha256=sha(a.ct10/'qualified-v1/decision.json'),old_capacity='fail',old_target_horizon=0.),retained_failed_build=dict(directory=str(a.root/'build-v1'),progress_sha256=sha(a.root/'build-v1/progress.json'),reason='warning-as-error misleading indentation in new digest loop; no native history generated'),source_scope='only read-only volume iteration,probe and qualification support changed; governing CT9 runtime never rebuilt; all physical fits,functions and signatures retained'))
    critical=['library/SolverUtils/Turbulence/RANSModels.cpp','library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPimple.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPicard.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionSparse.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallSerialSparse.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionOpenBoundary.cpp']
    import hashlib
    unchanged={}
    for name in critical:
        oldbytes=subprocess.check_output(['git','show','8b2d5661:'+name],cwd=a.root/'source');h=hashlib.sha256(oldbytes).hexdigest();require(h==sha(a.root/'source'/name),'governing source changed:'+name);unchanged[name]=h
    name='solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionDiagnostics.cpp';oldtext=subprocess.check_output(['git','show','8b2d5661:'+name],cwd=a.root/'source',text=True);newtext=(a.root/'source'/name).read_text();marker='std::string RANSWallEvolutionPatch::PreservationPacket(';require(oldtext[oldtext.index(marker):]==newtext[newtext.index(marker):],'existing diagnostic implementations changed')
    write_json(a.output/'source-scope.json',dict(starting_revision='8b2d56610ee12e8a2144a6becad2aa5150f263fb',controller_revision=revision,unchanged_governing_sources=unchanged,existing_PreservationPacket_and_later_diagnostics='byte-identical from starting head; new volume visitor inserted before them',changed_files=subprocess.check_output(['git','diff','--name-only','8b2d5661','HEAD'],cwd=a.root/'source',text=True).splitlines(),observer_storage='borrows existing governing quadrature rows; no second serialized/global observer volume catalog; intrinsic native row storage is unchanged'))
    if traction_catalog:
        import hashlib
        packed={}
        for digest,v in traction_catalog.items():
            q=copy.deepcopy(v)
            for name in ('ground','interface','interface_jump'):
                rows=q[name]
                if not rows:continue
                columns=list(rows[0]);require(all(set(r)==set(columns) for r in rows),'heterogeneous accepted traction record')
                q[name]=dict(schema='CT11ColumnarTractionRowsV1',columns=columns,rows=[[r[k] for k in columns] for r in rows])
            restored=copy.deepcopy(q)
            for name in ('ground','interface','interface_jump'):
                if isinstance(restored[name],dict):
                    table=restored[name];restored[name]=[dict(zip(table['columns'],r,strict=True)) for r in table['rows']]
            require(hashlib.sha256(json.dumps(restored,sort_keys=True,separators=(',',':')).encode()).hexdigest()==digest,'lossy traction publication')
            packed[digest]=q
        write_json(a.output/'traction-profiles.json',dict(schema='CT11LosslessTractionCatalogV1',profiles=packed,hash_format='SHA256 of restored original profile,sorted-key compact JSON; reconstruct each columnar row with its columns; all physical values,gradients,sides,geometry and order retained'))
    write_json(a.output/'loaded-libraries.json',catalog)
    files=[str(p.relative_to(a.output)) for p in a.output.iterdir() if p.is_file() and p.name not in ('request-manifest.json','completed-artifacts.json')];m.seal(a.output,files);print(json.dumps(decision,indent=2),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
