#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh selected-space pilot with unchanged CT11 physical classifiers and horizon."""
import argparse,json,math,os,resource,time
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted15_controls import ordering,COMPONENT
from ct_accepted13_orderings import native_run
from ct_accepted9_controls import enable
from ct_accepted9_pressure_run import sealed
from ct_accepted11_target import target_pair,bounded_budgets
from ct_accepted14_gate import cohort_pair,COHORTS
from ct_accepted10_run import comparison
from ct_accepted10_compare import gradients,outlet,regional_pressure
from ct_accepted3_metrics import inventory
from ct_accepted_pilot import replay,constitutive_comparison,artifact_hash
from ct_accepted5_traction import profiles
from ct_accepted_open import run,probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

PHASES={}
def measured(name,function):
    def call(*args,**kwargs):
        start=time.monotonic();cpu=time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
        try:return function(*args,**kwargs)
        finally:
            after=resource.getrusage(resource.RUSAGE_CHILDREN);r=PHASES.setdefault(name,dict(calls=0,self_CPU=0.,children_CPU=0.,wall=0.))
            r['calls']+=1;r['self_CPU']+=time.process_time()-cpu;r['children_CPU']+=after.ru_utime+after.ru_stime-child.ru_utime-child.ru_stime;r['wall']+=time.monotonic()-start
    return call

def numerical_receipts(path):
    rows=[json.loads(s) for s in path.read_text().splitlines()];factors={};solves=[]
    for row in rows:
        if row['event']=='factor':
            require(row['id'] not in factors and 0<row['factor_nnz']<=40000000 and 0<row['matrix_nnz']<=20000000,'actual produced factor guard');factors[row['id']]=row
        else:
            require(row['event']=='solve' and row['id'] in factors,'unbound linear solve receipt');rhs=row['rhs_maximum'];res=row['true_residual']
            require(all(math.isfinite(x) and x>=0 for x in (rhs,res,row['bound'])) and res<=2e-12*rhs+1e-25 and 0<=row['refinements']<=3,'actual produced original-system residual/refinement');solves.append(row)
    require(factors and solves,'missing numerical receipts')
    return dict(status='pass',path=str(path),sha256=sha(path),factor_count=len(factors),solve_count=len(solves),maximum_factor_entries=max(x['factor_nnz'] for x in factors.values()),minimum_fill_margin=min(40000000-x['factor_nnz'] for x in factors.values()),maximum_matrix_nonzeros=max(x['matrix_nnz'] for x in factors.values()),maximum_refinements=max(x['refinements'] for x in solves),maximum_residual_bound_fraction=max(x['true_residual']/x['bound'] for x in solves),factors=list(factors.values()))

def main():
    global bounded_budgets,probe,profiles,inventory
    n.physical_metrics=measured('independent_accepted_metrics',n.physical_metrics)
    bounded_budgets=measured('complete_endpoint_scalar_controls',bounded_budgets)
    probe=measured('endpoint_constitutive_probe',probe)
    profiles=measured('independent_traction_profiles',profiles)
    inventory=measured('independent_endpoint_inventory',inventory)
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','static','component','capacity','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--admission',type=Path)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1';os.environ.pop('NEKTAR_RANS_TIMING_FILE',None);n.reference=lru_cache(None)(n.reference)
    pre=sealed(a.static);component=sealed(a.component);capacity=sealed(a.capacity);require(pre['status']==component['status']==capacity['status']=='pass','pair prerequisites incomplete')
    co=n.read(a.contract);dts=co['dt'];require(tuple(dts) in COHORTS and co['T']==.002,'unadmitted short-pair contract')
    retry=dts!=[.001,.0005];require(co['schema']==('CT15AdmittedShortPairContractV2' if retry else 'CT15AdmittedShortPairContractV1'),'foreign cohort contract')
    classify=(lambda q:cohort_pair(q,dts)) if retry else target_pair
    require(co['static_sha256']==sha(a.static/'results.json') and co['component_sha256']==sha(a.component/'results.json') and co['capacity_sha256']==sha(a.capacity/'results.json'),'stale pair prerequisites')
    rt=runtime_identity(a.runtime/'IncNavierStokesSolver');pr=runtime_identity(a.runtime/'RANSPreservationProbe');request=n.read(a.static/'request-manifest.json')['request'];require(rt==request['runtime'] and pr==request['probe'],'foreign static runtime')
    if retry:
        require(a.admission is not None and co['timestep_admission_sha256']==sha(a.admission/'results.json'),'missing changed-alpha admission')
        admission=sealed(a.admission);ar=n.read(a.admission/'request-manifest.json')['request'];require(admission['status']=='pass' and ar['static']==sha(a.static/'results.json') and ar['runtime']==rt and ar['probe']==pr and admission['dt']==dts,'foreign/incomplete changed-alpha admission')
    c=pre['case'];ini=n.read(a.static/'target/frozen.json');require(sha(a.static/'target/frozen.json')==pre['native_packet_sha256'],'static observation rows changed');weights=n.inventory(c)[2];count=len(weights);first=None;initials={};reports={};zeros={}
    m=RequestManifest(a.output,dict(schema='CT15FreshSelectedPairV2' if retry else 'CT15FreshSelectedPairV1',source=sha(Path(__file__)),contract=sha(a.contract),static=sha(a.static/'results.json'),component=sha(a.component/'results.json'),capacity=sha(a.capacity/'results.json'),runtime=rt,probe=pr,case=c,dt=dts,T=.002,component_identity=COMPONENT,detailed_timing=False,physical_classifier_source=sha(Path(__file__).with_name('ct_accepted11_target.py')),request_specialization_source=sha(Path(__file__).with_name('ct_accepted14_gate.py')) if retry else None,timestep_admission=sha(a.admission/'results.json') if retry else None))
    out=dict(schema='CT11TargetPhysicalPairV1',CT15_schema='CT15FreshSelectedPairV2' if retry else 'CT15FreshSelectedPairV1',required_dt=dts,target=.002,case=c,prerequisites=dict(status='pass'),members={},initial={},detailed_timing=False)
    regions=np.asarray(ini['observation']['grid_regions'][:count]);wt=np.asarray(weights);col=n.col(c);reference_pressure=np.array([n.reference(col,q[2],bool(r))[3] for q,r in zip(ini['observation']['grid'][:count],regions,strict=True)]);initial_pressure=n.observations(ini,ini['observation']['state'],count)[:,4]
    for dt in dts:
        name=str(dt);dest=a.output/('zero-'+name);n.make(dest,a.runtime,c,dt,0);enable(dest);ordering(dest);os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(dest/'sparse-receipts.jsonl')
        e=native_run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log',300,600);z=n.read(dest/'case.evolution.json');require(e['exit_code']==0 and not z['failure'] and z['time']==0. and exact_vector_hash(z['state'])==pre['initial']['state_sha256'],'fresh selected actual initialization differs')
        require(z['state']==ini['observation']['state'] and (first is None or first==z['state']) and COMPONENT in z['signature'],'pair actual initialization/component');first=z['state'];initials[name]=ini;zeros[name]=(dest,z,e)
        out['initial'][name]=dict(actual_initial_sha256=exact_vector_hash(z['state']),history_sha256=artifact_hash(z['checkpoint']),native_report_sha256=sha(dest/'case.evolution.json'),configuration_sha256=sha(dest/'case.xml'),execution=e,numerical_receipts=numerical_receipts(dest/'sparse-receipts.jsonl'));write_json(a.output/'progress.json',out)
    out['same_actual_initial']=True
    for dt in dts:
        name=str(dt);zero,z,ze=zeros[name];folder=a.output/name;n.make(folder,a.runtime,c,dt,round(.002/dt),Path(z['checkpoint']),zero/'case.xml');ordering(folder);os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(folder/'sparse-receipts.jsonl')
        e=native_run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log',4000 if retry else 1800,4800 if retry else 2400);os.environ.pop('NEKTAR_RANS_SPARSE_RECEIPTS',None);r=n.read(folder/'case.evolution.json');require(r['accepted_states'][0]==z['state'],'new pair did not consume own genuine history')
        row=dict(dt=dt,time=r['time'],directory=str(folder),failure=r['failure'],signature=r['signature'],execution=e,zero_step_execution=ze,actual_initial_sha256=exact_vector_hash(z['state']),state_sha256=exact_vector_hash(r['state']),configuration_sha256=sha(folder/'case.xml'),native_report_sha256=sha(folder/'case.evolution.json'),history_sha256=artifact_hash(r['checkpoint']),initialization_receipt=r['initialization_receipt'])
        row['accepted_metrics']=[dict(time=t,**n.physical_metrics(ini,s,c,weights,count)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)];row['traction_history']=[dict(time=t,profile=profiles(ini,s,c)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)];row['actual_initial_traction_profile']=profiles(ini,z['state'],c)
        fixed=inventory(ini['test_basis'],ini['observation']['grid'],np.asarray(r['state']),ini['observation']['offsets'],'A95',ini['observation']['grid_sides']);row['domain']=dict(status='pass' if fixed['domain_failure_count']==0 and fixed['minima']['covariance']['value']>=-1e-10*c['K0'] else 'fail',fixed=fixed,native_audited_extrema=audited_extrema(r),native_evolution_failure=r['failure'])
        row['gradient_reference_errors']=gradients(ini,r['state'],c,weights,count);row['represented_outlet_traction']=outlet(ini,r['state'],c);row['regional_pressure']=[]
        for t,s in zip(r['audit_times'],r['accepted_states'],strict=True):
            pi=n.observations(ini,s,count)[:,4];row['regional_pressure'].append(dict(time=t,absolute_reference_error=regional_pressure(pi,reference_pressure,wt,regions,c['ustar']**2),drift=regional_pressure(pi,initial_pressure,wt,regions,c['ustar']**2)))
        row['budgets']=bounded_budgets(a.runtime,folder,ini,r,c) if r['steps'] else dict(gate=dict(status='missing'))
        grid,sides,_,_=n.inventory(c);request=replay(ini,r['state'],r['time'],grid,grid_sides=sides,export_constitutive=True);q,pe=probe(a.runtime,folder,request,'constitutive');row['constitutive']=constitutive_comparison(q,c) if pe['exit_code']==0 else dict(status='fail',failure=q.get('failure'));row['constitutive_execution']=pe
        row['corrections']=sum(x.get('bound_calls',0) for x in r['pimple_attempts']);row['extra_solves']=sum(x.get('additional_governing_equation_solves',0) for x in r['pimple_attempts']);row['actual_stage_checks']=[{k:v for k,v in x.items() if k in ('time','pressure_krylov_iterations','flow_matrix_nonzeros','predictor_factor_nonzeros','mass_factor_nonzeros','pressure_preconditioner_factor_nonzeros','predictor_linear_residual','pressure_linear_residual','native_full_G_action_defect','momentum_splitting_identity_defect','frozen_momentum_defect','corrected_continuity_gauge','bound_calls','additional_governing_equation_solves') or 'trace_defect' in k or 'domain' in k or k.endswith(('_ground','_outlet'))} for x in r['pimple_attempts']]
        row['numerical_receipts']=numerical_receipts(folder/'sparse-receipts.jsonl');out['members'][name]=row;reports[name]=r;write_json(a.output/'progress.json',out);print('CT15 fresh selected',name,r['time'],r['failure'],flush=True)
        if not r['failure'] and abs(r['time']-.002)<1e-13:
            require(sum(x['ordering']=='VectorCOLAMDV1' for x in row['numerical_receipts']['factors'])==r['steps'],'every actual predictor uses the admitted component exactly once per step')
        if e['exit_code']!=0 or r['failure'] or abs(r['time']-.002)>1e-13:
            out['qualification']=classify(out);out['first_failure']=dict(member=name,last_accepted_time=r['time'],failure=r['failure']);write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);return 1
    comparison(out,initials,reports,c,weights,count,classifier=classify)
    common=[{round(t,12):s for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)} for r in reports.values()];out['regional_temporal_pressure']=[dict(time=t,pressure=regional_pressure(*(n.observations(ini,x[t],count)[:,4] for x in common),wt,regions,c['ustar']**2)) for t in sorted(set(common[0])&set(common[1]))]
    out['clock_off_performance']={name:dict(total_native_CPU=row['execution']['process_cpu_seconds'],total_native_wall=row['execution']['wall_seconds'],accepted_steps=round(.002/row['dt']),zero_step_startup_CPU=row['zero_step_execution']['process_cpu_seconds'],scope='actual detailed-timing-off execution including startup/mandatory checks/output; no inferred exclusive phase speedup or long-run affordability') for name,row in out['members'].items()}
    out['workflow_phases']=PHASES
    out['workflow_phase_scope']='outer nonoverlapping wrappers; endpoint scalar wrapper contains its own native probes; do not add its children to native executions again; no detailed inner timers'
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT15 physical qualification',out['qualification'],flush=True);return 0 if out['qualification']['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
