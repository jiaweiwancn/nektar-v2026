#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only, identity-bound all-state WP2 physical-transit assessment.

Actual native coefficient fields are evaluated on immutable CG geometry rows.
This does not evolve, project, reset or repair any retained state.
"""
import argparse,math,time,resource,shutil
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,state_hash,sha
from wp2_all_state import Rows,US
from wp2_qualification import verify_report,quarter_budget_differences,homogeneity
from wp2_common_quadrature import compare_orders,verify_native_tuple_replay

TREQ=1.1942387693600194
K0={'CLS':.017598891430879275,'standard':.01838575534504143}


def analysis_identity():
    here=Path(__file__)
    return {name:sha(here.with_name(name)) for name in [here.name,'wp2_all_state.py','wp2_common_quadrature.py','wp2_qualification.py','cao_tamura_h2g_constitutive.py']}


def tensor_homogeneity(rows,current,cold):
    nz=(len(rows.meta['physical_z_partition'])-1)*8;weight=rows.weight.reshape(48,nz);wx=weight[:,0]/sum(weight[:,0]);out={}
    for name in ['R','T']:
        v=current[name].reshape(48,nz,3,3);mean=np.einsum('x,xzij->zij',wx,v)
        differences={'inlet_to_station':current[name]-cold[name],
                     'horizontal_variation':(v-mean[None,:,:,:]).reshape(-1,3,3)}
        out[name]={}
        for label,delta in differences.items():
            out[name][label]=[]
            for region in [0,1]:
                mask=rows.inner if region==0 else ~rows.inner;w=rows.weight[mask]
                out[name][label].append(float(np.sqrt(np.sum(w[:,None,None]*delta[mask]**2)/sum(w))/US**2))
    return dict(regional_relative_RMS=out,limit=.01,passed=all(v<=.01 for family in out.values() for row in family.values() for v in row))


def paths(root,model,phase,progress):
    prefix={'coarse':[f'trial-dt020-{model}/step2/case.evolution.json',f'runtime-equivalence-{model}/new/one/case.evolution.json'],
            'fine':[f'trial-dt010-{model}/step4/case.evolution.json'],
            'confirmation':[f'confirmation-trial-{model}/step2/case.evolution.json'],
            'recovery':[]}[phase]
    if progress['case_definition'].get('pressure_space','baseline-v1')!='baseline-v1':prefix=[]
    reports=[root/p for p in prefix]
    for obs in progress['observations']:
        p=Path(obs['report'])
        if p not in reports:reports.append(p)
    return reports


def assess(root,model,phase,rows,output,progress_path=None,independent_order=28):
    started=time.monotonic();progress_path=progress_path or root/f'preservation-{phase}-{model}/progress.json';progress_hash=sha(progress_path);p=load(progress_path)
    require(p['model']==model,'analysis trajectory model differs')
    require(p['status']=='completed_bounded_trajectory','trajectory is not complete: '+str(progress_path))
    require(p['actual_time']>=TREQ,'actual physical duration too short')
    recovery=phase=='recovery';history=[];last_hash=None;last_time=-1.;sources=[];replay=[];counts={};previous_values=None
    if recovery:
        # A zero-step seed report omits transient endpoint-sampler output.
        # The first evolved report includes that exact seed as accepted state
        # zero, together with the complete native monitoring inventory.
        seed_file=Path(p['prefix_observations'][0]['report']) if p.get('prefix_observations') and p['prefix_observations'][0]['time']==0. else root/f'recovery-seed-{model}/h0/case.evolution.json'
        first=load(paths(root,model,phase,p)[0])
        if p['case_definition'].get('pressure_space','baseline-v1')!='baseline-v1':
            require(first['wp2_per_state_monitoring'][0]['time']==0.,'pressure recovery must retain its fresh seed')
            seed=dict(time=0.,state=first['accepted_states'][0],fixed_reference=first['fixed_reference'],fixed_defect=first['fixed_defect'])
        else:seed=load(seed_file)
        require(first['wp2_per_state_monitoring'][0]['time']==seed['time']==0.,'recovery initial physical time differs')
        require(state_hash(first['accepted_states'][0])==state_hash(seed['state']),'recovery initial coefficients differ from saved seed')
        require(all(state_hash(first[k])==state_hash(seed[k]) for k in ['fixed_reference','fixed_defect']),'recovery seed reference/defect differs')
    final_start=p['actual_time']-.3;acc0=None;acc1=None;times=[];endpoint=None;start_report=None;initial_samples=None;cold_samples=None
    for file in paths(root,model,phase,p):
        d=load(file);checks=verify_report(d,US,K0[model],allow_recovery=recovery)
        require(len(d['accepted_states'])==len(d['wp2_per_state_monitoring']),'monitoring/state counts differ')
        if start_report is None:
            start_report=d;cold=rows.tuple(d['fixed_reference'],d)
            cold_samples=dict(d['wp2_common_samples'],values=cold['values'])
        require(state_hash(d['fixed_reference'])==state_hash(start_report['fixed_reference']) and
                state_hash(d['fixed_defect'])==state_hash(start_report['fixed_defect']),'reference/defect changed')
        sources.append(dict(path=str(file),sha256=sha(file),checkpoint=d['checkpoint'],checkpoint_sha256=sha(Path(d['checkpoint'])/'State.bin')))
        replay.append(dict(report=str(file),**rows.check_native(d)))
        counts[str(file)]=dict(attempts=d['implicit_callback_count'],accepted_executions=len(d['accepted_states'])-1)
        for index,(state,monitor,check) in enumerate(zip(d['accepted_states'],d['wp2_per_state_monitoring'],checks)):
            t=monitor['time'];h=state_hash(state)
            if abs(t-last_time)<2e-14:
                require(h==last_hash,'restart boundary changed accepted coefficients');continue
            require(t>last_time,'trajectory time order differs')
            v=rows.tuple(state,d);full=rows.full_errors(state,d,v)
            denser=rows.volume_audit(d,v,order=independent_order,allow_accuracy_failure=recovery)
            horizontal=homogeneity(dict(d['wp2_common_samples'],values=v['values']),cold_samples,US,K0[model])
            tensor_horizontal=tensor_homogeneity(rows,v,cold)
            if initial_samples is None:initial_samples=dict(d['wp2_common_samples'],values=v['values'].tolist())
            cov=min(x['covariance_minimum'] for x in full['regions']);require(cov>=1e-4*K0[model],'independent full-tensor common-point covariance failed')
            accuracy=all(max(x['mean_stress_tensor_relative_RMS'],x['reynolds_tensor_relative_RMS'])<=.01 for x in full['regions'])
            if not recovery:require(accuracy and horizontal['passed'] and tensor_horizontal['passed'],'full selected stress/homogeneity error budget failed')
            energy=dict(stored=float(np.dot(rows.weight,.5*(v['values'][:,0]**2+v['values'][:,1]**2)+v['values'][:,2])))
            if previous_values is not None:
                alpha=d['pimple_attempts'][index-1]['alpha'];du=v['values'][:,:2]-previous_values[:,:2];dk=v['values'][:,2]-previous_values[:,2]
                energy['physical_storage_rate']=float(alpha*np.dot(rows.weight,.5*np.sum((v['values'][:,:2]+previous_values[:,:2])*du,axis=1)+dk))
                energy['backward_euler_work']=float(alpha*np.dot(rows.weight,np.sum(v['values'][:,:2]*du,axis=1)+dk))
                energy['time_discretization_work']=float(.5*alpha*np.dot(rows.weight,np.sum(du**2,axis=1)))
                native_energy=d['state_audits'][index]['energy'];error=abs(energy['backward_euler_work']-native_energy[1])
                energy['native_mass_storage_absolute_error']=error;require(error<2e-10,'physical quadrature/native mass storage differs')
            history.append(dict(actual_time=t,state_sha256=h,physical_checks=check,full_tensor_checks=full,denser_volume_audit=denser,horizontal_homogeneity=horizontal,
                                full_tensor_horizontal_homogeneity=tensor_horizontal,full_tensor_accuracy_passed=accuracy,
                                energy=energy,native_energy_audit=d['state_audits'][index].get('energy'),
                                boundary_work=monitor['physical_checks']['boundary_status']['equilibrium_boundary_work'],
                                wall_mean_shear=monitor['physical_checks']['wall_mean_shear']))
            if t>=final_start-2e-14:
                if acc0 is None:acc0=np.zeros_like(v['values']);acc1=np.zeros_like(v['values'])
                acc0+=v['values'];acc1+=t*v['values'];times.append(t)
            last_time,last_hash=t,h;previous_values=v['values']
        endpoint=d
    require(abs(last_time-p['actual_time'])<2e-14,'terminal native time differs')
    require(all(o['physical_budgets_passed'] for o in p['observations']) or recovery,'per-state preservation failed')
    final_pass=history[-1]['physical_checks']['passed'] and history[-1]['denser_volume_audit']['passed'] and history[-1]['full_tensor_accuracy_passed'] and history[-1]['horizontal_homogeneity']['passed'] and history[-1]['full_tensor_horizontal_homogeneity']['passed']
    # Linear trend is measured on executed states only. The transit-scaled
    # slope is an accuracy-scale diagnostic, not a prediction of future flow.
    tt=np.asarray(times);den=float(np.sum((tt-np.mean(tt))**2));require(den>0,'unresolved final observation window')
    slope=(acc1-np.mean(tt)*acc0)/den
    sample=endpoint['wp2_common_samples'];zero=dict(sample,values=sample['reference_values']);moving=dict(sample,values=(np.asarray(sample['reference_values'])+TREQ*slope).tolist())
    trend=quarter_budget_differences(moving,zero,US,K0[model],.0065)
    result=dict(version='WP2AllStatePhysicalAssessmentV1',model=model,phase=phase,actual_time=last_time,Treq=TREQ,transit_units=last_time/TREQ,
        analysis_source_identity=analysis_identity(),independent_order=independent_order,
        unique_retained_states=len(history),sources=sources,counts=counts,native_tuple_replays=replay,history=history,
        final_physical_and_full_tensor_passed=final_pass,final_state_sha256=last_hash,final_window=dict(start=min(times),end=max(times),
            duration=max(times)-min(times),required=TREQ/4,actual_states=len(times),transit_scaled_linear_trend=trend,
            interpretation='least-squares slope of actual values on the final executed window; no infinite-time or future-state claim'),
        final_budgets=endpoint['wp2_final_turbulence_budget'],seconds=time.monotonic()-started,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        context_sha256=sha(rows.path),progress_sha256=progress_hash,governing_solves_in_this_analysis=0)
    require(sha(progress_path)==progress_hash,'trajectory changed during all-state assessment')
    require(result['final_window']['duration']>=TREQ/4-2e-14,'final observation window too short')
    atomic(output,result);return result,endpoint,cold_samples


def tensor_differences(left,right,rows_left,rows_right):
    a=rows_left.tuple(left['state'],left);b=rows_right.tuple(right['state'],right)
    require(np.array_equal(rows_left.physical,rows_right.physical) and np.array_equal(rows_left.weight,rows_right.weight),'tensor comparison geometry differs')
    out={}
    for region in [0,1]:
        mask=rows_left.inner if region==0 else ~rows_left.inner;w=rows_left.weight[mask];den=US**2*np.sqrt(sum(w))
        out[str(region)]={name:float(np.sqrt(np.sum(w[:,None,None]*(a[name][mask]-b[name][mask])**2))/den) for name in ['R','T']}
    return dict(passed=all(v<=.0025 for row in out.values() for v in row.values()),regional_full_tensor_RMS=out,limit=.0025)


def native_states(root,model,phase,progress_path=None):
    progress=load(progress_path or root/f'preservation-{phase}-{model}/progress.json');states={};context=None
    require(progress['model']==model,'native-state trajectory model differs')
    for file in paths(root,model,phase,progress):
        d=load(file)
        if context is None:context={k:d[k] for k in ['parameters','model_specification','wp2_common_samples']}
        for s,monitor in zip(d['accepted_states'],d['wp2_per_state_monitoring']):
            t=monitor['time'];key=round(t,12)
            if key in states:require(state_hash(states[key][1])==state_hash(s),'matched-time restart state changed')
            else:states[key]=(t,np.asarray(s,dtype=np.float64))
    return states,context


def all_common_time_differences(root,model,left,right,rows_left,rows_right):
    a,ca=native_states(root,model,left);b,cb=native_states(root,model,right);results=[]
    for key in sorted(set(a)&set(b)):
        ta,sa=a[key];tb,sb=b[key];require(abs(ta-tb)<2e-14,'common-time lookup exceeds native time allowance')
        va=rows_left.tuple(sa,ca);vb=rows_right.tuple(sb,cb)
        aa=dict(ca['wp2_common_samples'],values=va['values']);bb=dict(cb['wp2_common_samples'],values=vb['values'])
        check=quarter_budget_differences(aa,bb,US,K0[model],.0065);tensor={}
        for region in [0,1]:
            mask=rows_left.inner if region==0 else ~rows_left.inner;w=rows_left.weight[mask];den=US**2*np.sqrt(sum(w))
            tensor[str(region)]={name:float(np.sqrt(np.sum(w[:,None,None]*(va[name][mask]-vb[name][mask])**2))/den) for name in ['R','T']}
        results.append(dict(left_actual_time=ta,right_actual_time=tb,principal=check,full_tensor_RMS=tensor,
                            passed=check['passed'] and all(v<=.0025 for row in tensor.values() for v in row.values())))
    require(results and results[-1]['left_actual_time']>=TREQ,'matched-time comparison does not cover required horizon')
    return dict(passed=all(x['passed'] for x in results),native_common_times=len(results),results=results,
                scope='all executed common native times, including h0; no time interpolation or state relabeling; absolute pressures')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--phases',nargs='+',choices=['coarse','fine','confirmation','recovery'],default=['coarse','fine','confirmation','recovery'])
    p.add_argument('--reuse',nargs='*',type=Path,default=[])
    a=p.parse_args();a.output.mkdir(exist_ok=False);contexts={}
    for name in ['selected','confirmation']:contexts[name]=Rows(a.root/f'geometry-rows-{name}/result.json')
    summary=dict(version='WP2ApproachAnalysisV1',models={})
    for model in ['CLS','standard']:
        cases={};ends={};starts={}
        for phase in a.phases:
            rows=contexts['confirmation' if phase=='confirmation' else 'selected']
            reused=False
            for directory in a.reuse:
                prior=directory/f'{model}-{phase}-all-states.json'
                if not prior.exists():continue
                old=load(prior);progress=a.root/f'preservation-{phase}-{model}/progress.json'
                if old.get('analysis_source_identity')!=analysis_identity() or old['context_sha256']!=sha(rows.path) or old['progress_sha256']!=sha(progress):continue
                require(all(sha(Path(s['path']))==s['sha256'] and sha(Path(s['checkpoint'])/'State.bin')==s['checkpoint_sha256'] for s in old['sources']),'reused physical evidence changed')
                shutil.copyfile(prior,a.output/prior.name);cases[phase]=old;ends[phase]=load(Path(old['sources'][-1]['path']))
                starts[phase]=dict(ends[phase]['wp2_common_samples'],values=rows.tuple(ends[phase]['fixed_reference'],ends[phase])['values'])
                reused=True;break
            if not reused:cases[phase],ends[phase],starts[phase]=assess(a.root,model,phase,rows,a.output/f'{model}-{phase}-all-states.json')
        comparisons={}
        for label,x,y in [('timestep','coarse','fine'),('mesh','coarse','confirmation'),('final_window','fine','fine')]:
            if x not in cases or y not in cases:continue
            if label=='final_window':
                obs=load(a.root/f'preservation-{x}-{model}/progress.json')['observations'];left=load(Path(next(z['report'] for z in obs if abs(z['time']-.9)<2e-14)));pairs=[(.9,left,ends[x])]
            else:
                xp=load(a.root/f'preservation-{x}-{model}/progress.json')['observations'];yp=load(a.root/f'preservation-{y}-{model}/progress.json')['observations'];pairs=[]
                for t in [.3,.6,.9,1.2]:
                    l=load(Path(next(z['report'] for z in xp if abs(z['time']-t)<2e-14)));r=load(Path(next(z['report'] for z in yp if abs(z['time']-t)<2e-14)));pairs.append((t,l,r))
            comparisons[label]=[]
            for t,l,r in pairs:
                numerical=quarter_budget_differences(l['wp2_common_samples'],r['wp2_common_samples'],US,K0[model],.0065)
                tensor=tensor_differences(l,r,contexts['selected'],contexts['confirmation' if y=='confirmation' else 'selected'])
                quad=compare_orders(l,r);wall=abs(l['wp2_per_state_monitoring'][-1]['physical_checks']['wall_mean_shear']-r['wp2_per_state_monitoring'][-1]['physical_checks']['wall_mean_shear'])/US**2
                comparisons[label].append(dict(requested_observation=t,left_actual_time=l['time'],right_actual_time=r['time'],principal=numerical,full_tensors=tensor,
                    wall_traction_difference_over_ustar2=wall,wall_traction_quarter_budget_passed=wall<=.0025,quadrature_GL8_GL12=quad))
        homo={}
        for phase in a.phases:
            if phase=='recovery':continue
            # The independently checked radial seed is homogeneous; it is the
            # unchanged native inlet trace, not a regenerated reference field.
            homo[phase]=homogeneity(ends[phase]['wp2_common_samples'],starts[phase],US,K0[model])
        recovery={}
        if 'recovery' in cases and 'fine' in cases:
            for t in [.3,.6,.9,1.2]:
                rp=load(a.root/f'preservation-recovery-{model}/progress.json')['observations'];fp=load(a.root/f'preservation-fine-{model}/progress.json')['observations']
                r=load(Path(next(z['report'] for z in rp if abs(z['time']-t)<2e-14)));f=load(Path(next(z['report'] for z in fp if abs(z['time']-t)<2e-14)))
                recovery[str(t)]=quarter_budget_differences(r['wp2_common_samples'],f['wp2_common_samples'],US,K0[model],.0065)
        all_times={}
        if 'coarse' in cases and 'fine' in cases:
            all_times['timestep']=all_common_time_differences(a.root,model,'coarse','fine',contexts['selected'],contexts['selected'])
        if 'coarse' in cases and 'confirmation' in cases:
            all_times['mesh']=all_common_time_differences(a.root,model,'coarse','confirmation',contexts['selected'],contexts['confirmation'])
        summary['models'][model]=dict(cases={k:dict(path=str(a.output/f'{model}-{k}-all-states.json'),sha256=sha(a.output/f'{model}-{k}-all-states.json'),
              final_passed=v['final_physical_and_full_tensor_passed'],actual_time=v['actual_time']) for k,v in cases.items()},comparisons=comparisons,all_common_native_times=all_times,
              horizontal_homogeneity=homo,recovery_control_differences=recovery)
    atomic(a.output/'summary.json',summary)


if __name__=='__main__':main()
