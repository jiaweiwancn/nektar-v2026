#!/usr/bin/env python3
"""Linux-only M7 native accepted-state audits and bounded perturbations.

Physical failures are separate from verification failures. No stabilization,
floor, filtering, or change of constitutive coefficients is made here.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import itertools
import json
import math
import platform
import time
from pathlib import Path
import numpy as np
from cls_epawf_m7_cases import perturbation_session
from validate_cls_epawf_m4 import Runner, spaces, state_difference, native_history_difference
from validate_cls_epawf_m2 import Profiles

CHECKS = 0

def require(value, message):
    global CHECKS
    CHECKS += 1
    if not value:raise AssertionError(message)


def audits(raw, name):
    require(raw.get('diagnostics_schema') == 1 and not raw['failure'], name+': diagnostics missing/failed')
    require(len(raw['state_audits']) == len(raw['accepted_states']) == len(raw['stage_audits'])+1,
            name+': every accepted initial/stage state audited')
    summary = {'time':raw['time'], 'eta':raw['eta'], 'chi':raw['chi'],
               'states':len(raw['state_audits']), 'min_k':math.inf, 'min_outer_epsilon':math.inf,
               'min_covariance':math.inf, 'min_yx_tangent':math.inf,
               'maximum_stage_residual':max([r[4] for r in raw['stage_audits']],default=0.),
               'maximum_energy_closure':0., 'maximum_free_energy_balance':0.,
               'maximum_strong_divergence':0., 'maximum_normal_traction_jump':0.}
    for index,(state,audit) in enumerate(zip(raw['accepted_states'],raw['state_audits'])):
        require(audit['admissible'] == [1], name+': inadmissible accepted state')
        for key in ('source_volume_inner','source_volume_outer','physical_sample_K','mass_volume',
                    'interface_inner','interface_outer','reconstruction_native','reconstruction_dense'):
            row = audit[key]
            require(row[0]>0 and row[1]>=0 and row[5]==0 and np.all(np.isfinite(row)), name+': '+key)
            summary['min_k'] = min(summary['min_k'], row[1])
        if raw['upstream_transition']:
            require(audit['physical_sample_outer'][0]>0, name+': missing transition physical samplers')
        for key in ('source_volume_outer','interface_outer'):
            require(audit[key][2]>0, name+': outer epsilon positivity')
            summary['min_outer_epsilon'] = min(summary['min_outer_epsilon'],audit[key][2])
            summary['min_covariance'] = min(summary['min_covariance'],audit[key][3])
            summary['min_yx_tangent'] = min(summary['min_yx_tangent'],audit[key][4])
        e = audit['energy'];scale=max(1.,max(abs(v) for v in e[:7]))
        require(abs(e[7])/scale<1e-10, name+': true weak/direct energy accounting')
        summary['maximum_energy_closure']=max(summary['maximum_energy_closure'],abs(e[7])/scale)
        if index:
            require(max(audit['coupled_norms'])<1e-9, name+': unconverged coupled stage')
            require(abs(e[11])/scale<1e-8, name+': free-row energy balance including reaction and anchor')
            summary['maximum_free_energy_balance']=max(summary['maximum_free_energy_balance'],abs(e[11])/scale)
            dim=raw['dimension'];b=np.asarray(audit['momentum_k_load_and_boundary_reaction'])
            defect=raw['eta']*np.asarray(raw['defect_integrals'][:dim+1])
            require(np.max(np.abs(b[:dim+1]-b[dim+1:]-defect))<1e-8*scale,
                    name+': momentum/k balance including essential reactions')
        wall=audit['wall_sources']
        require(abs(wall[2])<1e-10*max(1.,abs(wall[0]),abs(wall[1])), name+': zero integral redistribution')
        div=audit['strong_continuity'];require(abs(div[2]+div[3]-div[4])<1e-9, name+': Gauss mass balance')
        summary['maximum_strong_divergence']=max(summary['maximum_strong_divergence'],div[0],div[1])
        summary['maximum_normal_traction_jump']=max(summary['maximum_normal_traction_jump'],audit['interface'][4])
    # Independently reconstruct every accepted native polynomial at all dense
    # audit locations. This uses the Jacobi basis oracle, not production rows.
    fields=spaces(raw);offsets=raw['offsets'];dim=raw['dimension'];independent_min=math.inf
    for state,audit in zip(raw['accepted_states'],raw['state_audits']):
        minimum=math.inf
        for el in fields[dim].elements:
            n=2*el['degree']+3
            axes=fields[dim].axes
            points=[]
            for index in itertools.product(range(n),repeat=dim):
                point=[0.,0.,0.]
                for a,j in zip(axes,index):point[a]=el['lower'][a]+(el['upper'][a]-el['lower'][a])*j/(n-1)
                points.append(point)
            ids,phi,_=fields[dim].rows(el,points)
            minimum=min(minimum,float(np.min(phi@np.asarray(state)[offsets[dim]+ids])))
        require(abs(minimum-audit['reconstruction_dense'][1])<1e-10, name+': independent dense polynomial minimum')
        independent_min=min(independent_min,minimum)
    summary['independent_dense_min_k']=independent_min
    # Trapezoidal diagnostic integrals: these are accounting defects sampled
    # in physical time, not a new definition of the IMEX stored-energy change.
    times=np.asarray(raw['audit_times']);rows=np.asarray([a['energy'] for a in raw['state_audits']])
    summary['cumulative_energy_accounting_defect']=float(np.trapz(rows[:,7],times))
    summary['final_energy']=rows[-1].tolist()
    summary['final_wall_sources']=raw['state_audits'][-1]['wall_sources']
    summary['final_interface']=raw['state_audits'][-1]['interface']
    if not raw['mapped']:
        energies=[]
        reference=np.asarray(raw['fixed_reference'])
        for state in (raw['accepted_states'][0],raw['accepted_states'][-1]):
            delta=np.asarray(state)-reference;energy=0.
            for a in range(dim):
                v=delta[offsets[a]:offsets[a+1]];energy+=.5*float(v@fields[a].mass()@v)
            energies.append(energy)
        summary['velocity_perturbation_energy']=energies
        summary['velocity_energy_gain']=energies[1]/energies[0] if energies[0]>1e-24 else None
    summary['turnover_time']=raw['parameters'][1]/raw['ustar']
    field=fields[0];z=.5*(raw['parameters'][1]+field.upper[2]);p=[.5*(field.lower[0]+field.upper[0]),0.,z]
    if dim==3:p[1]=.5*(field.lower[1]+field.upper[1])
    u=abs(field.evaluate(raw['fixed_reference'],offsets[0],p)[0])
    summary['flow_through_time']=(field.upper[0]-field.lower[0])/u if u>0 else None
    duration=raw['time']-raw['audit_times'][0]
    summary.update(duration=duration,turnover_units=duration/summary['turnover_time'],
                   flow_through_units=duration/summary['flow_through_time'] if u>0 else None)
    return summary


def campaign(run, mode):
    results={};raws={}
    def case(name, **kwargs):
        raw=run.evolution(name,generator=perturbation_session,**kwargs)
        raws[name]=raw;results[name]=audits(raw,name)
        return raw
    if mode=='pilot':
        case('joint',steps=1)
        return results
    # Freeze before execution; identity/control durations are explicitly
    # different from short nonlinear failure-diagnosis trajectories.
    matrix=[('velocity',dict(kind='velocity')),('k',dict(kind='k')),('epsilon',dict(kind='epsilon')),
            ('joint',dict(kind='joint')),('interface',dict(kind='interface')),
            ('joint-small',dict(kind='joint',amplitude=.0001)),
            ('joint-large',dict(kind='joint',amplitude=.01)),
            ('joint-wave2',dict(kind='joint',wave=2)),
            ('joint-dt2',dict(kind='joint',dt=.00125,steps=8)),
            ('joint-dt4',dict(kind='joint',dt=.000625,steps=16)),
            ('joint-p4',dict(kind='joint',degree=4)),
            ('joint-h2',dict(kind='joint',refine=2)),
            ('joint-unanchored',dict(kind='joint',eta=0)),
            ('control-unanchored',dict(kind='none',eta=0)),
            ('joint-3d',dict(kind='joint',dimension=3,steps=1)),
            ('transition-crossing',dict(kind='transition',transition=True,steps=2)),
            ('transition-control',dict(kind='none',transition=True,steps=2))]
    (run.output/'frozen-matrix.json').write_text(json.dumps({'schema':1,'matrix':matrix,
        'default':dict(amplitude=.001,dt=.0025,steps=4,degree=3,refine=1,eta=1,chi=1),
        'acceptance':{'coupled_residual':1e-9,'isolated_identity':1e-10,
            'MPI_physical_norm':1e-8,'no_everywhere_polynomial_positivity_claim':True,
            'robust_model_gate':'failed_as_specified; not a numerical pass criterion'}},indent=2)+'\n')
    for name,kwargs in matrix:case(name,**kwargs)
    # Actual native stationary benchmark lasting five flow-through times.
    control=case('five-flow-control-initial',kind='none',steps=0)
    T=results['five-flow-control-initial']['flow_through_time']
    long=case('five-flow-control',kind='none',steps=25,dt=T/5)
    require(max(long['step_audits'][-1][4:])<1e-9,'stationary five-flow benchmark drift')
    require(abs(results['five-flow-control']['flow_through_units']-5)<1e-12,'five flow-through duration')
    # Independent dt refinement measures the same finite-resolution nonlinear
    # trajectory. Spatial runs are reported separately; no conflated order.
    errors=[]
    for a,b in [('joint','joint-dt2'),('joint-dt2','joint-dt4')]:
        x,y=raws[a],raws[b];f=spaces(x);delta=np.asarray(x['state'])-np.asarray(y['state'])
        errors.append(sum(float(delta[x['offsets'][j]:x['offsets'][j+1]]@f[j].mass()@
                                delta[x['offsets'][j]:x['offsets'][j+1]]) for j in range(x['dimension']+2))**.5)
    require(errors[1]<errors[0], 'native temporal refinement must reduce trajectory difference')
    results['native_dt_refinement']={'differences':errors,'ratio':errors[0]/errors[1]}
    # Same-source MPI audit and complete physical field parity.
    parallel=run.evolution('joint-mpi2',generator=perturbation_session,ranks=2)
    results['joint-mpi2']=audits(parallel,'joint-mpi2')
    state_difference(raws['joint'],parallel,'M7 MPI all-field parity')
    for a,b in zip(raws['joint']['state_audits'],parallel['state_audits']):
        for name in a:require(np.allclose(a[name],b[name],rtol=1e-8,atol=1e-8),'MPI diagnostic parity '+name)
    first=case('joint-checkpoint',kind='joint',steps=2,checkpoint=run.output/'perturbed-history')
    resumed=case('joint-restart',kind='joint',steps=2,restart=first['checkpoint'])
    state_difference(raws['joint'],resumed,'M7 perturbed all-field restart')
    native_history_difference(raws['joint'],resumed,'M7 perturbed native histories')
    # h/p tests fix H,h,rK,rm. H study is a separately labelled model study.
    height_study=[]
    for H in (1.5,2.,2.5):
        key='height-'+str(H);r=case(key,kind='none',steps=0,height=H)
        nu,_,h,rK,rm,K,q,kappa=r['parameters']
        profile=Profiles(dict(geometry=[nu,H,h,rK,rm],k=K,q=q,kappa=kappa,chi=1,metric=[0,0]))
        require(profile.zv<0,'height reference uses fully turbulent auxiliary branch')
        zr=float(profile.d/profile.a)
        formula=H*((H+zr)/zr)**(-kappa/(.2295*math.sqrt(q)))
        UH=spaces(r)[0].evaluate(r['fixed_reference'],r['offsets'][0],[4.,0.,H])[0]
        native_z0=H*math.exp(-kappa*UH/r['ustar'])
        require(abs(native_z0-formula)<1e-10,'M7.matching_height_dependence continuum/native matching')
        height_study.append(dict(H=H,matching_velocity=UH,roughness_intercept=native_z0,formula=formula))
    require(abs(height_study[0]['roughness_intercept']-height_study[-1]['roughness_intercept'])>1e-6,
            'height study must retain model dependence rather than retuning roughness')
    results['matching_height_study']=height_study
    # A polynomial positive at native GLL/source GL8 abscissae but negative at
    # x=L/2 is rejected at the added reconstruction set, before any timestep.
    directory=run.output/'node-alias-reject';directory.mkdir()
    xml=perturbation_session(run.template,directory/'case.xml',kind='node_alias',chi=0,steps=0)
    run.command('node-alias-reject','IncNavierStokesSolver',[xml],expect_failure=True)
    rejected=json.loads((directory/'case.evolution.json').read_text())
    bad=rejected['state_audits'][-1]
    require(bad['reconstruction_native'][1]>0 and bad['source_volume_inner'][1]>0,
            'alias control must pass native reconstruction/source positivity subsets')
    require(bad['reconstruction_dense'][1]<0 and not rejected['accepted_states'],
            'dense audit must reject the between-node negative polynomial')
    results['node-alias-reject']={'native_min':bad['reconstruction_native'][1],
        'source_min':bad['source_volume_inner'][1],'dense_min':bad['reconstruction_dense'][1],
        'failure':rejected['failure'],'accepted_states':0}
    directory=run.output/'invalid-sampler-reject';directory.mkdir()
    xml=perturbation_session(run.template,directory/'case.xml',kind='invalid_sampler',chi=0,steps=0)
    run.command('invalid-sampler-reject','IncNavierStokesSolver',[xml],expect_failure=True)
    rejected=json.loads((directory/'case.evolution.json').read_text())
    require(not rejected['initialization_complete'] and not rejected['accepted_states']
            and 'sampled K' in rejected['failure'],'invalid initial sampler must preserve original failure before history exists')
    results['invalid-sampler-reject']={'failure':rejected['failure'],'accepted_states':0}
    return results


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--output',type=Path);parser.add_argument('--mode',choices=['pilot','all'],default='all')
    args=parser.parse_args();require(platform.system()=='Linux','scientific validation belongs on Linux')
    if args.output is None:args.output=Path.cwd()/('cls-epawf-m7-'+str(time.time_ns()))
    run=Runner(args.build,args.output)
    try:
        result=campaign(run,args.mode)
        receipt={'schema':1,'status':'passed','checks':CHECKS,'mode':args.mode,'results':result,
                 'model_stability':'failed_as_specified','positivity_invariance':'not_guaranteed',
                 'commands':run.commands,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (run.output/'receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    except Exception as error:
        (run.output/'failure.json').write_text(json.dumps({'status':'failed','checks':CHECKS,'error':str(error)},indent=2)+'\n')
        raise
    print(f'M7 native {args.mode}: {CHECKS} checks passed',flush=True)

if __name__=='__main__':main()
