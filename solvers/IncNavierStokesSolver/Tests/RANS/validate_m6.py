#!/usr/bin/env python3
"""Reproducible Linux neutral-ABL study; physical failures remain measured results.

The optional matrix-free backend evaluates the current weak operator. Results
are cached only when fixture, executable and audit hashes are identical.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import threading
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
THREADS={key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
os.environ.update(THREADS)
from generate_m6 import generate_abl
from generate_m4 import generate_cls, info, parameter
from generate_m3 import generate


def require(valid,message):
    if not valid:raise RuntimeError(message)


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rows(path):
    if not path.exists():return []
    return [{k:float(v) if k!='field' else v for k,v in r.items()} for r in csv.DictReader(path.open())]


def recorded_failure(directory, name, binaries, message):
    """Reduce an observed failed run without claiming a completed preservation test."""
    session=directory/(name+'.xml');log=directory/'solver.log'
    stop=directory/'engineering-stop.json'
    require(log.exists() and ('failed to converge' in log.read_text() or stop.exists()),
            'Missing native failure log or recorded engineering stop')
    parameters={}
    for node in ET.parse(session).getroot().findall('CONDITIONS/PARAMETERS/P'):
        key,value=node.text.split('=',1)
        try:parameters[key.strip()]=float(value)
        except ValueError:pass
    stages=read_rows(directory/(name+'.rans.csv'));errors=read_rows(directory/(name+'.rans-errors.csv'))
    require(stages and errors,'Failed run has no completed steps')
    fields={row['field'] for row in errors};last=max(row['time'] for row in errors)
    peak_l2={f:max(row['normalized_l2'] for row in errors if row['field']==f) for f in fields}
    peak_linf={f:max(row['pointwise_relative_linf'] for row in errors if row['field']==f) for f in fields}
    receipt=json.loads(stop.read_text()) if stop.exists() else None
    if receipt:
        require(receipt['binaries']==binaries and receipt['value']>=.001 and
                receipt['fixture_sha256']==digest(session),'Invalid engineering-stop receipt')
    return dict(status='profile_target_failure' if receipt else 'numerical_failure',message=message,
        fingerprint=dict(fixture=digest(session),binaries=binaries,ranks=4),
        planned_flow_times=parameters['NumSteps']*parameters['TimeStep']*parameters['Uref']/parameters['L'],
        completed_flow_times=last*parameters['Uref']/parameters['L'],completed_time=last,
        peak_l2=peak_l2,peak_linf=peak_linf,engineering_target_met=False,stages=stages,
        final_errors={row['field']:row for row in errors if row['time']==last},
        artifacts={p.name:digest(p) for p in directory.iterdir() if p.is_file() and p.suffix in ('.xml','.log','.csv')},
        engineering_stop=receipt)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--suite',choices=('pilots','study','restart','columns','robust','standard','all'),default='all')
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--study-output',type=Path)
    parser.add_argument('--extended-perturbation',action='store_true',help='Also repeat the five-flow 3D perturbation at dt=0.1')
    parser.add_argument('--record-existing-failure',action='store_true',help='Reduce an already stopped robust trial; never rerun its solver')
    args=parser.parse_args();require(sys.platform.startswith('linux'),'Use the authoritative Linux environment')
    require(not args.record_existing_failure or args.suite=='robust','Existing-failure reduction requires the robust suite')
    require(1<=args.jobs<=4,'At most four jobs with at most four ranks each (16 CPUs)')
    build=args.build.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    # Multiple drivers may share a cache, but only one may publish its files.
    cache_lock=(out/'.validation.lock').open('a')
    fcntl.flock(cache_lock,fcntl.LOCK_EX)
    solver=build/'IncNavierStokesSolver';audit=build/'RANSFieldAudit'
    binaries={p.name:digest(p) for p in (solver,audit)};records={};commands=[];lock=threading.Lock();verification=False

    def execute(cmd,directory,label,negative=None,timeout=14400):
        cmd=list(map(str,cmd))
        with lock:commands.append(dict(command=cmd,cwd=str(directory),timeout=timeout))
        result=subprocess.run(cmd,cwd=directory,env=os.environ.copy(),text=True,capture_output=True,timeout=timeout)
        log=result.stdout+result.stderr;(directory/(label+'.log')).write_text(log)
        require((result.returncode!=0 and negative in log) if negative else result.returncode==0,
                f'{label} failed: {directory/(label+".log")}')
        return log

    def run(name,ranks=1,modify=None,generator=generate_abl,negative=None,profiles=True,**options):
        directory=out/name;directory.mkdir(exist_ok=True);session=directory/(name+'.xml')
        spec=generator(session,**options)
        tree=ET.parse(session);root=tree.getroot()
        if generator==generate_abl:
            # Looser than verification cases, still far below the 0.1% target.
            parameter(root,'RANSReferenceMonitor',1)
            parameter(root,'RANSStageAudit',1 if name.startswith('backend_') else 0)
        if verification:
            parameter(root,'IterativeSolverTolerance',1e-15)
            parameter(root,'RANSMaxIterations',10000)
        if modify:modify(root)
        tree.write(session,encoding='utf-8',xml_declaration=True)
        fingerprint=dict(fixture=digest(session),binaries=binaries,ranks=ranks,negative=negative)
        cache=directory/'case.json'
        if cache.exists() and json.loads(cache.read_text()).get('fingerprint')==fingerprint:
            record=json.loads(cache.read_text());print('Reusing '+name,flush=True)
        else:
            print('Running '+name,flush=True)
            cmd=(['mpirun','--bind-to','none','--map-by','slot','-np',str(ranks)] if ranks>1 else [])+[solver,session.name]
            execute(cmd,directory,'solver',negative)
            record=dict(fingerprint=fingerprint,spec=spec,options=options,ranks=ranks)
            if negative:record['expected_rejection']=negative
            else:
                command=[audit,session.name,name+'.fld',name+'.fld']
                if profiles and generator==generate_abl:command+=['profiles.csv']
                raw=execute(command,directory,'audit')
                record['audit']=json.loads(raw[raw.index('{"time"'):])
                for key,suffix in (('stages','rans.csv'),('wall','awf.csv'),('errors','rans-errors.csv'),('residual','rans-residual.csv')):
                    record[key]=read_rows(directory/(name+'.'+suffix))
                require(len(record['stages'])==options.get('steps',spec['steps'] if spec else 4),'Incomplete timesteps')
                for row in record['stages']:
                    require(all(math.isfinite(v) for v in row.values()),'Nonfinite stage')
                    require(row['min_k']>0 and row['min_eps']>0 and row['guard_activations']==0,'Positivity guard')
                for row in record['wall']:
                    require(all(math.isfinite(v) for v in row.values()),'Nonfinite wall diagnostic')
                    require(row['guard_activations']==0 and row['min_kP']>0 and row['min_epsP']>0,'Strip guard')
                    require(row['cutoff_fallback_points']==row['regime_fallback_points']==row['range_fallback_points']==0,'Unexpected fallback')
                if name.startswith('backend_'):
                    residual=record['residual']
                    require(residual and all(math.isfinite(r['free_coefficient_l2']) for r in residual),'Missing/nonfinite native residual audit')
                    initial=max(r['free_coefficient_l2'] for r in residual if r['phase']==0)
                    final=max(r['free_coefficient_l2'] for r in residual if r['phase']==1)
                    require(final<max(1e-6,.01*initial),'Current weak stage residual did not converge')
                if record['errors']:
                    final=max(row['time'] for row in record['errors'])
                    record['final_errors']={row['field']:row for row in record['errors'] if row['time']==final}
                    record['peak_l2']={field:max(row['normalized_l2'] for row in record['errors'] if row['field']==field) for field in record['final_errors']}
                    primary=('u','k','eps' if options.get('model')=='standard' else 'eps_tilde','rans_nut','rans_R_xx','rans_R_yy','rans_R_zz','rans_R_xz')
                    record['peak_linf']={field:max(row['pointwise_relative_linf'] for row in record['errors'] if row['field']==field) for field in record['final_errors']}
                    record['engineering_target_met']=all(record['peak_l2'][f]<.001 and record['peak_linf'][f]<.001 for f in primary)
                record['artifacts']={p.name:digest(p) for p in directory.iterdir() if p.is_file() and p.suffix in ('.xml','.log','.csv','.his')}
            cache.write_text(json.dumps(record,indent=2)+'\n')
        with lock:
            records[name]=record
            (out/('progress-'+args.suite+'.json')).write_text(json.dumps(records,indent=2)+'\n')
        if 'final_errors' in record:print(name+' '+str({k:round(v['normalized_l2'],8) for k,v in record['final_errors'].items() if k in ('u','k','eps','eps_tilde','rans_nut')}),flush=True)
        return record

    def compare(a,b,label,tolerance=2e-8):
        raw=execute([audit,a[0],a[1],b],out,label)
        result=json.loads(raw[raw.index('{"time"'):])
        require(max(v['difference_linf'] for v in result['fields'].values())<tolerance,label+' field comparison failed')
        return result

    if args.suite in ('pilots','all'):
        for dim in (2,3):
            for backend in ('static','matrixfree'):
                run(f'backend_d{dim}_{backend}',ranks=1 if dim==2 else 4,dimension=dim,degree=4,layers=6,steps=2,dt=1,backend=backend,monitor=True)
            a=f'backend_d{dim}_static';b=f'backend_d{dim}_matrixfree'
            compare((out/a/(a+'.xml'),out/a/(a+'.fld')),out/b/(b+'.fld'),f'backend_compare_d{dim}',1e-6)
        run('outflow_pilot',dimension=2,layers=6,steps=4,dt=.5,backend='matrixfree',boundary='outflow',monitor=True)
        for model in ('standard','cls_switched'):
            run('tuned_'+model,model=model,layers=12,steps=10,dt=2,backend='matrixfree',monitor=True)

    if args.suite in ('study','all'):
        base=dict(degree=4,layers=12,nx=1,dt=.5,flow_times=5,backend='matrixfree',monitor=True,tolerance=1e-8)
        jobs=[]
        for dim in (2,3):
            for model in ('standard','cls_original','cls_fixed','cls_switched'):
                jobs.append((f'five_d{dim}_{model}',dict(base,model=model,dimension=dim,boundary='outflow',ranks=1 if dim==2 else 4)))
        # Refinement uses a fixed physical H/P/N; varying H is a distinct model test.
        for label,updates in (('h6',dict(layers=6)),('h24',dict(layers=24)),('p3',dict(degree=3)),('p5',dict(degree=5)),
                              ('dt025',dict(dt=.25)),('H1',dict(height=1)),('H4',dict(height=4)),
                              ('periodic',dict(boundary='periodic')),('molecular',dict(diffusion='Molecular')),
                              ('wrong_top',dict(wrong_top=True)),('wrong_kappa',dict(wrong_kappa=True)),
                              ('perturbed',dict(perturbation=.01))):
            jobs.append(('study_'+label,dict(base,model='cls_switched',dimension=2,boundary='outflow')|updates))
        for label,updates in (('h6',dict(layers=6)),('p3',dict(degree=3)),('dt025',dict(dt=.25)),('perturbed',dict(perturbation=.01))):
            jobs.append(('three_'+label,dict(base,model='cls_switched',dimension=3,boundary='outflow',ranks=4)|updates))
        def launch(job):
            name,opts=job
            try:return run(name,**opts)
            except RuntimeError as error:
                # A convergence/guard failure is a study outcome, never a pass.
                failure=dict(status='numerical_failure',message=str(error),options=opts,
                    stages=read_rows(out/name/(name+'.rans.csv')),wall=read_rows(out/name/(name+'.awf.csv')))
                (out/name/'failure.json').write_text(json.dumps(failure,indent=2)+'\n')
                with lock:records[name]=failure
                print(name+' NUMERICAL FAILURE: '+str(error),flush=True)
                return failure
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:list(pool.map(launch,jobs))

    if args.suite=='standard':
        for dim in (2,3):
            run('standard_perturbed_d'+str(dim),model='standard',dimension=dim,ranks=1 if dim==2 else 4,
                degree=4,layers=12,nx=1,dt=.5,flow_times=5,boundary='outflow',backend='matrixfree',
                perturbation=.01,tolerance=1e-8,monitor=True)

    if args.suite=='robust':
        def fine_linear(root):
            parameter(root,'IterativeSolverTolerance',1e-15)
            parameter(root,'RANSMaxIterations',10000)
        def robust(root):
            parameter(root,'IterativeSolverTolerance',1e-15)
            parameter(root,'RANSMaxIterations',10000)
            parameter(root,'RANSPicardRelaxation',.05)
            parameter(root,'RANSPressureRelaxation',.1)
        # The default study already tests the perturbed 3D problem for five
        # flows. Its smaller-step repetition is optional; unsteady time order
        # is established independently by the manufactured-solution sequence.
        if args.extended_perturbation:
            run('perturbed3d_small_dt',model='cls_switched',dimension=3,ranks=4,degree=4,layers=12,nx=1,dt=.1,flow_times=5,
                boundary='outflow',backend='matrixfree',perturbation=.01,tolerance=1e-7,monitor=True,modify=robust)
        if args.record_existing_failure:
            records['p5_3d']=recorded_failure(out/'p5_3d','p5_3d',binaries,'Recorded stopped refinement trial')
            if records['p5_3d']['engineering_stop']:
                commands.extend(records['p5_3d']['engineering_stop']['commands'])
        else:
            try:
                run('p5_3d',model='cls_switched',dimension=3,ranks=4,degree=5,layers=12,nx=1,dt=.25,flow_times=5,
                    boundary='outflow',backend='matrixfree',tolerance=1e-8,monitor=True,modify=fine_linear)
            except RuntimeError as error:
                records['p5_3d']=recorded_failure(out/'p5_3d','p5_3d',binaries,str(error))

    if args.suite=='columns':
        require(args.study_output is not None,'Columns require --study-output from the analytical study')
        def initial_file(root,file):
            initial=root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']");initial.clear();initial.attrib['NAME']='InitialConditions'
            for name in (v.text for v in root.findall('CONDITIONS/VARIABLES/V')):ET.SubElement(initial,'F',VAR=name,FILE=str(file))
            info(root,'RANSRestartMode','History')
        # Original-AWF periodic controls quantify its physical source imbalance.
        for dim in (2,3):
            try:
                run('original_periodic_d'+str(dim),model='cls_original',dimension=dim,ranks=1 if dim==2 else 4,
                    layers=12,nx=1,degree=4,dt=.25,flow_times=5,boundary='periodic',backend='matrixfree',tolerance=1e-8,monitor=True)
            except RuntimeError as error:
                records['original_periodic_d'+str(dim)]=dict(status='numerical_failure',message=str(error))
        study=args.study_output.resolve()
        sourcefile=study/'study_periodic/study_periodic.fld'
        run('column_settle',model='cls_switched',layers=12,nx=1,degree=4,dt=.5,flow_times=10,boundary='periodic',backend='matrixfree',tolerance=1e-8,
            monitor=True,modify=lambda r:initial_file(r,sourcefile))
        import numpy as np
        samples=list(csv.DictReader((out/'column_settle/profiles.csv').open()))
        reference={};fit_errors={};bounds=[2*15**(i/12) for i in range(13)]
        for field in {row['field'] for row in samples}:
            row=[r for r in samples if float(r['x'])==75 and r['field']==field]
            zz=np.array([float(r['height']) for r in row]);vv=np.array([float(r['value']) for r in row]);terms=[];fit_error=0
            for i,(a,b) in enumerate(zip(bounds,bounds[1:])):
                selected=((zz>=a-1e-9) if i==0 else (zz>a+1e-9)) & ((zz<=b+1e-9) if i==11 else (zz<b-1e-9))
                coeff=np.polynomial.polynomial.polyfit((2*zz[selected]-a-b)/(b-a),vv[selected],4)
                fitted=np.polynomial.polynomial.polyval((2*zz[selected]-a-b)/(b-a),coeff)
                fit_error=max(fit_error,float(np.max(abs(fitted-vv[selected]))))
                variable=f'((HEIGHT-({(a+b)/2:.17g}))/({(b-a)/2:.17g}))'
                polynomial='+'.join(f'({c:.17g})*{variable}^{j}' for j,c in enumerate(coeff))
                terms.append(f'(HEIGHT>=({a-1e-8 if i==0 else a:.17g}))*(HEIGHT'+('<=' if i==11 else '<')+f'({b+1e-8 if i==11 else b:.17g}))*({polynomial})')
            reference[field]='+'.join(terms);fit_errors[field]=fit_error
        require(max(fit_errors.values())<1e-9,'Numerical-column polynomial reconstruction is inaccurate')
        (out/'column-fit.json').write_text(json.dumps(dict(errors=fit_errors,expressions=reference),indent=2)+'\n')
        for dim in (2,3):
            def copy_column(root,dim=dim):
                functions={k:v.replace('HEIGHT','y' if dim==2 else 'z') for k,v in reference.items()}
                if dim==3:functions['w']=functions['v'];functions['v']='0'
                for label in ('InitialConditions','ExactSolution'):
                    for node in root.findall("CONDITIONS/FUNCTION[@NAME='"+label+"']/E"):
                        node.attrib['VALUE']=functions[node.attrib['VAR']]
                # Match the numerical column's weak traction at the outlet.
                region=root.find("CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF='1']")
                region.find("N[@VAR='u']").attrib['VALUE']='(2/3)*('+functions['k']+')-('+functions['rans_R_xx']+')-('+functions['p']+')'
                region.find("N[@VAR='"+('v' if dim==2 else 'w')+"']").attrib['VALUE']='-('+functions['rans_R_xz']+')'
                inlet=root.find("CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF='3']")
                for node in inlet.findall('D'):node.attrib['VALUE']=functions[node.attrib['VAR']]
                parameter(root,'IterativeSolverTolerance',1e-15)
            run('copied_column_d'+str(dim),model='cls_switched',dimension=dim,ranks=1 if dim==2 else 4,layers=12,nx=2,degree=4,dt=.1,
                flow_times=1,boundary='outflow',backend='matrixfree',tolerance=1e-8,monitor=True,modify=copy_column)

    if args.suite in ('restart','all'):
        verification=True
        def restart(root,file,mode='History'):
            initial=root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']");initial.clear();initial.attrib['NAME']='InitialConditions'
            for name in (v.text for v in root.findall('CONDITIONS/VARIABLES/V')):ET.SubElement(initial,'F',VAR=name,FILE=str(file))
            info(root,'RANSRestartMode',mode)
        for model in ('standard','cls_switched'):
            opts=dict(model=model,degree=4,layers=6,nx=1,dt=.05,backend='matrixfree',perturbation=0,monitor=True,tolerance=1e-8)
            a=run('restart_full_'+model,steps=8,**opts);prefix='restart_prefix_'+model
            run(prefix,steps=4,**opts);file=out/prefix/(prefix+'.fld')
            for ranks in (1,4):
                name=f'restart_r{ranks}_'+model
                run(name,steps=4,ranks=ranks,modify=lambda r:restart(r,file),**opts)
                compare((out/('restart_full_'+model)/('restart_full_'+model+'.xml'),out/('restart_full_'+model)/('restart_full_'+model+'.fld')),
                        out/name/(name+'.fld'),'compare_'+name,1e-7)
            checkpoint='checkpoint_prefix_'+model
            run(checkpoint,steps=4,checkpoint=1,**opts)
            first=out/checkpoint/(checkpoint+'_1.chk')
            name='checkpoint_restart_'+model
            run(name,steps=7,modify=lambda r:restart(r,first),**opts)
            compare((out/('restart_full_'+model)/('restart_full_'+model+'.xml'),out/('restart_full_'+model)/('restart_full_'+model+'.fld')),
                    out/name/(name+'.fld'),'compare_'+name,1e-7)
            run('restart_wrong_dt_'+model,steps=1,modify=lambda r:restart(r,file),negative='time or timestep mismatch',**dict(opts,dt=.1))
            run('restart_wrong_model_'+model,steps=1,modify=lambda r:(restart(r,file),parameter(r,'RANSCepsilon2',1.93)),negative='model/configuration mismatch',**opts)
        # Check old MPI payloads can be read by a different rank count in 3D.
        opts=dict(model='cls_switched',dimension=3,degree=3,layers=6,nx=1,dt=.02,backend='matrixfree',monitor=True)
        run('restart3_full',steps=6,ranks=4,**opts)
        run('restart3_prefix',steps=3,ranks=4,**opts)
        file3=out/'restart3_prefix/restart3_prefix.fld'
        run('restart3_r2',steps=3,ranks=2,modify=lambda r:restart(r,file3),**opts)
        compare((out/'restart3_full/restart3_full.xml',out/'restart3_full/restart3_full.fld'),out/'restart3_r2/restart3_r2.fld','compare_restart3',1e-6)
        opts=dict(model='cls_switched',degree=4,layers=6,nx=1,dt=.05,backend='matrixfree',perturbation=0,monitor=True,tolerance=1e-8)
        file=out/'restart_prefix_cls_switched/restart_prefix_cls_switched.fld'
        history=Path(records['restart_prefix_cls_switched']['audit'].get('history_directory',out/'restart_prefix_cls_switched/restart_prefix_cls_switched.rans-4'))
        for label,damage,expected in (
            ('checksum',lambda p:(p/'rank-0.bin').write_bytes((p/'rank-0.bin').read_bytes()[:-1]+b'X'),'payload checksum mismatch'),
            ('missing',lambda p:(p/'rank-0.bin').unlink(),'missing rank payload')):
            target=out/('damaged-'+label)
            if target.exists():shutil.rmtree(target)
            shutil.copytree(history,target);damage(target)
            run('restart_bad_'+label,steps=1,modify=lambda r:(restart(r,file),info(r,'RANSHistoryPath',str(target))),negative=expected,**opts)
        run('restart_wrong_order',steps=1,modify=lambda r:(restart(r,file),setattr(r.find('.//TIMEINTEGRATIONSCHEME/ORDER'),'text','1')),
            negative='model/configuration/order mismatch',**opts)
        repeated=run('history_repeat',model='standard',steps=2,layers=6,nx=1,degree=4,dt=.05,backend='matrixfree')
        old_history=Path(repeated['audit']['history_directory'])
        payload_before={p.name:digest(p) for p in old_history.iterdir() if p.is_file()}
        execute([solver,'history_repeat.xml'],out/'history_repeat','second_publication')
        raw=execute([audit,'history_repeat.xml','history_repeat.fld','history_repeat.fld'],out/'history_repeat','second_audit')
        new_history=Path(json.loads(raw[raw.index('{"time"'):])['history_directory'])
        require(new_history!=old_history and payload_before=={p.name:digest(p) for p in old_history.iterdir() if p.is_file()},
            'A repeated run overwrote history referenced by a backup field')
        repeated['backup_history_preserved']=dict(old=str(old_history),new=str(new_history),hashes=payload_before)
        for label,steps in (('full',6),('prefix',3)):
            run('imex1_'+label,generator=generate,kind='mms',order=1,dt=.025,steps=steps,modes=7,
                modify=lambda r:info(r,'RANSLinearSolver','MatrixFree'),profiles=False)
        run('imex1_restart',generator=generate,kind='mms',order=1,dt=.025,steps=3,modes=7,
            modify=lambda r:(info(r,'RANSLinearSolver','MatrixFree'),restart(r,out/'imex1_prefix/imex1_prefix.fld')),profiles=False)
        compare((out/'imex1_full/imex1_full.xml',out/'imex1_full/imex1_full.fld'),out/'imex1_restart/imex1_restart.fld','compare_imex1',1e-7)
        # The unsteady MMS independently tests multistep history and time order.
        for model,gen in (('standard',generate),('cls',generate_cls)):
            opts=dict(kind='mms',dt=.0125,modes=5 if model=='cls' else 7)
            select=lambda r:info(r,'RANSLinearSolver','MatrixFree')
            run('mms_full_'+model,generator=gen,steps=8,modify=select,profiles=False,**opts)
            prefix='mms_prefix_'+model
            run(prefix,generator=gen,steps=4,modify=select,profiles=False,**opts)
            run('mms_restart_'+model,generator=gen,steps=4,modify=lambda r:(select(r),restart(r,out/prefix/(prefix+'.fld'))),profiles=False,**opts)
            compare((out/('mms_full_'+model)/('mms_full_'+model+'.xml'),out/('mms_full_'+model)/('mms_full_'+model+'.fld')),
                    out/('mms_restart_'+model)/('mms_restart_'+model+'.fld'),'compare_mms_'+model,1e-7)
            for dt in (.025,.0125,.00625):
                run(f'unsteady_{model}_{dt:g}',generator=gen,kind='mms',dt=dt,steps=round(.1/dt),modes=5 if model=='cls' else 7,
                    modify=lambda r:info(r,'RANSLinearSolver','MatrixFree'),profiles=False)

    temporal_rates={}
    if args.suite in ('restart','all'):
        for model in ('standard','cls'):
            selected=[records[f'unsteady_{model}_{dt:g}'] for dt in (.025,.0125,.00625)]
            temporal_rates[model]={}
            for field in ('u','v','k','eps' if model=='standard' else 'eps_tilde','p'):
                errors=[r['audit']['fields'][field]['l2_error'] for r in selected]
                rates=[math.log(a/b,2) for a,b in zip(errors,errors[1:])]
                require(min(rates)>1.65,'Unsteady IMEX2 convergence failed: '+model+' '+field+' '+str(rates))
                temporal_rates[model][field]=rates
    report=dict(temporal_rates=temporal_rates,status='recorded',suite=args.suite,utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),threads=THREADS,
                maximum_workers=4*args.jobs,binaries=binaries,cases=records,commands=commands)
    (out/('results-'+args.suite+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],suite=args.suite,cases=len(records)),indent=2),flush=True)

if __name__=='__main__':main()
