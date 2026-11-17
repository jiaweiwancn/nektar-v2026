#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Gate-owned inner-P2 control, original cumulative ledger/deadline retained."""
import argparse,os,subprocess,time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import require,sha,THREADS
from run_cao_tamura_budget_wp2 import charge
from wp2_budget import wall_limit,pressure_scope
from wp2_all_state import Rows


def main():
    parser=argparse.ArgumentParser()
    for n in ['root','solver','probe','h1','preflight']:parser.add_argument('--'+n,type=Path,required=True)
    a=parser.parse_args();r=a.root;record=r/'pressure-waves-v2.json';require(not record.exists(),'pressure waves already owned')
    state=dict(version='WP2PressureWavesV2',status='waiting_initial_gates',pid=os.getpid(),source_sha256=sha(Path(__file__)),
               solver_sha256=sha(a.solver),probe_sha256=sha(a.probe),scope_sha256=sha(r/'pressure-control-scope.json'),events=[])
    atomic(record,state);env=dict(os.environ,**THREADS)
    def event(kind,**v):state['events'].append(dict(event_type=kind,epoch=time.time(),**v));atomic(record,state)
    def remaining():return wall_limit(r)-(time.time()-load(r/'origin.json')['started_epoch'])
    def launch(name,cmd):
        require(remaining()>120,'insufficient cumulative time for another bounded job')
        with (r/(name+'.log')).open('x') as log:j=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        event('launched',name=name,pid=j.pid,argv=cmd);return j
    try:
        require(pressure_scope(r) is not None,'explicit pressure scope required')
        while True:
            gates=[]
            for model in ['CLS','standard']:
                p=r/('pressure-initial-gates-'+model+'.json');d=load(p) if p.exists() else {}
                require(d.get('status')!='terminal_stop','pressure initial gate failed '+model+': '+str(d.get('failure')))
                gates.append(d)
            if all(d.get('status')=='passed_initial_pressure_gates' for d in gates):break
            require(remaining()>120,'cumulative deadline while waiting');time.sleep(5)
        require(all(d['solver_sha256']==sha(a.solver) for d in gates),'gate runtime differs')
        rss=max(d['peak_native_RSS_KiB'] for d in gates);require(rss<8*1024*1024,'pressure trial memory forecast exceeds scoped concurrency')
        available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))
        require(available>12*rss,'insufficient host memory for four bounded native cases')
        cost=max(d['continued_native_step_seconds_max'] for d in gates)
        # Balance the two 120-step legs per model in wave one, then 60-step
        # coarse/confirmation legs in wave two. Reserve 40% for the finer mesh.
        forecast=(120+60*1.4)*cost+1200
        require(forecast<remaining(),'pressure full comparison would exceed unchanged cumulative deadline')
        event('four_serial_cases_permitted',measured_trial_RSS_KiB=rss,available_memory_KiB=available,
              measured_continued_native_seconds_per_step=cost,remaining_seconds=remaining(),forecast_seconds=forecast,
              rank=1,math_threads=1,resource_snapshot=subprocess.getoutput('ps -eo pid,comm,pcpu,nlwp --sort=-pcpu | head -12'))
        while True:
            qp=r/'pressure-cold-order24/progress.json';quad=load(qp) if qp.exists() else {}
            require(quad.get('status')!='failed','pressure quadrature batch failed')
            if quad.get('status')=='completed_native_batch':break
            require(remaining()>120,'deadline while waiting for pressure quadrature gate');time.sleep(5)
        for model in ['CLS','standard']:
            reference=load(r/('pressure-trial-coarse-'+model)/'steps/case.evolution.json')
            data=load(r/'pressure-cold-order24'/(model+'-native/result.json'))
            row=Rows(r/'pressure-row-verification/mixed_p2/result.json')
            l=row.tuple(reference['fixed_reference'],reference);rr=row.tuple(data['spatial']['state'],reference)
            from wp2_qualification import quarter_budget_differences
            base=reference['wp2_common_samples'];check=quarter_budget_differences(dict(base,values=l['values']),dict(base,values=rr['values']),.07426793792419734,reference['parameters'][5],.0065)
            require(check['passed'],'higher production quadrature changes the principal fields')
            require(max(data['spatial']['unanchored_free_block_norms'][:4])<=1e-3,'higher-quadrature cold residual fails')
            event('production_quadrature_passed',model=model,changes=check,progress_sha256=sha(qp))
        # Freeze the finer native pressure row catalogue. No state is imported,
        # projected or advanced by this geometry-only export.
        xml=r/'pressure-cold-confirmation/standard-case/case.xml'
        common=load(r/'pressure-row-verification/mixed_p2/result.json')['physical_z_partition']
        charge(r,'pressure confirmation row catalogue',geometry_context_builds=1)
        packet=dict(operation='wp2_geometry_rows',inner_stress_model='Linear',expected_model='StandardKEpsilon',xml_sha256=sha(xml),
                    z_partition=common,pressure_rows='RegionalP2V1')
        data,receipt=run_probe(a.probe,xml,packet,r/'pressure-rows-confirmation')
        require(sum(data['sizes'])==11609,'unexpected confirmation native pressure dimensions')
        event('confirmation_rows_frozen',receipt=receipt,sizes=data['sizes'])
        common=['--root',str(r),'--solver',str(a.solver),'--h1',str(a.h1),'--preflight',str(a.preflight),'--pressure-space','inner-p2-v1']
        for wave in [[('recovery','recover',4,.01),('fine','preserve',4,.01)],
                     [('coarse','preserve',4,.02),('confirmation','preserve',5,.02)]]:
            jobs=[];state['status']='physical_wave';atomic(record,state)
            for phase,mode,level,dt in wave:
                for model in ['CLS','standard']:
                    name='pressure-'+phase+'-'+model;output=r/name;require(not output.exists(),'pressure output already exists')
                    cmd=['/usr/bin/python3',str(Path(__file__).with_name('run_wp2_preservation.py')),*common,'--output',str(output),
                         '--model',model,'--mode',mode,'--level',str(level),'--dt',str(dt),'--endpoint','1.2']
                    jobs.append((name,launch(name,cmd),output,phase))
            while any(j.poll() is None for _,j,_,_ in jobs):time.sleep(5)
            failure=[]
            for name,j,out,phase in jobs:
                p=out/'progress.json';d=load(p) if p.exists() else {}
                good=j.returncode==0 and d.get('status')=='completed_bounded_trajectory' and d.get('actual_time',0)>=1.1942387693600194
                if phase!='recovery':good=good and d.get('final_physical_budgets_passed',False)
                event('native_terminal',name=name,exit_code=j.returncode,status=d.get('status'),actual_time=d.get('actual_time'),accepted_steps=d.get('accepted_steps'),
                      progress_sha256=sha(p) if p.exists() else None,passed_required_native_stage=good)
                if not good:failure.append(name+': '+str(d.get('failure')))
            require(not failure,'pressure wave stopped before successors: '+'; '.join(failure))
        state['status']='readback';atomic(record,state)
        # Analyses remain read-only and cannot relaunch native evolution.
        jobs=[]
        for model in ['CLS','standard']:
            cmd=['/usr/bin/python3',str(Path(__file__).with_name('analyze_wp2_pressure_control.py')),'--root',str(r),
                 '--selected-rows',str(r/'pressure-row-verification/mixed_p2/result.json'),'--confirmation-rows',str(r/'pressure-rows-confirmation/result.json'),
                 '--model',model,'--output',str(r/('pressure-assessment-'+model))]
            jobs.append((model,launch('pressure-analysis-'+model,cmd)))
        while any(j.poll() is None for _,j in jobs):time.sleep(5)
        results={}
        for model,j in jobs:
            require(j.returncode==0,'pressure readback failed '+model)
            p=r/('pressure-assessment-'+model)/'result.json';d=load(p);results[model]=dict(passed=d['passed'],path=str(p),sha256=sha(p))
        state.update(status='completed_bounded_pressure_control',results=results,both_models_passed=all(x['passed'] for x in results.values()))
    except Exception as e:state.update(status='terminal_stop',failure=str(e));raise
    finally:atomic(record,state)


if __name__=='__main__':main()
