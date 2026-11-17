#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Cold/trial/transaction gates for the explicit inner-P2 pressure control."""
import argparse,os,resource,time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from wp2_native import Run
from wp2_all_state import Rows,US
from wp2_qualification import physical_gates,quarter_budget_differences,verify_report
from run_cls_epawf_picard import require,sha,THREADS
from run_cao_tamura_budget_wp2 import charge
from diagnose_cls_accepted_covariance import run_probe


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','probe','h1','preflight']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True)
    p.add_argument('--resume-qualified-trials',action='store_true')
    a=p.parse_args();os.environ.update(THREADS);r=a.root;record=r/('pressure-initial-gates-'+a.model+'.json')
    if record.exists():
        require(a.resume_qualified_trials,'pressure gate already has an owner')
        old=load(record);require(old['status']=='terminal_stop' and old['failure'] in ("main.<locals>.event() got multiple values for argument 'name'", "dict() got multiple values for keyword argument 'name'"),'only the identified bookkeeping stop may resume')
        require(old['solver_sha256']==sha(a.solver),'resumed gate runtime changed')
        archive=record.with_name(record.stem+'.before-bookkeeping-fix-'+sha(record)[:12]+'.json')
        require(not archive.exists(),'this exact bookkeeping stop was already resumed')
        atomic(archive,old)
    d=dict(version='WP2PressureInitialQualificationV1',status='waiting_cold_gates',model=a.model,solver_sha256=sha(a.solver),events=[])
    def save():atomic(record,d)
    def event(kind,**values):d['events'].append(dict(event_type=kind,epoch=time.time(),**values));save()
    save()
    try:
        from wp2_budget import wall_limit
        for batch in ['pressure-cold-primary','pressure-cold-confirmation']:
            while True:
                require(time.time()-load(r/'origin.json')['started_epoch']<wall_limit(r),'cumulative WP2 deadline')
                path=r/batch/'progress.json';s=load(path) if path.exists() else {}
                require(s.get('status')!='failed','cold native gate failed '+batch)
                if s.get('status')=='completed_native_batch':break
                time.sleep(5)
            for model in ['CLS','standard']:
                native=load(r/batch/(model+'-native/result.json'))
                checks=physical_gates(native['spatial'],US,.017598891430879275 if model=='CLS' else .01838575534504143,complete=True)
                require(native['boundary']['passed'] and checks['passed'],'cold boundary or full physical gate failed '+batch+' '+model)
                require(max(native['spatial']['analytical_free_block_norms'][:4])<1e-9,'analytical free weak loads differ')
            event('cold_batch_passed',batch=batch,progress_sha256=sha(path))
        rows=Rows(r/'pressure-row-verification/mixed_p2/result.json')
        K=.017598891430879275 if a.model=='CLS' else .01838575534504143
        reports={};cost=[]
        for name,dt,n in [('coarse',.02,2),('fine',.01,4)]:
            output=r/('pressure-trial-'+name+'-'+a.model)
            if output.exists():
                require(a.resume_qualified_trials,'trial output already exists')
                saved=load(output/'progress.json')
                require(saved['status']=='passed_pressure_trial' and saved['solver_sha256']==sha(a.solver) and
                        saved['case_definition']==dict(level=4,order=16,projection='RadialEnergyV1',pressure_space='inner-p2-v1') and
                        saved['dt']==dt and saved['accepted_steps']==n,'retained qualified trial identity differs')
                item=saved['results']['steps'];require(sha(Path(item['report']))==item['report_sha256'],'retained trial report changed')
                report=load(Path(item['report']));x=Run.__new__(Run);x.result=saved;x.output=output
                event('reused_qualified_trial',name=name,report_sha256=item['report_sha256'])
            else:
                x=Run(r,a.solver,a.h1,a.preflight,output,a.model,4,16,'RadialEnergyV1',False,'inner-p2-v1')
                charge(r,'pressure trial '+a.model+' '+name,timestep_trials=1)
                report=x.native('steps',n,dt)
            verify_report(report,US,K)
            replay=rows.check_native(report)
            extra=[rows.volume_audit(report,rows.tuple(s,report),order=36) for s in report['accepted_states']]
            require(all(z['passed'] for z in extra),'predeclared GL36 trial accuracy failed')
            require(all(z['native_full_G_action_defect']<1e-10 for z in report['pimple_attempts']),'full pressure action failed')
            x.result.update(status='passed_pressure_trial',native_tuple_replay=replay,independent_GL36=extra,accepted_steps=n,dt=dt)
            x.save();reports[name]=report;cost.append(x.result['results']['steps'])
            event('trial_passed',name=name,progress_sha256=sha(x.output/'progress.json'),seconds=cost[-1]['seconds'])
        difference=quarter_budget_differences(reports['coarse']['wp2_common_samples'],reports['fine']['wp2_common_samples'],US,K,.0065)
        require(difference['passed'],'pressure trial adjacent timestep quarter-budget failed')
        x=Run(r,a.solver,a.h1,a.preflight,r/('pressure-transaction-'+a.model),a.model,4,16,'RadialEnergyV1',False,'inner-p2-v1')
        x.fixture(.02)
        continued=[]
        for c in cost:
            for line in (Path(c['report']).parent/'timing.jsonl').read_text().splitlines():
                try:row=__import__('json').loads(line)
                except ValueError:continue
                if row['name']=='native.step':continued.append(row['wall_seconds'])
        require(continued,'missing measured native step timing')
        d.update(continued_native_step_seconds_max=max(continued),native_step_seconds=continued,
            status='passed_initial_pressure_gates',quarter_timestep=difference,cost=cost,peak_native_RSS_KiB=max(int(c['usage'].split()[-1]) for c in cost),
            seconds_per_physical_step_including_cold=[c['seconds']/n for c,n in zip(cost,[2,4])],transaction_progress_sha256=sha(x.output/'progress.json'),
            scope='Complete physical cold and trial checks, nonsingular native pressure Schur applications and exact full pressure action, fresh GL36, actual dt=.02 split/rollback. No long recovery inferred.')
    except Exception as e:d.update(status='terminal_stop',failure=str(e));raise
    finally:save()


if __name__=='__main__':main()
