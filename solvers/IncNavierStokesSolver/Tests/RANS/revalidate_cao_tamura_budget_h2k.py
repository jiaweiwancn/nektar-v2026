#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Requalify a validation-only repair; reuse accepted native trajectories.

The certified seed, physical executable, libraries, probe and producer sources
must remain unchanged. No ordinary trajectory is rerun. One focused deliberate
post-assignment attempt verifies the corrected accepted-parent rollback check.
"""
import argparse,os
from pathlib import Path
from cao_tamura_h2d_state import load
from compare_cls_epawf_pimple import same_history
from run_cao_tamura_budget_h2f_startup import StartupRun
from cao_tamura_budget_h2g import case
from run_cls_epawf_picard import require,sha,state_hash,THREADS


def main():
    p=argparse.ArgumentParser()
    for n in ('root','solver','certificate','preflight','h1','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);cert=load(a.certificate);src=Path(__file__).resolve().parents[4]
    for name,h in cert['source_hashes'].items():require(sha(src/name)==h,'certified producer/probe source changed '+name)
    old=load(a.root/'pilot/progress.json');require(old['status']=='terminal_stop' and old.get('failure')=='accepted checkpoint missing','only identified validation-only stop may be requalified')
    integrity=load(a.root/'rollback-native-integrity.json');require(integrity['passed'] and all(integrity['checks'].values()),'original native rollback integrity failed')
    run=StartupRun(a.solver,a.certificate,a.preflight,a.h1,a.output,builder=lambda folder,h1,pre,flat=False:case(folder,h1,pre,flat=flat));reports={}
    try:
        for name,steps,dt,local in [('coarse-h0',0,2e-6,0),('coarse8',8,2e-6,8),('fine-h0',0,1e-6,0),('fine16',16,1e-6,16),('coarse4',4,2e-6,4),('coarse4-restart',8,2e-6,4)]:
            record=old['results'][name];require(sha(Path(record['report']))==record['report_sha256'],'accepted report changed '+name);d=load(Path(record['report']));run.out['results'][name]=record;run.validate_result(name,d,steps,dt);require(d['implicit_callback_count']==local and len(d['pimple_attempts'])==local,'callback count differs '+name);reports[name]=d
        for name in ('coarse-h0','fine-h0'):require(state_hash(reports[name]['state'])==cert['state_sha256'] and reports[name]['initialization_lift_calls_iterations']==[0,0],'reused import changed state or used Lift')
        require(same_history(reports['coarse8'],reports['coarse4-restart']),'reused direct/split bits differ')
        record=old['results']['post-assignment-reject'];require(sha(Path(record['report']))==record['report_sha256'],'original deliberate report changed');run.out['results']['original-rejection-revalidated']=record;run.validate_result('original-rejection-revalidated',load(Path(record['report'])),0,2e-6,True,Path(reports['coarse-h0']['checkpoint']))
        run.out.update(ordinary_trajectories_reused=True,ordinary_accepted_step_executions_reused=32,original_deliberate_rejection_revalidated=True,certified_native_sources_unchanged=True,validation_driver_sha256=sha(src/'solvers/IncNavierStokesSolver/Tests/RANS/run_cao_tamura_budget_h2f_startup.py'))
        run.save();run.rollback(reports['coarse-h0'],2e-6)
        run.out.update(status='completed_bounded_pair_revalidated',new_ordinary_physical_steps=0,new_deliberate_physical_attempts=1,scope='validation-only repair; original accepted trajectories retained, native and certificate closure unchanged; startup evidence only')
    except Exception as error:run.out.update(status='terminal_stop',failure=str(error));raise
    finally:run.save()


if __name__=='__main__':main()
