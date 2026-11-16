# SPDX-License-Identifier: MIT
"""Fail-fast new-policy native transactions. No historical certificate reuse."""
import os,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
from cao_tamura_h2d_state import load
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint
from run_cls_epawf_picard import require,sha,state_hash,set_parameter,write,THREADS
from run_cao_tamura_budget_wp2 import charge
from wp2_neutral_cases import case
from run_cao_tamura_budget_h2f_startup import verify_retained_rollback


class Run:
    def __init__(self,root,solver,h1,preflight,output,model,level=0,order=16,projection='RadialL2V1',resume=False):
        self.root,self.solver,self.h1,self.pre,self.output=root,solver,h1,load(preflight),output
        self.model,self.level,self.order,self.projection=model,level,order,projection
        if resume:
            self.result=load(output/'progress.json')
            require(self.result['status']=='completed_bounded_trajectory' and self.result['model']==model and
                    self.result['solver_sha256']==sha(solver) and self.result['case_definition']==dict(level=level,order=order,projection=projection),
                    'resume requires the same completed bounded trajectory, case and immutable runtime')
            return
        output.mkdir(exist_ok=False);self.result=dict(version='WP2NativeProtocolV1',status='running',model=model,
            solver=str(solver),solver_sha256=sha(solver),results={},case_definition=dict(level=level,order=order,projection=projection))
        self.save()
    def save(self):atomic(self.output/'progress.json',self.result)
    def native(self,label,steps=0,dt=2e-6,parent=None,reject=False,pout=0.,mismatch=False,perturb=False):
        charge(self.root,label,physical_attempts=steps)
        xml=case(self.output/label,self.h1,self.pre,self.model,self.level,self.order,steps,dt,parent,pout,self.projection)
        if reject:
            t=ET.parse(xml).getroot();set_parameter(t,'RANSWallPimplePostAssignRejectStep',(checkpoint(parent)['steps'] if parent else 0)+1);write(t,xml)
        if perturb:
            require(parent is None,'recovery perturbation is fresh setup, not restart-time repair')
            from run_cls_epawf_picard import set_info
            t=ET.parse(xml).getroot();set_info(t,'RANSWallWP2InitialPerturbation','OuterKBumpV1');write(t,xml)
        env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(self.solver.parent/'lib'),NEKTAR_RANS_TIMING_FILE=str(xml.parent/'timing.jsonl'))
        started=time.monotonic()
        remaining=43200-(time.time()-load(self.root/'origin.json')['started_epoch'])
        require(remaining>0,'WP2 aggregate wall envelope exhausted')
        with (xml.parent/'solver.log').open('w') as log:
            p=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(self.solver),str(xml)],cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=min(1800,remaining))
        report=xml.with_suffix('.evolution.json');d=load(report) if report.exists() else {}
        self.result['results'][label]=dict(code=p.returncode,seconds=time.monotonic()-started,xml_sha256=sha(xml),
            report=str(report),report_sha256=sha(report) if report.exists() else None,failure=d.get('failure'),usage=(xml.parent/'usage.txt').read_text(),
            requested_attempt_work=steps,observed_physical_attempts=d.get('implicit_callback_count',0),
            accepted_step_executions=max(0,d.get('steps',0)-(checkpoint(parent)['steps'] if parent else 0)));self.save()
        require(int(self.result['results'][label]['usage'].split()[-1])<=24*1024*1024,'WP2 per-process 24GiB RSS envelope exceeded')
        if mismatch:
            require(p.returncode!=0 and 'signature' in (xml.parent/'solver.log').read_text().lower(),'incompatible boundary restart was not rejected')
            return d
        require(p.returncode==(1 if reject else 0),'unexpected native process exit '+label)
        require(d and d.get('failure')==('deliberate post-assignment explicit rejection' if reject else ''),'unexpected native failure '+label+': '+str(d.get('failure')))
        expected=(checkpoint(parent)['steps'] if parent else 0)+(0 if reject else steps)
        require(d['steps']==expected and d['dt']==dt,'accepted count or timestep differs')
        cp=Path(parent) if reject else Path(d['checkpoint']);c=checkpoint(cp)
        require(c['signature']==d['signature'] and c['time']==d['time'] and c['steps']==d['steps'],'checkpoint identity/time/count differs')
        require(state_hash(c['coupled'])==state_hash(d['state']),'checkpoint/report coupled bits differ')
        for k in ('reference','defect'):require(state_hash(c[k])==state_hash(d['fixed_'+k]),'checkpoint/report '+k)
        require([state_hash(v) for v in c['history']]==[state_hash(v) for v in d['native_history']] and state_hash(c['history_times'])==state_hash(d['native_history_times']),'checkpoint/report histories differ')
        require(d['implicit_callback_count']==steps and len(d['pimple_attempts'])==steps,'reported attempt/solve inventory differs from requested native work')
        for row in d['pimple_attempts']:
            for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(row[k]==1,'one-sweep count '+k)
            require(row['solved_epsilon_retained_bitwise'],'solved epsilon overwritten')
        if not reject:require(all(x['accepted_awf_hard_admissible']==[1] for x in d['state_audits']),'native hard audit failed')
        return d
    def fixture(self):
        try:
            h0=self.native('h0');cp=Path(h0['checkpoint']);before=sha(cp/'State.bin')
            direct=self.native('direct2',2,parent=cp)
            one=self.native('split1',1,parent=cp)
            two=self.native('split1-restart',1,parent=Path(one['checkpoint']))
            require(same_history(direct,two),'new-policy split histories differ')
            bad=self.native('post-assignment-reject',1,parent=cp,reject=True)
            verify_retained_rollback(h0,bad,cp,before)
            self.native('changed-boundary-restart',0,parent=cp,pout=.001,mismatch=True)
            require(sha(cp/'State.bin')==before,'negative restart changed checkpoint')
            self.result.update(status='passed_fixture',split_bitwise=True,post_assignment_rollback=True,
                incompatible_boundary_restart_rejected=True,accepted_step_executions=4,deliberate_rejections=1,
                parent_checkpoint_sha256=before)
        except Exception as e:self.result.update(status='terminal_stop',failure=str(e));raise
        finally:self.save()
