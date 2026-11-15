#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""V2 exact-runtime startup successor: changed affine fixture or bounded hill pair."""
import argparse,os,subprocess,time,math
from pathlib import Path
import xml.etree.ElementTree as ET
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic,identity
from run_cao_tamura_budget_h2d import case
from run_cls_epawf_picard import THREADS,sha,state_hash,require,set_info,set_parameter,write
from compare_cls_epawf_pimple import same_history

class StartupRun:
    def __init__(self,solver,certificate,preflight,h1,output,flat=False,builder=None):
        self.solver=solver.resolve();self.certpath=certificate.resolve();self.cert=load(certificate);self.pre=load(preflight);self.h1=h1;self.root=output;self.flat=flat;self.builder=builder or case
        c=self.cert;require(c['version']=='CertifiedNeutralStartupV2' and c['complete_native_admissible'],'genuine V2 certificate required')
        require(state_hash(c['state'])==c['state_sha256'],'certificate state checksum')
        require(sha(self.solver)==c['runtime']['sha256'] and str(self.solver.parent/'lib')==c['runtime']['library_path'],'select certified runtime; do not refreeze it after issuance')
        for path,h in c['runtime']['libraries'].items():require(sha(Path(path))==h,'certified library changed')
        require(sha(Path(c['native_receipt_file']))==c['native_receipt_sha256'],'native receipt changed')
        receipt=load(Path(c['native_receipt_file']));require(receipt['passed'],'native receipt rejected')
        if not flat:require(receipt['gates']['optimizer_interior_margin'] and min(receipt['audit']['accepted_awf_covariance'])>=1e-4*receipt['K0']-1e-12 and receipt['extra_sampling']['minimum_covariance']>=1e-4*receipt['K0']-1e-12,'positive native/extra hill startup margin required')
        output.mkdir(exist_ok=False);self.out=dict(version='H2fV2StartupV1',status='running',certificate_sha256=sha(certificate),seed_sha256=c['state_sha256'],base_sha256=c['immutable_base_sha256'],runtime=c['runtime'],results={},flat_fixture=flat);self.save()
    def save(self):atomic(self.root/'progress.json',self.out)
    def native(self,name,dt,steps=0,checkpoint=None,reject=False):
        xml=self.builder(self.root/name,self.h1,self.pre,flat=self.flat);tree=ET.parse(xml).getroot();set_parameter(tree,'TimeStep',dt);set_parameter(tree,'NumSteps',steps);set_info(tree,'RANSWallCertifiedStartupFile',str(self.certpath));set_info(tree,'RANSWallCertifiedStartupSHA256',sha(self.certpath));set_info(tree,'RANSWallCertifiedStartupVersion','CertifiedNeutralStartupV2')
        if checkpoint:set_info(tree,'RANSWallRestart',str(checkpoint))
        if reject:set_parameter(tree,'RANSWallPimplePostAssignRejectStep',1)
        write(tree,xml);env=dict(os.environ,**THREADS);env['LD_LIBRARY_PATH']=self.cert['runtime']['library_path'];env['NEKTAR_RANS_TIMING_FILE']=str(xml.parent/'timing.jsonl');started=time.monotonic()
        with (xml.parent/'solver.log').open('w') as log:p=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(self.solver),str(xml)],cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        report=xml.with_suffix('.evolution.json');parse_error=None
        try:d=load(report) if report.exists() else {}
        except (ValueError,OSError) as error:d={};parse_error=str(error)
        self.out['results'][name]=dict(code=p.returncode,failure=d.get('failure'),report_parse_error=parse_error,report=str(report),report_sha256=sha(report) if report.exists() else None,seconds=time.monotonic()-started,usage=(xml.parent/'usage.txt').read_text(),accepted_steps=d.get('steps'));self.save();return d
    def checkpoint_steps(self,path):
        from diagnose_cls_accepted_covariance import checkpoint
        return checkpoint(path)['steps']
    def validate_checkpoint(self,d,expected,dt):
        from diagnose_cls_accepted_covariance import checkpoint
        cp=Path(d['checkpoint']);require((cp/'State.bin').is_file(),'accepted checkpoint missing')
        c=checkpoint(cp);require(c['steps']==expected and c['dt']==dt and c['order']==1,'checkpoint step/dt/order differs')
        require(state_hash(c['coupled'])==state_hash(d['state']),'checkpoint differs from reported accepted state')
        require(math.isfinite(c['time']) and abs(c['time']-d['time'])<1e-12,'checkpoint physical time differs')
        for name in ('reference','defect'):require(state_hash(c[name])==state_hash(d['fixed_'+name]),'checkpoint changed '+name)
    def validate_result(self,name,d,expected,dt,deliberate=False):
        record=self.out['results'][name];require(record['code']==(1 if deliberate else 0),'unexpected native process status '+name)
        require(Path(record['report']).is_file() and not record.get('report_parse_error') and isinstance(d,dict) and 'failure' in d,'missing or invalid native report '+name)
        success=d.get('failure')=='deliberate post-assignment explicit rejection' if deliberate else d.get('failure') in (None,'')
        require(success,'unexpected native rejection '+name+': '+str(d.get('failure')))
        require(d.get('steps')==expected and 'state' in d and 'checkpoint' in d,'unexpected accepted-step count/report '+name)
        self.validate_checkpoint(d,expected,dt)
        if not deliberate:self.accepted(d)
    def checked_native(self,name,dt,steps=0,checkpoint=None):
        expected=(self.checkpoint_steps(checkpoint) if checkpoint else 0)+steps
        d=self.native(name,dt,steps,checkpoint);self.validate_result(name,d,expected,dt);require(len(d.get('pimple_attempts',[]))==steps and d.get('implicit_callback_count')==steps,'missing/extra native callback or attempt '+name);return d
    def accepted(self,d):
        require(d and not d.get('failure'),'native rejection is terminal; inspect preserved attempt')
        for v in d.get('pimple_attempts',[]):
            require(v.get('scheduled_complete') is True,'incomplete governing attempt accepted')
            for key in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(v[key]==1,'one-sweep solve count changed')
            require(v['solved_epsilon_retained_bitwise'],'solved epsilon overwritten')
    def rollback(self,parent,dt):
        cp=Path(parent['checkpoint']);h=sha(cp/'State.bin');bad=self.native('post-assignment-reject',dt,1,cp,True);self.validate_result('post-assignment-reject',bad,parent['steps'],dt,True);require(bad.get('failure')=='deliberate post-assignment explicit rejection','did not reach deliberate post-assignment rejection');require(same_history(parent,bad) and sha(cp/'State.bin')==h,'rollback changed parent state/history/bytes')
        for k in ('fixed_reference','fixed_defect'):require(state_hash(parent[k])==state_hash(bad[k]),'rollback changed '+k)
        require(parent['physical_samples']==bad['physical_samples'],'rollback changed physical samples')
    def pair(self):
        try:
            coarse=self.checked_native('coarse-h0',2e-6)
            require(state_hash(coarse['state'])==self.cert['state_sha256'] and coarse['initialization_lift_calls_iterations']==[0,0],'coarse seed import changed or used Lift')
            q=self.checked_native('coarse8',2e-6,8,Path(coarse['checkpoint']))
            fine=self.checked_native('fine-h0',1e-6)
            require(state_hash(fine['state'])==self.cert['state_sha256'] and fine['initialization_lift_calls_iterations']==[0,0],'fine seed import changed or used Lift')
            f=self.checked_native('fine16',1e-6,16,Path(fine['checkpoint']))
            first=self.checked_native('coarse4',2e-6,4,Path(coarse['checkpoint']))
            second=self.checked_native('coarse4-restart',2e-6,4,Path(first['checkpoint']))
            require(same_history(q,second),'curved 4+4 restart differs')
            self.rollback(coarse,2e-6)
            self.out.update(status='completed_bounded_pair',interpretation='startup/adjacent timestep evidence only; no temporal order/developed-flow/experimental claim')
        except Exception as e:
            self.out.update(status='terminal_stop',failure=str(e));raise
        finally:self.save()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','certificate','preflight','h1','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);StartupRun(a.solver,a.certificate,a.preflight,a.h1,a.output).pair()
if __name__=='__main__':main()
