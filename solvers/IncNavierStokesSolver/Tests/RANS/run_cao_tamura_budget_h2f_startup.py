#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""V2 exact-runtime startup successor: changed affine fixture or bounded hill pair."""
import argparse,os,subprocess,time
from pathlib import Path
import xml.etree.ElementTree as ET
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic,identity
from run_cao_tamura_budget_h2d import case
from run_cls_epawf_picard import THREADS,sha,state_hash,require,set_info,set_parameter,write
from compare_cls_epawf_pimple import same_history

class StartupRun:
    def __init__(self,solver,certificate,preflight,h1,output,flat=False):
        self.solver=solver.resolve();self.certpath=certificate.resolve();self.cert=load(certificate);self.pre=load(preflight);self.h1=h1;self.root=output;self.flat=flat
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
        xml=case(self.root/name,self.h1,self.pre,flat=self.flat);tree=ET.parse(xml).getroot();set_parameter(tree,'TimeStep',dt);set_parameter(tree,'NumSteps',steps);set_info(tree,'RANSWallCertifiedStartupFile',str(self.certpath));set_info(tree,'RANSWallCertifiedStartupSHA256',sha(self.certpath));set_info(tree,'RANSWallCertifiedStartupVersion','CertifiedNeutralStartupV2')
        if checkpoint:set_info(tree,'RANSWallRestart',str(checkpoint))
        if reject:set_parameter(tree,'RANSWallPimplePostAssignRejectStep',1)
        write(tree,xml);env=dict(os.environ,**THREADS);env['LD_LIBRARY_PATH']=self.cert['runtime']['library_path'];env['NEKTAR_RANS_TIMING_FILE']=str(xml.parent/'timing.jsonl');started=time.monotonic()
        with (xml.parent/'solver.log').open('w') as log:p=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(self.solver),str(xml)],cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        report=xml.with_suffix('.evolution.json');d=load(report) if report.exists() else {};self.out['results'][name]=dict(code=p.returncode,failure=d.get('failure'),report=str(report),report_sha256=sha(report) if report.exists() else None,seconds=time.monotonic()-started,usage=(xml.parent/'usage.txt').read_text(),accepted_steps=d.get('steps'));self.save();return d
    def accepted(self,d):
        require(d and not d.get('failure'),'native rejection is terminal; inspect preserved attempt')
        for v in d.get('pimple_attempts',[]):
            for key in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(v[key]==1,'one-sweep solve count changed')
            require(v['solved_epsilon_retained_bitwise'],'solved epsilon overwritten')
    def rollback(self,parent,dt):
        cp=Path(parent['checkpoint']);h=sha(cp/'State.bin');bad=self.native('post-assignment-reject',dt,1,cp,True);require(bad.get('failure')=='deliberate post-assignment explicit rejection','did not reach deliberate post-assignment rejection');require(same_history(parent,bad) and sha(cp/'State.bin')==h,'rollback changed parent state/history/bytes')
        for k in ('fixed_reference','fixed_defect'):require(state_hash(parent[k])==state_hash(bad[k]),'rollback changed '+k)
        require(parent['physical_samples']==bad['physical_samples'],'rollback changed physical samples')
    def pair(self):
        try:
            coarse=self.native('coarse-h0',2e-6);fine=self.native('fine-h0',1e-6);self.accepted(coarse);self.accepted(fine)
            for d in (coarse,fine):require(state_hash(d['state'])==self.cert['state_sha256'] and d['initialization_lift_calls_iterations']==[0,0],'seed import changed or used Lift')
            q=self.native('coarse8',2e-6,8,Path(coarse['checkpoint']));f=self.native('fine16',1e-6,16,Path(fine['checkpoint']))
            if not q.get('failure'):
                self.accepted(q);first=self.native('coarse4',2e-6,4,Path(coarse['checkpoint']));self.accepted(first);second=self.native('coarse4-restart',2e-6,4,Path(first['checkpoint']));self.accepted(second);require(same_history(q,second),'curved 4+4 restart differs');self.rollback(coarse,2e-6)
            if not f.get('failure'):self.accepted(f)
            self.out.update(status='completed_bounded_pair' if not q.get('failure') and not f.get('failure') else 'terminal_numerical_rejection',interpretation='startup/adjacent timestep evidence only; no temporal order/developed-flow/experimental claim')
        except Exception as e:self.out.update(status='terminal_stop',failure=str(e));raise
        finally:self.save()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','certificate','preflight','h1','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);StartupRun(a.solver,a.certificate,a.preflight,a.h1,a.output).pair()
if __name__=='__main__':main()
