#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Conditional 8/16 hill pair bound to one exact independently certified seed."""
import argparse,os,subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic,freeze,CampaignLock
from run_cao_tamura_budget_h2d import case
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,set_parameter,write
from compare_cls_epawf_pimple import same_history

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('solver','preflight','h1','certificate','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();c=load(a.certificate)
    require(c['version']=='CertifiedNeutralStartupV1' and c['complete_native_admissible'],'exact complete native certificate required')
    require(c['state_sha256']==state_hash(c['state']) and c['complete_audit']['accepted_awf_hard_admissible']==[1],'certificate state/gate mismatch')
    require(c.get('extra_sampling') and c['extra_sampling']['minimum_covariance']>0,'independent extra-point covariance margin required')
    require(sha(Path(c['validation_report']))==c['validation_sha256'],'certificate native report changed')
    # No output experiment is created until all external certificate gates pass.
    a.output.mkdir(exist_ok=False);lock=CampaignLock(a.output);os.environ.update(THREADS)
    solver,runtime=freeze(a.solver,a.output/'runtime');pre=load(a.preflight);out=dict(status='started',certificate_sha256=sha(a.certificate),seed_sha256=c['state_sha256'],runtime=runtime,results={})
    def native(name,dt,steps=0,cp=None,reject=False):
        xml=case(a.output/name,a.h1,pre);t=ET.parse(xml).getroot();set_parameter(t,'TimeStep',dt);set_parameter(t,'NumSteps',steps)
        set_info(t,'RANSWallCertifiedStartupFile',str(a.certificate.resolve()));set_info(t,'RANSWallCertifiedStartupSHA256',sha(a.certificate))
        if cp:set_info(t,'RANSWallRestart',str(cp))
        if reject:set_parameter(t,'RANSWallPimplePostAssignRejectStep',1)
        write(t,xml);env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),NEKTAR_RANS_TIMING_FILE=str(xml.parent/'timing.jsonl'))
        with (xml.parent/'solver.log').open('w') as f:r=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(solver),str(xml)],cwd=xml.parent,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
        report=xml.with_suffix('.evolution.json');d=load(report) if report.exists() else {}
        out['results'][name]=dict(code=r.returncode,report=str(report),sha256=sha(report) if report.exists() else None,usage=(xml.parent/'usage.txt').read_text(),failure=d.get('failure'));atomic(a.output/'progress.json',out)
        require(d and (bool(d.get('failure'))==reject),'native rejection is terminal for this leg; inspect retained attempt')
        for v in d.get('pimple_attempts',[]):
            if not v.get('rejected',False):
                for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(v[k]==1,'schedule changed')
                require(v['solved_epsilon_retained_bitwise'],'epsilon overwritten')
        return d
    try:
        coarse=native('coarse-h0',2e-6);fine=native('fine-h0',1e-6)
        require(state_hash(coarse['state'])==state_hash(fine['state'])==c['state_sha256'],'import changed exact certified state')
        q=native('coarse8',2e-6,8,Path(coarse['checkpoint']));f=native('fine16',1e-6,16,Path(fine['checkpoint']))
        first=native('coarse4',2e-6,4,Path(coarse['checkpoint']));second=native('coarse4-restart',2e-6,4,Path(first['checkpoint']))
        require(same_history(q,second),'curved split restart differs')
        parent=Path(coarse['checkpoint']);old=sha(parent/'State.bin');bad=native('post-assignment-reject',2e-6,1,parent,True)
        require(same_history(coarse,bad) and sha(parent/'State.bin')==old,'curved rejection changed parent/history')
        out.update(status='completed_bounded_pair',interpretation='startup and adjacent timestep sensitivity only; no order or experimental validation claim')
    except Exception as e:out.update(status='terminal_stop',failure=str(e));raise
    finally:atomic(a.output/'progress.json',out);lock.stream.close()
if __name__=='__main__':main()
