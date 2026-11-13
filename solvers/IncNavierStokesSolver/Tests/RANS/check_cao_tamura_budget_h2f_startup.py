#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Minimal unchanged restart-negative and V1 checks beside the nonzero V2 fixture."""
import argparse,os,subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import THREADS,sha,state_hash,require,set_info,set_parameter,write
from run_cao_tamura_budget_h2f_startup import StartupRun
from run_cao_tamura_budget_h2d import case
from compare_cls_epawf_pimple import same_history

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('fixture','old-v1','old-h0','preflight','h1','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);cert=a.fixture/'certificate.json';c=load(cert);solver=a.fixture/'runtime/IncNavierStokesSolver';seed=load(a.fixture/'evolution/h0/case.evolution.json');cp=Path(seed['checkpoint']);h=sha(cp/'State.bin');run=StartupRun(solver,cert,a.preflight,a.h1,a.output,flat=True)
    bad=run.native('ordinary-changed-dt',1e-6,1,cp);require(bad.get('failure') and 'checkpoint signature mismatch' in bad['failure'],'ordinary changed-dt restart accepted');require(sha(cp/'State.bin')==h,'negative restart changed parent bytes')
    corrupt=load(cert);row=load(a.fixture/'verification.json')['tests']['permitted_rows'][0];corrupt['state'][row]+=1e-7;corrupt['state_sha256']=state_hash(corrupt['state']);badcert=a.output/'altered-after-certification.json';atomic(badcert,corrupt);r=StartupRun(solver,badcert,a.preflight,a.h1,a.output/'corrupt',flat=True);bad=r.native('reject-state-receipt',2e-6);require(bad.get('failure') and 'certificate rejected' in bad['failure'],'changed state was accepted under old receipt');require(sha(cp/'State.bin')==h,'rejected certificate changed original parent')
    xml=case(a.output/'V1-regression',a.h1,load(a.preflight),flat=True);tree=ET.parse(xml).getroot();set_parameter(tree,'TimeStep',2e-6);set_parameter(tree,'NumSteps',0);set_info(tree,'RANSWallCertifiedStartupFile',str(a.old_v1));set_info(tree,'RANSWallCertifiedStartupSHA256',sha(a.old_v1));set_info(tree,'RANSWallCertifiedStartupVersion','CertifiedNeutralStartupV1');write(tree,xml)
    with (xml.parent/'solver.log').open('w') as log:p=subprocess.run([str(solver),str(xml)],cwd=xml.parent,env=dict(os.environ,LD_LIBRARY_PATH=c['runtime']['library_path']),stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    d=load(xml.with_suffix('.evolution.json'));require(p.returncode==0 and not d['failure'] and same_history(d,load(a.old_h0)),'V1 state/history regression changed')
    atomic(a.output/'checks.json',dict(passed=True,ordinary_changed_dt_rejected=True,changed_state_old_receipt_rejected=True,V1_native_h0_history_bitwise=True,parent_bytes_preserved=True,reports=run.out['results'],V1_report=str(xml.with_suffix('.evolution.json'))))
if __name__=='__main__':main()
