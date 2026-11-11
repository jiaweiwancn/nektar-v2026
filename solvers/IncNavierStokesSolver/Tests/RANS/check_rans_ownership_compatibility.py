#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native fresh-history restart, rollback, cache and shared CLS controls."""
import argparse
import os
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_pimple import prepare, check_report
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info, set_parameter, write
from rans_duration_contract import require
from rans_request_manifest import RequestManifest, sha, runtime_identity
from rans_qualification import THREADS, write_json
from qualify_rans_preservation import load
from rans_preservation_oracle import Reference

OPERATOR='GeneratedVolumeOwnerV1'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--previous-root',type=Path,required=True);a=p.parse_args();r=a.root.resolve()
    out=r/'compatibility-v1';runtime=r/'runtime';solver=runtime/'IncNavierStokesSolver'
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
    old=a.previous_root/'campaigns/A95-N-Accepted-duration-0.002-B5/T-0.0160'
    require(old.exists(),'old_operator_fixture_missing')
    request=dict(schema='NativeOwnershipCompatibilityV1',operator=OPERATOR,runtime=runtime_identity(solver),
                 old_configuration_sha256=sha(old/'case.xml'),old_native_report_sha256=sha(old/'case.evolution.json'))
    seal=RequestManifest(out,request,False);events=[];checks=[]
    def execute(name,xml,expected=True):
        folder=xml.parent;before=os.times();start=time.monotonic();mapped={}
        with (folder/'solver.log').open('x') as stream:
            proc=subprocess.Popen([str(solver),'case.xml'],cwd=folder,stdout=stream,stderr=subprocess.STDOUT,env=dict(os.environ))
            while proc.poll() is None:
                maps=Path('/proc')/str(proc.pid)/'maps'
                if maps.exists():
                    for line in maps.read_text().splitlines():
                        words=line.split()
                        if len(words)>=6 and words[-1].startswith('/'):
                            path=Path(words[-1])
                            if path.is_file() and str(path) not in mapped:mapped[str(path)]=sha(path)
                time.sleep(.03)
        after=os.times();code=proc.returncode
        events.append(dict(name=name,exit_code=code,configuration_sha256=sha(xml),
            actual_mapped_files=mapped,process_cpu_seconds=(after.children_user-before.children_user)+(after.children_system-before.children_system),
            wall_seconds=time.monotonic()-start,expected_execution='success' if expected else 'rejection'))
        write_json(out/'progress.json',dict(execution=events,checks=checks))
        require(bool(mapped) and str(solver) in mapped,'actual_loaded_identity_missing')
        if expected:
            require(code==0,'native_execution_failed:'+name)
            report=load(folder/'case.evolution.json');require(not report['failure'],'native_physical_failure:'+name)
            require(report['numerical_operator']==OPERATOR,'wrong_native_operator')
            return report
        return code,(folder/'solver.log').read_text()
    def newer(name,steps,restart=None,cache='ExactStepV1',reject=0):
        folder=out/name;ref=Reference('E94E95','Natural','Accepted')
        xml=case_xml(folder,'E94E95','Natural','Accepted',steps=steps,dt=1e-4,restart=restart,cache=cache,
            descriptor_probe=runtime/'RANSNewModelsProbe',radial_knots=ref.knots(),degree=3,
            quadrature=(32,8,8,12,8),checkpoint_interval=1,reject=reject)
        return xml
    try:
        full=execute('new-full',newer('new-full',4));first=execute('new-first',newer('new-first',2))
        restart=Path(first['checkpoint']);continued=execute('new-restart',newer('new-restart',2,restart))
        checks.append(dict(name='same_version_restart',status='pass' if same_history(full,continued) else 'fail'))
        off=execute('new-cache-off',newer('new-cache-off',4,cache='Off'))
        checks.append(dict(name='cache_state_and_history',status='pass' if same_history(full,off) else 'fail'))
        reject_xml=newer('new-reject',3,reject=3)
        code,log=execute('new-reject',reject_xml,False);rejected=load(reject_xml.parent/'case.evolution.json')
        require('deliberate PIMPLE rejection' in rejected['failure'],'rollback_wrong_failure_reason')
        checks.append(dict(name='complete_rejected_step_rollback',status='pass' if same_history(first,rejected) else 'fail',exit_code=code,reason=rejected['failure']))
        legacy=out/'old-restart';legacy.mkdir();tree=ET.parse(old/'case.xml');root=tree.getroot()
        old_report=load(old/'case.evolution.json');old_checkpoint=Path(old_report['checkpoint'])
        set_info(root,'RANSWallRestart',str(old_checkpoint));set_info(root,'RANSWallCheckpoint',str(legacy/'checkpoint'))
        set_parameter(root,'NumSteps',0);write(root,legacy/'case.xml')
        code,log=execute('old-restart',legacy/'case.xml',False)
        require(code!=0 and ('signature' in log.lower() or 'fingerprint' in log.lower()),'old_operator_rejection_reason_missing')
        checks.append(dict(name='old_operator_native_restart_rejected',status='pass',reason=log[-2400:],old_signature=old_report.get('signature')))
        def cls(name,steps,restart=None,cache='ExactStepV1',reject=0):
            xml=prepare(out/name,1e-4,steps,diagnostic=True,restart=restart,checkpoint=True,reject_step=reject)
            root=ET.parse(xml).getroot();set_info(root,'RANSWallEvaluationReuse',cache);write(root,xml);return xml
        cfull=execute('cls-full',cls('cls-full',2));check_report(cfull)
        cfirst=execute('cls-first',cls('cls-first',1));check_report(cfirst)
        cresume=execute('cls-restart',cls('cls-restart',1,Path(cfirst['checkpoint'])));check_report(cresume,True)
        coff=execute('cls-cache-off',cls('cls-cache-off',2,cache='Off'));check_report(coff)
        checks.extend([dict(name='ordinary_CLS_restart',status='pass' if same_history(cfull,cresume) else 'fail'),
                       dict(name='ordinary_CLS_cache',status='pass' if same_history(cfull,coff) else 'fail')])
        rejected_xml=cls('cls-reject',2,reject=2);code,log=execute('cls-reject',rejected_xml,False)
        crejected=load(rejected_xml.parent/'case.evolution.json')
        require('deliberate PIMPLE rejection' in crejected['failure'],'CLS_rollback_wrong_reason')
        checks.append(dict(name='ordinary_CLS_rollback',status='pass' if same_history(cfirst,crejected) else 'fail',reason=crejected['failure']))
        result=dict(status='pass',checks=checks,execution=events,numerical_operator=OPERATOR,
            scope='native same-version fresh history and actual old-history rejection; CLS diagnostic infrastructure control, not renewed CLS physical qualification')
    except Exception as error:result=dict(status='fail',failure=repr(error),checks=checks,execution=events,numerical_operator=OPERATOR)
    write_json(out/'results.json',result)
    seal.seal(out,['results.json',*[str(p) for p in out.glob('*/solver.log')],*[str(p) for p in out.glob('*/case.xml')]])
    print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
