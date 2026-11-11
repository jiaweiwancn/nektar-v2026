#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Preserve the MatrixFree startup failure; complete scoped native history tests."""
import argparse,os,subprocess,time,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_duration_contract import require
from qualify_rans_preservation import load
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_pimple import prepare,check_report
from run_cls_epawf_picard import set_info,set_parameter,write


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--old-root',type=Path,required=True)
    p.add_argument('--output-label',default='directional-history-v3');a=p.parse_args()
    r=a.root.resolve();runtime=r/'runtime';source=r/'directional-history-v1';out=r/a.output_label
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
    q=load(source/'request-manifest.json')['request'];require(RequestManifest(source,q,True).completed(source),'history_source_not_sealed')
    old=load(source/'results.json');require(old['status']=='fail' and old['failure']=="ValueError('directional_axis_points_history_not_rejected')",'history_supplement_wrong_failure')
    require(q['runtime']==runtime_identity(runtime/'IncNavierStokesSolver'),'history_supplement_runtime_changed')
    seal=RequestManifest(out,dict(schema='DirectionalNativeHistorySupplementV2',source_request_sha256=sha(source/'request-manifest.json'),
        source_failed_result_sha256=sha(source/'results.json'),runtime=q['runtime'],mesh_sha256=q['mesh_sha256'],
        negative_fixture_collection='NoCollection, before inherited driver allocation; same production checkpoint reader and strict per-axis signature',
        startup_policy='explicit case.opt and native-clone .opt select NoCollection and prevent missing-opt automatic timing; input files are hashed, not checkpoint modifications',
        scope='reuse only recomputed passing target-layout histories, then independently execute both axis-signature rejections and ordinary CLS controls'))
    events=[];checks=[]
    def run(name,xml,expected=True):
        before=os.times();start=time.monotonic();maps={};folder=xml.parent
        with (folder/'solver.log').open('x') as f:
            proc=subprocess.Popen([str(runtime/'IncNavierStokesSolver'),'case.xml'],cwd=folder,stdout=f,stderr=subprocess.STDOUT)
            while proc.poll() is None:
                path=Path('/proc')/str(proc.pid)/'maps'
                if path.exists():
                    for line in path.read_text().splitlines():
                        w=line.split()
                        if len(w)>=6 and w[-1].startswith('/'):
                            lib=Path(w[-1])
                            if lib.is_file() and str(lib) not in maps:maps[str(lib)]=sha(lib)
                time.sleep(.03)
        after=os.times();events.append(dict(name=name,exit_code=proc.returncode,actual_loaded_maps=maps,
            configuration_sha256=sha(xml),wall_seconds=time.monotonic()-start,
            process_cpu_seconds=after.children_user-before.children_user+after.children_system-before.children_system))
        write_json(out/'progress.json',dict(events=events,checks=checks))
        report=load(folder/'case.evolution.json') if (folder/'case.evolution.json').exists() else {}
        require(str(runtime/'IncNavierStokesSolver') in maps,'history_loaded_identity_missing')
        if expected:require(proc.returncode==0 and not report.get('failure'),'history_native_execution:'+name)
        return report,proc.returncode,(folder/'solver.log').read_text()
    try:
        reports={n:load(source/n/'case.evolution.json') for n in ('directional-full','directional-first','directional-restart','directional-cache-off','directional-reject')}
        full=reports['directional-full'];first=reports['directional-first'];require(len(full['state'])==2769,'actual_directional_count_changed')
        require(same_history(full,reports['directional-restart']) and same_history(full,reports['directional-cache-off']) and
            same_history(first,reports['directional-reject']),'reused_history_primitive_mismatch')
        require('deliberate PIMPLE rejection' in reports['directional-reject']['failure'],'reused_rollback_wrong_reason')
        for name in reports:
            checks.append(dict(name='reuse-'+name,status='pass',report_sha256=sha(source/name/'case.evolution.json'),scope='sealed target 8x8-point history; recomputed raw state and complete native history equality'))
        for axis in (0,1):
            name='changed-axis-points-'+str(axis);folder=out/name;folder.mkdir();root=ET.parse(source/name/'case.xml').getroot()
            for e in root.findall('COLLECTIONS'):root.remove(e)
            ET.SubElement(root,'COLLECTIONS',{'DEFAULT':'NoCollection'})
            # SessionReader enables timing when an .opt sidecar is absent even
            # with an explicit XML default. Pin the existing NoCollection path
            # for both original and native-clone sessions before construction.
            for sidecar in ('case.opt','case-m4-native-rank0.opt'):
                (folder/sidecar).write_text('<NEKTAR><COLLECTIONS DEFAULT="NoCollection"/></NEKTAR>\n')
            set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'));write(root,folder/'case.xml')
            _,code,log=run(name,folder/'case.xml',False)
            require(code!=0 and 'signature mismatch' in log,'axis_signature_rejection_wrong_reason:'+str(axis))
            checks.append(dict(name=name,status='pass',reason='checkpoint signature mismatch',collection='NoCollection'))
        def cls(name,steps,restart=None,cache='ExactStepV1',reject=0):
            xml=prepare(out/name,1e-4,steps,diagnostic=True,restart=restart,checkpoint=True,reject_step=reject)
            root=ET.parse(xml).getroot();set_info(root,'RANSWallEvaluationReuse',cache);write(root,xml);return xml
        cf,_,_=run('cls-full',cls('cls-full',2));check_report(cf)
        ci,_,_=run('cls-first',cls('cls-first',1));check_report(ci)
        cr,_,_=run('cls-restart',cls('cls-restart',1,Path(ci['checkpoint'])));check_report(cr,True)
        co,_,_=run('cls-cache-off',cls('cls-cache-off',2,cache='Off'));check_report(co)
        cj,code,log=run('cls-reject',cls('cls-reject',2,reject=2),False)
        require(code!=0 and 'deliberate PIMPLE rejection' in cj.get('failure',''),'CLS_rollback_wrong_reason')
        require(same_history(cf,cr) and same_history(cf,co) and same_history(ci,cj),'CLS_history_or_rollback_mismatch')
        checks.append(dict(name='ordinary_CLS_restart_cache_rollback',status='pass'))
        oldfolder=a.old_root/'campaigns/A95-N-Accepted-duration-0.002-B5/T-0.0160';old_report=load(oldfolder/'case.evolution.json')
        folder=out/'old-isotropic-zero-step';folder.mkdir();root=ET.parse(oldfolder/'case.xml').getroot()
        set_info(root,'RANSWallRestart',old_report['checkpoint']);set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'));set_parameter(root,'NumSteps',0);write(root,folder/'case.xml')
        continued,_,_=run('old-isotropic-zero-step',folder/'case.xml');require(same_history(old_report,continued),'pinned_isotropic_history_changed')
        checks.append(dict(name='same_operator_pinned_isotropic_history',status='pass'))
        result=dict(status='pass',actual_unknowns=2769,checks=checks,execution=events,
            retained_first_failure_sha256=sha(source/'results.json'),scope='qualified target layout and strict axis signatures with NoCollection negative fixtures; unequal-point MatrixFree startup is not qualified')
    except Exception as e:result=dict(status='fail',failure=repr(e),checks=checks,execution=events)
    write_json(out/'results.json',result);seal.seal(out,[p for p in out.rglob('*') if p.is_file() and p.name not in ('request-manifest.json','completed-artifacts.json')])
    print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
