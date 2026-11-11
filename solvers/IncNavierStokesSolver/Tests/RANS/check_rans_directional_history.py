#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh directional cache/restart/rollback and ordinary CLS infrastructure checks."""
import argparse,os,subprocess,time,xml.etree.ElementTree as ET
from pathlib import Path
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_pimple import prepare,check_report
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_preservation_oracle import Reference
from rans_request_manifest import RequestManifest,runtime_identity,sha
from rans_qualification import THREADS,write_json
from rans_duration_contract import require
from qualify_rans_preservation import load


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--old-root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'directional-history-v1';runtime=r/'runtime'
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
    seal=RequestManifest(out,dict(schema='DirectionalNativeHistoryV1',runtime=runtime_identity(runtime/'IncNavierStokesSolver'),
        mesh_sha256=sha(a.mesh),dt=1e-4,steps=2,layout='6,5/4,3',unknowns=2769,
        old_corrected_isotropic_report_sha256=sha(a.old_root/'campaigns/A95-N-Accepted-duration-0.002-B5/T-0.0160/case.evolution.json')))
    events=[];checks=[]
    def run(name,xml,expected=True):
        before=os.times();started=time.monotonic();maps={};folder=xml.parent
        with (folder/'solver.log').open('x') as log:
            process=subprocess.Popen([str(runtime/'IncNavierStokesSolver'),'case.xml'],cwd=folder,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                path=Path('/proc')/str(process.pid)/'maps'
                if path.exists():
                    for line in path.read_text().splitlines():
                        words=line.split()
                        if len(words)>=6 and words[-1].startswith('/'):
                            lib=Path(words[-1])
                            if lib.is_file() and str(lib) not in maps:maps[str(lib)]=sha(lib)
                time.sleep(.03)
        after=os.times();events.append(dict(name=name,exit_code=process.returncode,wall_seconds=time.monotonic()-started,
            process_cpu_seconds=after.children_user-before.children_user+after.children_system-before.children_system,
            actual_loaded_maps=maps,configuration_sha256=sha(xml)))
        write_json(out/'progress.json',dict(events=events,checks=checks))
        require(str(runtime/'IncNavierStokesSolver') in maps,'directional_loaded_solver_identity_missing')
        report=load(folder/'case.evolution.json') if (folder/'case.evolution.json').exists() else {}
        if expected:require(process.returncode==0 and not report.get('failure'),'directional_native_history_execution:'+name)
        return report,process.returncode,(folder/'solver.log').read_text()
    def native(name,steps,restart=None,cache='ExactStepV1',reject=0):
        xml=case_xml(out/name,'A95','Natural','Accepted',steps=steps,dt=1e-4,restart=restart,cache=cache,reject=reject,
            descriptor_probe=runtime/'RANSNewModelsProbe',radial_knots=load(a.mesh)['knots'],degree=5,
            quadrature=(48,12,12,12,12),boundary_policy='FrozenReferenceV1')
        root=ET.parse(xml).getroot()
        for e in root.findall('./EXPANSIONS/E'):e.set('NUMMODES','5,4' if e.get('FIELDS')=='p' else '7,6')
        write(root,xml);return xml
    try:
        full,_,_=run('directional-full',native('directional-full',2))
        require(len(full['state'])==2769,'directional_actual_native_count_not_2769')
        first,_,_=run('directional-first',native('directional-first',1))
        restarted,_,_=run('directional-restart',native('directional-restart',1,Path(first['checkpoint'])))
        off,_,_=run('directional-cache-off',native('directional-cache-off',2,cache='Off'))
        rejected,code,log=run('directional-reject',native('directional-reject',2,reject=2),False)
        require(code!=0 and 'deliberate PIMPLE rejection' in rejected.get('failure',''),'directional_rollback_wrong_reason')
        checks += [dict(name=name,status='pass' if value else 'fail') for name,value in
            [('same_version_restart',same_history(full,restarted)),('cache_state_and_history',same_history(full,off)),
             ('complete_rollback',same_history(first,rejected))]]
        require(all(c['status']=='pass' for c in checks),'directional_state_or_history_mismatch')
        for axis in (0,1):
            name='changed-axis-points-'+str(axis);xml=native(name,0,Path(first['checkpoint']));root=ET.parse(xml).getroot()
            for e in root.findall('./EXPANSIONS/E'):
                points=e.get('NUMPOINTS').split(',');points[axis]=str(int(points[axis])+1);e.set('NUMPOINTS',','.join(points))
            write(root,xml);_,code,log=run(name,xml,False)
            require(code!=0 and 'signature mismatch' in log,'directional_axis_points_history_not_rejected')
            checks.append(dict(name=name,status='pass',reason='checkpoint signature mismatch'))
        def cls(name,steps,restart=None,cache='ExactStepV1',reject=0):
            xml=prepare(out/name,1e-4,steps,diagnostic=True,restart=restart,checkpoint=True,reject_step=reject)
            root=ET.parse(xml).getroot();set_info(root,'RANSWallEvaluationReuse',cache);write(root,xml);return xml
        cfull,_,_=run('cls-full',cls('cls-full',2));check_report(cfull)
        cfirst,_,_=run('cls-first',cls('cls-first',1));check_report(cfirst)
        cresume,_,_=run('cls-restart',cls('cls-restart',1,Path(cfirst['checkpoint'])));check_report(cresume,True)
        coff,_,_=run('cls-cache-off',cls('cls-cache-off',2,cache='Off'));check_report(coff)
        creject,code,log=run('cls-reject',cls('cls-reject',2,reject=2),False)
        require(code!=0 and 'deliberate PIMPLE rejection' in creject.get('failure',''),'CLS_rollback_wrong_reason')
        checks += [dict(name=name,status='pass' if value else 'fail') for name,value in
            [('CLS_restart',same_history(cfull,cresume)),('CLS_cache',same_history(cfull,coff)),('CLS_rollback',same_history(cfirst,creject))]]
        # Same operator and isotropic layout keep ordinary checkpoint compatibility.
        old=a.old_root/'campaigns/A95-N-Accepted-duration-0.002-B5/T-0.0160';old_report=load(old/'case.evolution.json')
        folder=out/'old-isotropic-zero-step';folder.mkdir();root=ET.parse(old/'case.xml').getroot()
        set_info(root,'RANSWallRestart',old_report['checkpoint']);set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'))
        set_parameter(root,'NumSteps',0);write(root,folder/'case.xml')
        old_resume,_,_=run('old-isotropic-zero-step',folder/'case.xml')
        require(old_resume['state']==old_report['state'] and same_history(old_report,old_resume),'ordinary_corrected_isotropic_history_changed')
        checks.append(dict(name='pinned_isotropic_checkpoint_compatibility',status='pass',scope='zero-step same-operator compatibility only'))
        require(all(c['status']=='pass' for c in checks),'directional_history_checks_failed')
        result=dict(status='pass',actual_unknowns=len(full['state']),checks=checks,execution=events,scope='fresh two-step histories, axis-specific signature rejection and focused ordinary CLS control; no duration qualification')
    except Exception as e:result=dict(status='fail',failure=repr(e),checks=checks,execution=events)
    write_json(out/'results.json',result)
    seal.seal(out,[p for p in out.rglob('*') if p.is_file() and p.name not in ('request-manifest.json','completed-artifacts.json')])
    print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
