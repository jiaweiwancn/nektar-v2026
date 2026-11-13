#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Final-runtime regression and independent import/gradient/tensor receipts."""
import argparse,os,json,subprocess,math
from pathlib import Path
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic,freeze
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,set_parameter,write
from compare_cls_epawf_pimple import same_history
from run_cao_tamura_budget_h2d import case
from diagnose_cls_accepted_covariance import run_probe
from validate_cao_tamura_budget_h2d import certificate

from cao_tamura_h2d_state import load

def stages(a,b):
    count=0
    for x,y in zip(a['pimple_attempts'],b['pimple_attempts']):
        require(x['post_flow_sha256']==y['post_flow_sha256'] and x['pre_k_sha256']==y['pre_k_sha256'],'native flow/epsilon stages differ')
        if x.get('k_snapshot') and Path(x['k_snapshot']).exists():
            px,py=load(Path(x['k_snapshot'])),load(Path(y['k_snapshot']))
            for k in ('committed','post_flow','pre_k','raw','corrected','history','fixed_reference'):require(state_hash(px[k])==state_hash(py[k]),'packet vector differs '+k)
            count+=1
    require(count>0,'no complete retained stage packets');return count

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','solver','probe','preflight','h1','old_h2c'):p.add_argument('--'+k.replace('_','-'),type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);r=a.root;out=r/'final-checks-v3';out.mkdir(exist_ok=False);solver,runtime=freeze(a.solver,out/'runtime');probe,pimage=freeze(a.probe,out/'probe-runtime');pre=load(a.preflight);summary=dict(status='running',runtime=runtime,probe=pimage,tests={});atomic(out/'progress.json',summary)
    v=r/'import-validation-v1';old=load(v/'progress.json');require(old['status']=='passed','prior import suite incomplete')
    ordinary=load(v/'ordinary2/case.evolution.json');imported=load(v/'import2/case.evolution.json');summary['tests']['import_stage_packets_bitwise']=stages(ordinary,imported)
    parent=load(v/'import-h0/case.evolution.json');bad=load(v/'post-assignment-reject/case.evolution.json')
    for k in ('state','fixed_reference','fixed_defect'):require(state_hash(parent[k])==state_hash(bad[k]),'rollback changed '+k)
    require(parent['physical_samples']==bad['physical_samples'],'rollback changed physical samples')
    summary['tests']['rollback_fields_reference_defect_samples_history']=True
    cert=v/'certificate.json'
    def native(name,steps=0,cp=None,mutate=None):
        xml=case(out/name,a.h1,pre,flat=True);t=ET.parse(xml).getroot();set_parameter(t,'NumSteps',steps);set_info(t,'RANSWallCertifiedStartupFile',str(cert));set_info(t,'RANSWallCertifiedStartupSHA256',sha(cert))
        if cp:set_info(t,'RANSWallRestart',str(cp))
        if mutate:mutate(t)
        write(t,xml)
        with (xml.parent/'solver.log').open('w') as f:code=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(solver),str(xml)],cwd=xml.parent,env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib')),stdout=f,stderr=subprocess.STDOUT,timeout=1800).returncode
        d=load(xml.with_suffix('.evolution.json'));summary.setdefault('native',{})[name]=dict(code=code,failure=d.get('failure'),report=str(xml.with_suffix('.evolution.json')),usage=(xml.parent/'usage.txt').read_text());atomic(out/'progress.json',summary);return d
    fresh=native('final-fresh-h0');require(not fresh['failure'] and same_history(fresh,parent),'fresh exact seed/history differs after literal-zero bypass')
    require(fresh['initialization_lift_calls_iterations']==[0,0],'certified fresh seed called Lift')
    q=native('final4',4,Path(fresh['checkpoint']));require(q.get('initialization_lift_calls_iterations')==[0,0],'certified restart called Lift');baseline=load(a.old_h2c/'flat-v1/flat4/case.evolution.json');require(not q['failure'] and same_history(q,baseline),'final runtime differs from frozen H2c affine path')
    summary['tests']['final_runtime_H2c_4_step_state_history_bitwise']=True;summary['tests']['final_runtime_H2c_stage_packets_bitwise']=stages(q,baseline)
    rollback=native('final-post-assignment-reject',1,Path(fresh['checkpoint']),lambda t:set_parameter(t,'RANSWallPimplePostAssignRejectStep',1))
    require(rollback['failure'] and same_history(fresh,rollback),'final no-Lift rollback changed accepted state/history')
    for k in ('fixed_reference','fixed_defect'):require(state_hash(fresh[k])==state_hash(rollback[k]),'final rollback changed '+k)
    summary['tests']['final_runtime_post_assignment_rollback']=True
    bad=native('final-changed-dt',1,Path(fresh['checkpoint']),lambda t:set_parameter(t,'TimeStep',1e-6));require(bad['failure'],'ordinary changed-dt restart accepted')
    summary['tests']['ordinary_changed_dt_rejected']=True
    def zero_perturbation(t):
        f=t.find('./CONDITIONS/FUNCTION[@NAME="InitialPerturbationInner"]')
        if f is None:f=ET.SubElement(t.find('CONDITIONS'),'FUNCTION',NAME='InitialPerturbationInner')
        for e in list(f):f.remove(e)
        ET.SubElement(f,'E',VAR='u',VALUE='1e-9')
    bad=native('reject-post-import-perturbation',mutate=zero_perturbation);require(bad['failure'] and not bad.get('initialization_complete',True),'post-import perturbation reached histories')
    def geometry_change(t):
        for v in t.findall('./GEOMETRY/VERTEX/V'):
            x=list(map(float,v.text.split()))
            if abs(x[1]-.00005)<1e-14:x[1]=.000055;v.text=' '.join(map(str,x))
    bad=native('reject-actual-geometry',mutate=geometry_change);require(bad['failure'] and not bad.get('initialization_complete',True),'changed geometry accepted')
    summary['tests']['fresh_history_import_guards']=True
    # Frozen candidate cannot be promoted by supplying only a boolean gate.
    best=load(r/'best-final-v1/progress.json');d=load(Path(best['results']['65-hill']['result_file']))
    try:certificate(d,Path(best['results']['65-hill']['result_file']));raise RuntimeError('failed seed certified')
    except (RuntimeError,ValueError) as e:require('cannot certify failed native inventory' in str(e),'unexpected certificate gate')
    summary['tests']['rejected_seed_certificate_refused']=True
    xml=case(out/'uncached',a.h1,pre);t=ET.parse(xml).getroot();set_info(t,'RANSWallAuditReuse','Off');write(t,xml)
    off,receipt=run_probe(probe,xml,dict(operation='h2c_seed',inner_stress_model='Linear',xml_sha256=sha(xml)),out/'uncached-probe')
    on=load(r/'replay-final-v1/PhysicalHeightInterfaceV2-probe/result.json')
    for k in ('state_sha256','audit','scalar_inventory','seed_constitutive','seed_target_witnesses','signed_derivative_checks','moving_foot_fixture'):require(off[k]==on[k],'cache changed '+k)
    summary['tests']['corrected_seed_cache_on_off_bitwise']=True;summary['uncached_receipt']=receipt
    summary['status']='passed';atomic(out/'progress.json',summary)
if __name__=='__main__':
    try:main()
    except Exception as error:
        import sys
        if '--root' in sys.argv:
            path=Path(sys.argv[sys.argv.index('--root')+1])/'final-checks-v3/progress.json'
            if path.exists():
                d=load(path);d.update(status='failed',failure=str(error));atomic(path,d)
        raise
