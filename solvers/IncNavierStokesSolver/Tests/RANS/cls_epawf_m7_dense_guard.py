"""A dense-reconstruction rejection with all other admissibility sets valid."""
import sys
sys.dont_write_bytecode=True
import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def verify(run):
    from cls_epawf_m7_cases import perturbation_session
    from cls_epawf_m3_cases import write
    directory=run.output/'dense-only-reject';directory.mkdir()
    path=perturbation_session(run.template,directory/'case.xml',kind='none',chi=0,steps=0)
    root=ET.parse(path).getroot()
    for node in root.findall('CONDITIONS/FUNCTION'):
        if node.get('NAME','').startswith('InitialPerturbation'):
            for field in node.findall('E'):
                if field.get('VAR')=='k':
                    # This initial-condition localization is zero above the
                    # first element. It is not a limiter or closure clipping.
                    field.set('VALUE','-1.001*Qstar*0.16*4*(x/8)*(1-x/8)*(max(0*t,0.1-y)/0.1)^2'
                              if node.get('NAME').endswith('Inner') else '0')
    write(root,path)
    run.command('dense-only-reject','IncNavierStokesSolver',[path],expect_failure=True)
    report=directory/'case.evolution.json';raw=json.loads(report.read_text());audit=raw['state_audits'][-1]
    if raw['accepted_states'] or not raw['rejected_state'] or audit['reconstruction_dense'][1]>=0:
        raise AssertionError('Dense-only negative initial polynomial was accepted or not recorded')
    checked={}
    for name in ('source_volume_inner','source_volume_outer','physical_sample_K','prescribed_source_inner',
                 'mass_volume','interface_inner','interface_outer','reconstruction_native'):
        row=audit[name]
        if row[0]<=0 or row[1]<=0 or row[5]!=0:
            raise AssertionError('Another admissibility set failed before dense reconstruction: '+name+': '+str(row))
        checked[name]=row
    if audit['reconstruction_dense'][5]<=0 or 'sampled admissibility' not in raw['failure']:
        raise AssertionError('Wrong rejection mechanism')
    return dict(success=True,expected_rejection=True,sole_failed_admissibility_set='reconstruction_dense',
                native_report=str(report),native_report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
                xml=str(path),xml_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                positive_other_sets=checked,dense=audit['reconstruction_dense'],failure=raw['failure'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--source',type=Path)
    args=parser.parse_args()
    if platform.system()!='Linux':raise RuntimeError('Scientific guard validation belongs on Linux')
    source=(args.source or Path(__file__).resolve().parents[4]).resolve()
    sys.path.insert(0,str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
    from validate_cls_epawf_m4 import Runner
    run=Runner(args.build,args.output)
    before=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    result=dict(schema=1,source_commit=before,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),status='failed')
    try:result.update(result=verify(run),status='passed')
    except Exception as error:result['error']=str(error)
    result['commands']=run.commands
    assert subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()==before
    (run.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','error') if k in result}),flush=True)
    return 0 if result['status']=='passed' else 1

if __name__=='__main__':raise SystemExit(main())
