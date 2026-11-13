#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Sequential native CT7 checks after the isolated build preparation failure.

Run through ct_accepted3_command.py. New immutable output destinations are
mandatory. A completed command never grants physical/reference admission.
"""
import argparse,json,os,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted_open import run
from rans_request_manifest import RequestManifest,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','old','ct4','ct5','output','build-output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--suffix',required=True)
    p.add_argument('--checks',nargs='+',choices=('build','ctest','top','history','domain','regression'),default=['build','ctest','top','history','domain','regression'])
    a=p.parse_args();t=Path(__file__).parent
    m=RequestManifest(a.output,dict(schema='CT7SequentialNativeChecksV1',source=sha(Path(__file__)),checks=a.checks,
        build_destination=str(a.build_output),suffix=a.suffix,contract=sha(a.root/'source/docs/plans/cao-tamura-accepted-7-evidence/contract.json')))
    rows={}
    def execute(label,cmd):
        e=run([str(x) for x in cmd],a.output,label+'.log');rows[label]=e
        write_json(a.output/'progress.json',rows);print('CT7 native',label,'exit',e['exit_code'],flush=True)
        return e['exit_code']==0
    built=True
    if 'build' in a.checks:built=execute('build',[sys.executable,'-B',t/'ct_accepted6_build.py','--root',a.root,'--previous',a.ct4,'--output',a.build_output])
    rt=a.build_output/'runtime'
    if built:
        os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
        if 'ctest' in a.checks:execute('ctest',['ctest','--test-dir',a.root/'build','-R','^CT[67]','--parallel','1','--output-on-failure'])
        top=True
        if 'top' in a.checks:top=execute('top',[sys.executable,'-B',t/'ct_accepted7_top.py','--runtime',rt,'--output',a.root/('top-'+a.suffix)])
        if 'history' in a.checks and top:execute('history',[sys.executable,'-B',t/'ct_accepted7_history.py','--runtime',rt,'--output',a.root/('top-history-'+a.suffix)])
        if 'domain' in a.checks:execute('domain',[sys.executable,'-B',t/'ct_accepted7_domain.py','--runtime',rt,'--reference',a.root/'reference-v1/results.json','--output',a.root/('domain-'+a.suffix)])
        if 'regression' in a.checks:execute('regression',[sys.executable,'-B',t/'ct_accepted6_regression.py','--runtime',rt,'--ct5',a.ct5,'--member',a.old/'evolution-v2/A95-p4-dt002','--output',a.root/('regression-'+a.suffix)])
    result=dict(schema='CT7SequentialNativeChecksV1',executions=rows,execution='complete',
        status='pass' if rows and all(v['exit_code']==0 for v in rows.values()) else 'failed_required_check',physical_admission=False)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1

if __name__=='__main__':raise SystemExit(main())
