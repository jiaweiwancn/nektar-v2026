#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Sequential bounded CT7 verification; execution never substitutes admission.

Use one owned resource launcher. Every output destination must be new.
The separately launched fresh P4 finer trajectory must already be closed.
"""
import argparse,importlib.metadata,json,os,platform,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted_open import run,require
from rans_request_manifest import RequestManifest,sha,runtime_identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','old','ct4','ct5','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();r=a.root;s=r/'source';t=Path(__file__).parent;oldrt=a.old/'build-v12/runtime'
    source=subprocess.check_output(['git','-c','core.preloadIndex=false','-c','index.threads=1','rev-parse','HEAD'],cwd=s,text=True).strip()
    manifest=RequestManifest(a.output,dict(schema='CT7BoundedVerificationDriverV1',source_revision=source,
        source_sha256=sha(Path(__file__)),contract=sha(s/'docs/plans/cao-tamura-accepted-7-evidence/contract.json'),baseline_runtime=runtime_identity(oldrt/'IncNavierStokesSolver'),
        environment=dict(python=sys.version,platform=platform.platform(),packages={n:importlib.metadata.version(n) for n in ('numpy','scipy','mpmath')},threads={k:os.environ[k] for k in THREADS})))
    rows={}
    def execute(label,command):
        e=run([str(x) for x in command],a.output,label+'.log');rows[label]=e
        write_json(a.output/'progress.json',rows);print('CT7',label,'exit',e['exit_code'],flush=True)
        return e['exit_code']
    data=s/'docs/plans/cao-tamura-validation-data';contract=s/'docs/plans/cao-tamura-accepted-7-evidence/contract.json'
    audit=execute('frozen-data',[sys.executable,'-B',t/'cao_tamura_data.py','--data',data,'--report',a.output/'frozen-data.json'])
    math_ok=audit==0 and execute('displaced-reference',[sys.executable,'-B',t/'ct_accepted7_reference.py','--data',data,'--contract',contract,'--output',r/'reference-v1'])==0
    built=execute('build',[sys.executable,'-B',t/'ct_accepted6_build.py','--root',r,'--previous',a.ct4,'--output',r/'build-v1'])==0
    if built:
        rt=r/'build-v1/runtime';os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
        execute('ctest',['ctest','--test-dir',r/'build','-R','^CT6','--parallel','1','--output-on-failure'])
        top=execute('top-operator',[sys.executable,'-B',t/'ct_accepted7_top.py','--runtime',rt,'--output',r/'top-v1'])
        if top==0:execute('top-history-budgets',[sys.executable,'-B',t/'ct_accepted7_history.py','--runtime',rt,'--output',r/'top-history-v1'])
        if math_ok:execute('domain-preflight',[sys.executable,'-B',t/'ct_accepted7_domain.py','--runtime',rt,'--reference',r/'reference-v1/results.json','--output',r/'domain-v1'])
        execute('flat-curved-regressions',[sys.executable,'-B',t/'ct_accepted6_regression.py','--runtime',rt,'--ct5',a.ct5,'--member',a.old/'evolution-v2/A95-p4-dt002','--output',r/'regression-v1'])
    # Diagnostic scope is independent of the physical fit/root outcome.
    os.environ['LD_LIBRARY_PATH']=str(oldrt/'lib')
    finer=execute('finer-p4-full-audit',[sys.executable,'-B',t/'ct_accepted6_audit.py','--runtime',oldrt,'--member',r/'ct6-finer-p4','--output',r/'ct6-finer-p4-audit'])
    if finer==0:
        execute('temporal',[sys.executable,'-B',t/'ct_accepted7_temporal.py','--old',a.old,'--fine',r/'ct6-finer-p4','--audit',r/'ct6-finer-p4-audit','--output',r/'temporal-v1'])
        decision=json.loads((r/'temporal-v1/results.json').read_text())
        if decision['matched_p3_eligible']:
            code=execute('fresh-p3-finer',[sys.executable,'-B',t/'ct_accepted6_evolve.py','--runtime',oldrt,'--contract',contract,'--admission',a.old/'admission-v5.json','--model','A95','--degree','3','--dt','.0005','--target','.016','--output',r/'ct6-finer-p3'])
            if code==0:
                ac=execute('finer-p3-full-audit',[sys.executable,'-B',t/'ct_accepted6_audit.py','--runtime',oldrt,'--member',r/'ct6-finer-p3','--output',r/'ct6-finer-p3-audit'])
                if ac==0:execute('matched-degree',[sys.executable,'-B',t/'ct_accepted7_temporal.py','--old',a.old,'--fine',r/'ct6-finer-p4','--audit',r/'ct6-finer-p4-audit','--fine-p3',r/'ct6-finer-p3','--audit-p3',r/'ct6-finer-p3-audit','--output',r/'temporal-degree-v1'])
    result=dict(schema='CT7BoundedVerificationDriverV1',execution='complete',executions=rows,
        required_numerical_process_status='pass' if all(v['exit_code']==0 for v in rows.values()) else 'failed_or_unresolved_required_check',
        physical_admission=False,scope='do not infer physical, full-domain or model-transfer admission from process completion; final primitive aggregation is separate')
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    return 0 if all(v['exit_code']==0 for v in rows.values()) else 1

if __name__=='__main__':raise SystemExit(main())
