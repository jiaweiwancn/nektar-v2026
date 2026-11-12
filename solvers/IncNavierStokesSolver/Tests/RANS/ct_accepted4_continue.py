#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finish only primitive-admitted CT4 six-millimetre pairs, never replace a run."""
import argparse,json,os,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted4_jobs import CASES,label
from ct_accepted4_run import qualify
from ct_accepted4_reference import write_json
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','runtime','initial-root','evolution-root','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();s=Path(__file__).parent;cases=[x for x in CASES if x[1]==6.]
    contract=a.root/'source/docs/plans/cao-tamura-accepted-4-evidence/continuation-contract.json'
    m=RequestManifest(a.output,dict(schema='CT4ConditionalSixV1',runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),
        contract_sha256=sha(contract),source_sha256=sha(Path(__file__)),
        short_receipts={label(*x):sha(a.evolution_root/label(*x)/'T0.008/results.json') for x in cases}))
    result=dict(diagnostic_execution='completed_selections',cases={})
    for model,h,degree in cases:
        key=label(model,h,degree);parent=a.evolution_root/key/'T0.008';q=json.loads((parent/'results.json').read_text())
        initial=a.initial_root/key;ini=json.loads((initial/'results.json').read_text())
        require(q['case']==ini['case'] and q['qualification']==qualify(q),'foreign_or_forged_short_gate')
        if q['qualification']['status']!='pass':
            result['cases'][key]=dict(status='not_admitted',short_sha256=sha(parent/'results.json'),reason=q['qualification']['reasons']);continue
        dest=a.evolution_root/key/'T0.02'
        if not dest.exists():
            cmd=[sys.executable,'-B',str(s/'ct_accepted4_run.py'),'--mode','advance','--runtime',str(a.runtime),
                '--references',str(a.root/'reference-v1/results.json'),'--contract',str(contract),'--model',model,'--height',str(h),'--degree',str(degree),
                '--output',str(dest),'--initial',str(initial),'--parent',str(parent),'--target','.02']
            with (a.output/(key+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
        else:code=None
        qq=json.loads((dest/'results.json').read_text());seal=json.loads((dest/'completed-artifacts.json').read_text())
        require(all(sha(dest/f)==v for f,v in seal['files'].items()) and qq['case']==ini['case'] and qq['qualification']==qualify(qq),'changed_continuation')
        result['cases'][key]=dict(status=qq['qualification']['status'],exit_code=code,receipt_sha256=sha(dest/'results.json'),body=qq['qualification'])
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT4 conditional selections', {k:v['status'] for k,v in result['cases'].items()})
    return 1 if any(q['status']=='fail' for q in result['cases'].values()) else 0
if __name__=='__main__':raise SystemExit(main())
