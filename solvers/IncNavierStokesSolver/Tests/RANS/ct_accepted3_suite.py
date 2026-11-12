#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Execute the three independent remaining fixed cases serially, retain failures."""
import argparse, json, os, subprocess, sys
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from rans_request_manifest import RequestManifest,sha,runtime_identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','ct1','ct2','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();cases=['A95-XZminus-P4','A93-Zminus-P3','E94E95-Zminus-P3'];s=Path(__file__).parent
    manifest=RequestManifest(a.output,dict(schema='CT3IndependentFixedCasesV1',cases=cases,
        runtime=runtime_identity(a.ct2/'runtime/IncNavierStokesSolver'),contract_sha256=sha(a.root/'source/docs/plans/cao-tamura-accepted-3-evidence/contract.json'),
        source_sha256={f:sha(s/f) for f in ('ct_accepted3_suite.py','ct_accepted3_jobs.py','ct_accepted3_run.py','ct_accepted3_metrics.py','ct_accepted3_pressure.py','ct_accepted_budget_gate.py')}))
    result=dict(diagnostic_execution='incomplete',cases={})
    for case in cases:
        cmd=[sys.executable,'-B',str(s/'ct_accepted3_jobs.py'),'--root',str(a.root),'--ct1',str(a.ct1),'--ct2',str(a.ct2),
             '--case',case,'--output',str(a.root/case/'campaign-v1')]
        with (a.output/(case+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
        file=a.root/case/'campaign-v1/results.json'
        result['cases'][case]=dict(exit_code=code,command=cmd,result_sha256=sha(file) if file.exists() else None,
            body=json.loads(file.read_text()) if file.exists() else None)
        write_json(a.output/'progress.json',result)
    result['diagnostic_execution']='completed_all_independent_cases';write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    print('all fixed cases attempted; classifications retained')
    return 0 if all(q['exit_code']==0 for q in result['cases'].values()) else 1
if __name__=='__main__':raise SystemExit(main())
