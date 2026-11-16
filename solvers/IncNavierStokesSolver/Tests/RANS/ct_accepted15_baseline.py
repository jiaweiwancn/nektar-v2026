#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read complete CT14 prerequisite bodies and verify actual loader bindings."""
import argparse,json,os,resource,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from rans_request_manifest import sha,runtime_identity
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import require
from ct_accepted_feasibility import write_json

def main():
    p=argparse.ArgumentParser();p.add_argument('--ct14',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(exist_ok=False);rt=a.ct14/'build-runtime-v1/runtime';os.environ['LD_LIBRARY_PATH']=str(rt/'lib');receipts={}
    for name in ('build-runtime-v1','capacity-v4','component-v1','static-v1','timestep-v1','pair-v2','collect-final-v1'):
        body=sealed(a.ct14/name);req=json.loads((a.ct14/name/'request-manifest.json').read_text())
        receipts[name]=dict(result_sha256=sha(a.ct14/name/'results.json'),request_sha256=sha(a.ct14/name/'request-manifest.json'),completion_sha256=sha(a.ct14/name/'completed-artifacts.json'),keys=sorted(body))
    runtimes={n:runtime_identity(rt/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe')}
    q=json.loads((a.ct14/'static-v1/request-manifest.json').read_text())['request']
    require(runtimes['IncNavierStokesSolver']==q['runtime'] and runtimes['RANSPreservationProbe']==q['probe'],'CT14 actual runtime/loaded-library bytes differ from sealed admission')
    pair=sealed(a.ct14/'pair-v2');native={}
    for dt,row in pair['members'].items():
        report=Path(row['directory'])/'case.evolution.json';require(sha(report)==row['native_report_sha256'],'changed oracle trajectory')
        r=json.loads(report.read_text());snapshots={}
        for x in r['pimple_attempts']:
            for key in ('flow_snapshot','k_snapshot'):
                if x.get(key):
                    f=Path(x[key]);require(sha(f)==x[key+'_sha256'],'changed actual stage witness');snapshots[str(f)]=dict(sha256=sha(f),bytes=f.stat().st_size)
        native[dt]=dict(report_sha256=sha(report),steps=r['steps'],times=r['audit_times'],snapshots=snapshots,execution=row['execution'],checkpoint=r['checkpoint'])
    resources={str(f.parent.name):json.loads(f.read_text()) for f in a.ct14.glob('resources-*/execution.json')}
    result=dict(schema='CT15BaselineVerificationV1',status='pass',runtimes=runtimes,complete_receipts=receipts,native=native,historical_CT14_CPU=sum(r['process_cpu_seconds'] for r in resources.values()),historical_family_count=len(resources),historical_scope='separate sunk ledger; CT11-13 remain in their immutable reports; not future production cost')
    write_json(a.output/'results.json',result);print(json.dumps(dict(status='pass',receipts=len(receipts),members=len(native),historical_CPU=result['historical_CT14_CPU'])),flush=True)
if __name__=='__main__':main()
