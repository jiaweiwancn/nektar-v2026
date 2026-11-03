#!/usr/bin/env python3
"""Read-only requalification of stored unforced M7 accepted prefixes."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
sys.dont_write_bytecode=True
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
assert platform.system()=='Linux'
source=a.source.resolve();a.output.mkdir(parents=True,exist_ok=False)
sys.path.insert(0,str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
import validate_cls_epawf_m7_transition_diagnosis as diagnosis
import validate_cls_epawf_m7_transition_flow as flow
import validate_cls_epawf_m7 as native

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
original=json.loads(a.input.read_text());assert original['status']=='diagnosis_collected' and not original['collection_errors']
revision=git('rev-parse','HEAD');assert not git('status','--porcelain')
changed=git('diff','--name-only','bf202be','--','library','solvers/IncNavierStokesSolver').splitlines()
assert not [x for x in changed if Path(x).suffix in ('.cpp','.h','.hpp','.cc','.cxx','.c')]
for path,expected in native.FROZEN.items():assert sha(source/path)==expected,path
receipt=copy.deepcopy(original)
receipt.update(replay_schema=1,audit_source_commit=revision,audit_source_clean=True,
               original_diagnosis_path=str(a.input),original_diagnosis_sha256=sha(a.input),
               replay_script_sha256=sha(__file__),started_replay_utc=datetime.now(timezone.utc).isoformat())
audit=flow.Audit();views={};artifacts=[]
try:
 assert sha(original['program']['path'])==original['program']['sha256']==original['program']['sha256_end']
 for record in receipt['cases']:
  report=Path(record['native_report']);xml=Path(record['xml']['path'])
  assert sha(report)==record['native_report_sha256'] and sha(xml)==record['xml']['sha256']
  command=record.get('command') or record.get('original_command')
  assert sha(command['log'])==command['log_sha256']
  artifacts.extend([dict(path=str(report),sha256=sha(report)),dict(path=str(xml),sha256=sha(xml)),dict(path=command['log'],sha256=sha(command['log']))])
  raw=json.loads(report.read_text());views[record['name']]=diagnosis.analyze(raw,record,audit)
 receipt['pairs'],_=diagnosis.compare_prefixes(views,audit)
 receipt.update(replay_status='passed',replay_checks=audit.checks,replay_errors=audit.errors,
                artifact_pins=artifacts,passage_qualified=original['passage_qualified'],
                scope='original solver execution/source provenance retained; independent current covariance, energy and accepted-prefix audit only')
 assert all(r['observation_verified'] for r in receipt['cases'])
except Exception as error:receipt.update(replay_status='failed',replay_failure=str(error))
assert git('rev-parse','HEAD')==revision and not git('status','--porcelain')
assert sha(original['program']['path'])==original['program']['sha256']
receipt['finished_replay_utc']=datetime.now(timezone.utc).isoformat()
(a.output/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt.get(k) for k in ('replay_status','replay_checks','passage_qualified','replay_failure')}),flush=True)
raise SystemExit(0 if receipt['replay_status']=='passed' else 1)
