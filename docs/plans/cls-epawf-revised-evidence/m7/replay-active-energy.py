import hashlib,json,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
sys.dont_write_bytecode=True
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
source=root/'m7-audit-source';tests=source/'solvers/IncNavierStokesSolver/Tests/RANS'
sys.path.insert(0,str(tests))
from cls_epawf_m7_energy_oracle import verify
from validate_cls_epawf_m7 import audits
out=root/'m7-energy-preflight-003';out.mkdir(exist_ok=False)
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
assert not git('status','--porcelain')
receipt=dict(schema=1,source_commit=git('rev-parse','HEAD'),started_utc=datetime.now(timezone.utc).isoformat(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=[],status='failed')
try:
 for name,path in [('flat-joint',root/'m7-attempt-004/native/joint/case.evolution.json'),('unforced-control',root/'m7-flow-pilot-001/preflight-control/case.evolution.json'),('unforced-pulse',root/'m7-flow-pilot-001/preflight-pulse/case.evolution.json')]:
  raw=json.loads(path.read_text());sha=hashlib.sha256(path.read_bytes()).hexdigest()
  record=dict(name=name,path=str(path),sha256=sha,native_audit=audits(raw,name),independent_energy=verify(raw))
  assert hashlib.sha256(path.read_bytes()).hexdigest()==sha
  receipt['cases'].append(record);print(name,record['independent_energy']['scope'],flush=True)
 receipt['status']='passed'
except Exception as error:receipt['error']=str(error)
assert git('rev-parse','HEAD')==receipt['source_commit'] and not git('status','--porcelain')
receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'cases':len(receipt['cases']),'error':receipt.get('error')}),flush=True)
raise SystemExit(0 if receipt['status']=='passed' else 1)
