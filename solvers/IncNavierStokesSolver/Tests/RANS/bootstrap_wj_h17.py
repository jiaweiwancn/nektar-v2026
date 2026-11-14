#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One new H17 root/ledger; selected immutable inputs, no native execution."""
import json,os,shutil,subprocess,time,hashlib
from pathlib import Path
ROOT=Path('/home/jay/build/wj-h17-exact-row-reuse-20261114')
OLD=Path('/home/jay/build/wj-h16-repaired-monitor-20261114')
def pin(p):
 p=Path(p);h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest())
def save(p,v):
 p=Path(p);t=p.with_suffix('.tmp');t.write_text(json.dumps(v,indent=2)+'\n')
 with t.open('rb')as f:os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ROOT.mkdir()
 for d in ('operations','science','inputs','tmp'): (ROOT/d).mkdir()
 b=dict(schema='WJH17LedgerV1',root=str(ROOT),closed=False,ceiling_seconds=7200,setup_charged_seconds=60,closure_reserve_seconds=300,starting_revision='066a1224cc1ed627bf2eaecdc780668a1cf68c11',small_charged_attempts=0,target_reserved_attempts=0,full_audit_reserved=0,full_audit_policies=[],entries={'selected-bootstrap-v1':dict(status='RESERVED',phase='operation',reserved_seconds=120)})
 save(ROOT/'accounting.json',b);start=time.monotonic()
 try:
  fs=subprocess.check_output(['findmnt','-n','-o','FSTYPE,SOURCE,TARGET','-T',str(ROOT)],text=True).strip();assert fs.startswith('ext4 ')
  assert shutil.disk_usage(ROOT).free>12*1024**3
  p=ROOT/'tmp/write-readback.bin';data=bytes(range(256))*4096
  with p.open('xb')as f:f.write(data);f.flush();os.fsync(f.fileno())
  assert p.read_bytes()==data
  index=json.loads((OLD/'authoritative-index-v1.json').read_text());selected=[]
  for key in ('publication','closure','accounting'):
   v=index[key];assert pin(v['path'])==v;selected.append(v)
  oldbook=json.loads(Path(index['accounting']['path']).read_text());assert oldbook['closed']and oldbook['total_charged_seconds']==1265.2252584086964
  for v in json.loads((OLD/'publication-v1/controls/runtime-build-v1.json').read_text())['executables']:
   assert pin(v['ELF']['path'])==v['ELF'];assert pin(v['wrapper']['path'])==v['wrapper'];selected.extend((v['ELF'],v['wrapper']))
  # Selected original session; neither an input regeneration nor a root traversal.
  xml=Path('/home/jay/build/wj-h13-fullheight-paired-domain-20261113/publication-v1/paired-case/flat/case.xml');v=pin(xml);assert v['sha256']=='137566f8704fed63d881dfbb5b804a02776a09c2081b9f29c42cadd91e794511';selected.append(v)
  cache=pin('/home/jay/build/wj-h15-persistent-neutral-20261114/science/observer-cache/immutable-cache.pickle');assert cache['bytes']==211582940 and cache['sha256']=='800ac28c73ecb6498e6504b4e14e774db2cff893f2c0aaf418cb7036d50d00b0';selected.append(cache)
  shutil.copytree(OLD/'source',ROOT/'source')
  save(ROOT/'inputs/selected-bootstrap.json',dict(status='PASS',filesystem=fs,free_bytes=shutil.disk_usage(ROOT).free,selected=selected,H16_ledger_unchanged=index['accounting'],historical_source=pin(OLD/'source-intent.json'),native_work=0))
  e=dict(status='PASS',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start,native_work=0)
 except Exception as error:
  e=dict(status='OBSTRUCTED',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start,reason=repr(error),native_work=0);b['entries']['selected-bootstrap-v1']=e;save(ROOT/'accounting.json',b);raise
 b['entries']['selected-bootstrap-v1']=e;save(ROOT/'accounting.json',b);save(ROOT/'operations/bootstrap-receipt.json',e)
 print(json.dumps(dict(status='PASS',root=str(ROOT),seconds=e['seconds'],selected=pin(ROOT/'inputs/selected-bootstrap.json'))),flush=True)
if __name__=='__main__':main()
