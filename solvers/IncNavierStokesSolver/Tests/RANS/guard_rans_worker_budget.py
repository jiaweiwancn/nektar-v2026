#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Pause only this experiment's families when unrelated jobs consume headroom."""
from pathlib import Path
import argparse,os,subprocess,json,signal,time
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);args=parser.parse_args()
r=args.root.resolve();started=time.time();paused=set();maximum=0
if (r/'worker-guard-status.json').exists():maximum=json.loads((r/'worker-guard-status.json').read_text()).get('maximum_observed_active_before_actions',0)
initial=r/'resource-paused.json'
if initial.exists():paused.update(json.loads(initial.read_text())['paused_pids'])
while not (r/'worker-guard.stop').exists():
 rows={}
 for line in subprocess.check_output(['ps','-eo','pid,ppid,nlwp,stat,comm,args'],text=True).splitlines()[1:]:
  p=line.split(None,5)
  if len(p)==6:rows[int(p[0])]=dict(parent=int(p[1]),threads=int(p[2]),state=p[3],name=p[4],args=p[5])
 relevant=lambda d:d['name'].startswith(('IncNavier','RANS','Nek','cc1plus','f951','orted','mpiexec')) or (d['name'].startswith('python') and 'nektar' in d['args'])
 try:leases=json.loads((r/'worker-leases.json').read_text()) if (r/'worker-leases.json').exists() else []
 except (json.JSONDecodeError,FileNotFoundError):time.sleep(.05);continue
 families=[];covered=set()
 for root in leases:
  if root not in rows or str(r) not in rows[root]['args']:continue
  fam={root}
  while True:
   more={pid for pid,d in rows.items() if d['parent'] in fam}-fam
   if not more:break
   fam.update(more)
  families.append((root,fam));covered.update(fam)
 external=sum(d['threads'] for pid,d in rows.items() if pid not in covered and relevant(d) and not d['state'].startswith(('T','Z')))
 allowed=4 if external<=16 else 3 if external<=24 else 2 if external<=40 else 1 if external<=48 else 0
 chosen=set()
 for root,fam in families[:allowed]:chosen.update(fam)
 events=[]
 for root,fam in families:
  for pid in sorted(fam,reverse=True):
   if pid not in rows:continue
   if pid not in chosen and not rows[pid]['state'].startswith(('T','Z')):
    try:os.kill(pid,signal.SIGSTOP);paused.add(pid);events.append(['pause',root,pid])
    except ProcessLookupError:pass
   elif pid in chosen and (pid in paused or rows[pid]['state'].startswith('T')):
    try:os.kill(pid,signal.SIGCONT);paused.discard(pid);events.append(['resume',root,pid])
    except ProcessLookupError:pass
 active=sum(d['threads'] for pid,d in rows.items() if relevant(d) and not d['state'].startswith(('T','Z')));maximum=max(maximum,active)
 if events:
  with (r/'worker-guard-events.jsonl').open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),external=external,active_before_actions=active,allowed_families=allowed,events=events))+'\n')
 if int(time.time()-started)%10==0:(r/'worker-guard-status.json').write_text(json.dumps(dict(external=external,active_before_actions=active,allowed_families=allowed,paused_pids=sorted(paused),maximum_observed_active_before_actions=maximum,time_unix=time.time()))+'\n')
 time.sleep(.5)
for pid in paused:
 try:os.kill(pid,signal.SIGCONT)
 except ProcessLookupError:pass
