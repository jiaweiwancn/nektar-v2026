#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Pause only this experiment's families when unrelated jobs consume headroom."""
from pathlib import Path
import argparse,os,subprocess,json,signal,time
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--project-limit',type=int);args=parser.parse_args()
r=args.root.resolve();started=time.time();paused=set();maximum=0;active_seconds={};last_clock=time.monotonic()
projectMaximum=0
if (r/'worker-guard-status.json').exists():maximum=json.loads((r/'worker-guard-status.json').read_text()).get('maximum_observed_active_before_actions',0)
initial=r/'resource-paused.json'
if initial.exists():paused.update(json.loads(initial.read_text())['paused_pids'])
while not (r/'worker-guard.stop').exists():
 rows={}
 for line in subprocess.check_output(['ps','-eo','pid,ppid,nlwp,stat,comm,args'],text=True).splitlines()[1:]:
  p=line.split(None,5)
  if len(p)==6:rows[int(p[0])]=dict(parent=int(p[1]),threads=int(p[2]),state=p[3],name=p[4],args=p[5])
 relevant=lambda d:d['name'].startswith(('IncNavier','RANS','Nek','cc1plus','f951','orted','mpiexec','gdb','lldb')) or (d['name'].startswith('python') and 'nektar' in d['args'])
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
 # Each complete serial family reserves eight threads, including ORTE and
 # its waiting Python parent. Keep sixteen of the 72 threads unallocated.
 allowed=max(0,min(4,(56-external)//8))
 if args.project_limit is not None:
  # This opt-in package permits one family; both waiting control helpers count.
  allowed=min(allowed,1)
  project=2+sum(d['threads'] for pid,d in rows.items() if pid in covered and not d['state'].startswith(('T','Z')))
  if project>projectMaximum:
   with (r/'project-highwater.jsonl').open('a') as f:
    f.write(json.dumps(dict(time_unix=time.time(),project_threads=project,control_helpers=2,
     processes=[dict(pid=pid,parent=d['parent'],threads=d['threads'],name=d['name'],state=d['state']) for pid,d in rows.items() if pid in covered and not d['state'].startswith(('T','Z'))]))+'\n')
  projectMaximum=max(projectMaximum,project)
  if project>args.project_limit:allowed=0
 chosen=set()
 for root,fam in families[:allowed]:chosen.update(fam)
 elapsed=time.monotonic()-last_clock;last_clock=time.monotonic()
 held_roots=set()
 for root,fam in families[:allowed]:
  live_children=[pid for pid in fam if pid!=root and pid in rows and not rows[pid]['state'].startswith('Z')]
  for pid in live_children:active_seconds[pid]=active_seconds.get(pid,0.)+(elapsed if not rows[pid]['state'].startswith('T') else 0.)
  # Let an interrupted native child finish before waking its waiting Python
  # parent: suspension must not consume a subprocess wall-time deadline.
  # A genuinely hung child still releases the parent's timeout after its
  # generous active-time allowance.
  limit=14400 if any(rows[p]['name'].startswith('IncNavier') for p in live_children) else 1800
  batch_children=[p for p in live_children if rows[p]['name'].startswith(('IncNavier','RANSPreserv','RANSWallResol'))]
  # Interactive coefficient probes and loader queries need their Python parent
  # to read/write pipes; never wait for those services to exit first.
  if rows[root]['state'].startswith('T') and batch_children and max(active_seconds.get(p,0.) for p in batch_children)<limit:held_roots.add(root)
 events=[]
 for root,fam in families:
  for pid in sorted(fam,reverse=True):
   if pid not in rows:continue
   if pid not in chosen and not rows[pid]['state'].startswith(('T','Z')):
    try:os.kill(pid,signal.SIGSTOP);paused.add(pid);events.append(['pause',root,pid])
    except ProcessLookupError:pass
   elif pid in chosen and pid not in held_roots and (pid in paused or rows[pid]['state'].startswith('T')):
    try:os.kill(pid,signal.SIGCONT);paused.discard(pid);events.append(['resume',root,pid])
    except ProcessLookupError:pass
 # Count every owned helper thread, irrespective of its executable's name.
 active=sum(d['threads'] for pid,d in rows.items() if (pid in covered or relevant(d)) and not d['state'].startswith(('T','Z')))
 if active>maximum:
  # Retain attribution even when no owned family is running and therefore no
  # pause event is emitted. Do not write unrelated command-line arguments.
  with (r/'worker-guard-highwater.jsonl').open('a') as f:
   f.write(json.dumps(dict(time_unix=time.time(),active_before_actions=active,
    external=external,owned_roots=[root for root,_ in families],
    processes=[dict(pid=pid,parent=d['parent'],threads=d['threads'],name=d['name'],state=d['state'],owned=pid in covered)
     for pid,d in rows.items() if (pid in covered or relevant(d)) and not d['state'].startswith(('T','Z'))]))+'\n')
 maximum=max(maximum,active)
 if events:
  with (r/'worker-guard-events.jsonl').open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),external=external,active_before_actions=active,allowed_families=allowed,events=events))+'\n')
 if int(time.time()-started)%10==0:(r/'worker-guard-status.json').write_text(json.dumps(dict(external=external,active_before_actions=active,allowed_families=allowed,paused_pids=sorted(paused),maximum_observed_active_before_actions=maximum,time_unix=time.time()))+'\n')
 time.sleep(.5)
for pid in paused:
 try:os.kill(pid,signal.SIGCONT)
 except ProcessLookupError:pass
if args.project_limit is not None:
 (r/'resource-final.json').write_text(json.dumps(dict(maximum_aggregate=maximum,maximum_project=projectMaximum,aggregate_limit=72,project_limit=args.project_limit,
  aggregate_status='pass' if maximum<=72 else 'fail',project_status='pass' if projectMaximum<=args.project_limit else 'fail',control_helpers_counted=2))+'\n')
