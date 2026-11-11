from pathlib import Path
import os,time,json,hashlib
r=Path(__file__).parent;leader=717466;proc=Path('/proc');records={};library_hashes={}
leader_start=(proc/str(leader)/'stat').read_text().rsplit(')',1)[1].split()[19]
while (proc/str(leader)/'stat').exists():
 try:
  if (proc/str(leader)/'stat').read_text().rsplit(')',1)[1].split()[19]!=leader_start:break
 except FileNotFoundError:break
 parents={}
 for folder in proc.iterdir():
  if not folder.name.isdigit():continue
  try:
   stat=(folder/'stat').read_text().rsplit(')',1)[1].split();parents[int(folder.name)]=int(stat[1])
  except (OSError,IndexError):pass
 family={leader}
 while True:
  more={pid for pid,parent in parents.items() if parent in family}-family
  if not more:break
  family.update(more)
 for pid in family:
  folder=proc/str(pid)
  try:
   command=(folder/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
   if not any(n in command.split(' ',1)[0] for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')):continue
   status=dict(line.split(':',1) for line in (folder/'status').read_text().splitlines() if ':' in line)
   row=records.setdefault(str(pid),dict(command=command,first_seen=time.time(),maximum_RSS_kB=0,maximum_HWM_kB=0,maximum_threads=0,actual_maps={}))
   row['last_seen']=time.time();row['maximum_RSS_kB']=max(row['maximum_RSS_kB'],int(status.get('VmRSS','0').split()[0]))
   row['maximum_HWM_kB']=max(row['maximum_HWM_kB'],int(status.get('VmHWM','0').split()[0]));row['maximum_threads']=max(row['maximum_threads'],int(status['Threads']))
   for line in (folder/'maps').read_text().splitlines():
    words=line.split()
    if len(words)<6 or not words[-1].startswith('/'):continue
    path=Path(words[-1])
    if path.is_file():
     digest=library_hashes.setdefault(str(path),None)
     if digest is None:library_hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
     row['actual_maps'][str(path)]=library_hashes[str(path)]
  except (OSError,KeyError,ValueError):pass
 (r/'actual-owned-runtime-maps-faces.json').write_text(json.dumps(dict(schema='OwnedRuntimeObservationV1',leader=leader,leader_start_ticks=leader_start,records=records,status='observing',scope='one-second read-only samples; RSS/HWM are observed maxima; process CPU is recorded separately per numerical job'),indent=2)+'\n')
 time.sleep(1)
(r/'actual-owned-runtime-maps-faces.json').write_text(json.dumps(dict(schema='OwnedRuntimeObservationV1',leader=leader,leader_start_ticks=leader_start,records=records,status='completed',scope='read-only samples, no debugger or process signals; CPU costs in job receipts'),indent=2)+'\n')
