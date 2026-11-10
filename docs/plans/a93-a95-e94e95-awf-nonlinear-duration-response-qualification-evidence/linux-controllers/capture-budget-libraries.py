from pathlib import Path
import os,time,json,hashlib
r=Path(__file__).resolve().parent;remaining={'runtime-budget-v4'}
while remaining:
 for folder in list(Path('/proc').iterdir()):
  if not folder.name.isdigit():continue
  try:exe=Path(os.readlink(folder/'exe'))
  except (OSError,FileNotFoundError):continue
  name=next((n for n in remaining if exe==r/n/('IncNavierStokesSolver' if n=='runtime' else 'RANSPreservationProbe')),None)
  if name is None:continue
  try:
   paths=sorted({x.split()[-1] for x in (folder/'maps').read_text().splitlines() if len(x.split())>=6 and x.split()[-1].startswith('/') and '.so' in x.split()[-1]})
   libs={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths if Path(p).is_file()}
   d=dict(schema='ActualLoadedLibrariesV1',pid=int(folder.name),executable=str(exe),executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),libraries=libs,time_unix=time.time(),environment={a.split('=',1)[0]:a.split('=',1)[1] for a in (folder/'environ').read_bytes().decode().split('\x00') if '=' in a and a.split('=',1)[0] in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')})
   dest=r/name/'actual-loaded-libraries-current.json';assert not dest.exists();dest.write_text(json.dumps(d,indent=2)+'\n');remaining.remove(name);print(name,len(libs),flush=True)
  except (OSError,ProcessLookupError):pass
 time.sleep(1)
