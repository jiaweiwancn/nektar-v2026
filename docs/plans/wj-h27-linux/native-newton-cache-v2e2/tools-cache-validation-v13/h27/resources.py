# SPDX-License-Identifier: MIT
"""Small standalone Linux resource/identity helpers; no old campaign imports."""
import fcntl,hashlib,json,os,subprocess
from pathlib import Path

THREADS={k:'1'for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}
def identity(path):
    p=Path(path);h=hashlib.sha256()
    with p.open('rb')as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return dict(path=str(p.resolve()),bytes=p.stat().st_size,sha256=h.hexdigest())
def atomic(path,value,replace=False):
    path=Path(path)
    if path.exists()and not replace:raise FileExistsError(path)
    tmp=path.with_name(path.name+f'.tmp-{os.getpid()}')
    with tmp.open('x')as f:json.dump(value,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path);fd=os.open(path.parent,os.O_RDONLY);os.fsync(fd);os.close(fd)
def update(path,change):
    path=Path(path)
    with path.with_name(path.name+'.lock').open('a')as f:
        fcntl.flock(f,fcntl.LOCK_EX);value=json.loads(path.read_text());change(value);atomic(path,value,True)
    return value
def workers():
    rows=[]
    for line in subprocess.check_output(['ps','-eo','pid=,ppid=,nlwp=,args='],text=True).splitlines():
        parts=line.strip().split(None,3)
        if len(parts)!=4:continue
        pid,parent,n,args=parts;command=Path(args.split()[0]).name
        if command.startswith(('IncNavierStokes','RANS','WJ','python','cc1','c++','g++','gcc','clang','gfortran','f951','lto1','orted','mpirun','mpiexec','ninja','make','cmake'))or('ld-linux'in command and'.elf'in args):rows.append(dict(pid=int(pid),parent=int(parent),threads=int(n),args=args))
    return dict(workers=sum(r['threads']for r in rows),processes=rows)
def tree(pid):
    todo=[pid];seen=set();rss=0;groups=set()
    while todo:
        p=todo.pop()
        if p in seen:continue
        seen.add(p)
        try:
            data=dict(line.split(':',1)for line in(Path('/proc')/str(p)/'status').read_text().splitlines()if':'in line);rss+=int(data.get('VmRSS','0 kB').split()[0])*1024;groups.add(os.getpgid(p));todo.extend(map(int,(Path('/proc')/str(p)/'task'/str(p)/'children').read_text().split()))
        except (OSError,ProcessLookupError):pass
    return rss,groups
def size(path):
    p=Path(path)
    if p.is_file():return p.stat().st_size
    total=0
    for f in p.rglob('*'):
        try:
            if f.is_file():total+=f.stat().st_size
        except FileNotFoundError:pass
    return total
