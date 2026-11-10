# SPDX-License-Identifier: MIT
"""Write-once intent and package-wide reservations for bounded native studies.

A stale artifact is an error, never permission to rerun/overwrite an attempt.
Unreconciled reservations remain charged after crashes. Measurements have their
own versioned lineage and charge wall time, but reserve zero scheduled steps.
"""
import contextlib
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
import xml.etree.ElementTree as ET


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':')).encode()).hexdigest()

def file_hash(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def read(path):return json.loads(Path(path).read_text())

def atomic(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.tmp-'+str(os.getpid()))
    with tmp.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

@contextlib.contextmanager
def locked(path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX);yield

def write_once(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with locked(path.with_suffix(path.suffix+'.lock')):
        if path.exists():
            if read(path)!=value:raise ValueError('immutable intent/design mismatch: '+str(path)+'; use a new experiment root')
        else:atomic(path,value)
    return value

def runtime(path):
    path=Path(path).resolve();libs={}
    if __import__('platform').system()=='Linux':
        text=subprocess.check_output(['ldd',str(path)],text=True)
        for name in re.findall(r'(?:=>\s+)?(/[^\s()]+)',text):
            p=Path(name).resolve()
            if p.is_file():libs[str(p)]=file_hash(p)
    return dict(path=str(path),sha256=file_hash(path),libraries=libs)

def xml_identity(path):
    tree=ET.parse(path).getroot()
    # Explicit semantic tree in addition to exact original XML bytes.
    def node(e):return [e.tag,dict(sorted(e.attrib.items())),(e.text or '').strip(),[node(c) for c in e]]
    return dict(sha256=file_hash(path),semantic_sha256=digest(node(tree)))

def file_identity(path):
    p=Path(path).resolve()
    return dict(path=str(p),bytes=p.stat().st_size,sha256=file_hash(p))

def artifacts(paths):return {str(Path(p).resolve()):file_identity(p) for p in paths if Path(p).is_file()}

def verify_artifacts(index):
    for p,identity in index.items():
        if not Path(p).is_file() or file_identity(p)!=identity:raise ValueError('immutable artifact mismatch: '+p)

def execution_intent(argv,inputs,code,runtime_identity,settings):
    return dict(argv=list(map(str,argv)),inputs=artifacts(inputs),code=code,runtime=runtime_identity,settings=settings)

class Ledger:
    def __init__(self,root,policy,cap=1500,wall_cap=10800):
        self.file=Path(root)/('ledger-'+policy+'.json');self.cap=cap;self.wall_cap=wall_cap
        self.file.parent.mkdir(parents=True,exist_ok=True)
        with self.lock():
            if not self.file.exists():atomic(self.file,dict(schema='WJWorkLedgerV1',policy=policy,attempt_cap=cap,wall_cap_seconds=wall_cap,entries={}))
            d=read(self.file)
            if d['attempt_cap']!=cap or d['wall_cap_seconds']!=wall_cap:raise ValueError('ledger limit mismatch')
    def lock(self):return locked(self.file.with_suffix('.lock'))
    def begin_session(self,label,started_epoch=None):
        with self.lock():
            d=read(self.file);sessions=d.setdefault('sessions',{})
            if label in sessions:raise ValueError('session already exists')
            sessions[label]=dict(started_epoch=time.time() if started_epoch is None else started_epoch,status='RUNNING')
            atomic(self.file,d)
        self.session=label
    def finish_session(self,label):
        with self.lock():
            d=read(self.file);e=d['sessions'][label]
            if e['status']!='RUNNING':raise ValueError('session already closed')
            e.update(seconds=time.time()-e['started_epoch'],status='FINISHED');atomic(self.file,d)
    def totals(self,d=None):
        d=d or read(self.file);entries=list(d['entries'].values())
        command_seconds=sum(e.get('seconds',max(e['reserved_seconds'],time.time()-e['started_epoch'])) for e in entries)
        # Process-session wall time includes Python algebra, quadrature, SVD,
        # report/checkpoint I/O and cache verification. Command times within a
        # session are not counted twice. Crashed sessions remain conservatively
        # charged until explicitly reconciled with an immutable process receipt.
        sessions=list(d.get('sessions',{}).values())
        seconds=sum(e.get('seconds',max(0,time.time()-e['started_epoch'])) for e in sessions)
        seconds+=sum(e.get('seconds',max(e['reserved_seconds'],time.time()-e['started_epoch'])) for e in entries if not e.get('session'))
        return dict(attempts=sum(e.get('actual_attempts',e['reserved_attempts']) for e in entries),
            seconds=seconds,command_seconds=command_seconds)
    def reserve(self,label,intent,attempts,seconds):
        if attempts<0 or int(attempts)!=attempts or seconds<=0:raise ValueError('invalid reservation')
        with self.lock():
            d=read(self.file)
            if label in d['entries']:raise ValueError('immutable attempt already reserved: '+label)
            totals=self.totals(d)
            if totals['attempts']+attempts>self.cap or totals['seconds']+seconds>self.wall_cap:raise RuntimeError('RESOURCE_LIMIT: package work reservation exceeds cap')
            d['entries'][label]=dict(intent_sha256=digest(intent),reserved_attempts=attempts,reserved_seconds=seconds,started_epoch=time.time(),status='RUNNING',session=getattr(self,'session',None))
            atomic(self.file,d)
    def finish(self,label,attempts,seconds,receipt):
        with self.lock():
            d=read(self.file);e=d['entries'][label]
            if e['status']!='RUNNING':raise ValueError('immutable attempt already reconciled')
            if attempts is not None and attempts>e['reserved_attempts']:raise ValueError('native attempts exceeded reservation')
            if attempts is not None:e['actual_attempts']=attempts
            e.update(observed_attempts=attempts,seconds=seconds,status='FINISHED' if attempts is not None else 'FINISHED_UNKNOWN_ATTEMPTS',receipt=receipt);atomic(self.file,d)
