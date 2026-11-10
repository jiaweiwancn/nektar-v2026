#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable request admission and hash-verified explicit reuse."""
import hashlib
import json
from pathlib import Path
from rans_qualification import write_json

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class RequestManifest:
    def __init__(self,directory,request,resume=False):
        self.directory=Path(directory);self.path=self.directory/'request-manifest.json'
        payload=dict(schema='RANS-ImmutableRequestV1',request=request)
        payload['sha256']=hashlib.sha256(canonical(payload)).hexdigest()
        if self.directory.exists():
            if not resume:raise FileExistsError('destination exists: explicit identical --resume required')
            if not self.path.exists():raise ValueError('existing destination has no immutable request manifest')
            stored=json.loads(self.path.read_text())
            if canonical(stored)!=canonical(payload):raise ValueError('immutable request mismatch; no receipts changed')
        else:
            if resume:raise ValueError('resume destination does not exist')
            self.directory.mkdir(parents=True)
            with self.path.open('x') as f:f.write(json.dumps(payload,indent=2)+'\n')
        self.payload=payload
    def completed(self,folder,receipt='completed-artifacts.json'):
        path=Path(folder)/receipt
        if not path.exists():return False
        record=json.loads(path.read_text())
        if record['request_sha256']!=self.payload['sha256']:raise ValueError('completed artifact request differs')
        for name,digest in record['files'].items():
            if sha(Path(folder)/name)!=digest:raise ValueError('completed numerical evidence changed: '+name)
        return True
    def seal(self,folder,files,receipt='completed-artifacts.json'):
        path=Path(folder)/receipt
        record=dict(request_sha256=self.payload['sha256'],files={str(p):sha(Path(folder)/p) for p in files})
        if path.exists():
            if json.loads(path.read_text())!=record:raise ValueError('refusing to replace completed numerical evidence')
        else:write_json(path,record)
