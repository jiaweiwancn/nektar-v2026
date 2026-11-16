#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bound only owned stdout copies; scientific reports/checkpoints stay separate."""
import hashlib,os,threading
from pathlib import Path
from wj_h15_observer import atomic,identity


class BoundedLog:
    def __init__(self,path,chunk=262144,keep=2,prefix=16384):
        self.path=Path(path);self.chunk=chunk;self.keep=keep;self.prefix=prefix
        assert chunk>0 and keep>=2 and 0<=prefix<=chunk
        assert not any(self.path.parent.glob(self.path.name+'*'))
        self.file=self.path.open('xb');self.used=0;self.total=0;self.rotations=0
        self.hash=hashlib.sha256();self.first=bytearray();self.error=None;self.thread=None

    def attach(self,pipe):
        assert self.thread is None
        def consume():
            try:
                while True:
                    data=pipe.read(65536)
                    if not data:break
                    self.hash.update(data);self.total+=len(data)
                    self.first.extend(data[:max(0,self.prefix-len(self.first))])
                    while data:
                        if self.used==self.chunk:
                            self.file.close()
                            for n in range(self.keep,0,-1):
                                old=self.path if n==1 else self.path.with_name(self.path.name+'.'+str(n-1))
                                new=self.path.with_name(self.path.name+'.'+str(n))
                                if old.exists():old.replace(new)
                            self.file=self.path.open('xb');self.used=0;self.rotations+=1
                        part=data[:self.chunk-self.used];self.file.write(part);self.used+=len(part);data=data[len(part):]
            except Exception as e:self.error=repr(e)
            finally:pipe.close()
        self.thread=threading.Thread(target=consume,name='bounded-native-stdout',daemon=True);self.thread.start()

    def finish(self):
        if self.thread:self.thread.join(timeout=10);assert not self.thread.is_alive(),'stdout pump did not finish'
        self.file.flush();os.fsync(self.file.fileno());self.file.close()
        first=self.path.with_name(self.path.name+'.first');first.write_bytes(self.first)
        assert self.error is None,self.error
        retained=[identity(p)for p in sorted(self.path.parent.glob(self.path.name+'*'))if p.is_file()]
        assert sum(p['bytes']for p in retained)<=self.chunk*(self.keep+1)+self.prefix
        packet=dict(status='PASS_BOUNDED_STDOUT',chunk_bytes=self.chunk,previous_chunks_retained=self.keep,
            first_prefix_bytes=self.prefix,total_emitted_bytes=self.total,full_stream_sha256=self.hash.hexdigest(),
            rotations=self.rotations,retained=retained,
            scope='stdout trace only; complete native/independent scientific records, failure reports and recovery payloads have their own retention')
        atomic(self.path.with_name(self.path.name+'.retention.json'),packet);return packet
