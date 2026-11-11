#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independently verify transported bytes, receipt links and Git source blobs.

Linux heavy fields/binaries stay remote. This verifier distinguishes those
recorded/remote-verified identities from bytes actually transported and read.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import posixpath
import subprocess
import time


def require(ok, message):
    if not ok:raise ValueError(message)


class GitBlobs:
    def __init__(self, repository):
        self.repo=repository
        self.process=subprocess.Popen(['git','-C',str(repository),'cat-file','--batch'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        self.cache={}
    def data(self,revision,name):
        key=revision+':'+name
        if key not in self.cache:
            self.process.stdin.write(key.encode()+b'\n');self.process.stdin.flush()
            header=self.process.stdout.readline().split()
            require(len(header)==3 and header[1]==b'blob','missing Git blob: '+key)
            size=int(header[2]);data=self.process.stdout.read(size)
            require(len(data)==size and self.process.stdout.read(1)==b'\n','truncated Git blob')
            # Cache digests separately below; release large source payloads.
            return data
        return self.cache[key]
    def hash(self,revision,name):
        return hashlib.sha256(self.data(revision,name)).hexdigest()
    def close(self):
        self.process.stdin.close();require(self.process.wait(timeout=10)==0,'Git reader failed')


def identities(value):
    if isinstance(value,dict):
        if isinstance(value.get('path'),str) and isinstance(value.get('sha256'),str):yield value
        for child in value.values():yield from identities(child)
    elif isinstance(value,list):
        for child in value:yield from identities(child)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('directory','repository','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();start=time.time();root=a.directory.resolve()
    transport=json.loads((root/'transport.json').read_bytes());decoded={};packets={};encoded_bytes=raw_bytes=0
    for name,row in transport['files'].items():
        path=root/name;encoded=path.read_bytes()
        require(len(encoded)==row['transport_bytes'] and hashlib.sha256(encoded).hexdigest()==row['transport_sha256'],'transport size/hash differs: '+name)
        data=gzip.decompress(encoded) if row['encoding']=='gzip mtime0' else encoded
        require(len(data)==row['original']['bytes'] and hashlib.sha256(data).hexdigest()==row['uncompressed_sha256']==row['original']['sha256'],'decoded bytes differ: '+name)
        decoded[row['original']['path']]=dict(sha256=row['uncompressed_sha256'],bytes=len(data),local=name)
        if name.removesuffix('.gz').endswith('.json'):packets[name.removesuffix('.gz')]=json.loads(data)
        encoded_bytes+=len(encoded);raw_bytes+=len(data)
    checked=set();remote=set()
    for packet in packets.values():
        for identity in identities(packet):
            path=identity['path']
            if path in decoded:
                require(identity['sha256']==decoded[path]['sha256'] and ('bytes' not in identity or identity['bytes']==decoded[path]['bytes']),
                        'transported producer/observer/receipt link differs: '+path)
                checked.add(path)
            else:remote.add((path,identity['sha256']))
    manifest=packets['manifest.json'];blobs=GitBlobs(a.repository);source_count=archive_count=0
    try:
        source=packets['identity/producer-source-content.json'];revision=manifest['native_source_revision']
        raw=subprocess.check_output(['git','-C',str(a.repository),'ls-tree','-r','-z',revision])
        modes={item.split(b'\t',1)[1].decode():item.split(None,1)[0].decode() for item in raw.split(b'\0') if item}
        for name,wanted in source.items():
            actual=blobs.hash(revision,name)
            if actual!=wanted and modes.get(name)=='120000':
                target=blobs.data(revision,name).decode();target=posixpath.normpath(posixpath.join(posixpath.dirname(name),target))
                require(not target.startswith(('/','../')),'external symlink source cannot be verified locally')
                actual=blobs.hash(revision,target)
            require(actual==wanted,'producer source/Git blob differs: '+name);source_count+=1
        for key,packet in packets.items():
            if key.startswith('identity/observer-') and key.endswith('-source.json'):
                for name,wanted in packet['members'].items():
                    require(blobs.hash(packet['revision'],name)==wanted,'observer code/Git blob differs: '+name);archive_count+=1
    finally:blobs.close()
    result=dict(schema='WJ6IndependentLocalVerificationV1',status='PASS',transport_sha256=hashlib.sha256((root/'transport.json').read_bytes()).hexdigest(),
        transported_files=len(transport['files']),encoded_bytes=encoded_bytes,decoded_bytes=raw_bytes,
        verified_receipt_links=len(checked),producer_git_blobs=source_count,observer_git_blobs=archive_count,
        remote_only_identity_count=len(remote),seconds=time.time()-start,
        scope='all transported encoded/decoded bytes and available receipt links read locally; producer/observer source checked against local Git objects; heavy fields/binaries remain remotely verified recorded identities')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
