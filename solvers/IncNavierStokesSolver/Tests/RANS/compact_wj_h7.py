#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Byte-verify a Linux publication, then omit bulk boundary/station row copies.

Metadata-only export: numerical results and heavy native files are unchanged.
The original seal/transport and each omitted array's containing packet remain
bound explicitly. This command performs no closure evaluation or native solve.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha(b):return hashlib.sha256(b).hexdigest()
def load(p):return json.loads(p.read_text(),parse_int=lambda v:-0.0 if v=='-0'else int(v))
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def compact(source,out):
    seal=load(source/'seal.json');out.mkdir(parents=True,exist_ok=False);changed=[]
    for item in seal['files']:
        src=source/item['name'];raw=src.read_bytes();assert len(raw)==item['bytes']and sha(raw)==item['sha256']
        dst=out/item['name'];dst.parent.mkdir(parents=True,exist_ok=True)
        if src.name.startswith('endpoint-H'):
            d=load(src);omitted={}
            for key in ('boundary_rows','stations'):
                value=d['open_boundary'].pop(key);omitted[key]=dict(count=len(value),original_sealed_packet=item,native_source=d['native'])
            d['bulk_rows_retained_remotely']=omitted;write(dst,d);changed.append(dict(original=item,compact=dict(name=item['name'],bytes=dst.stat().st_size,sha256=sha(dst.read_bytes()))))
        else:dst.write_bytes(raw)
    shutil.copyfile(source/'seal.json',out/'parent-publication-seal.json')
    write(out/'compact-bridge.json',dict(parent_seal_sha256=sha((source/'seal.json').read_bytes()),changed=changed,
        reason='keep complete boundary/station inventories in native heavy packets, preserve all numerical decisions, weak loads, interface rows and metrics',numerical_reruns=0,metadata_source_sha256=sha(Path(__file__).read_bytes())))
    write(out/'seal.json',dict(files=[dict(name=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=sha(p.read_bytes()))for p in sorted(out.rglob('*'))if p.is_file()and p.name!='seal.json']))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path);a=p.parse_args();compact(a.source,a.output)
