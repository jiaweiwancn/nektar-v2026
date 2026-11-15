#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Publish compact evidence while retaining complete sealed numerical bodies externally."""
import argparse,hashlib,json
from pathlib import Path

def compact(body):
    omitted=[]
    def visit(value,path):
        if isinstance(value,dict):return {k:visit(v,path+'/'+k) for k,v in value.items()}
        if isinstance(value,list):
            if len(value)>=10000 and all(type(x) in (int,float) for x in value):
                encoded=json.dumps(value,separators=(',',':'),allow_nan=False).encode();r=dict(entries=len(value),JSON_values_sha256=hashlib.sha256(encoded).hexdigest(),scope='complete numerical vector retained only in the hash-bound external body; publication view is not a checkpoint/replay input')
                omitted.append(dict(path=path,**r));return dict(external_numerical_vector=r)
            return [visit(v,path+'/'+str(i)) for i,v in enumerate(value)]
        return value
    return visit(body,''),omitted

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--external',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    raw=a.input.read_bytes();body,omitted=compact(json.loads(raw));body['compact_publication']=dict(schema='CT14CompactPublicationV1',complete_external_body=a.external,complete_body_sha256=hashlib.sha256(raw).hexdigest(),omitted_vectors=omitted,rule='no numerical primitive/classification changed; only large numerical vectors replaced by count/hash/external scope')
    a.output.write_text(json.dumps(body,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(output=str(a.output),bytes=a.output.stat().st_size,omitted=[dict(path=r['path'],entries=r['entries']) for r in omitted])))
if __name__=='__main__':main()
