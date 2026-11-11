#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent local Git-blob/tree proof for transported WJ-7 source records."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
from verify_wj6_transport import GitBlobs


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--repository',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    def read(name):
        f=a.directory/name
        return json.loads(f.read_bytes() if f.exists() else gzip.decompress(f.with_name(f.name+'.gz').read_bytes()))
    sources=read('source-test-lineage.json')['sources'];eq=read('native-equivalence.json')
    blobs=GitBlobs(a.repository);trees={};count=0
    try:
        for root,info in sources.items():
            tree=subprocess.check_output(['git','-C',str(a.repository),'rev-parse',info['revision']+'^{tree}'],text=True).strip()
            assert tree==info['tree'],'source Git tree differs';trees[info['revision']]=tree
            for path,identity in info['python_sources'].items():
                name=path[len(root)+1:]
                assert blobs.hash(info['revision'],name)==identity['sha256'],'source Python/Git bytes differ: '+name
                count+=1
        revision=sources['/media/jay/新加卷/nektar-v2026-runs/wj7-20261111/source-v8']['revision']
        for name,wanted in eq['native_inputs'].items():assert blobs.hash(revision,name)==wanted,'native source/Git bytes differ: '+name
        for packet,name in ((read('source-test-lineage.json'),'collect_wj7_lineage.py'),(read('final-seal.json'),'seal_export_wj7.py')):
            expected=packet['code']['sha256']
            assert blobs.hash('1a7af9e','solvers/IncNavierStokesSolver/Tests/RANS/'+name)==expected,'closing source/Git blob differs'
    finally:blobs.close()
    result=dict(schema='WJ7IndependentGitSourcesV1',status='PASS',verified_python_source_blobs=count,
        verified_native_source_blobs=len(eq['native_inputs']),verified_closing_source_blobs=2,trees=trees,
        transport_sha256=hashlib.sha256((a.directory/'transport.json').read_bytes()).hexdigest(),
        scope='local Git objects read independently; ELF/dependency/checkpoint bytes remain remote verification')
    assert not a.output.exists();a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
