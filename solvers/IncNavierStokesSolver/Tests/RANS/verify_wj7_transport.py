#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent byte and semantic-parent verification; standard library only."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time


def require(ok,message):
    if not ok:raise ValueError(message)


def parse(raw):return json.loads(raw,parse_int=float)
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':')).encode()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True);p.add_argument('--parent',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    start=time.time();root=a.directory;transport=parse((root/'transport.json').read_bytes());packets={};encoded=decoded=0;links={}
    for name,info in transport['files'].items():
        raw=(root/name).read_bytes();require(len(raw)==info['transport_bytes'] and hashlib.sha256(raw).hexdigest()==info['transport_sha256'],'encoded transport differs: '+name)
        encoded+=len(raw);original=gzip.decompress(raw) if info['encoding'].startswith('gzip') else raw
        require(len(original)==info['original']['bytes'] and hashlib.sha256(original).hexdigest()==info['uncompressed_sha256']==info['original']['sha256'],'decoded transport differs: '+name)
        decoded+=len(original);links[info['original']['path']]=info['original']
        logical=name[:-3] if name.endswith('.gz') else name
        if logical.endswith('.json'):packets[logical]=parse(original)
    correction=packets['correction.json'];oldmap=parse(gzip.decompress((a.parent/'accepted-stage-lineage.json.gz').read_bytes()))
    expected={}
    for policy in ('EP','Budget'):
        prior=parse(gzip.decompress((a.parent/(policy+'.json.gz')).read_bytes()))
        for name,e in prior['execution'].items():
            if e.get('report'):expected[e['report']['path']]=(policy,name,e['report'])
    require(set(expected)==set(correction['entries']) and len(expected)==36,'semantic correction report coverage differs')
    changed=0
    for path,c in correction['entries'].items():
        policy,name,report=expected[path]
        require(c['report']==report and c['policy']==policy and c['case']==name,'correction does not bind sealed execution producer')
        old=oldmap[c['original_key']]
        require(digest(old)==c['original_stage_packet_sha256'],'original packet byte-decoded semantics differ')
        for key,wanted in c['stage_fields_sha256'].items():require(digest(old[key])==wanted,'retained callback/stage field differs: '+key)
        changed+=c['corrected']
    require(changed==18 and correction['corrected_count']==18,'corrected matrix/report count differs')
    stage=packets['accepted-stage-lineage.json'];native_reports={}
    for policy in ('EP','Budget'):
        for name,e in packets[policy+'-executions.json'].items():
            if e.get('report'):native_reports[e['report']['path']]=e['report']
    require(set(stage)==set(native_reports),'new report-role coverage differs')
    for path,s in stage.items():require(s['report']==native_reports[path],'new stage identity is not its execution report')
    def identities(value):
        if isinstance(value,dict):
            if all(k in value for k in ('path','sha256','bytes')):yield value
            for v in value.values():yield from identities(v)
        elif isinstance(value,list):
            for v in value:yield from identities(v)
    verified=set();remote=set()
    for packet in packets.values():
        for identity in identities(packet):
            if identity['path'] in links:
                require(identity==links[identity['path']],'transported receipt/source link differs')
                verified.add(identity['path'])
            else:remote.add((identity['path'],identity['sha256']))
    require(not a.output.exists(),'independent verification output must be new')
    result=dict(schema='WJ7IndependentLocalVerificationV1',status='PASS',transport_sha256=hashlib.sha256((root/'transport.json').read_bytes()).hexdigest(),
        files=len(transport['files']),encoded_bytes=encoded,decoded_bytes=decoded,verified_receipt_links=len(verified),
        independently_bound_original_reports=36,corrected_report_roles=18,new_execution_reports=len(stage),remote_only_identities=len(remote),seconds=time.time()-start,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='all transported bytes read locally; corrected roles independently compared with sealed local WJ6 execution records and retained stage fields; raw native reports/checkpoints/matrices remain recorded remote verification')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
