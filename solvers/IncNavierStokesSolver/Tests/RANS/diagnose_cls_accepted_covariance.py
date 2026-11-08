#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict zero-PDE archived PIMPLE stages, native identities and tensor oracle."""
import argparse
from copy import deepcopy
from datetime import datetime,timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import struct
import subprocess
import time

from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,write_json
from continue_cls_epawf_pimple import atomic

_created_output=None
STAGES=('committed','post_flow','pre_k','raw','corrected')
CASES={'refined':'selected-refined-2G','coarse':'faithful-2G-campaign',
       'one':'selected-regularized-saturated-05O'}


def checkpoint(path):
    """Decode only; the native reader independently checks checksum/signature."""
    payload=(path/'State.bin').read_bytes();cursor=24
    def value(fmt):
        nonlocal cursor
        v=struct.unpack_from('<'+fmt,payload,cursor);cursor+=struct.calcsize('<'+fmt);return v[0] if len(v)==1 else list(v)
    def vector():
        nonlocal cursor
        n=value('Q');result=list(struct.unpack_from('<'+'d'*n,payload,cursor));cursor+=8*n;return result
    n=value('Q');signature=payload[cursor:cursor+n].decode();cursor+=n
    order,steps,t,dt,eta=[value(fmt) for fmt in ('Q','Q','d','d','d')]
    state=vector();slots=value('Q');history=[vector() for _ in range(slots)];times=vector();coupled=vector();reference=vector();defect=vector()
    require(cursor==len(payload)-8,'unexpected native checkpoint layout')
    return dict(signature=signature,order=order,steps=steps,time=t,dt=dt,eta=eta,differential_state=state,
        history=history,history_times=times,coupled=coupled,reference=reference,defect=defect)


def read_report(path):
    if path.suffix=='.gz':
        with gzip.open(path,'rt') as f:return json.load(f,parse_int=float)
    return load(path)


def common_packet(cp):
    c=checkpoint(cp)
    return dict(native_checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),native_signature=c['signature'],
        native_differential_state=c['differential_state'],native_history=c['history'],native_history_times=c['history_times'])


def run_probe(probe,xml,packet,folder,success=True):
    folder.mkdir(exist_ok=False);write_json(folder/'input.json',packet)
    env=os.environ.copy();env.update(THREADS);env['LD_LIBRARY_PATH']=str(probe.parent/'lib')
    command=[str(probe),str(xml),'input.json','result.json'];start=time.monotonic()
    with (folder/'probe.log').open('w') as log:r=subprocess.run(command,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    receipt=dict(command=command,exit_code=r.returncode,seconds=time.monotonic()-start,
        input_sha256=sha(folder/'input.json'),log_sha256=sha(folder/'probe.log'),output=str(folder/'result.json'))
    if (folder/'result.json').exists():receipt['output_sha256']=sha(folder/'result.json')
    write_json(folder/'receipt.json',receipt)
    require((r.returncode==0)==success,'unexpected diagnostic outcome '+str(folder)+': '+(folder/'probe.log').read_text()[-1500:])
    return (load(folder/'result.json') if success else (folder/'probe.log').read_text()),receipt


def value(row,state,offset=0,gradient=False):
    import numpy as np
    ids=np.array(row['ids'],dtype=int)+offset;x=np.array(state)[ids];constant=np.array(row['constant'])
    # Independently center on the element constant by least squares (the
    # native evaluator uses its first vertex). These evaluate the same CG
    # polynomial; the different contraction supplies an arithmetic cross-check.
    base=float(x@constant/(constant@constant));centered=x-base*constant
    if gradient:return np.array(row['derivative'])@centered
    return float(np.array(row['phi'])@centered+base)


def cls_tensor(k,epsilon,g):
    import numpy as np
    S=(g+g.T)*k/epsilon;W=(g-g.T)*k/epsilon
    ss=float(np.sum(S*S));oo=float(np.sum(W*W));q=max((ss/2)**.5,(oo/2)**.5)
    C=.3*(1. if q>8 else -__import__('math').expm1(-.36*__import__('math').exp(.75*q)))/(1+.35*q**1.5)
    I=np.eye(3);quad=k*C*(-.1*(S@S-ss*I/3)+.1*(W@S-S@W)+.26*(-W@W-oo*I/3))
    cubic=k*C*(-10*C*C*((S@S)@W-W@(S@S))+(-5*C*C*ss+5*C*C*oo)*S)
    return 2*k/3*I-k*C*S+quad+cubic


def independent_witnesses(result,archived):
    import numpy as np
    checks=[]
    for stage,d in result['diagnostics'].items():
        state=archived[stage];offset=[int(x) for x in d['offsets']]
        points=[d['minimum_outer'],*d['one_sided_witnesses']]
        for p in points:
            if p['inner']:continue
            g=np.zeros((3,3));g[0]=value(p['velocity_evaluation_row'],state,offset[0],True);g[2]=value(p['velocity_evaluation_row'],state,offset[1],True)
            k=value(p['velocity_evaluation_row'],state,offset[2]);eps=value(p['epsilon_evaluation_row'],state,offset[3])
            require(abs(k-p['k'])<2e-12 and np.max(abs(g-np.array(p['velocity_gradient'])))<2e-11,'independent native-row reconstruction differs')
            # Current dependent target is an explicit audit quantity, rather
            # than the retained trace; its semantics remain visible.
            e=p['epsilon'] if p['epsilon_semantics'].startswith('refreshed') else eps
            R=cls_tensor(k,e,g);native=np.array(p['same_state_Original']['tensor'])
            error=float(np.max(abs(R-native)));require(error<2e-11,'independent Original tensor differs')
            eigen=np.linalg.eigvalsh(R);ev=p['same_state_Original']['eigenvalues'];require(np.max(abs(eigen-ev))<2e-11,'independent tensor spectrum differs')
            selected=np.array(p['selected']['tensor']);direction=np.array(p['selected']['critical_eigenvector'])
            require(np.max(abs(selected@direction-p['selected']['eigenvalues'][0]*direction))<2e-11,'critical eigenvector residual')
            checks.append(dict(stage=stage,geometry=p['geometry'],coordinate=p['coordinate'],tensor_error=error,eigenvalue=float(eigen[0])))
    return checks


def main():
    global _created_output
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',required=True,type=Path)
    p.add_argument('--old-root',type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--native-report',type=Path);p.add_argument('--native-checkpoint',type=Path)
    p.add_argument('--case-label',default='focused-control')
    p.add_argument('--cases',nargs='+',choices=CASES,default=list(CASES));p.add_argument('--contract-tests',action='store_true')
    a=p.parse_args();require(platform.system()=='Linux','authoritative diagnostics run on Linux');os.environ.update(THREADS)
    custom=a.native_report is not None
    require((custom and a.native_checkpoint is not None and a.old_root is None) or
        (not custom and a.native_checkpoint is None and a.old_root is not None),
        'select either archived campaign root or explicit native report plus compatible checkpoint')
    require(Path(a.case_label).name==a.case_label and a.case_label not in ('.','..'),'invalid diagnostic label')
    a.probe=a.probe.resolve();a.output=a.output.resolve();a.output.mkdir(exist_ok=False);_created_output=a.output
    result=dict(schema=1,started_UTC=datetime.now(timezone.utc).isoformat(),status='started',probe=str(a.probe),probe_sha256=sha(a.probe),
        production_governing_solves=0,diagnostic_governing_solves=0,physical_commits=0,threads=THREADS,cases={})
    atomic(a.output/'summary.json',result)
    for label in ([a.case_label] if custom else a.cases):
        root=a.output/label;root.mkdir()
        if custom:
            reportfile=a.native_report.resolve();segment=reportfile.parent;accepted=a.native_checkpoint.resolve()
            accepted_hash=sha(accepted/'State.bin')
        else:
            progress=load(a.old_root/CASES[label]/'trajectory/progress.json');segment=Path(progress['segment'])
            reportfile=segment/'case.evolution.json'
            if not reportfile.exists():reportfile=reportfile.with_suffix('.json.gz')
            accepted=Path(progress['checkpoint']);accepted_hash=progress['checkpoint_sha256']
        native=read_report(reportfile);attempt=native['pimple_attempts'][-1];rawfile=Path(attempt['k_snapshot'])
        require(native['failure'],'archived terminal-stage diagnostic needs the actual failed native report')
        require(sha(rawfile)==attempt['k_snapshot_sha256'],'original terminal packet hash differs')
        cp=root/'checkpoint';cp.mkdir();shutil.copy2(accepted/'State.bin',cp/'State.bin')
        require(sha(cp/'State.bin')==accepted_hash,'accepted checkpoint identity differs')
        shutil.copy2(rawfile,root/'archived.json');archived=load(root/'archived.json');c=checkpoint(cp)
        require(state_hash(c['coupled'])==state_hash(archived['committed']),'packet committed state differs from accepted native checkpoint')
        require(state_hash(c['coupled'])==state_hash(native['state']) and c['steps']==native['steps'] and c['time']==native['time'],
            'rolled-back native report does not identify the selected accepted checkpoint')
        require([state_hash(h) for h in c['history']]==[state_hash(h) for h in native['native_history']],
            'rolled-back native report history differs from accepted checkpoint')
        witness=native['physical_range_events'][-1]['audit']['accepted_awf_outer_minimum_point'][:3]
        packet=dict(common_packet(cp),operation='archived_stages',archived_packet=str(root/'archived.json'),
            archived_packet_sha256=sha(root/'archived.json'),stage_hashes={n:state_hash(archived[n]) for n in STAGES},
            expected_history_sha256=state_hash(archived['history']),state=archived['committed'],committed=archived['committed'],
            history=archived['history'],fixed_reference=archived['fixed_reference'],fixed_external_coefficients=archived['fixed_external_coefficients'],
            alpha=archived['alpha'],time=archived['time'],inner_stress_model=native['inner_stress_model'],witnesses=[witness])
        xml=segment/'case.xml';data,receipt=run_probe(a.probe,xml,packet,root/'stages')
        require(data['governing_solves']==data['physical_commits']==0 and not data['native_history_advanced'] and data['initialization_work']==[0,0],'forensic operation advanced or lifted state')
        checks=independent_witnesses(data,archived)
        rows={name:dict(inner_minimum=d['minimum_inner']['selected']['eigenvalues'][0],outer_minimum=d['minimum_outer']['selected']['eigenvalues'][0],
            outer_coordinate=d['minimum_outer']['coordinate'],witnesses=d['one_sided_witnesses']) for name,d in data['diagnostics'].items()}
        entry=dict(status='passed',checkpoint_sha256=sha(cp/'State.bin'),source_report=str(reportfile),source_report_sha256=sha(reportfile),
            archived_packet=str(rawfile),archived_packet_sha256=sha(rawfile),stages=rows,independent_tensor_checks=checks,
            receipt=receipt,stage_hashes=packet['stage_hashes'],history_sha256=packet['expected_history_sha256'],contract_rejections={})
        result['cases'][label]=entry;atomic(a.output/'summary.json',result)
        if a.contract_tests and label=='refined':
            mutations={'packet_hash':lambda q:q.update(archived_packet_sha256='0'*64),
                'stage_hash':lambda q:q['stage_hashes'].update(raw='0'*64),
                'state_placeholder':lambda q:q['state'].__setitem__(0,q['state'][0]+.01),
                'callback_history':lambda q:q['history'].__setitem__(0,q['history'][0]+.01),
                'native_history':lambda q:q['native_history'][1].__setitem__(0,q['native_history'][1][0]+.01),
                'signature':lambda q:q.update(native_signature='incompatible'),
                'time':lambda q:q.update(time=q['time']+.0001),
                'reference':lambda q:q['fixed_reference'].__setitem__(0,q['fixed_reference'][0]+.01)}
            for name,mutate in mutations.items():
                bad=deepcopy(packet);mutate(bad);error,rc=run_probe(a.probe,xml,bad,root/('reject-'+name),False)
                entry['contract_rejections'][name]=dict(receipt=rc,reason=error[-2000:]);atomic(a.output/'summary.json',result)
        require(sha(cp/'State.bin')==entry['checkpoint_sha256'] and sha(rawfile)==entry['archived_packet_sha256'],'forensic replay changed archived files')
    result.update(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat());atomic(a.output/'summary.json',result);print(a.output/'summary.json')


if __name__=='__main__':
    try:main()
    except Exception as error:
        if _created_output is not None:
            path=_created_output/'summary.json'
            previous=load(path) if path.exists() else {}
            previous.update(status='diagnostic_failed',failure=str(error))
            atomic(path,previous)
        raise
