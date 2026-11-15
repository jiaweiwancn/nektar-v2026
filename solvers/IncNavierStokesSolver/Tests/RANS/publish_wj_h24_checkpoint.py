#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small H24 publication checkpoints A/B, without scientific reevaluation."""
import argparse,gzip,json,shutil,tarfile
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,verify,atomic
from run_wj_h24 import usage,frozen

def publish(root,stage):
    cfg=frozen(root);book=load(root/'accounting.json');out=root/f'checkpoint-{stage.lower()}-v1';out.mkdir();copies=[]
    def copy(src,name):
        src=Path(src)
        if not src.exists():return
        p=out/name;p.parent.mkdir(parents=True,exist_ok=True)
        if src.suffix=='.json'and src.stat().st_size>100000:
            p=Path(str(p)+'.gz')
            with gzip.GzipFile(filename=str(p),mode='wb',mtime=0)as f:f.write(src.read_bytes())
        else:shutil.copyfile(src,p)
        copies.append(dict(original=identity(src),stored=identity(p)))
    for n in ('execution-config.json','source-manifest.json','inputs/selected.json','inputs/h23-closed-link.json','inputs/historical-indices.json',
        'science/actual-launch-binding-proof.json','science/runtime-build-v1.json','science/runtime-dependencies-v1.json',
        'science/target-cache-compatibility.json','science/preflight/result.json','preflight-host.json'):
        copy(root/n,n)
    for z in ('M','C','F'):copy(root/f'science/{z}/execution-diff.json',f'{z}/execution-diff.json')
    for label in ('runtime-input-setup-v1','actual-launch-bindings-v1','focused-controls-v1','retained-triplet-preflight-v1','freeze-admission-v1'):
        for n in ('request.json','receipt.json','command.log'):copy(root/f'operations/{label}/{n}',f'operations/{label}/{n}')
    status=dict(schema='WJH24PublicationCheckpointV1',stage=stage,coherent_bindings='PASS',native_capsule=cfg['native_capsule'],
        Python_sources=cfg['source_manifest'],frozen_configuration=identity(root/'execution-config.json'),historical_H23='CLOSED_UNCHANGED',
        native_startup='NOT_YET_RUN',M='PENDING',C='PENDING',F='PENDING',physical_triplet_preflight='PASS_EXACT_RETAINED_02',
        expected_remaining_advances=28,budget=usage(book),scope='publication point, continue automatically; no scientific PASS inferred')
    if stage=='B':
        q=load(root/'science/M/qualification.json');strict=load(root/'science/M/final-strict/result.json')
        assert q['status']==strict['status']=='PASS'and strict['selected_absolute_step']==16
        status.update(native_startup='M PASS, C/F pending',M='QUALIFIED_M16_T.04_STRICT_PASS',expected_remaining_advances=20)
        for n in ('qualification.json','outcome.json','final-strict/result.json','final-strict/request.json','monitor/startup-restored-parity.json',
            'monitor/startup.verification-commit.json','monitor/step-16.qualification-commit.json','monitor/latest-fully-verified-checkpoint.json'):
            copy(root/'science/M'/n,'M/'+n)
        for n in ('measured-forecast-M.json',):copy(root/'science'/n,n)
        for step in (10,12,14,16):copy(root/f'science/M/monitor/step-{step}.result.json',f'M/step-{step}.result.json')
        for label in ('governor-M-v1','strict-M-v1'):
            for n in ('request.json','receipt.json','attempt-counts.json'):copy(root/f'operations/{label}/{n}',f'operations/{label}/{n}')
        from wj_h23_continuation import branch_forecast
        from study_wj_h23 import forecast_costs
        costs=forecast_costs(root);estimate=sum(branch_forecast(z,**costs[z])for z in ('C','F'))
        u=usage(book);assert u['native']+1.25*estimate<=4200 and u['total']+1.25*estimate<=6000
        status['remaining_cost']=dict(estimate_seconds=estimate,headroom=1.25,used=u,selected_forecasts=costs,protected_strict=240,closure=300)
    atomic(out/'status.json',status);atomic(out/'accounting-snapshot.json',book);atomic(out/'copy-bindings.json',copies)
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='immutable publication snapshot; campaign ledger remains OPEN'))
    archive=root/f'h24-checkpoint-{stage.lower()}-v1.tar.gz'
    with tarfile.open(archive,'w:gz')as t:
        for p in sorted(out.rglob('*')):
            if p.is_file():t.add(p,arcname=str(p.relative_to(root)))
    atomic(root/f'checkpoint-{stage.lower()}-transport.json',dict(archive=identity(archive),manifest=identity(out/'manifest.json')))
    print(json.dumps(dict(status=status,archive=identity(archive),manifest=identity(out/'manifest.json'))),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--stage',choices=('A','B'),required=True)
    a=p.parse_args();publish(a.root,a.stage)
