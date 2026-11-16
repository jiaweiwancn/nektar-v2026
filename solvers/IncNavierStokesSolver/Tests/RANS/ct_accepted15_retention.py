#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""SelectiveReplayV1: sealed, owned scratch retirement after verified replay.

This does not preserve every original file losslessly. Native restart states,
runtime/configuration, commands, scientific receipts and replay witnesses are
retained. Callers must actually execute and compare the fresh-process replay.
"""
import hashlib,json,os
from pathlib import Path
from rans_request_manifest import canonical,sha
from ct_accepted_open import require

POLICY='SelectiveReplayV1'
def atomic(path,value):
    path=Path(path);temporary=path.with_name(path.name+'.pending')
    require(not path.exists() and not temporary.exists(),'existing or interrupted publication')
    with temporary.open('xb') as stream:
        stream.write(canonical(value)+b'\n');stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)

def owned(root,name):
    path=root/name
    require(not Path(name).is_absolute() and '..' not in Path(name).parts and path.resolve().is_relative_to(root.resolve()),'foreign retirement path')
    require(not path.is_symlink(),'symlink retirement path')
    return path

def verify(root,publication,binding):
    root=Path(root);publication=Path(publication)
    require((publication/'seal.json').is_file() and not (publication/'manifest.json.pending').exists(),'incomplete publication seal')
    seal=json.loads((publication/'seal.json').read_text());manifest=json.loads((publication/'manifest.json').read_text())
    require(seal['manifest_sha256']==sha(publication/'manifest.json'),'corrupt publication manifest')
    require(manifest['policy']==POLICY and manifest['binding']==binding,'stale/foreign native identity')
    require(manifest['verified_replay']['status']=='pass' and manifest['verified_replay']['decisions_exact_equal'] is True,'unverified replay')
    require(manifest['verified_replay']['maximum_scaled_array_difference']<=2e-10 and manifest['verified_replay']['maximum_scaled_state_pressure_traction_difference']<=1e-8,'replay exceeds prospective tolerance')
    for name,row in manifest['required_files'].items():
        f=owned(root,name);require(f.is_file() and sha(f)==row['sha256'],'missing/corrupt required restart or audit:'+name)
    require(len(manifest['readable_restart_states'])>=2,'sole recoverable state may not be pruned')
    for name in manifest['readable_restart_states']:
        require(any(k==name or k.startswith(name+'/') for k in manifest['required_files']),'unbound readable restart')
    for name,row in manifest['retired_payloads'].items():
        f=owned(root,name)
        if f.exists():require(sha(f)==row['sha256'],'changed disposable payload:'+name)
    return manifest

def publish(root,publication,binding,required,restart_states,eligible,replay):
    root=Path(root).resolve();publication=Path(publication);publication.mkdir(exist_ok=False)
    required={str(Path(f).resolve().relative_to(root)):dict(sha256=sha(f),bytes=Path(f).stat().st_size) for f in required}
    payloads={str(Path(f).resolve().relative_to(root)):dict(sha256=sha(f),bytes=Path(f).stat().st_size) for f in eligible}
    require(not set(required)&set(payloads),'required recovery evidence marked disposable')
    manifest=dict(policy=POLICY,binding=binding,required_files=required,readable_restart_states=[str(Path(f).resolve().relative_to(root)) for f in restart_states],retired_payloads=payloads,verified_replay=replay)
    atomic(publication/'manifest.json',manifest)
    atomic(publication/'seal.json',dict(manifest_sha256=sha(publication/'manifest.json')))
    verify(root,publication,binding);return manifest

def retire(root,publication,binding):
    root=Path(root);manifest=verify(root,publication,binding)
    # The complete group is validated before its first unlink. Only disposable
    # files from this newly owned CT15 root appear in retired_payloads.
    for name in manifest['retired_payloads']:
        f=owned(root,name)
        if f.exists():f.unlink()
    return sum(v['bytes'] for v in manifest['retired_payloads'].values())
