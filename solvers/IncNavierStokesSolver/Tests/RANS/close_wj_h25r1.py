#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify retained provenance, then faithfully close both independent accounts."""
import argparse,datetime,hashlib,json,os,subprocess
from pathlib import Path
from run_wj_h4 import load,save
from run_wj_h25 import usage
from wj_h15_observer import identity,verify


def pins(value):
    if isinstance(value,dict):
        if {'path','sha256','bytes'}<=value.keys():yield value
        for v in value.values():yield from pins(v)
    elif isinstance(value,list):
        for v in value:yield from pins(v)


def validate(root):
    base=root/'resumption-r1';records={};rot=load(base/'science/C-sparse-restart-v1/rotation.json')
    assert load(base/'science/cadence-comparison-v1.json')['status']=='PASS_UNPERTURBED_CADENCE_AND_RESTART_ONLY'
    assert load(base/'science/native-independent-replays-v1.json')['status']=='PASS'
    assert load(base/'science/C-sparse-restart-v1/strict-after-rotation/qualification.json')['status']=='PASS_POST_ROTATION_TERMINAL_READ'
    removed={p['path']for r in rot['records']if r['id']in rot['removed']for p in r['payloads']}
    selected=load(base/'science/implementation-selected-v1.json')
    sources=[base/f'source-manifest-v{i}.json'for i in (1,2,3)]
    sources += [base/f'science/runtime-build-{n}-v1.json'for n in ('seed','scatter','memo')]
    sources += [root/'inputs/baseline.json',base/'science/implementation-selected-v1.json']
    sources += sorted(base.glob('tools-manifest-v*.json'))
    for d in ('fresh-C-v1','C-wide-v1','C-sparse-restart-v1','C-dense-restart-v1'):
        sources.append(base/'science'/d/'qualification.json')
    sources += [base/'science/C-sparse-restart-v1/strict/qualification.json']
    sources += [base/'science/C-sparse-restart-v1/strict-after-rotation/qualification.json']
    for source in sources:
        for pin in [identity(source),*pins(load(source))]:
            p=pin['path']
            if p in records:
                assert records[p].get('sha256')==pin['sha256'],'conflicting retained pin'
                continue
            if p in removed:
                assert not Path(p).exists();records[p]=dict(**pin,status='EXPECTED_OWNED_SUPERSEDED_RECOVERY_ROTATION');continue
            verify(pin);records[p]=dict(**pin,status='PASS')
    for folder in (root,base):
        a=load(folder/'accounting.json')
        current=Path(os.environ['WJ_H25_OPERATION']).name
        assert all((folder==base and label==current and v['status']=='RESERVED'and v['phase']=='passive')or
            (v['status']not in ('RESERVED','RUNNING')and not v.get('unknown_attempts'))for label,v in a['entries'].items())
        for label in a['entries']:
            p=folder/'operations'/label/'native-exec-maps.json'
            if p.exists():assert all(v['status']=='PASS'for v in load(p))
    out=base/'final-validation.json';assert not out.exists()
    save(out,dict(status='PASS_RETAINED_IDENTITIES',records=list(records.values()),
        selected_runtime=selected['runtime'],expected_rotated_recoveries=sorted(removed),
        scope='current bytes of source/runtime/library/input/selected trajectory and strict-read pins; numerical qualification remains exactly its recorded scope'))
    print('verified identities',len(records),'expected rotated payloads',len(removed),flush=True)


def close(root):
    base=root/'resumption-r1';assert load(base/'final-validation.json')['status']=='PASS_RETAINED_IDENTITIES'
    assert load(base/'science/complete-cost-v1.json')['status']=='TARGET_NOT_MET'
    # Only the present Python control process may refer to H25 now. Unrelated
    # scientific jobs are preserved; no process receives a signal here.
    ps=subprocess.check_output(['ps','-eo','pid=,args='],text=True).splitlines()
    live=[s for s in ps if str(root)in s and any(k in s for k in ('.elf','cmake --build','run_wj_h25r1.py --root'))]
    assert not live,live
    accounts={};before={}
    for name,folder in (('original',root),('continuation',base)):
        p=folder/'accounting.json';a=load(p);assert not a['closed']
        assert all(v['status']not in ('RESERVED','RUNNING')and not v.get('unknown_attempts')for v in a['entries'].values())
        for label,v in a['entries'].items():
            receipt=folder/'operations'/label/'receipt.json'
            if receipt.exists():
                r=load(receipt)
                for k in ('status','phase','seconds'):assert r[k]==v[k],(label,k)
        u=usage(a);assert all(u[k]<=a['ceilings'][k]for k in u)
        assert a['reserved_attempts']<=a['maximum_attempts']and a['governor_blocks']<=a['maximum_governor_blocks']
        before[name]=identity(p);snapshot=base/(name+'-account-before-close.json');assert not snapshot.exists();snapshot.write_bytes(p.read_bytes())
        ceilings=dict(a['ceilings']);charges=json.dumps(a['entries'],sort_keys=True)
        a.update(closed=True,closure_status='TARGET_NOT_MET',closed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            closure_scope='WJ-H25R1 checkpoint C; original limits and all historical charges preserved')
        save(p,a);after=load(p)
        assert after['ceilings']==ceilings and json.dumps(after['entries'],sort_keys=True)==charges and usage(after)==u
        accounts[name]=dict(before=before[name],closed=identity(p),usage=u,ceilings=ceilings,
            attempts=a['reserved_attempts'],advancing_blocks=a['governor_blocks'],
            ledger_only_control_entries=[k for k in a['entries']if not(folder/'operations'/k/'receipt.json').exists()])
    cumulative={k:sum(v['usage'][k]for v in accounts.values())for k in ('native','build','passive','total')}
    save(base/'closure-accounts.json',dict(status='CLOSED_TARGET_NOT_MET',accounts=accounts,cumulative=cumulative,
        fixed_reserves_fully_charged=dict(original_setup=120,original_closure=600,continuation_setup=300,continuation_closure=900),
        scope='charged Linux campaign operation wall time plus fully charged fixed control/publication reserves; Mac source work, reading and unavailable-host waiting are separate elapsed work, not one-case cost',
        live_owned_native_build_processes=live))
    print('closed campaign',json.dumps(cumulative),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('validate','close'),required=True);a=p.parse_args()
    (validate if a.mode=='validate'else close)(a.root)
