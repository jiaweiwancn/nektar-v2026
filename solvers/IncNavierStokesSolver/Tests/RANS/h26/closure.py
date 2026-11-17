#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only final provenance/resource audit, then close only the new ledger."""
import argparse,datetime,hashlib,json,os,platform,subprocess
from pathlib import Path
from run import scientific_workers,size,update_journal

def load(p):return json.loads(p.read_text())
def identity(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def write(p,value):
    assert not p.exists(),str(p)
    p.write_text(json.dumps(value,indent=2)+'\n')
def checked(item):
    p=Path(item['path']);actual=identity(p);assert actual['sha256']==item['sha256']and actual['bytes']==item['bytes'],str(p);return actual
def audit(root,version,output):
    source=load(root/f'source-v{version}.json');files=[checked(f)for f in source['files']]
    runtime=load(root/'science'/f'runtime-v{version}.json')
    for row in runtime['executables']:
        for key in ('executable','wrapper','link'):checked(row[key])
    checked(runtime['loader'])
    index=load(output/'linux-index.json');checked_index=[]
    for item in index:
        if Path(item['path'])==root/'accounting.json':continue # current mutable journal, closed separately
        checked_index.append(checked(item))
    libraries=set();maps=[];nonfiles=set()
    for p in sorted((root/'operations').glob(f'*-v{version}.native-maps')):
        maps.append(identity(p))
        for line in p.read_text().splitlines():
            parts=line.split(maxsplit=5)
            if len(parts)!=6 or not parts[5].startswith('/'):continue
            path=Path(parts[5])
            if path.is_file():libraries.add(path)
            else:nonfiles.add(str(path))
    assert maps
    loaded=[identity(p)for p in sorted(libraries)]
    assert any('WJNativeFlowDriver.elf'in p['path']for p in loaded)and any('WJStockFlowDriver.elf'in p['path']for p in loaded)
    old=root.parent/'wj-h25-production-20261116';originals=[]
    for p,expected in ((old/'accounting.json','1d11e9cfb1871c3ca2dd723578d26479ba7b0871e555e009bb4324208565ef1a'),(old/'resumption-r1/accounting.json','c5903152bb6521ebc45dd9636f3878622d378ff959623274264239529d0b24a7')):
        actual=identity(p);assert actual['sha256']==expected and load(p)['closed'];originals.append(actual)
    old_manifest=load(old/'resumption-r1/source-manifest-v3.json');old_files=[checked(f)for f in old_manifest['files']]
    # Report auxiliary bytecode separately from the pinned scientific capsule.
    caches=[dict(path=str(p),mtime=p.stat().st_mtime)for p in (old/'resumption-r1/source-v3').rglob('*.pyc')if p.stat().st_mtime>=1794886106]
    a=load(root/'accounting.json');receipts=[];failures=[];peak_build=peak_science=host=0;admission=0
    for p in sorted((root/'operations').glob('*/receipt.json')):
        row=load(p);receipts.append(identity(p));host=max(host,row.get('peak_host_workers',0));admission=max(admission,row.get('allocation',{}).get('workers',0)+row.get('reservation',0));rss=row.get('peak_tree_rss_bytes',0)
        if row.get('phase')=='build':peak_build=max(peak_build,rss)
        else:peak_science=max(peak_science,rss)
        if row['status']!='PASS':failures.append(dict(label=p.parent.name,status=row['status'],phase=row.get('phase'),seconds=row.get('seconds'),reason=row.get('reason'),exit=row.get('exit'),receipt=identity(p),expected_negative_control='restart-reject-height'in p.parent.name))
    code_bytes=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','build','runtime')));science_bytes=size(root/'science');free=os.statvfs(root).f_bavail*os.statvfs(root).f_frsize
    assert code_bytes<=8*1024**3 and science_bytes<=8*1024**3 and free>=15*1024**3 and peak_build<=16*1024**3 and peak_science<=8*1024**3 and host<=72 and admission<=72
    checkout=Path('/media/jay/新加卷1/nektar-v2026');environment=dict(os.environ,GIT_OPTIONAL_LOCKS='0')
    revision=subprocess.check_output(['git','-C',str(checkout),'rev-parse','HEAD'],text=True,env=environment).strip();status=subprocess.check_output(['git','-C',str(checkout),'status','--short'],text=True,env=environment)
    models=sorted({line.split(':',1)[1].strip()for line in Path('/proc/cpuinfo').read_text().splitlines()if line.startswith('model name')})
    summary=load(output/'summary.json');assert summary['H26_INTEGRATION_READINESS']=='READY'and summary['H26_NUMERICS']=='PASS'and summary['H26_COST']=='PASS'
    result=dict(status='PASS',tested_source=source['revision'],source_files_checked=len(files),scientific_old_source_files_checked=len(old_files),indexed_immutable_files_checked=len(checked_index),runtime_manifest=identity(root/'science'/f'runtime-v{version}.json'),native_map_snapshots=maps,actual_loaded_files=loaded,mapped_virtual_or_deleted_files=sorted(nonfiles),old_closed_accounts_unchanged=originals,old_auxiliary_bytecode_created_since_takeover=caches,source_build_runtime_bytes=code_bytes,science_bytes=science_bytes,free_bytes=free,peak_build_tree_rss_bytes=peak_build,peak_science_or_analysis_tree_rss_bytes=peak_science,peak_recorded_actual_host_workers=host,maximum_start_admission_bound=admission,live_workers=scientific_workers(),filesystem=subprocess.check_output(['findmnt','-n','-o','FSTYPE,SOURCE,TARGET','-T',str(root)],text=True).strip(),host=platform.uname()._asdict(),cpu_models=models,controller_cpu_affinity=sorted(os.sched_getaffinity(0)),general_checkout=dict(path=str(checkout),revision=revision,status=status,note='read-only snapshot; no intended H26 source was synchronized into this unrelated working tree'),account_before_close=identity(root/'accounting.json'),operation_receipts=receipts,failed_operations=failures,bytes_check_note='heavy H26 files are confined to the new ext4 namespace; old pinned source/data/account files were checked without writes')
    write(output/'runtime-resource-audit.json',result);print('PASS runtime/source/resources',host,peak_science,code_bytes,science_bytes)
def close(root,version,output):
    summary=load(output/'summary.json');audit=load(output/'runtime-resource-audit.json');assert summary['H26_INTEGRATION_READINESS']=='READY'and audit['status']=='PASS'
    def finish(book):
        assert not book['closed']and not any(v['status']=='RESERVED'for v in book['entries'].values())
        label='final-publication-epilogue';assert label not in book['entries'];book['entries'][label]=dict(status='FULLY_CHARGED',phase='evidence',seconds=120.,scope='final ledger close, manifest/transport checks and read-only publication verification; conservatively prepaid')
        used=sum(v['seconds']for v in book['entries'].values());assert used<=21600
        book['closed']=True;book['closed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();book['closure']=dict(status='ALL_H26_GATES_PASSED',tested_source=load(root/f'source-v{version}.json')['revision'],used_seconds=used,remaining_seconds=21600-used,summary=identity(output/'summary.json'),runtime_resource_audit=identity(output/'runtime-resource-audit.json'))
    update_journal(root/'accounting.json',finish);book=load(root/'accounting.json');write(output/'accounting-final.json',book)
    manifest=[identity(p)for p in sorted(output.iterdir())if p.is_file()];receipts=[identity(p)for p in sorted((root/'operations').glob('*/receipt.json'))];requests=[identity(p)for p in sorted((root/'operations').glob('*/request.json'))];write(output/'publication-manifest.json',dict(files=manifest,operation_receipts=receipts,operation_requests=requests,closed_account=identity(root/'accounting.json'),closure=book['closure']));print(json.dumps(book['closure'],indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--version',type=int,required=True);p.add_argument('--mode',choices=('audit','close'),required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert a.root.name=='wj-h26-native-flow-20261117';globals()[a.mode](a.root,a.version,a.output)
