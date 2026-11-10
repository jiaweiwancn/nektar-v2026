#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal completed WJ-4 measurements, then bind final accounting without cycles.

The inner manifest/receipt is an immutable publication-time snapshot. The outer
publication.json binds that receipt, completed command logs and final ledger
bytes after the publication work itself has been reconciled. It never rewrites
either original studies or an existing publication.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from wj_immutable import Ledger,artifacts,read,write_once,file_hash
from run_cls_epawf_picard import THREADS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','observer-root','study','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();observer=a.observer_root.resolve();study=a.study.resolve();out=a.output.resolve()
    if platform.system()!='Linux':raise ValueError('Linux is authoritative')
    if out.exists() or (root/'figures').exists():raise ValueError('publication output is write-once; choose a new root/lineage')
    ledgers={p:Ledger(root,p) for p in ('EP','Budget')}
    for policy,l in ledgers.items():
        ledger=read(l.file)
        if any(s['status']=='RUNNING' for s in ledger.get('sessions',{}).values()) or any(e['status']=='RUNNING' for e in ledger['entries'].values()):
            raise ValueError('finish/reconcile active work before final publication: '+policy)
        s=read(study/policy/'summary.json')
        if not s.get('qualification') or s.get('read_only_reprocessing',{}).get('governing_attempts')!=0:
            raise ValueError('final immutable observer summary is incomplete: '+policy)
        if read(root/('all-step-transactions-'+policy+'.json')).get('status')!='PASS':
            raise ValueError('all-step transaction replay is incomplete: '+policy)
    scripts=Path(__file__).resolve().parent;env=os.environ.copy();env.update(THREADS)
    allocation=root/'publication-allocation.json'
    subprocess.run([sys.executable,'-B',str(scripts/'wj4_resources.py'),'--reserve','8','--output',str(allocation)],env=env,check=True)
    label=str(out)+'/publication';intent=dict(study=artifacts([study/p/'summary.json' for p in ('EP','Budget')]),
        code=artifacts([scripts/p for p in ('finalize_wj4.py','publish_wj4.py','plot_wj4.py')]))
    reserved=[]
    try:
        for l in ledgers.values():l.reserve(label,intent,0,300);reserved.append(l)
    except Exception:
        for l in reserved:l.finish(label,0,0,'publication not launched; another policy could not reserve work')
        raise
    start=time.time();failure=None
    try:
        for name,argv in (
            ('final-figures',[sys.executable,'-B',str(scripts/'plot_wj4.py'),'--study',str(study),'--output',str(root/'figures')]),
            ('final-publication',[sys.executable,'-B',str(scripts/'publish_wj4.py'),'--root',str(root),'--observer-root',str(observer),'--study',str(study),'--output',str(out)])):
            subprocess.run([sys.executable,'-B',str(scripts/'wj3_receipts.py'),'--root',str(observer),'--label',name,'--',*argv],env=env,check=True)
    except Exception as error:
        failure=str(error);raise
    finally:
        elapsed=time.time()-start
        for l in ledgers.values():l.finish(label,0,elapsed,str(out/'receipt.json') if not failure else failure)
    for policy,l in ledgers.items():shutil.copyfile(l.file,out/('final-ledger-'+policy+'.json'))
    shutil.copyfile(observer/'command-receipts.json',out/'final-observer-command-receipts.json')
    seal=dict(schema='WJ4FinalPublicationV1',inner_receipt_sha256=file_hash(out/'receipt.json'),
        complete_publication_seconds=elapsed,governing_attempts=0,
        final_accounting=artifacts([out/'final-ledger-EP.json',out/'final-ledger-Budget.json',out/'final-observer-command-receipts.json']),
        totals={p:l.totals() for p,l in ledgers.items()},inner_snapshot='manifest/receipt preserve the exact before-completion accounting snapshot; this outer seal binds reconciled final publication work')
    write_once(out/'publication.json',seal)
    for f in out.rglob('*'):
        if f.is_file():f.chmod(0o644)
    print(json.dumps(seal))


if __name__=='__main__':main()
