#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Resume only a proven zero-launch resource denial in a fresh process.

No rejected physical step is retried in process. The immutable CP24 transaction
itself remains the frozen engine's one rejection followed by native4+4 restart.
"""
import argparse
import json
from pathlib import Path
import sys
import time


def prelaunch_alias(receipt, name):
    if receipt and receipt.get('exit_code') == 125 and receipt.get('actual_attempts') == 0 and receipt.get('resource', {}).get('not_launched') is True:
        return name+'-worker-resume1'
    return name


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('root','study','driver-scripts','solver','probe'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--policy',choices=('EP','Budget'),required=True)
    a=parser.parse_args();sys.path.insert(0,str(a.driver_scripts.resolve()))
    from run_wj6_native import Study
    from wj6_resources import scientific_workers
    from wj_immutable import artifacts,file_identity,write_once,read
    from wj_qualification import component
    class Resumed(Study):
        def native(self,name,*args,**kwargs):
            original=self.root/name/'receipt.json'
            receipt=json.loads(original.read_text()) if original.exists() else None
            alias=prelaunch_alias(receipt,name)
            if alias!=name:
                write_once(self.root/(alias+'-lineage.json'),dict(original=file_identity(original),
                    actual_attempts=0,reason='known pre-Popen worker denial; fresh external invocation, unchanged immutable CP24',code=artifacts([Path(__file__)])))
            # This is bounded waiting BEFORE launch, never a governing retry.
            # Leave room for the native/helper group and a peer phase change.
            samples=[];start=time.monotonic()
            while True:
                snapshot=scientific_workers();samples.append(dict(elapsed=time.monotonic()-start,**snapshot))
                if snapshot['workers']<=56:break
                if time.monotonic()-start>=60:
                    raise RuntimeError('RESOURCE_LIMIT: bounded transaction prelaunch allocation wait')
                time.sleep(2)
            path=self.root/(alias+'-allocation-wait.json')
            if not path.exists():write_once(path,dict(samples=samples,reserved_workers=16,limit=72,code=artifacts([Path(__file__)])))
            result=super().native(alias,*args,**kwargs)
            if result[1].get('failure'):
                raise RuntimeError(result[1]['failure'])
            return result
    s=Resumed(argparse.Namespace(output=a.study,policy=a.policy,solver=a.solver.resolve(),probe=a.probe.resolve(),mode='transactions'))
    try:s.transactions()
    except Exception as error:
        s.s['components']['transactions']=component('RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(error) else 'OBSTRUCTED',str(error))
        s.s.setdefault('failures',[]).append(dict(mode='transaction_external_resume',reason=str(error),code=artifacts([Path(__file__)])))
        s.save();raise
    finally:s.ledger.finish_session(s.session);s.save()
    print(s.summary)


if __name__=='__main__':main()
