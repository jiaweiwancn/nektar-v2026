#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Sequential, owned CT2 post-reboot jobs; independent failures do not cancel others.
Invoke through ct_accepted_command.py. This does not admit a physical continuation.
"""
import argparse,json,os,sys
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_open import run
from ct_accepted_feasibility import write_json
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();s=Path(__file__).parent;rt=r/'runtime';os.environ.update(LD_LIBRARY_PATH=str(rt/'lib'),TMPDIR=str(r/'tmp'))
    (r/'logs-post-reboot').mkdir(exist_ok=False);jobs=[]
    def add(name,script,*args):jobs.append((name,[sys.executable,'-B',str(s/script),*[str(x) for x in args]]))
    base=['--runtime',rt,'--previous',a.previous]
    add('recover-primary','ct_accepted2_diagnosis.py',*base,'--mode','recover','--trajectories',r/'matrix-v1','--preparation',r/'initial-v2','--output',r/'recovered-primary-v1')
    add('initial-companion','ct_accepted2_diagnosis.py',*base,'--mode','initial','--layouts','Zminus','--output',r/'companion-initial-v1')
    add('advance-companion','ct_accepted2_diagnosis.py',*base,'--mode','advance','--layouts','Zminus','--preparation',r/'companion-initial-v1','--output',r/'companion-v1')
    add('boundary-verification','ct_accepted_open.py','--runtime',rt,'--output',r/'boundary-verification-v1')
    add('history-closed-controls','ct_accepted_open_controls.py','--runtime',rt,'--old-runtime',a.previous/'runtime-v2','--output',r/'history-controls-v1')
    add('audit-primary','ct_accepted2_audit.py','--runtime',rt,'--preparation',r/'initial-v2','--trajectories',r/'recovered-primary-v1','--output',r/'primary-audits-v1')
    add('audit-companion','ct_accepted2_audit.py','--runtime',rt,'--preparation',r/'companion-initial-v1','--trajectories',r/'companion-v1','--output',r/'companion-audits-v1')
    add('shared-model-smokes','ct_accepted2_diagnosis.py',*base,'--mode','smoke','--layouts','Zminus','--output',r/'shared-smokes-v1')
    add('dataset-audit','cao_tamura_data.py','--data',s.parents[3]/'docs/plans/cao-tamura-validation-data','--report',r/'dataset-audit.json')
    add('wall-heights','ct_accepted2_wall.py','--runtime',rt,'--previous',a.previous,'--data',s.parents[3]/'docs/plans/cao-tamura-validation-data','--output',r/'wall-heights-v1')
    executions=[]
    for name,cmd in jobs:
        print('starting',name,flush=True)
        e=run(cmd,r/'logs-post-reboot',name+'.log');executions.append(dict(name=name,**e));write_json(r/'post-reboot-execution.json',executions)
        print('finished',name,e['exit_code'],flush=True)
    return 0 if all(x['exit_code']==0 for x in executions) else 1
if __name__=='__main__':raise SystemExit(main())
