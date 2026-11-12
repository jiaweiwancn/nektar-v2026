#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small frozen-runtime regressions for the H2c provider integration."""
import argparse,os,subprocess,time
from pathlib import Path
from continue_cls_epawf_pimple import atomic,freeze,CampaignLock
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import THREADS,load,require,sha
from run_cao_tamura_budget_h2a import prepare
from run_cls_accepted_awf import case_xml

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','baseline','h1','parent','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output);os.environ.update(THREADS)
    new,ni=freeze(a.solver.resolve(),a.output/'new-runtime');old,oi=freeze(a.baseline.resolve(),a.output/'baseline-runtime')
    out=dict(status='running',runtime=ni,baseline=oi,tests={},runs={})
    def save():atomic(a.output/'progress.json',out)
    def run(name,solver,finite=False,steps=2,cp=None):
        folder=a.output/name
        xml=prepare(folder,a.h1,steps,5e-6,cp or a.parent,True) if finite else case_xml('Accepted-AWF-05G',folder,steps,dt=.000625,restart=cp,companion='regularized',companion_reason='neutral frozen-runtime regression only',checkpoint_every=1)
        env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),TMPDIR=str(folder));start=time.monotonic()
        cmd=[str(solver),str(xml)]
        with (folder/'solver.log').open('w') as f:r=subprocess.run(cmd,cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=900)
        file=xml.with_suffix('.evolution.json');q=load(file) if file.exists() else None
        require(r.returncode==0 and q and not q['failure'],'regression native rejection '+name+': '+(folder/'solver.log').read_text()[-1500:])
        out['runs'][name]=dict(command=cmd,seconds=time.monotonic()-start,xml_sha256=sha(xml),report=str(file),report_sha256=sha(file));save();return q
    save()
    try:
        for finite,mode in [(False,'NeutralLogV1'),(True,'FiniteInflowV1')]:
            b=run(mode+'-baseline2',old,finite);n=run(mode+'-new2',new,finite)
            require(same_history(b,n),'old '+mode+' state/history changed')
            for before,after in zip(b['pimple_attempts'],n['pimple_attempts']):
                for key in ('history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256'):
                    require(before[key]==after[key],mode+' governing stage changed '+key)
            out['tests'][mode+'-state-history-stages-bitwise']=True;save()
        first=run('neutral-new1',new,False,1);second=run('neutral-new1-restart',new,False,1,Path(first['checkpoint']))
        require(same_history(load(a.output/'NeutralLogV1-new2/case.evolution.json'),second),'neutral split restart differs')
        out['tests']['neutral-split-restart-bitwise']=True
        out['status']='passed';save()
    except Exception as e:out.update(status='failed',failure=str(e));save();raise
    finally:lock.stream.close()
if __name__=='__main__':main()
