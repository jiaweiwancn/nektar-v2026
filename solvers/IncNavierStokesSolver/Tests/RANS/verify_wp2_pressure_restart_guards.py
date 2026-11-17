#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict native changed-dt and cross-pressure checkpoint rejection."""
import argparse
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha
from wp2_native import Run
from run_cao_tamura_budget_wp2 import charge


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','h1','preflight']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True);p.add_argument('--replay-rejections',action='store_true');a=p.parse_args()
    h0=load(a.root/('pressure-transaction-'+a.model)/'h0/case.evolution.json');parent=Path(h0['checkpoint']);before=sha(parent/'State.bin')
    results={}
    for label,space,dt in [('changed-dt','inner-p2-v1',.01),('changed-pressure','baseline-v1',.02)]:
        output=a.root/('pressure-negative-'+label+'-'+a.model)
        reason='signature' if label=='changed-dt' else 'invalid checkpoint payload size'
        if a.replay_rejections:
            saved=load(output/'progress.json');item=saved['results']['reject']
            require(saved['solver_sha256']==sha(a.solver) and item['code']!=0 and item['observed_physical_attempts']==0,
                    'retained rejection runtime/count/exit differs')
            require(sha(Path(item['report']))==item['report_sha256'],'retained rejection report changed')
            native=load(Path(item['report']));expected='Wall IMEX: checkpoint signature mismatch' if label=='changed-dt' else 'Wall IMEX: invalid checkpoint payload size'
            require(native['failure']==expected,'retained rejection cause differs')
            archive=output/'before-rejection-message-recheck.json';require(not archive.exists(),'rejection already replayed');atomic(archive,saved)
            require(reason in (output/'reject/solver.log').read_text().lower(),'retained native log rejection differs')
            x=Run.__new__(Run);x.result=saved;x.output=output
        else:
            charge(a.root,'native pressure checkpoint negative '+a.model+' '+label,readback_negative_guards=1)
            x=Run(a.root,a.solver,a.h1,a.preflight,output,a.model,4,16,'RadialEnergyV1',False,space)
            x.native('reject',0,dt,parent=parent,mismatch=True,mismatch_reason=reason)
        require(sha(parent/'State.bin')==before,'negative restart altered parent bytes')
        x.result.update(status='passed_native_incompatibility_rejection',parent_checkpoint_sha256=before,governing_solves=0)
        x.save();results[label]=dict(path=str(x.output/'progress.json'),sha256=sha(x.output/'progress.json'))
    atomic(a.root/('pressure-negative-guards-'+a.model+'.json'),dict(version='WP2PressureNativeRestartGuardsV1',passed=True,model=a.model,
        source_sha256=sha(Path(__file__)),solver_sha256=sha(a.solver),parent_checkpoint_sha256=before,results=results,governing_solves=0))


if __name__=='__main__':main()
