#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H2k exact seen/fresh/import ledger; evolution is a separate post-push action."""
import argparse,os
from pathlib import Path
from cao_tamura_h2d_state import load
from run_cao_tamura_budget_h2k import charge
from run_cao_tamura_budget_h2i import repeatable,fresh,pilot,COLD_HASH
from run_cao_tamura_budget_h2f_startup import StartupRun
from cao_tamura_budget_h2g import case
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,state_hash,sha,THREADS


def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=('seen','fresh-import','pilot'))
    for n in ('root','context','candidate','probe','solver','preflight','h1','seen','certificate','output'):p.add_argument('--'+n,type=Path)
    p.add_argument('--source-revision');a=p.parse_args();os.environ.update(THREADS)
    if a.phase=='seen':
        charge(a,'repeatable native qualification',native_checks=1);repeatable(a);return
    if a.phase=='fresh-import':
        require(load(a.seen)['passed'],'native/seen failure cannot consume fresh grid')
        charge(a,'one frozen hill GL28 audit',fresh_audits=1);fresh(a)
        receipt=a.output/'receipt.json';r=load(receipt)
        if not r['passed']:return
        cf=a.output/'certificate.json';c=load(cf);src=Path(__file__).resolve().parents[4]
        names=['solvers/IncNavierStokesSolver/Tests/RANS/run_cao_tamura_budget_h2k.py','solvers/IncNavierStokesSolver/Tests/RANS/qualify_cao_tamura_budget_h2k.py','solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_startup_receipt.py']
        require(a.source_revision and len(a.source_revision)==40,'exact H2k source revision required')
        c['source_revision']=a.source_revision;c['source_hashes'].update({n:sha(src/n) for n in names});c['H2k']=dict(version='H2kSingleFreshEpochV1',source_revision=a.source_revision,context_sha256=sha(a.context),fresh_receipt_sha256=sha(receipt),aggregate_ledger_file=str(a.root/'aggregate.json'),aggregate_ledger_sha256=sha(a.root/'aggregate.json'),physical_startup_goal='open until bounded pair/restart/rollback pass');atomic(cf,c)
        run=StartupRun(a.solver,cf,a.preflight,a.h1,a.output/'zero-import',builder=lambda folder,h1,pre,flat=False:case(folder,h1,pre,flat=flat))
        try:
            zero=run.checked_native('hill-zero',2e-6)
            require(state_hash(zero['state'])==c['state_sha256'] and state_hash(zero['state'])!=COLD_HASH,'hill candidate/reference distinction or exact state failed')
            require(state_hash(zero['fixed_reference'])==c['fixed_reference_sha256']==COLD_HASH and state_hash(zero['fixed_defect'])==c['fixed_defect_sha256'],'hill immutable reference/defect failed')
            require(zero['initialization_lift_calls_iterations']==[0,0] and zero['time']==0 and zero['steps']==0,'hill zero import used setup projection or advance')
            atomic(a.output/'hill-import.json',dict(version='H2kExactHillImportV1',passed=True,state_sha256=c['state_sha256'],base_sha256=COLD_HASH,certificate_sha256=sha(cf),fresh_receipt_sha256=sha(receipt),reference_sha256=c['fixed_reference_sha256'],defect_sha256=c['fixed_defect_sha256'],zero_initialization_lift_projection=True,physical_steps=0,native_import=run.out))
        except Exception as error:
            atomic(a.output/'hill-import.json',dict(version='H2kExactHillImportV1',passed=False,error=str(error),native_import=run.out));raise
        return
    c=load(a.certificate);b=c.get('H2k');require(b and b['version']=='H2kSingleFreshEpochV1','H2k epoch binding required')
    imported=load(a.certificate.parent/'hill-import.json');require(imported['passed'] and imported['certificate_sha256']==sha(a.certificate) and imported['state_sha256']==c['state_sha256'],'exact production hill import prerequisite failed/changed')
    require(sha(a.context)==b['context_sha256'],'H2k selected context changed')
    require(sha(Path(b['aggregate_ledger_file']))==b['aggregate_ledger_sha256'],'H2k certified epoch ledger changed');pilot(a)


if __name__=='__main__':main()
