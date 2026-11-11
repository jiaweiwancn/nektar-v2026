#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Final H2a pressure/projection fixture and unchanged native-step proof."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import platform
import subprocess

from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import atomic,freeze
from diagnose_cls_accepted_covariance import run_probe
from run_cao_tamura_budget_h2a import prepare
from run_cls_epawf_picard import THREADS,load,require,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','probe','focused','h1','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();require(platform.system()=='Linux','authoritative Linux validation required')
    a.output.mkdir(parents=True,exist_ok=False)
    solver,identity=freeze(a.solver.resolve(),a.output/'runtime')
    probe,pi=freeze(a.probe.resolve(),a.output/'probe-runtime')
    zero=load(a.focused/'open-h0/case.evolution.json');parent=Path(zero['checkpoint']);before=sha(parent/'State.bin')
    packet=load(a.focused/'open-h0-readback/input.json')
    result,receipt=run_probe(probe,a.focused/'open-h0/case.xml',packet,a.output/'fixture')
    require(result['boundary_fixture']['status']=='passed' and result['initialization_work']==[0,0],
        'final absolute-pressure/active-projection/flux fixture failed')
    def one(label,image):
        xml=prepare(a.output/label,a.h1,1,restart=parent,packets='EveryAttempt')
        env={**os.environ,**THREADS,'VECLIB_MAXIMUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1',
            'LD_LIBRARY_PATH':str(image.parent/'lib'),'TMPDIR':str(xml.parent)}
        with (xml.parent/'solver.log').open('w') as log:
            code=subprocess.run([str(image),str(xml)],cwd=xml.parent,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300).returncode
        q=load(xml.with_suffix('.evolution.json'));require(code==0 and not q['failure'],'final native single-step comparison rejected')
        return q
    with ThreadPoolExecutor(max_workers=2) as pool:
        old=pool.submit(one,'verified-v1',a.focused/'runtime/IncNavierStokesSolver')
        new=pool.submit(one,'final',solver);old,new=old.result(),new.result()
    require(same_history(old,new) and old['state_audits']==new['state_audits'],'final projection helper changed finite governing path or audit')
    for k in ('history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256'):
        require(old['pimple_attempts'][0][k]==new['pimple_attempts'][0][k],'final helper changed stage '+k)
    require(sha(parent/'State.bin')==before,'immutable parent checkpoint changed')
    atomic(a.output/'summary.json',dict(status='passed',runtime=identity,probe=pi,fixture_receipt=receipt,
        fixture=result['boundary_fixture'],parent_sha256=before,complete_stage_state_history_audit_bitwise=True,
        physical_steps_per_independent_verification_leg=1,diagnostic_projection_bound_calls=1,
        diagnostic_governing_solves=0,native_verification_governing_solves_per_equation=2))
    print(a.output/'summary.json')


if __name__=='__main__':main()
