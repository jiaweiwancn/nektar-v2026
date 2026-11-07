#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native attempt-65 identity, active restart/rollback and earlier-step checks."""
import argparse
import os
import shutil
from pathlib import Path
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import header,atomic
from run_cls_epawf_height_comparison import case_xml
from run_cls_epawf_picard import load,write_json,run,require,sha,state_hash,THREADS
from run_cls_epawf_pimple import check_report,identity


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--old-root',type=Path,required=True);p.add_argument('--replay',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--phase',choices=('first','transactions','earlier'),default='first')
    a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=False)
    summary=dict(phase=a.phase,status='started',solver=identity(a.solver.resolve()),production_governing_solves=0,
        note='independent verification trajectories; scheduled equation counts checked in every native report')
    def save(**kw):summary.update(kw);atomic(a.output/'summary.json',summary)
    archive=a.old_root/'terminal-U05-graded';old=load(archive/'attempt.json')
    cp64=a.output/'checkpoint64';shutil.copytree(archive/'checkpoint',cp64)
    def native(name,cp,steps,backend='RetainedPool',hard_reject=0,executable=None):
        folder=a.output/name;xml=case_xml('U05-graded',folder,steps,restart=cp,inventory='ConstitutiveV3',
            bound_solver=backend,hard_reject=hard_reject)
        rc=run([str(executable or a.solver.resolve()),str(xml)],folder,'solver.log',7200);write_json(folder/'receipt.json',rc)
        d=load(folder/'case.evolution.json');check_report(d,True)
        require(state_hash(d['accepted_states'][0])==state_hash(header(cp,True)['coupled_state']),'native startup changed accepted coupled state')
        require(d['initial_k_bound_check']['activations']==0,'startup correction changed checkpoint')
        require(d['initial_k_bound_check']['work_controls'].get('finite_qp_solver','WorkingSet')==backend,'native backend mismatch')
        return d,rc
    if a.phase=='first':
        replay=load(a.replay/'result.json');require(replay['k_bound']['success'] and replay['complete_hard_admissible'],'saved-candidate repair is not hard-audited')
        start,rc=native('startup64',cp64,0);require(rc['exit_code']==0 and not start['failure'],'compatible checkpoint startup failed')
        accepted,rc=native('step65',cp64,1);save(native_receipt=rc)
        require(rc['exit_code']==0 and not accepted['failure'] and accepted['steps']==65,'native attempt65 did not accept')
        attempt=accepted['pimple_attempts'][-1];q=load(Path(attempt['k_snapshot']))
        for key in ('committed','post_flow','pre_k','raw','history','transport_rows','transport_residual','raw_increment'):
            require(state_hash(q[key])==state_hash(old[key]),'native first-step identity differs: '+key)
        require(q['matrix_sha256']==old['matrix_sha256'] and q['alpha']==old['alpha'] and q['time']==old['time'],'native frozen k system/time differs')
        for key in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):
            require(attempt[key]==1,'wrong production schedule: '+key)
        require(accepted['implicit_callback_count']==1 and attempt['additional_governing_equation_solves']==0,'hidden governing solve/callback')
        import numpy as np
        mass=np.fromfile(a.replay/'pool.json.Mff.bin',dtype='<f8');rows=np.array(q['transport_rows'],dtype=int);mass=mass.reshape((len(rows),len(rows)))
        difference=(np.array(q['corrected'])-np.array(replay['corrected']))[rows]
        norm=float(np.sqrt(max(0.,difference@mass@difference)))
        require(norm<1e-11,'native and exact-candidate diagnostic projections disagree materially')
        save(status='passed',startup_unchanged=True,raw_frozen_system_history_exact=True,one_sweep_verified=True,
            native_step65_checkpoint=header(Path(accepted['checkpoint'])),native_projection_mass_difference=norm,
            native_projection_exact=state_hash(q['corrected'])==state_hash(replay['corrected']),native_attempt=attempt,
            test_governing_solves_per_equation=1)
    elif a.phase=='transactions':
        whole,rc=native('whole2',cp64,2);require(rc['exit_code']==0 and not whole['failure'],'two-step repaired trajectory failed')
        half,rc=native('half1',cp64,1);require(rc['exit_code']==0 and not half['failure'],'split first step failed')
        cp65=Path(half['checkpoint']);split,rc=native('split1',cp65,1);require(rc['exit_code']==0 and not split['failure'],'split restart failed')
        save(restart_exact=same_history(whole,split))
        zero,rc=native('startup65',cp65,0);require(rc['exit_code']==0 and not zero['failure'],'accepted repaired checkpoint startup failed')
        rejected,rc=native('forced-rejection66',cp65,1,hard_reject=66)
        require('deliberate' in rejected['failure'].lower() and rejected['pimple_attempts'][-1]['bound_activations']==1,'forced rejection did not follow an active correction')
        save(status='passed',rollback_exact=same_history(zero,rejected),rejection_time_history_preserved=True,
            test_governing_solves_per_equation=5,accepted_checkpoint=header(cp65))
    else:
        matches=[q.parent for q in (a.old_root/'screening/U05-graded/trajectory').glob('checkpoint*/State.bin') if header(q.parent)['steps']==8]
        require(len(matches)==1,'immutable earlier active checkpoint8 missing');cp8=a.output/'checkpoint8';shutil.copytree(matches[0],cp8)
        current,rc=native('earlier-working',cp8,1,backend='WorkingSet');require(rc['exit_code']==0 and not current['failure'],'legacy backend regression failed')
        pinned,rc=native('earlier-pinned',cp8,1,backend='WorkingSet',executable=a.old_root/'screening/runtime/IncNavierStokesSolver')
        require(rc['exit_code']==0 and not pinned['failure'],'pinned earlier active fixture failed');save(legacy_exact=same_history(current,pinned))
        pooled,rc=native('earlier-pool',cp8,1);require(rc['exit_code']==0 and not pooled['failure'],'earlier retained-pool correction failed')
        import numpy as np
        rows=np.array(old['transport_rows'],dtype=int);mass=np.fromfile(a.replay/'pool.json.Mff.bin',dtype='<f8').reshape((len(rows),len(rows)))
        diff=(np.array(pooled['state'])-np.array(pinned['state']))[rows]
        save(status='passed',earlier_pool_mass_difference=float(np.sqrt(max(0.,diff@mass@diff))),
            earlier_pool_exact=state_hash(pooled['state'])==state_hash(pinned['state']),test_governing_solves_per_equation=3)
    print(a.output/'summary.json')


if __name__=='__main__':main()
