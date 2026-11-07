#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Archive a terminal attempt and prepare, but never execute, a zero-PDE replay."""
import argparse
import shutil
from pathlib import Path
from continue_cls_epawf_pimple import load_progress,header
from run_cls_epawf_picard import load,write_json,require,sha,state_hash

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    progress=load_progress(a.campaign/'progress.json');require(progress['status']=='numerically_rejected','only terminal numerical attempts may be archived')
    segment=Path(progress['segment']);report=load(segment/'case.evolution.json');attempt=report['pimple_attempts'][-1]
    packet=Path(attempt['k_snapshot']);require(sha(packet)==attempt['k_snapshot_sha256'],'terminal raw packet changed')
    q=load(packet);matrix=Path(q['matrix_file']);require(sha(matrix)==q['matrix_sha256'],'saved transport A_k changed')
    checkpoint=Path(progress['checkpoint']);require(sha(checkpoint/'State.bin')==progress['checkpoint_sha256'],'accepted checkpoint changed')
    receipt=header(checkpoint,True);require(state_hash(receipt['coupled_state'])==state_hash(q['committed']),'packet is not based on the accepted checkpoint')
    shutil.copytree(checkpoint,a.output/'checkpoint');shutil.copy2(segment/'case.evolution.json',a.output/'terminal.evolution.json')
    shutil.copy2(segment/'case.xml',a.output/'case.xml');shutil.copy2(packet,a.output/'attempt.json');shutil.copy2(matrix,a.output/'transport_A.bin')
    request=dict(operation='k_bound',state=q['raw'],pre_k=q['pre_k'],history=q['history'],committed=q['committed'],
        time=q['time'],alpha=q['alpha'],transport_rows=q['transport_rows'],raw_sha256=state_hash(q['raw']),
        native_checkpoint=str((a.output/'checkpoint').resolve()),checkpoint_sha256=progress['checkpoint_sha256'],native_signature=report['signature'],
        fixed_external_coefficients=q['fixed_external_coefficients'],projection_only=False,require_history_h0=False,inner_stress_model='SaturatedR2')
    write_json(a.output/'replay-input.json',request)
    write_json(a.output/'manifest.json',dict(status='terminal_preserved_no_retry',accepted_steps=progress['accepted_steps'],accepted_time=progress['accepted_time'],
        attempted_time=q['time'],failure=attempt['failure'],checkpoint_sha256=progress['checkpoint_sha256'],original_segment=str(segment),
        raw_sha256=state_hash(q['raw']),latest_repair_sha256=state_hash(q['corrected']),history_sha256=state_hash(q['history']),
        packet_sha256=sha(packet),transport_A_sha256=q['matrix_sha256'],metric='native PicardKMass Mff; transport_A.bin is not the metric',
        bound=q['k_bound'],production_counts={k:attempt[k] for k in ['outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls']},
        replay_executed=False,accepted_state_or_history_modified=False,
        files={str(f.relative_to(a.output)):sha(f) for f in a.output.rglob('*') if f.is_file()}))
    print(a.output/'manifest.json')

if __name__=='__main__':main()
