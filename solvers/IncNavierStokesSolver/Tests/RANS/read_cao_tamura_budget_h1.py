#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict zero-Lift readback of actual H1 native h0 and accepted endpoint."""
import argparse
from copy import deepcopy
import math
from pathlib import Path
import platform
from continue_cls_epawf_pimple import atomic,freeze
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from run_cls_epawf_picard import load,require,sha,state_hash


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('probe','campaign','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();require(platform.system()=='Linux','authoritative native Linux reader required')
    a.output.mkdir(parents=True,exist_ok=False)
    probe,identity=freeze(a.probe.resolve(),a.output/'runtime')
    summary=dict(status='started',probe=identity,campaign=str(a.campaign),observations={},contract_rejections={},governing_solves=0,physical_commits=0)
    atomic(a.output/'summary.json',summary)
    for name in ('initial','whole8'):
        report=load(a.campaign/name/'case.evolution.json');cp=Path(report['checkpoint']);saved=checkpoint(cp)
        packet=dict(common_packet(cp),operation='finite_inflow_h1',state=saved['coupled'],history=saved['coupled'],
            time=saved['time'],alpha=0.,inner_stress_model='Linear',
            permitted_queries_file=str(a.campaign/'input/permitted-queries.txt'),permitted_queries_sha256=sha(a.campaign/'input/permitted-queries.txt'))
        q,receipt=run_probe(probe,a.campaign/name/'case.xml',packet,a.output/name)
        require(q['governing_solves']==q['physical_commits']==0 and not q['native_history_advanced'] and q['initialization_work']==[0,0],
            'H1 readback altered native state/history or called Lift')
        require(q['state_sha256']==state_hash(saved['coupled']) and q['reference_sha256']==state_hash(saved['reference']),'readback state/reference differs')
        fits={}
        for quantity in ('U','minusRxz','sigmaU'):
            rows=[r for r in q['permitted_flat_native_queries'] if r['quantity']==quantity]
            weight=math.fsum(r['weight'] for r in rows);valid=math.fsum(r['weight'] for r in rows if r['valid'])
            fits[quantity]=dict(weight=weight,valid_weight=valid,omitted_weight=weight-valid,
                weighted_RMSE=math.sqrt(math.fsum(r['weight']*(r['native_prediction']-r['target'])**2 for r in rows if r['valid'])/weight) if valid==weight else None,
                invalid_queries=[r for r in rows if not r['valid']])
        summary['observations'][name]=dict(receipt=receipt,checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),
            time=saved['time'],step=saved['steps'],state_sha256=q['state_sha256'],reference_sha256=q['reference_sha256'],
            complete_audit=q['complete_audit'],metrics=q['physical_projection_metrics'],
            native_selected_spatial_blocks=q['selected_spatial_block_norms'],resolution_error=q['residual_resolution_native_max_error'],
            provider=q['provider'],physical_profiles=q['physical_profiles'],profile_columns=q['physical_profile_columns'],
            actual_native_flat_data_comparison=fits)
        atomic(a.output/'summary.json',summary)
        if name=='whole8':
            for label,mutate in [
                ('checkpoint_hash',lambda v:v.update(checkpoint_sha256='0'*64)),
                ('state',lambda v:v['state'].__setitem__(0,v['state'][0]+1e-6)),
                ('native_history',lambda v:v['native_history'][0].__setitem__(0,v['native_history'][0][0]+1e-6)),
                ('signature',lambda v:v.update(native_signature='incompatible'))]:
                bad=deepcopy(packet);mutate(bad)
                error,rr=run_probe(probe,a.campaign/name/'case.xml',bad,a.output/('reject-'+label),False)
                summary['contract_rejections'][label]=dict(receipt=rr,reason=error[-1800:])
            require(sha(cp/'State.bin')==packet['checkpoint_sha256'],'immutable parent checkpoint bytes changed')
    summary['status']='passed';atomic(a.output/'summary.json',summary);print(a.output/'summary.json')


if __name__=='__main__':main()
