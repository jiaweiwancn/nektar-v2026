#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete a native report's steady screen using actual full epsilon mass.

This read-only confirmation works for the legacy flat cases too. The strict
native reader binds every selected checkpoint/history/reference; no production
audit upgrade, constitutive substitution, field correction or PDE is invoked.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone, timedelta
import json
import math
import os
from pathlib import Path
import shutil

from analyse_cls_accepted_awf import report
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint, common_packet, run_probe
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--common-panel',type=Path,required=True);parser.add_argument('--case',required=True)
    parser.add_argument('--probe',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--contract-tests',action='store_true');a=parser.parse_args();os.environ.update(THREADS)
    import numpy as np
    a.output.mkdir(parents=True,exist_ok=False);source=load(a.common_panel/'summary.json');case=source['cases'][a.case]
    policy=source['policy'];steps=case['selected_steps'];observations={};mass=None;field=None;packets={};inputs={}
    out=dict(status='running',case=a.case,probe_sha256=sha(a.probe),script_sha256=sha(Path(__file__)),
        common_panel_summary_sha256=sha(a.common_panel/'summary.json'),policy=policy,observations={},windows=[],
        governing_solves=0,physical_commits=0,native_history_advanced=False,
        metric='actual full physical consistent outer-epsilon mass; retain prescribed and matching coefficients; fixed common E0*sqrt(12) normalization',
        started_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat())
    atomic(a.output/'summary.json',out)
    for step in steps:
        row=case['observations'][str(int(step))];cp=Path(row['receipt']['checkpoint']);saved=checkpoint(cp)
        native=report(Path(row['receipt']['source_report']));xml=Path(row['receipt']['source_report']).parent/'case.xml'
        require(sha(cp/'State.bin')==row['receipt']['checkpoint_sha256'] and state_hash(saved['coupled'])==row['receipt']['state_sha256'], 'immutable source checkpoint differs')
        packet=dict(common_packet(cp),operation='native_epsilon_metric',state=saved['coupled'],history=saved['coupled'],
            alpha=0.,time=saved['time'],inner_stress_model=native['inner_stress_model'])
        result,receipt=run_probe(a.probe,xml,packet,a.output/f'step-{int(step):06d}')
        require(result['governing_solves']==result['physical_commits']==0 and not result['native_history_advanced'] and
            result['selected_state_unchanged'] and result['initialization_work']==[0,0], 'metric observation changed the selected state')
        offset=int(result['offset']);coeff=result['coefficients'];require(state_hash(coeff)==state_hash(saved['coupled'][offset:offset+len(coeff)]),'native epsilon coefficient bits differ')
        current=np.array(result['matrix']);require(current.shape==(len(coeff),len(coeff)) and np.max(abs(current-current.T))<1e-14, 'mass shape/symmetry differs')
        np.linalg.cholesky(current)
        constant=np.array(result['constant_coefficients']);volume=float(result['physical_volume'])
        require(abs(constant@current@constant-volume)<1e-10 and volume>0, 'physical constant norm differs')
        if mass is None:mass=current;field=result['field']
        else:require(np.array_equal(mass,current) and field==result['field'], 'native epsilon geometry/mass changed between windows')
        observations[int(step)]=result;packets[int(step)]=packet;inputs[int(step)]=xml
        out['observations'][str(int(step))]=dict(time=saved['time'],checkpoint_sha256=sha(cp/'State.bin'),
            state_sha256=state_hash(saved['coupled']),native_history_sha256=[state_hash(v) for v in saved['history']],
            native_history_times_sha256=state_hash(saved['history_times']),coefficients_sha256=state_hash(coeff),
            physical_volume=volume,constant_mass_error=abs(constant@current@constant-volume),receipt=receipt)
        atomic(a.output/'summary.json',out)
    scale=policy['common_normalization']['epsilon']*math.sqrt(policy['common_normalization']['outer_epsilon_volume'])
    for first,last,original in zip(steps,steps[1:],case['windows']):
        left=observations[int(first)];right=observations[int(last)];delta=np.array(right['coefficients'])-np.array(left['coefficients'])
        norm=math.sqrt(max(0.,float(delta@mass@delta)));terms=[]
        require(left['quadrature_rows']==right['quadrature_rows'],'element/side quadrature provenance differs')
        for q in left['quadrature_rows']:
            ids=np.array(q['ids'],dtype=int);value=float(np.array(q['phi'])@delta[ids]);terms.append(q['weight']*value*value)
        quadrature=math.sqrt(math.fsum(terms));require(abs(norm-quadrature)<1e-12,'consistent mass/independent GL4 epsilon difference disagree')
        normalized=norm/scale;passed=normalized<=policy['normalized_window_L2_change_maximum']['epsilon']
        out['windows'].append(dict(start=left['physical_time'],end=right['physical_time'],full_outer_epsilon_L2_change=norm,
            full_outer_epsilon_normalized_change=normalized,mass_quadrature_norm_error=abs(norm-quadrature),
            full_outer_epsilon_pass=passed,common_panel_screen_pass=original['common_screen_pass'],
            complete_screen_pass=passed and original['common_screen_pass']))
    out['contract_rejections']={}
    if a.contract_tests:
        first=int(steps[0]);base=packets[first];variants={}
        p=deepcopy(base);p['state'][0]+=.125;variants['state']=(p,'state differs from checkpoint')
        p=deepcopy(base);p['native_history'][0][0]+=.125;variants['history']=(p,'native history slot differs')
        p=deepcopy(base);p['alpha']=1.;variants['alpha']=(p,'requires alpha=0')
        p=deepcopy(base);p['native_signature']='incompatible';variants['signature']=(p,'physical signature differs')
        corrupt=a.output/'invalid-checksum-fixture';corrupt.mkdir();binary=bytearray((Path(base['native_checkpoint'])/'State.bin').read_bytes());binary[-9]^=1
        (corrupt/'State.bin').write_bytes(binary);p=deepcopy(base);p['native_checkpoint']=str(corrupt);p['checkpoint_sha256']=sha(corrupt/'State.bin');variants['checksum']=(p,'checksum')
        for name,(p,reason) in variants.items():
            text,receipt=run_probe(a.probe,inputs[first],p,a.output/('reject-'+name),False)
            require(reason in text,'unexpected contract rejection: '+name+': '+text[-500:]);out['contract_rejections'][name]=receipt
    count=sum(w['complete_screen_pass'] for w in out['windows'])
    out.update(status='passed',passing_windows=count,classification='candidate_steady_selected_operator_at_declared_screen' if count==3 else 'not_steady_by_declared_screen',
        completed_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat())
    atomic(a.output/'summary.json',out);print(a.output/'summary.json')


if __name__=='__main__':main()
