#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small immutable H5 publication. Heavy startup fields stay in the run root."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from run_wj_h5 import load,save,sha,THREADS
from study_wj_h5 import identity,write_csv


def main(root):
    out=root/'publication';out.mkdir(exist_ok=False);cal=load(root/'reference/calibration.json')
    results={};witnesses={};table=[];interfaces=[];heavy=[]
    for H in ('5','5.5','6'):
        folder=root/f'observations/H{H}mm';v=load(folder/'result.json');results[H]=v
        gates=dict(v['gates']);bridge=None
        if H=='5':
            bridge=load(folder/'independent-weak-v2.json');gates['independent_weak']=bridge['status']=='PASS'
        status='OBSTRUCTED'if not gates['native_initialization']else 'PASS'if all(gates.values())else 'UNRESOLVED'
        table.append(dict(H_w_mm=float(H),state_size=v['state_size'],native_initialization='PASS'if gates['native_initialization']else 'OBSTRUCTED',
            final_passive_status=status,raw_assessment_status=v['status'],weak='PASS'if gates['independent_weak']else 'NOT_RUN',
            principal_min_over_nu_eff=v['represented']['min_planar_principal_over_nu_eff'],
            min_covariance=v['represented']['min_covariance'],max_raw_divergence=v['represented']['max_raw_divergence'],
            max_U_error=v['represented']['max_errors']['U'],max_k_error=v['represented']['max_errors']['k'],
            max_epsilon_error=v['represented']['max_errors']['epsilon'],
            quadrature='PASS'if gates['quadrature']else 'UNRESOLVED',governing_timesteps=0))
        witnesses[H]=v['principal_witness_replays']
        for row in v['interface']:interfaces.append(dict(H_w_mm=float(H),**row))
        for path in (folder/'result.json',folder/'physical-rows.json',folder/'outer-loads.json',root/f'cases/H{H}mm/case.evolution.json',root/f'cases/H{H}mm/case.xml',root/f'cases/H{H}mm/case-design.json'):
            heavy.append(identity(path))
        if v['checkpoint']:heavy.append(identity(v['checkpoint']['path']))
        if v['native']:heavy.append(identity(v['native']['path']))
        results[H]=dict(status=status,gates=gates,raw_assessment=identity(folder/'result.json'),
            represented=v['represented'],interface=v['interface'],spatial_norms=v['spatial_norms'],
            native_audit=v['native_audit'],quadrature=v['quadrature'],
            weak_bridge=identity(folder/'independent-weak-v2.json')if bridge else None,
            weak=bridge['independent_weak']if bridge else v['independent_weak'],
            producer=v['producer'],xml=v['xml'],checkpoint=v['checkpoint'],native=v['native'],
            phase='fresh CP0'if gates['native_initialization']else 'uncommitted startup candidate; no accepted checkpoint')
    write_csv(out/'initializations.csv',table)
    save(out/'initializations.json',results);save(out/'principal-witnesses.json',witnesses);save(out/'interface-traces.json',interfaces)
    for name in ('calibration.json','candidates.csv','observations.csv','fit-residuals.csv'):shutil.copyfile(root/'reference'/name,out/name)
    for path in (root/'figures').iterdir():
        if path.is_file():shutil.copyfile(path,out/path.name)
    eq=load(root/'runtime-equivalence.json');actual=[]
    for label in ('CP0-H5mm','CP0-H5.5mm','CP0-H6mm','passive-H5mm'):
        records=load(root/'operations'/label/'loaded-maps.json');spatial=[v for v in records if Path(v['path']).name=='libSpatialDomains.so.5.9.0']
        assert len(spatial)==1 and spatial[0]['sha256']=='0c23e0de180a233bd88f70b265ae9865aeac456e6482f8ff0fa1b4be84d69dab'
        actual.append(dict(operation=label,actual_mappings=identity(root/'operations'/label/'loaded-maps.json'),SpatialDomains=spatial[0]))
    save(out/'runtime.json',dict(governing_producer=eq['governing_producer'],observer_producer=eq['observer_producer'],
        compiled_source_count=len(eq['compiled_inputs']),verified_ELF_count=len(eq['qualified_ELFs']),executables=eq['executables'],
        original_equivalence=identity(root/'runtime-equivalence.json'),actual_mappings=actual,threads=THREADS,
        arithmetic_policy=eq['reused_arithmetic_policy'],native_rebuild=False))
    receipts={}
    logs=out/'logs';logs.mkdir()
    for path in sorted((root/'operations').glob('*/receipt.json')):
        r=load(path);receipts[path.parent.name]=dict(r,receipt=identity(path))
        shutil.copyfile(path.with_name('command.log'),logs/(path.parent.name+'.log'))
    save(out/'command-receipts.json',receipts)
    gates=dict(exact_neutral_outer_and_finite_nu_wall='PASS',mean_only_calibration='PASS',
        fig3_continuum_candidates='PASS',all_7mm_outside_height_interval='EMPTY',
        native_5mm=results['5']['status'],native_5p5mm=results['5.5']['status'],native_6mm=results['6']['status'],
        governing_timesteps=0,governing_solves=0,planned_initializations=3,executed_initializations=3,extra_initialization='NOT_RUN',
        full_hill_matching_height='NOT_SELECTED',finite_depth_physical_inlet='NOT_IMPLEMENTED_IN_H5',
        roughness_sublayer_validity='UNRESOLVED: physical planes 1.7/1.8/1.9 cube heights; homogenization assumed',
        old_H2_H4='PAUSED; original outcomes/denials unchanged',EP_and_physical_hill='NOT_RUN',
        native_limitation='5.5/6mm startup Q3 projection has principal-negative face states; no field/mesh/constant repair or retry')
    save(out/'component-gates.json',gates)
    # Retain complete heavy identities, without copying checkpoints or weak rows.
    revision=lambda v:subprocess.check_output(['git','rev-parse',v],cwd=Path(__file__).parent,text=True).strip()
    shutil.copyfile(root/'accounting.json',out/'accounting-reserved.json')
    save(out/'manifest.json',dict(root=str(root),generating_revision=eq['generating_revision'],
        passive_revision=revision('5e19336'),corrected_weak_revision=revision('6499527'),publication_revision=subprocess.check_output(
            ['git','rev-parse','HEAD'],cwd=Path(__file__).parent,text=True).strip(),
        specifications=eq['specifications'],heavy=heavy,
        accounting_snapshot=identity(out/'accounting-reserved.json'),setup=identity(root/'setup.json'),
        figures_scope='continuum fit/reference and native CP0 representation, no time evolution',
        original_weak_assessment='UNRESOLVED retained; independent-weak-v2 corrects the exterior projected-test contraction only'))
    book=load(root/'accounting.json');assert book['cp0_reservations']==3 and book['governing_attempts']==0
    known=sum(v.get('seconds',v['reserved_seconds'])for v in book['entries'].values())
    save(out/'accounting-at-publication.json',dict(book,known_operation_seconds=known,
        conservative_total_including_full_setup_and_closure=known+book['setup_seconds']+book['closure_reserve'],
        note='publication reservation may still be open; final closure receipt is separate'))
    # Byte identity is sufficient only together with the explicit producer roles
    # above. There is no completed-stage packet in a zero-step H5 experiment.
    save(out/'seal.json',dict(files=[dict(name=str(p.relative_to(out)),**identity(p))for p in sorted(out.rglob('*'))if p.is_file()],
        governing_timesteps=0,source='H5 source-built runtime reuse, three fresh startup attempts and passive observations'))
    print(json.dumps(dict(component_gates=gates,initializations=table)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();main(a.root)
