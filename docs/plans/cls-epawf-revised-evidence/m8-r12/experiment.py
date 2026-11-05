#!/usr/bin/env python3
"""Prepare immutable R12 frozen-QP or matched native inputs; run only on Linux."""
import argparse
import importlib.util
import json
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--phase',choices=('Replay','Complete','Verify','Compare'),required=True)
    p.add_argument('--cuts',type=int,default=96)
    p.add_argument('--direction',type=Path)
    p.add_argument('--reference',type=Path)
    p.add_argument('--direction-seconds',type=float,default=0.)
    a=p.parse_args()
    assert platform.system()=='Linux' and 1<=a.cuts<=512
    e=a.epoch.resolve();evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r12_prepare')
    old=h.load(evidence/'m8-r11/decision-ledger.json')
    prior=module(evidence/'m8-r10/experiment.py','r10_prepare')
    prior.prepare(e,a.out,h,evidence)
    manifest=h.load(a.out/'predictor.json')
    tree=ET.parse(a.out/'case.xml');info=tree.find('.//SOLVERINFO')
    for item in list(info):
        if item.get('PROPERTY')=='RANSWallR10Experiment':info.remove(item)
    ET.SubElement(info,'I',PROPERTY='RANSWallR12Experiment',VALUE=a.phase)
    params=tree.find('.//PARAMETERS');ET.SubElement(params,'P').text='RANSWallR12CutLimit = '+str(a.cuts)
    tree.write(a.out/'case.xml',encoding='unicode')
    assert h.predictor.protected_xml(a.out/'case.xml')==manifest['protected_physical_xml_sha256']
    control=next(x for x in old['policies'] if x['name']=='L1')
    h.bound(control['direction'])
    manifest['r12']=dict(seed=control['direction'],seed_cost=control['direction_seconds'],direction_cost=a.direction_seconds)
    if a.direction:manifest['r12']['direction']=h.pin(a.direction)
    if a.reference:manifest['r12']['reference']=h.pin(a.reference)
    manifest['current_target_xml']=h.pin(a.out/'case.xml')
    (a.out/'predictor.json').write_text(json.dumps(manifest,indent=2)+'\n')
    note=dict(frozen_investigation_seconds=900,matched_native_seconds=600,phase=a.phase,cut_admission_limit=a.cuts,
        strategies=['corrected dual cut exchange with evidence-based bounded extension','feasible-start method only if evidence requires'],
        native_trials_per_policy=9,Nf_tangents=3,Nf_armijo=10,Nf_source_values=50,
        root_gate=2e-12,action_gate=1e-7,original_row_gate=5e-12,
        KKT_stationarity_gate=1e-9,KKT_active_equality_gate=1e-9,KKT_complementarity_gate=1e-9,
        normalized_objective_upper_bound_allowance=2e-9,reference_objective_agreement=2e-8,primal_dual_gap=2e-8,
        exact_frozen_upper_Psi=272.9590401766774,
        source_compatibility=h.load(a.out/'budget.json')['physical_value_and_source_compatibility'])
    (a.out/'budget.json').write_text(json.dumps(note,indent=2)+'\n')

if __name__=='__main__':main()
