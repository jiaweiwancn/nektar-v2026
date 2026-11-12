#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Lightweight complete CT2 preservation and stale/forged primitive controls."""
import argparse,copy,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted_budget_gate import physical_budgets,required_times
from ct_accepted3_controls import negatives
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha
from ct_accepted3_run import make_case
from ct_accepted2_diagnosis import function_identity


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ct2',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    manifest=RequestManifest(a.output,dict(schema='CT3CompleteBudgetNegativesV1',source_sha256=sha(Path(__file__))))
    result=dict(status='pass',existing_negatives=negatives(a.ct2),historical_CT2={},additional_negatives={})
    require(required_times([.008,.02,.2])['status']=='pass','complete_milestone_positive')
    missing=required_times([.02,.2]);require(missing['status']=='fail' and missing['missing_required_times']==[.008],'missing_time_pass_label_not_rejected')
    result['additional_negatives']['missing_required_milestone_with_pass_label']=dict(status='pass',enclosing_status='pass',observed=missing)
    for folder in ('companion-audits-v1','continue02-audits-v1','continue2-audits-v1'):
        data=json.loads((a.ct2/folder/'results.json').read_text());rows={}
        for dt in ('0.002','0.001'):
            name='Zminus-'+dt;q=physical_budgets(data['cases'][name]);require(q['status']=='pass','complete_historical_CT2_changed:'+folder+name);rows[name]=q
        result['historical_CT2'][folder]=dict(receipt_sha256=sha(a.ct2/folder/'results.json'),members=rows)
    good=json.loads((a.ct2/'continue2-audits-v1/results.json').read_text())['cases']['Zminus-0.002']
    for b in good['budgets'].values():
        for q in b['reconstruction'].values():q['normalization']=q['maximum_defect']/q['scaled_defect'] if q['scaled_defect'] else 1.
    good['budget_binding']=dict(state='actual',configuration='same',time=.2,model='A95')
    for b in good['budgets'].values():b['binding']=good['budget_binding']
    require(physical_budgets(good)['status']=='pass','complete_bound_budget_control')
    for name,change,reason in [
        ('foreign_budget_state',lambda q:q['budgets']['native'].update(binding=dict(q['budget_binding'],state='foreign')),'state_configuration_budget_identity_mismatch'),
        ('forged_reconstruction_label',lambda q:q['budgets']['q64f16']['reconstruction']['k'].update(maximum_defect=1.),'inconsistent_reconstruction_primitives'),
        ('nonfinite_other_required_source',lambda q:q['budgets']['native']['equations']['whole_k'].update(instantaneous_source=float('inf')),'nonfinite_equation_primitive')]:
        q=copy.deepcopy(good);change(q);got=physical_budgets(q)
        require(got['status']=='fail' and reason in json.dumps(got),'wrong_negative_reason:'+name)
        result['additional_negatives'][name]=dict(status='pass',intended_reason=reason,observed=got)
    source=a.ct2/'companion-initial-v1/work/Zminus-0.002'
    c=json.loads((source/'construction.json').read_text())
    out=a.output/'metadata-only-no-evolution';make_case(out,source/'case.xml',c,3,.001,12)
    actual=json.loads((out/'construction.json').read_text());root=ET.parse(out/'case.xml').getroot()
    require(actual['dt']==.001 and actual['steps']==12 and actual['segment_duration']==.012 and actual['xml_sha256']==sha(out/'case.xml'),'construction_execution_metadata_not_actual')
    require(function_identity(root)==function_identity(ET.parse(source/'case.xml').getroot()),'metadata_fix_changed_physical_functions')
    result['construction_metadata_control']=dict(status='pass',actual_dt=actual['dt'],segment_steps=actual['steps'],segment_duration=actual['segment_duration'],
        physical_functions='identical',role='writer-only regression; no native solve, new candidate, or trajectory rerun')
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']);print('complete CT2 records pass; missing/stale/forged primitives reject');return 0
if __name__=='__main__':raise SystemExit(main())
