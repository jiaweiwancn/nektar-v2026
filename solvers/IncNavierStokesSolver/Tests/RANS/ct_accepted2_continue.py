#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict matched Zminus continuation, only after actual endpoint/domain/budget gates."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_feasibility import write_json
from ct_accepted_open import require,run
from ct_accepted2_diagnosis import inventory,make_case,analyze_attempt,analytic,function_identity
from ct_accepted2_audit import values
from ct_accepted_pilot import compare
from rans_request_manifest import RequestManifest,sha,runtime_identity

def endpoint_domain(initial,report,c):
    o=initial['observation'];s=np.array(report['state']);off=np.array(o['offsets'],int)
    out=dict(native=inventory(o['volume'],[q['coordinate'] for q in o['volume']],s,off),fixed=inventory(initial['test_basis'],o['grid'],s,off))
    out['status']='pass' if all(x['domain_failure_count']==0 and x['covariance_minimum'] is not None and x['covariance_minimum']['value']>=-1e-10*c['K0'] for x in (out['native'],out['fixed'])) else 'fail'
    return out

def physical_budgets(row):
    checks={};all_values={}
    for label,b in row['budgets'].items():
        ok=b.get('reconstruction_status')=='pass' and b.get('accuracy_status')=='pass'
        for name,q in b.get('equations',{}).items():
            scale=q['physical_scale'];defect=q['storage_backward_difference']-q['stage_source']-q['stage_physical_boundary_flux']
            ok=ok and abs(defect-q['frozen_stage_physical_defect'])<=2e-10*max(1.,scale) and abs(defect)/scale<=.02
            all_values.setdefault(name,[]).append(defect/scale)
        for q in b.get('reconstruction',{}).values():ok=ok and q['scaled_defect']<=2e-10
        checks[label]=dict(status='pass' if ok else 'fail')
    require(set(checks)=={'native','q48f16','q64f16','q64f32'},'missing_frozen_budget_control')
    uncertainty={n:max(v)-min(v) for n,v in all_values.items()}
    return dict(status='pass' if all(c['status']=='pass' for c in checks.values()) and max(uncertainty.values())<=.02 else 'fail',controls=checks,scaled_signed_defect_spreads=uncertainty,uncertainty_limit=.02)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--parent',type=Path,required=True);p.add_argument('--audits',type=Path,required=True);p.add_argument('--target',type=float,choices=[.02,.2],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();rt=r/'runtime';os.environ['LD_LIBRARY_PATH']=str(rt/'lib');prep=r/'companion-initial-v1';pd=json.loads((a.parent/'results.json').read_text());parent=Path(pd.get('trajectory_directory',a.parent/'work'));ad=json.loads((a.audits/'results.json').read_text());pr=json.loads((prep/'results.json').read_text());weights=pr['global_weights'];n=pr['global_count'];start=.008 if a.target==.02 else .02
    manifest=RequestManifest(a.output,dict(schema='CT2ConditionalSameMethodContinuationV1',runtime=runtime_identity(rt/'IncNavierStokesSolver'),parent_sha256=sha(a.parent/'results.json'),audits_sha256=sha(a.audits/'results.json'),preparation_sha256=sha(prep/'results.json'),source_sha256=sha(Path(__file__)),start=start,target=a.target,dt=[.002,.001],limits=dict(temporal=1e-3,mass=1e-8,scalar=.02,covariance_relative_negative_tolerance=1e-10)))
    admission=dict(status='fail',members={});initials={};parents={};constructions={};sampled={}
    for dt in (.002,.001):
        name='Zminus-'+str(dt);folder=parent/name;report=json.loads((folder/'case.evolution.json').read_text());ini=json.loads((prep/'work'/name/'initial.json').read_text());c=json.loads((folder/'construction.json').read_text())
        domain=endpoint_domain(ini,report,c);budget=physical_budgets(ad['cases'][name])
        ok=not report['failure'] and abs(report['time']-start)<1e-13 and abs(report['dt']-dt)<1e-15 and domain['status']=='pass' and budget['status']=='pass' and ad['cases'][name]['mass']['scaled_defect']<=1e-8
        admission['members'][name]=dict(status='pass' if ok else 'fail',actual_domain=domain,budgets=budget,parent_report_sha256=sha(folder/'case.evolution.json'))
        q=values(ini,report['state'],n);q['observation']['time']=start;sampled[name]=q;initials[name]=ini;parents[name]=report;constructions[name]=c
    temporal=compare(sampled['Zminus-0.002'],sampled['Zminus-0.001'],weights,c);admission['temporal']=temporal
    # Independent analytic pressure representation on a fixed physical band.
    fine=json.loads((prep/'work/Zminus-0.002/initial.json').read_text());coarse=json.loads((r/'initial-v2/work/B-0.002/initial.json').read_text());finep=values(fine,fine['observation']['state'],n);coarsep=values(coarse,coarse['observation']['state'],n)
    coords=np.array(fine['observation']['grid'][:n]);mask=(coords[:,2]>=4.)&(coords[:,2]<=10.);w=np.array(weights)[mask];root=ET.parse(prep/'work/Zminus-0.002/case.xml').getroot();p0=analytic(root,'Outer','p',coords[mask,2])
    errors=[float(np.sqrt(w@((np.array(q['observation']['grid_values'])[mask,4]-p0)**2)/w.sum())) for q in (coarsep,finep)]
    b=json.loads((r/'initial-v2/results.json').read_text())['layouts']['B']['0.002']['pressure_projection'];f=pr['layouts']['Zminus']['0.002']['pressure_projection']
    representation=dict(critical_band_analytic_Pi_RMS=dict(B=errors[0],Zminus=errors[1]),pressure_gradient_ratio=abs(f['native_gradient_at_witness']/b['native_gradient_at_witness']),signed_common_load_ratio=abs(f['signed_common_test']['sum']/b['signed_common_test']['sum']),scope='fixed test/normalization and frozen physical functions; local representation evidence, not a general h-convergence claim')
    representation['status']='pass' if errors[1]<errors[0] and representation['pressure_gradient_ratio']<1 and representation['signed_common_load_ratio']<1 else 'fail'
    admission['representation']=representation
    admission['status']='pass' if all(x['status']=='pass' for x in admission['members'].values()) and temporal['status']=='pass' and representation['status']=='pass' else 'fail'
    write_json(a.output/'admission.json',admission)
    if admission['status']!='pass':
        result=dict(diagnostic_execution='pass',physical_completion='unattempted_failed_admission',admission=admission,cases={});write_json(a.output/'results.json',result);manifest.seal(a.output,['admission.json','results.json']);print('not admitted');return 1
    work=a.output/'work';work.mkdir();result=dict(schema='CT2MatchedContinuationV1',diagnostic_execution='pass',physical_completion='fail',admission=admission,cases={},comparisons={},trajectory_directory=str(work))
    observed={}
    for dt in (.002,.001):
        name='Zminus-'+str(dt);p=parents[name];ini=initials[name];c=constructions[name];folder=work/name
        make_case(folder,parent/name/'case.xml',c,'B',dt,round((a.target-start)/dt),Path(p['checkpoint']))
        require(function_identity(ET.parse(folder/'case.xml').getroot())==function_identity(ET.parse(parent/name/'case.xml').getroot()),'continuation_changed_physical_functions')
        e=run([rt/'IncNavierStokesSolver','case.xml'],folder,'solver.log');report=json.loads((folder/'case.evolution.json').read_text())
        require(report['accepted_states'][0]==p['state'],'continuation_initial_state_not_parent')
        domain=endpoint_domain(ini,report,c)
        row=dict(status='pass' if not report['failure'] and abs(report['time']-a.target)<1e-13 and domain['status']=='pass' else 'physical_failure',time=report['time'],steps=report['steps'],failure=report['failure'],requested_time=a.target,actual_endpoint_domain=domain,execution=e,native_report_sha256=sha(folder/'case.evolution.json'),state_sha256=exact_vector_hash(report['state']),attempts=analyze_attempt(ini,report,c),checkpoint=report['checkpoint'],history=report['native_history_times'])
        result['cases'][name]=row;observed[name]={}
        for t,s in zip(report['audit_times'],report['accepted_states'],strict=True):
            q=values(ini,s,n);q['observation']['time']=t;observed[name][round(t,12)]=q
        write_json(a.output/'progress.json',result)
    aa,bb=observed['Zminus-0.002'],observed['Zminus-0.001'];result['comparisons']['Zminus']=[compare(aa[t],bb[t],weights,c) for t in sorted(set(aa)&set(bb))]
    result['physical_completion']='pass' if all(x['status']=='pass' for x in result['cases'].values()) else 'fail'
    write_json(a.output/'results.json',result);manifest.seal(a.output,['admission.json','results.json']);print('continuation',result['physical_completion']);return 0 if result['physical_completion']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
