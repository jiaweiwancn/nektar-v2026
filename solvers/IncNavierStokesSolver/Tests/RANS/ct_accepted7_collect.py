#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT7 evidence, keeping blocked physical scopes and failed MMS gates.

Use the owned launcher after numerical families have closed. Attach this
collection family's own closed resource receipt after it returns. No history,
physical input, historical receipt or qualification threshold is modified.
"""
import argparse,json,os,subprocess
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted_open import require
from ct_accepted7_gate import top,history,resources,temporal,finite
from ct_accepted_budget_gate import physical_budgets
from rans_request_manifest import RequestManifest,sha


def read(folder):
    folder=Path(folder);seal=json.loads((folder/'completed-artifacts.json').read_text())
    require(all(sha(folder/n)==v for n,v in seal['files'].items()),'changed sealed receipt:'+str(folder))
    return json.loads((folder/'results.json').read_text())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','old','top','history','domain','temporal','regression','build','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--resource',type=Path,action='append',required=True)
    a=p.parse_args();s=a.root/'source';refs=read(a.root/'reference-v1');tp=read(a.top);hist=read(a.history);dom=read(a.domain);time=read(a.temporal);reg=read(a.regression);build=read(a.build)
    quad=read(a.root/'top-frozen-quadrature-v3');fine=read(a.root/'ct6-finer-p4');fine_audit=read(a.root/'ct6-finer-p4-audit')
    primitives={model:json.loads((a.history/('full' if model=='A95' else model)/'physical-budgets.json').read_text()) for model in ('A93','A95','E94E95')}
    closed={str(folder):json.loads((folder/'resource-final.json').read_text()) for folder in a.resource}
    m=RequestManifest(a.output,dict(schema='CT7CompactEvidenceV1',source=sha(Path(__file__)),receipts={str(f):sha(f/'results.json') for f in [a.root/'reference-v1',a.top,a.history,a.domain,a.temporal,a.regression,a.build]},
        resources={str(f):sha(f/'resource-final.json') for f in a.resource},contract=sha(s/'docs/plans/cao-tamura-accepted-7-evidence/contract.json')))
    tops=top(tp);histories=history(hist,primitives);times=temporal(time);res=resources(closed)
    qerrors=[]
    if set(quad.get('cases',{}))!={'A95-p3','A93-p3','E94E95-p3','A95-p4'}:qerrors.append('missing_frozen_case')
    for name,orders in quad.get('cases',{}).items():
        if set(orders)!={'48','64'}:qerrors.append(name+':missing_frozen_order')
        for order,row in orders.items():
            if not row.get('frozen_state_and_boundary_identical'):qerrors.append(name+':frozen_input_changed')
            equations=row.get('analytical_weak',{}).get('equations',{})
            if set(equations)!={'momentum_x','momentum_z','k_shared_once','outer_epsilon','inner_continuity','outer_continuity'}:qerrors.append(name+':missing_weak_equation')
            for equation,value in equations.items():
                v,scale=value.get('free_maximum'),value.get('scale')
                if not finite(v) or v<0 or not finite(scale) or scale<=0 or v/scale>2e-10:qerrors.append(name+':'+order+':'+equation)
    boundary_checks={}
    for name in ('full','first','continued','shifted','A93','E94E95'):
        path=a.history/name/'case.evolution.json';report=json.loads(path.read_text());receipt=report.get('initialization_receipt',{})
        ok=receipt.get('open_explicit_implicit_agreement') and receipt.get('explicit_boundary_sha256') and receipt.get('explicit_boundary_sha256')==receipt.get('implicit_boundary_sha256')
        boundary_checks[name]=dict(status='pass' if ok else 'fail',receipt=receipt,report_sha256=sha(path))
    if any(v['status']!='pass' for v in boundary_checks.values()):histories.update(status='fail',consumed_boundary_failure=True)
    domain_errors=[]
    if set(dom.get('cases',{}))!={'flat-p3','flat-p4','terrain-p3','terrain-p4'}:domain_errors.append('domain_pair_coverage')
    for label,row in dom.get('cases',{}).items():
        count=row.get('native_count',{}).get('stored_slots')
        if label.startswith('flat') and (not finite(count) or count<=3000):domain_errors.append('missing_or_unexpected_capacity:'+label)
        if set(row.get('native_geometry',{}))!={'16','32'}:domain_errors.append('missing_geometry_orders:'+label)
        for q,g in row.get('native_geometry',{}).items():
            keys=('minimum_sampled_native_Jacobian','volume_defect','expected_physical_volume')
            if any(not finite(g.get(k)) for k in keys) or g['minimum_sampled_native_Jacobian']<=0 or g['expected_physical_volume']<=0 or g['volume_defect']>2e-11*g['expected_physical_volume']:domain_errors.append('native_geometry:'+label+':'+q)
    cert=dom.get('analytic_geometry_certificate',{})
    if any(not finite(cert.get(k)) or cert[k]<=0 for k in ('minimum_interval_Jacobian','minimum_interval_xq','collar_sufficient_Jacobian_bound','outer_vertical_Jacobian_bound')):domain_errors.append('analytic_geometry_certificate')
    required_reg={'fresh_CT5_A95_two_steps','unchanged_restart','deliberate_rejection_rollback','incompatible_geometry_history'}
    regression=bool(set(reg.get('checks',{}))==required_reg and all(v['status']=='pass' for v in reg['checks'].values()))
    failed_roots={model:[h for h in (5.5,6.,6.5,6.99) if refs['cases'][model+'-'+str(h)]['roughness_match']['status']!='root'] for model in ('A93','A95','E94E95')}
    require(failed_roots['A95']==[5.5,6.,6.5,6.99],'unexpected central fit result: a new physical admission contract is required')
    decision=dict(schema='CT-Accepted-7-DecisionV1',scopes=dict(
        A_CT6_temporal=times,
        B_experimental_reference=dict(status='blocked',reason='central displaced A95 resistance has no original-interval root below nominal7mm; original Accepted epsilon trace also lacks displacement',failed_roots_mm=failed_roots,
            proposed_shifted_sink='independent candidate study only; not an enabled native policy',physical_evolution_admitted=False),
        C_Richards_Norris=dict(operator=tops,history_and_budgets=histories,regression='pass' if regression else 'fail',
            frozen_reference_integration=dict(status='pass' if not qerrors else 'fail',reasons=qerrors),scope='d0 synthetic flat implementation; not the displaced experimental inlet'),
        D_geometry_capacity=dict(preflight='pass' if not domain_errors else 'fail',reasons=domain_errors,evolution_capacity='blocked_3000_guard',
            counts={k:v['native_count'] for k,v in dom['cases'].items() if k.startswith('flat')},coupled_collar_to_roof='not_implemented_or_qualified'),
        E_full_domain_preservation=dict(status='not_executed_blocked',required=['complete physical reference','supported scalable assembly/solver'],admitted=False),
        F_terrain=dict(status='not_admitted',scope='no experimental hill evolution, developed precursor, or physical validation'),
        resources=res),historical_CT4_CT6_decisions='unchanged',curved_model_transfer_eligible=False)
    # A process completion is distinct from each mathematical/numerical gate.
    fine_budget=physical_budgets(fine_audit)
    decision['fresh_CT6_finer_member']=dict(physical_completion=fine['physical_completion'],diagnostic_status=fine['diagnostic_status'],
        budget_primitives=fine_budget,captured_k_matrix=fine_audit['captured_k_matrix'],mass=fine_audit['mass'])
    decision['required_implemented_checks']='pass' if tops['status']==histories['status']==res['status']==fine_budget['status']=='pass' and not domain_errors and not qerrors and regression else 'failed_or_incomplete'
    decision['CT6_new_temporal_status']=times['status']
    for name,body in [('decision.json',decision),('fit-and-reference.json',refs),('top-operator.json',tp),('top-history.json',hist),('domain-preflight.json',dom),('CT6-temporal.json',time),('regression.json',reg),
        ('frozen-top-integration.json',quad),('consumed-history-boundaries.json',boundary_checks),('CT6-finer-physical-budgets.json',fine_audit)]:write_json(a.output/name,body)
    write_json(a.output/'CT6-finer-member.json',{k:v for k,v in fine.items() if k not in ('trajectory','zero')})
    write_json(a.output/'CT6-finer-native-identity.json',{k:fine['trajectory'].get(k) for k in ('report_sha256','state_sha256','initial_sha256','checkpoint','history_sha256','stored_slots','execution')})
    for model,body in primitives.items():write_json(a.output/('physical-budgets-'+model+'.json'),body)
    g=['git','-c','core.preloadIndex=false','-c','index.threads=1','-c','core.filemode=false']
    changes=subprocess.check_output(g+['diff','--name-only','4e253fe03f4c6460040f9359e57f578426f2a7b8..HEAD'],cwd=s,text=True).splitlines()
    provenance=dict(schema='CT7SourceRuntimeArchiveV1',reviewed_start='4e253fe03f4c6460040f9359e57f578426f2a7b8',
        source_revision=subprocess.check_output(g+['rev-parse','HEAD'],cwd=s,text=True).strip(),source_dirty=subprocess.check_output(g+['status','--porcelain'],cwd=s,text=True),
        changed_source_sha256={name:sha(s/name) for name in changes if (s/name).is_file()},build=build,
        governing_CT6_runtime=json.loads((a.root/'ct6-finer-p4/request-manifest.json').read_text())['request']['runtime'],
        archives={str(f):dict(results_sha256=sha(f/'results.json'),request_sha256=sha(f/'request-manifest.json')) for f in [a.root/'reference-v1',a.top,a.history,a.domain,a.temporal,a.regression,a.build,a.root/'ct6-finer-p4',a.root/'ct6-finer-p4-audit']},
        specification_sha256={n:sha(s/'docs/literature'/n) for n in ('A93_A95_E94_E95_neutral_ABL_EP_AWF_combined.tex','accepted_awf_appendix_A93_A95_E94_E95.tex')},
        threads=THREADS,resource_receipts=closed,collection_resource='attach closed resource-final.json and execution.json from this family before final report/qualification')
    provenance['host_storage']=subprocess.check_output(['findmnt','-n','-o','SOURCE,FSTYPE,TARGET','/media/jay/新加卷1'],text=True).strip()
    provenance['host_kernel']=subprocess.check_output(['uname','-a'],text=True).strip()
    failures={}
    for name in ('build-v1/build-1.log','domain-v2/flat-p3/count.log','top-v2/results.json','top-frozen-quadrature-v1/results.json',
        'top-frozen-quadrature-v2/A95-p3-f96/frozen.json','top-v6/results.json'):
        path=a.root/name;failures[name]=dict(path=str(path),sha256=sha(path),tail=path.read_text()[-3000:],historical_status='retained_failed_preparation; no physical trajectory was admitted')
    write_json(a.output/'failed-preparations.json',failures)
    write_json(a.output/'provenance.json',provenance)
    m.seal(a.output,[str(f.relative_to(a.output)) for f in a.output.glob('*.json') if f.name not in ('request-manifest.json','completed-artifacts.json')])
    print('CT7 implemented checks',decision['required_implemented_checks'],'CT6 temporal',times['status'],'physical/full-domain admission blocked',flush=True)
    return 0 if decision['required_implemented_checks']==times['status']=='pass' else 1

if __name__=='__main__':raise SystemExit(main())
