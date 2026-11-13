# SPDX-License-Identifier: MIT
"""Recompute CT7 scopes without turning a diagnostic into physical admission."""
import math
from ct_accepted_budget_gate import physical_budgets

TOP_CASES={'A95-p3','A93-p3','E94E95-p3','A95-p4'}
TOP_EXTRA={'pressure_offset','state_independent_epsilon_flux','displaced_original_wall','unsupported_cache'}

def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)

def top(receipt):
    reasons=[];checks=receipt.get('checks',{})
    if set(checks)!=TOP_CASES|TOP_EXTRA:reasons.append('missing_or_foreign_top_control')
    for name in TOP_CASES:
        row=checks.get(name,{})
        if row.get('status')!='pass':reasons.append(name+':native_operator_not_completed')
        defects=row.get('defects',{})
        expected={'all_row','boundary_advection','boundary_traction','boundary_diffusion',
            *('frozen_'+str(i) for i in (0,2,3)),*('residual_'+str(i) for i in (0,2,3))}
        if set(defects)!=expected:reasons.append(name+':missing_operator_primitive')
        if any(not finite(x) or x<0 or x>2e-10 for x in defects.values()):reasons.append(name+':operator_accuracy')
        weak=row.get('continuous_native_test_maximum')
        if not finite(weak) or weak<0 or weak>2e-10:reasons.append(name+':continuous_weak_accuracy')
        if row.get('additional_pressure_null_modes')!=0 or not finite(row.get('pressure_minimum_singular_ratio')) or row['pressure_minimum_singular_ratio']<=1e-10:
            reasons.append(name+':pressure_rank')
        if not finite(row.get('constant_pressure_action')) or row['constant_pressure_action']<=1e-8:reasons.append(name+':pressure_level_not_fixed')
        for key in ('zero_epsilon_flux_negative_signal','reversed_shear_negative_signal'):
            if not finite(row.get(key)) or row[key]<=1e-8:reasons.append(name+':nondiscriminating_negative:'+key)
        ini=row.get('initialization_receipt',{})
        if not ini.get('open_explicit_implicit_agreement') or not ini.get('explicit_boundary_sha256') or ini.get('explicit_boundary_sha256')!=ini.get('implicit_boundary_sha256'):
            reasons.append(name+':consumed_boundary_identity')
    offset=checks.get('pressure_offset',{})
    for key in ('coefficient_defect','free_residual_defect'):
        if not finite(offset.get(key)) or offset[key]<0 or offset[key]>2e-10:reasons.append('pressure_offset:'+key)
    flux=checks.get('state_independent_epsilon_flux',{})
    if not finite(flux.get('datum')) or flux['datum']>=0 or not finite(flux.get('wrong_evolving_diffusivity_signal')) or flux['wrong_evolving_diffusivity_signal']<=1e-8:
        reasons.append('prescribed_physical_epsilon_flux')
    for name,word in [('displaced_original_wall','displacement'),('unsupported_cache','caching')]:
        if word not in checks.get(name,{}).get('reason','').lower():reasons.append('wrong_expected_rejection:'+name)
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,scope='small d0 synthetic top operator; not experimental reference or full-domain preservation')

def history(receipt,primitives):
    reasons=[]
    expected={'restart','rollback','old_boundary_history','pressure_level_evolution'}
    if set(receipt.get('checks',{}))!=expected:reasons.append('missing_or_foreign_history_control')
    for name in expected:
        if receipt.get('checks',{}).get(name,{}).get('status')!='pass':reasons.append('history_control:'+name)
    if set(receipt.get('models',{}))!={'A93','A95','E94E95'}:reasons.append('missing_or_foreign_model_smoke')
    budgets={}
    for model in ('A93','A95','E94E95'):
        row=receipt.get('models',{}).get(model,{})
        binding=primitives.get(model,{}).get('budget_binding',{})
        if row.get('model')!=model or not row.get('resolved_descriptor','').startswith(model+':OriginalHighRe-v1') or ':resolved='+model+':OriginalHighRe-v1' not in row.get('native_signature',''):
            reasons.append(model+':foreign_native_model_configuration')
        if binding.get('model')!=model or binding.get('resolved_descriptor')!=row.get('resolved_descriptor') or any(binding.get(k)!=row.get(v) or not binding.get(k)
            for k,v in [('native_report_sha256','report_sha256'),('state_sha256','state_sha256'),('configuration_sha256','configuration_sha256')]):
            reasons.append(model+':budget_native_identity_mismatch')
        if row.get('completion')!='pass' or not finite(row.get('time')) or abs(row['time']-2e-6)>1e-15:reasons.append(model+':incomplete_history')
        budgets[model]=physical_budgets(primitives.get(model,{}))
        if budgets[model]['status']!='pass':reasons.append(model+':physical_budget_primitives')
        mass=row.get('mass',{});net=mass.get('outward_net');scale=mass.get('scale')
        if not finite(net) or not finite(scale) or scale<=0 or abs(net)/scale>1e-8:reasons.append(model+':mass')
        matrix=row.get('captured_k_matrix',{});error=matrix.get('maximum_defect');scale=matrix.get('scale')
        if not finite(error) or error<0 or not finite(scale) or scale<=0 or error/scale>2e-10:reasons.append(model+':actual_frozen_matrix')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,budgets=budgets,scope='two tiny unforced synthetic steps, restart/rollback/pressure invariance; not a physical preservation campaign')

def resources(receipts):
    reasons=[]
    if not receipts:reasons.append('missing_closed_resource_receipts')
    for name,r in receipts.items():
        if any(not finite(r.get(key)) or r[key]>limit for key,limit in [('maximum_project',16),('maximum_aggregate',72)]):reasons.append(name+':worker_limit')
        if r.get('project_status')!='pass' or r.get('aggregate_status')!='pass' or r.get('control_helpers_counted')!=2:reasons.append(name+':resource_receipt_incomplete')
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,receipts=receipts)

def temporal(receipt):
    rows={};reasons=[]
    if set(receipt.get('orders',{}))!={'U','w','k','outer_epsilon','Pi'}:reasons.append('missing_or_foreign_temporal_field')
    for name in ('U','w','k','outer_epsilon','Pi'):
        r=receipt.get('orders',{}).get(name,{});x=r.get('D002_001');y=r.get('D001_0005')
        if not finite(x) or not finite(y) or min(x,y)<0:rows[name]=dict(status='fail',reason='missing_or_nonfinite_primitive');reasons.append(name);continue
        order=math.log2(x/y) if min(x,y)>1e-6 else None
        status=('pass' if .8<=order<=1.2 else 'fail') if order is not None else 'unresolved_floor'
        rows[name]=dict(D002_001=x,D001_0005=y,order=order,status=status,floor=1e-6,band=[.8,1.2])
        if status=='fail' or (name=='Pi' and status!='pass'):reasons.append(name)
    return dict(status='pass' if not reasons else 'failed_or_unresolved',reasons=reasons,orders=rows,absolute_pressure=True)
