# SPDX-License-Identifier: MIT
"""Complete V2 producer contract; the production importer remains strict."""
from pathlib import Path
import math
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require,sha,state_hash

NATIVE_DIV_TOL=1e-12

def divergence_gate(native,limit):
    explicit=native['native_divergence_GL16'];regional=native['audit'].get('strong_continuity',[])
    require(explicit['order']==16 and len(explicit['regional_L2'])==2,'missing explicit GL16 divergence')
    value=math.hypot(*explicit['regional_L2'])
    require(abs(value-explicit['L2'])<=NATIVE_DIV_TOL,'GL16 regional aggregate differs')
    identity=len(regional)>=2 and abs(math.hypot(*regional[:2])-value)<=NATIVE_DIV_TOL
    return dict(value=value,regional=explicit['regional_L2'],audit_identity=identity,passed=identity and value<=limit+NATIVE_DIV_TOL,tolerance=NATIVE_DIV_TOL)

def complete_receipt(report,model,limit,probe_runtime,solver_runtime,extra=None,extra_binding=None):
    report=Path(report);native=load(report);a=native['audit'];v=native['extra_sampling'] if extra is None else extra
    K0=model['K0'];base=native['outer_initializer_contract']['base_sha256'];state=native['state_sha256'];div=divergence_gate(native,limit)
    require(state_hash(native['state'])==state and state_hash(native['fixed_reference'])==base==native['reference_sha256'],'receipt state/base/reference binding')
    require(state_hash(native['fixed_defect'])==native['defect_sha256'],'receipt defect binding')
    require(v['state_sha256']==state and v['cold_base_sha256']==base and v['model_specification']==model,'receipt extra/model identity differs')
    gates=dict(native_admissible=a['admissible']==[1] and a['accepted_awf_hard_admissible']==[1],budget=a['accepted_awf_scalar_domain'][8]<=1e-10,algebraic=max(a['coupled_norms'][4:])<1e-9,optimizer_interior_margin=min(a['accepted_awf_covariance'])>=1e-4*K0-1e-12,extra_covariance=min(v['inner_minimum_covariance'],v['outer_minimum_covariance'])>=1e-4*K0-1e-12,native_GL16=div['passed'],extra_divergence=v['strong_divergence_L2']<=limit+1e-8,velocity_L2=v['velocity_change_L2']<=.12*(1+1e-10),sampled_velocity_increment=v['maximum_sampled_velocity_increment']<=.5+1e-10,scalars=bool(v['scalar_positive_finite'] and v['scalar_envelopes_passed'] and v['integrals_complete']))
    gates['strong_divergence']=gates['native_GL16'] and gates['extra_divergence']
    return dict(version='CompleteNativeV2ReceiptV1',passed=all(gates.values()),gates=gates,K0=K0,divergence_limit=limit,native_GL16=div,state_sha256=state,base_sha256=base,probe_output_file=str(report),probe_output_sha256=sha(report),audit=a,extra_sampling=v,extra_binding=extra_binding,model_specification=model,reference_sha256=native['reference_sha256'],defect_sha256=native['defect_sha256'],geometry_identity=native['geometry_identity'],inventory_identity=native['inventory_identity'],patch_signature=native['patch_signature'],validation_probe=probe_runtime,solver_runtime=solver_runtime)

def validate_receipt(receipt,require_pass=True):
    require(receipt.get('version')=='CompleteNativeV2ReceiptV1','incompatible V2 receipt interpretation')
    for key in ('passed','gates','state_sha256','base_sha256','probe_output_file','probe_output_sha256','model_specification','divergence_limit','validation_probe','solver_runtime','extra_sampling'):
        require(key in receipt,'incomplete V2 receipt: '+key)
    for key in ('extra_covariance','strong_divergence','native_GL16','optimizer_interior_margin'):
        require(key in receipt['gates'],'incomplete V2 receipt gate: '+key)
    binding=receipt.get('extra_binding')
    if binding:
        f=Path(binding['file']);require(sha(f)==binding['sha256'],'V2 extra receipt changed');require(load(f).get('fresh_sampling')==receipt['extra_sampling'],'V2 extra sampling differs')
    report=Path(receipt['probe_output_file']);require(sha(report)==receipt['probe_output_sha256'],'V2 native report changed')
    new=complete_receipt(report,receipt['model_specification'],receipt['divergence_limit'],receipt['validation_probe'],receipt['solver_runtime'],receipt['extra_sampling'],receipt.get('extra_binding'))
    for key in ('state_sha256','base_sha256','reference_sha256','defect_sha256','geometry_identity','inventory_identity','patch_signature'):
        require(receipt.get(key)==new[key],'V2 receipt binding differs: '+key)
    require(receipt['gates']==new['gates'] and receipt['passed']==new['passed'],'V2 receipt gates do not reflect native results')
    if require_pass:require(new['passed'],'V2 receipt is not qualified')
    return new
