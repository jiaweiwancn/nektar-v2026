# SPDX-License-Identifier: MIT
"""WJ-6 two-trajectory uncertainty and complete, evidence-bound phase gates."""
import copy
import math
import numpy as np
from wj3_measure import FIELDS
from wj5_measure import boundary_envelope
from wj_qualification import aggregate, component, finite
from wj_immutable import verify_artifacts
from run_cls_epawf_picard import require


def individual_arithmetic(packet):
    """Retain100x; account each trajectory separately before adding the pair.

    WJ-3's stored comparison estimate has an internal factor2 for a pair.
    Here the native propagation/formation estimates give an individual bound,
    with its own evaluation sensitivity; no asymmetric trajectory is replaced
    by the other trajectory's terminal estimate.
    """
    require(packet['safety_factor']==100,'arithmetic reserve changed')
    raw=np.array([packet['raw_one_step_physical_error_estimate'][f] for f in FIELDS])
    growth=np.array(packet['accumulated_propagation_factors'])
    evaluation=np.array(packet['evaluation_sensitivity'])
    values=100*(growth*raw+evaluation)
    values[5]+=2*values[2]/3
    factors=np.array(packet['sampled_conversion_factors'])
    sampled=values*factors;sampled[5]=sampled[4]+2*sampled[2]/3
    require(finite(values.tolist()+sampled.tolist()) and min(values)>=0 and min(sampled)>=0,'invalid individual arithmetic envelope')
    result=copy.deepcopy(packet)
    result['estimated_comparison_floor']=dict(zip(FIELDS,map(float,values)))
    result['estimated_sampled_comparison_floor']=dict(zip(FIELDS,map(float,sampled)))
    for key in ('estimated_comparison_floor','estimated_sampled_comparison_floor'):
        result[key]['velocity']=math.hypot(result[key]['u'],result[key]['v'])
    result['comparison_accounting']='individual100x*(propagation*raw+evaluation); pair adds both separate trajectory envelopes'
    return result


def sum_envelopes(a,b):
    require(set(a)==set(b) and all(set(a[f])==set(b[f]) for f in a),'incomplete paired uncertainty')
    value={f:{norm:a[f][norm]+b[f][norm] for norm in a[f]} for f in a}
    require(finite(value) and all(x>=0 for r in value.values() for x in r.values()),'invalid paired uncertainty')
    return value


def paired_envelope(ta,tb):
    floors=[];boundary=[];coverage=[]
    for tag,t in (('a',ta),('b',tb)):
        record=dict(initial=t['initial'],trajectory=[])
        for c in t.get('trajectory',[]):
            q=dict(c)
            if c.get('arithmetic'):q['arithmetic']=individual_arithmetic(c['arithmetic'])
            record['trajectory'].append(q)
        f,b,p=boundary_envelope({tag:record})
        floors.append(f);boundary.append(b);coverage.append(p)
    return sum_envelopes(*floors),sum_envelopes(*boundary),dict(
        contributions=coverage,scope='sum of separately estimated trajectories; all actual row maps, own units and full-domain pressure gauge; factor100 retained')


def required_accuracy(records,expected):
    parts={'coverage':component('PASS' if set(records)==set(expected) else 'NOT_RUN','every declared observation required')}
    for key,c in records.items():
        parts[str(key)]=c.get('errors',{}).get('accuracy_gate',component('NOT_RUN','missing accuracy packet'))
    return aggregate(parts)


def fresh_startup_gate(packet):
    if not isinstance(packet,dict) or any(k not in packet for k in ('exact_physical_primary_gauge_history','native_labels_differ','checkpoints')):
        return component('NOT_RUN','missing fresh-startup compatibility packet')
    if packet['exact_physical_primary_gauge_history'] is not True or packet['native_labels_differ'] is not True or len(packet['checkpoints'])!=2:
        return component('FAIL','incompatible physical startup or native dt labels')
    verify_artifacts(packet['checkpoints'])
    return component('PASS','equal physical primary/gauge/history, distinct native dt labels; exact saved bytes verified')


def common_extension_authorized(policies):
    return set(policies)=={'EP','Budget'} and all(
        p.get('primary',{}).get('status')=='PASS' and
        p.get('transactions',{}).get('status')=='PASS' and
        p.get('extension_resource',{}).get('status')=='PASS' for p in policies.values())
