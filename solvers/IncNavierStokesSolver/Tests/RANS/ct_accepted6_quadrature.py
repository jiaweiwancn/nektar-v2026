# SPDX-License-Identifier: MIT
"""Frozen-state CT6 integration uncertainty, separate from reconstruction.

The current-state residual norms cover momentum, scalars, continuity and
constraints. Complete frozen scalar load vectors provide an additional
direction-sensitive check. No changed load is used to advance a state.
"""
import math
from ct_accepted6_gate import finite
from ct_accepted_feasibility import parameters

CONTROLS = {'native', 'q48f16', 'q64f16', 'q64f32'}
BLOCKS = ('U', 'w', 'k', 'outer_epsilon', 'continuity_inner',
          'continuity_outer', 'unused_gauge', 'constraints')


def extract(packet, budget):
    o = packet['observation']
    return dict(norms=o['audit_alpha_zero']['coupled_norms'],
        volume=math.fsum(r['weight'] for r in o['volume']),
        frozen_loads={k:packet['physical_scalar_budget']['reconciliation']['assembled'][k]['native_raw_frozen_load']
                      for k in ('k', 'outer_epsilon')},
        budget_defects={k:r['frozen_stage_physical_defect']/r['physical_scale']
                        for k,r in budget['equations'].items()},
        actual_orders=o['actual_quadrature'])


def classify(rows, comparison, model):
    reasons, checks = [], {}
    if set(rows) != CONTROLS:
        return dict(status='incomplete', reasons=['missing_or_foreign_frozen_quadrature_control'], checks={})
    p = parameters(model);time_scale=.5/.2
    scales = dict(k=p['q']*.16, outer_epsilon=.4**3/(p['kappa']*.5))

    def limit(name):
        row=comparison.get('degree',{}).get(name,{})
        difference=row.get('spatial_difference',{}).get('scaled_weighted_RMS')
        if not finite(difference) or difference<0:
            reasons.append('missing_or_nonfinite_degree_difference:'+name)
            return 2e-8
        return max(2e-8,.1*difference)

    def offer(key, value, bound):
        checks[key]=dict(change=value,limit=bound,
                        status='pass' if finite(value) and 0<=value<=bound else 'fail')
        if checks[key]['status']!='pass':reasons.append('quadrature_uncertainty:'+key)

    base=rows['native']
    if not finite(base['volume']) or base['volume']<=0:
        return dict(status='fail',reasons=['invalid_physical_volume'],checks={})
    for label,row in rows.items():
        if len(row['norms'])!=len(BLOCKS):
            reasons.append('missing_residual_blocks:'+label);continue
        if label=='native':continue
        offer(label+':physical_volume',abs(row['volume']-base['volume'])/base['volume'],2e-11)
        # Native Norms uses u*=.4 and H=.5. Convert its rate scales to
        # the declared velocity .2 and reference time H/.2=2.5.
        factors=(4.,4.,2.,2.,2.,2.,1.,1.)
        for j,name in enumerate(BLOCKS):
            value=abs(row['norms'][j]-base['norms'][j])*factors[j]
            offer(label+':alpha_zero_mass_dual_norm:'+name,value,limit(name) if j<4 else 2e-8)
        for name,scale in scales.items():
            a,b=base['frozen_loads'][name],row['frozen_loads'][name]
            if len(a)!=len(b) or not a:
                reasons.append('frozen_scalar_row_identity:'+label+':'+name);continue
            maximum=max(abs(x-y) if finite(x) and finite(y) else math.inf for x,y in zip(a,b,strict=True))
            offer(label+':frozen_all_row_load:'+name,maximum*time_scale/(scale*base['volume']),limit(name))
        for name,value in row['budget_defects'].items():
            field='outer_epsilon' if name=='outer_epsilon' else 'k'
            offer(label+':physical_budget:'+name,abs(value-base['budget_defects'][name])*time_scale,limit(field))
    return dict(schema='CT6FrozenIntegrationUncertaintyV1',status='pass' if not reasons else 'fail',
        reasons=sorted(set(reasons)),checks=checks,actual_orders={k:r['actual_orders'] for k,r in rows.items()},
        reference_time=time_scale,floor=2e-8,fraction_of_measured_degree_difference=.1,
        scope='same actual coefficients/reference/live essential boundaries; current unanchored native mass-dual norm changes, complete frozen scalar row-load changes and physical source/storage/face-budget changes; not temporal or spatial convergence')
