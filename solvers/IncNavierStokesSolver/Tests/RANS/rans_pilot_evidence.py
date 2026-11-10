# SPDX-License-Identifier: MIT
"""Identity-checked supplements; never rewrite historical numerical receipts."""
import hashlib,json,math
from pathlib import Path
from rans_request_manifest import RequestManifest
from rans_preservation_oracle import Reference

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

EPSILON_BUDGET_SCALED_TOLERANCE=1e-10

def epsilon_budget_gate(budgets,ref,width=2.):
    scale=width*(ref.Z-ref.H)*(ref.A/ref.H)*ref.us/ref.H
    checks=[]
    for b in budgets:
        value=b['free_constant_test_defect'];scaled=abs(value)/scale
        checks.append(dict(time=b['time'],alpha=b['alpha'],absolute_defect=value,scale=scale,scaled_defect=scaled,
            scaled_limit=EPSILON_BUDGET_SCALED_TOLERANCE,status='pass' if math.isfinite(scaled) and scaled<=EPSILON_BUDGET_SCALED_TOLERANCE else 'fail'))
    return dict(status='pass' if checks and all(x['status']=='pass' for x in checks) else 'fail',checks=checks,
        normalization='outer physical volume * epsilon(H) * u_star/H; actual time-discrete load minus prescribed-row reaction',
        scope='constant test only; distinct from the alpha-zero spatial residual and refreshed dependent trace lag')

def comparison_body(root,path):
    d=load(path)
    required=('status','time','coarse','fine','temporal_control','fields','actual_temporal_initial_and_boundary_inputs_identical')
    if any(k not in d for k in required):raise ValueError('incomplete required comparison body: '+str(path))
    if set(d['fields'])!={'u','v','k','epsilon_outer','Pi'}:raise ValueError('missing fieldwise comparison')
    cases=[];identities=[]
    for role in ('coarse','fine','temporal_control'):
        folder=root/'campaigns'/d[role];manifest=load(folder/'request-manifest.json');req=manifest['request'];admission=RequestManifest(folder,req,True)
        endpoint=folder/('T-'+format(d['time'],'.4f'))
        if not admission.completed(endpoint):raise ValueError('unsealed comparison endpoint: '+str(endpoint))
        metric=load(endpoint/'metrics.json');initial=load(folder/'T-0.0016/metrics.json')['initialization_receipt']
        if abs(metric['time']-d['time'])>1e-13 or metric['status']!='pass':raise ValueError('comparison terminal time/status mismatch')
        keys=('raw_initial_sha256','effective_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256')
        identities.append({k:initial[k] for k in keys})
        cases.append(dict(role=role,label=d[role],request_sha256=manifest['sha256'],request_manifest_file_sha256=sha(folder/'request-manifest.json'),
            model=req['model'],policy=req['policy'],wall=req['wall'],descriptor=req['descriptor'],
            space={k:req[k] for k in ('knots','degree','pressure_degree','columns','quadrature')},dt=req['dt'],endpoint_time=metric['time'],
            actual_consumed_initial_reference_boundary={k:initial[k] for k in keys},raw_effective_max_difference=initial['raw_effective_max_difference'],
            target_generated_boundary_max_difference=initial.get('source_target_boundary_max_difference'),
            runtime=req['runtime'],endpoint=dict(source=str(endpoint),metrics_sha256=sha(endpoint/'metrics.json'),samples_sha256=sha(endpoint/'endpoint-samples.json'),
                native_report_sha256=sha(endpoint/'case.evolution.json'),state_sha256=metric['final_state_sha256'],history=metric['history_identity']),
            analytical_reference_errors=metric['reference_error'],change_from_actual_initial=metric['change_from_initial']))
    if identities[0]!=identities[2] or not d['actual_temporal_initial_and_boundary_inputs_identical']:raise ValueError('actual consumed temporal inputs differ')
    if any(c['raw_effective_max_difference']!=0 for c in cases):raise ValueError('effective initialization changed')
    return dict(d,source=str(path),source_sha256=sha(path),cases=cases,
        limits=dict(spatial_error_ratio=.85,temporal_fraction=.1,temporal_scaled_absolute_floor=1e-7,spatial_error_floor=1e-10),
        identity_checks=dict(actual_temporal_initial_reference_and_boundary=True,coarse_fine_same_continuum_reference=True,
            coefficient_byte_identity_required_only_for_same_space=True))

def required_comparisons(root,directory=None):
    directory=directory or root/'natural-schedule';records={}
    for model in ('E','A95','A93'):
        for wall in ('EP','Accepted'):
            for t in (.0016,.016):
                name=f'{model}-N-{wall}-space-{t}';path=directory/(name+'.json')
                try:records[name]=comparison_body(root,path)
                except Exception as e:records[name]=dict(status='fail',failure=repr(e),source=str(path))
    return records
