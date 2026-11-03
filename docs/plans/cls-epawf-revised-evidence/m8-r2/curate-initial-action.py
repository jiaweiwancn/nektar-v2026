"""Select compact fields from retained Linux initial-action evidence only."""
import argparse
import hashlib
import json
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--epoch',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
def pin(path):
    return dict(path=str(path.resolve()),bytes=path.stat().st_size,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())
def load(path):
    return json.loads(path.read_text())
first=a.epoch/'initial-jacobian-001'
second=a.epoch/'initial-k-diagnosis-001'
paths=[first/'receipt.json',first/'validation.json',second/'receipt.json',second/'validation.json']
before=[pin(path) for path in paths]
r,v,dr,dv=map(load,paths)
assert r['success'] and dr['success'] and v['analysis_success']
assert before[1]['sha256']=='52491ec3974757cf4b177b0502360171e5715448c68ab1ceb19b8bf564e38977'
assert before[3]['sha256']=='2224a4e387a92d0a8a160e817bb9265534a3d7291384faf025d23ddea8b6ecf0'
columns={key:value for key,value in v['continuity_columns'].items()
    if key!='all_column_absolute_physical_errors'}
pressures=[{key:item[key] for key in ('kind','absolute_physical_error',
    'fixed_1e_minus7_gate_met','nonlinear_source_and_trace_norm')} for item in v['pressure_actions']]
matrix=[dict(path=path,**value) for path,value in r['inputs'].items() if path.endswith('.solve1.J0.bin')]
assert len(matrix)==1
result=dict(schema=1,milestone='M8-R2 initial auxiliary J0 action diagnosis',
    method='Metadata and existing result-field selection only; no scientific/native calls',
    curator=pin(Path(__file__)),authoritative_evidence=before,
    producing_source_commit=r['producing_source_commit'],producing_solver=r['producing_solver'],
    producing_runtime_sections=r['producing_runtime_sections'],
    initial_analysis_source=r['source_commit'],component_analysis_source=dr['source_commit'],
    selected_input_pins_unchanged_during_analyses=True,native_runtime_after_check_pending=True,
    decoded_initial_state_sha256=v['initial_h0_float64_LE_sha256'],
    source_bound_h0_scope=v['initial_h0_binding'],matrix=matrix[0],
    alpha=v['auxiliary_alpha'],time=v['physical_stage_time'],
    numerical_threads=r['numerical_threads'],analysis_CPU=r['cpu_affinity'],
    checks=dict(initial=v['checks'],component=dv['checks']),
    initial_derivative_gates_passed=r['initial_derivative_gates_passed'],
    original_independent_FD_gate=v['outer_k_action']['gate'],
    exact_affine_continuity_columns=columns,exact_analytic_pressure_actions=pressures,
    analytic_component_assumptions=dv['assumptions'],analytic_controls=dv['analytic_controls'],
    analytic_controls_met=dv['analytic_controls_met'],
    stored_J_minus_analytic=dv['stored_J_minus_analytic'],
    stored_J_roughness_support_error=dv['stored_J_roughness_support_error'],
    stored_J_outer_support_error=dv['stored_J_outer_support_error'],
    retained_FD_decomposition=dv['retained_FD_decomposition'],
    fixed_absolute_action_tolerance=dv['fixed_absolute_action_tolerance'],
    fixed_absolute_gate_met=dv['fixed_absolute_gate_met'],
    localization_classification=dv['localization_classification'],
    original_failed_absolute_gate_preserved=dv['original_failed_absolute_gate_preserved'],
    rigorous_derivative_defect_proven=dv['rigorous_derivative_defect_proven'],
    numerical_causality_claim=False,target_stage_status_assessed=False,
    completed_native_nonlinear_record_used=False,physical_trajectory_qualification=False,M8_closed=False,
    scope='Only the initial Linear-inner-law first auxiliary matrix. Later iterates and original target require separate records. Large vectors/matrices remain on Linux.')
assert before==[pin(path) for path in paths]
a.output.parent.mkdir(parents=True,exist_ok=True)
with a.output.open('x') as stream:
    json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
print(json.dumps(pin(a.output)))
