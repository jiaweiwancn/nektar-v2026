"""Select compact fields from the retained Linux original-target J0 check."""
import argparse
import hashlib
import json
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--run',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
def pin(path):
    return dict(path=str(path.resolve()),bytes=path.stat().st_size,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())
paths=[a.run/'receipt.json',a.run/'validation.json']
before=[pin(path) for path in paths]
r,v=[json.loads(path.read_text()) for path in paths]
assert r['success'] and r['native_calls']==0 and v['native_calls']==0
assert before[1]['sha256']=='29396c6623849f0257869e77624eafb6cae6637ab9c36f0709c1e90b98522d76'
keys=('schema','scope','target_alpha','auxiliary_alpha','physical_stage_time',
    'initial_h0_float64_LE_sha256','native_field_descriptors_sha256',
    'protected_physics_sha256','target_initial_guess_binding','actual_log_iteration0',
    'checks','mass_change_control_met','target_minus_auxiliary_minus_exact_mass_change',
    'target_J_minus_transferred_independent_analytic','target_support_errors',
    'target_initial_k_action_gate_met','original_initial_derivative_gates_passed',
    'original_component_absolute_gate_met','original_failed_gates_preserved',
    'pressure_action_transfer','dense_matrix_files','fixed_absolute_action_tolerance',
    'fresh_residual_endpoints','no_rigorous_floating_point_error_proof',
    'rigorous_derivative_defect_proven','nonlinear_terminal_cause_assessed',
    'target_stage_status_assessed','final_native_history_slots_assessed',
    'coupled_retry_requested','physical_trajectory_qualification','M8_closed')
result={key:v[key] for key in keys}
result.update(milestone='M8-R2 actual original-target initial action check',
    method='Metadata/result-field selection only; no scientific/native calls',
    curator=pin(Path(__file__)),authoritative_evidence=before,
    producing_source_commit=r['producing_source_commit'],producing_solver=r['producing_solver'],
    producing_runtime_sections=r['producing_runtime_sections'],
    analysis_source_commit=r['source_commit'],numerical_threads=r['numerical_threads'],
    cpu_affinity=r['cpu_affinity'],native_calls=0,
    exact_affine_continuity_columns={key:value for key,value in
        v['exact_affine_continuity_columns'].items() if key!='all_column_absolute_physical_errors'})
assert before==[pin(path) for path in paths]
a.output.parent.mkdir(parents=True,exist_ok=True)
with a.output.open('x') as stream:
    json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
print(json.dumps(pin(a.output)))
