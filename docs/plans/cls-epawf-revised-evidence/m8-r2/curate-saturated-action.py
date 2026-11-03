"""Select compact fields from the retained Linux SaturatedR2 initial J0 check."""
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
def pick(value,keys):
    return {key:value[key] for key in keys}
paths=[a.run/'receipt.json',a.run/'validation.json']
before=[pin(path) for path in paths]
r,v=[json.loads(path.read_text()) for path in paths]
assert r['success'] and v['analysis_success'] and r['native_calls']==0
assert before[1]['sha256']=='5f5df2c4254e5a305880d8a61cefe38db21462da8aacaf8daa86c136286268ab'
action=v['outer_k_action']
levels=[]
for level in action['levels']:
    row={key:value for key,value in level.items() if key not in ('plus','minus')}
    for side in ('plus','minus'):
        row[side]=pick(level[side],('branches','endpoint_source','evaluated',
            'spatial_residual_float64_LE_sha256','state_float64_LE_sha256'))
    levels.append(row)
result=dict(schema=1,milestone='M8-R2 SaturatedR2 initial auxiliary J0 action check',
    method='Metadata/result-field selection only; no scientific/native calls',
    curator=pin(Path(__file__)),authoritative_evidence=before,
    producing_source_commit=r['producing_source_commit'],producing_solver=r['producing_solver'],
    producing_runtime_sections=r['producing_runtime_sections'],analysis_source_commit=r['source_commit'],
    selected_cpu=r['selected_cpu'],cpu_affinity=r['cpu_affinity'],numerical_threads=r['numerical_threads'],
    checks=v['checks'],native_calls=0,scope=v['scope'],alpha=v['auxiliary_alpha'],time=v['physical_stage_time'],
    initial_h0_float64_LE_sha256=v['initial_h0_float64_LE_sha256'],matrix=v['dense_matrix_files'],
    initial_derivative_gates_passed=v['initial_derivative_gates_passed'],
    outer_k_action_gate=action['gate'],outer_k_action_levels=levels,
    exact_affine_continuity_columns={key:value for key,value in
        v['continuity_columns'].items() if key!='all_column_absolute_physical_errors'},
    pressure_actions=[pick(row,('kind','absolute_physical_error','fixed_1e_minus7_gate_met',
        'nonlinear_source_and_trace_norm')) for row in v['pressure_actions']],
    solved_path_assessed=False,target_stage_status_assessed=False,
    physical_trajectory_qualification=False,M8_closed=False,
    interpretation='The original absolute finite-difference action gate remains failed. Finer residual subtraction has no rigorous error bound; this is not a new independent analytic component diagnosis or a nonlinear-causality claim.')
assert before==[pin(path) for path in paths]
a.output.parent.mkdir(parents=True,exist_ok=True)
with a.output.open('x') as stream:
    json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
print(json.dumps(pin(a.output)))
