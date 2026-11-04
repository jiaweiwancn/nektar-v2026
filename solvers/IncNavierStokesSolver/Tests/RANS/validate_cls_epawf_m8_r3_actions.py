#!/usr/bin/env python3
"""Linux-only independent complete VALUE-action checks for R3 fixed states."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode = True
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
for _name in THREADS: os.environ[_name] = '1'
import numpy as np
from cls_epawf_m8_r3_action_oracle import ValueOracle, CD, LD
from cls_epawf_m8_r1_nonlinear import Audit

ABSOLUTE_ACTION_GATE = 1e-7
COMPLEX_STEPS = ('1e-20', '1e-30')


def require(value, message):
    if not value:
        raise ValueError(message)


def identity(path):
    path = Path(path).resolve(); result = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''): result.update(part)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=result.hexdigest())


def load(path):
    return json.loads(Path(path).read_text(), parse_int=lambda x: -0.0 if x == '-0' else int(x),
        parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))


def receipt_binds(value, pin):
    if isinstance(value, dict):
        if value.get('path') == pin['path'] and value.get('sha256') == pin['sha256']: return True
        for name,item in value.items():
            if name.endswith('_path') and item == pin['path'] and value.get(name[:-5]+'_sha256') == pin['sha256']:
                return True
        return any(receipt_binds(item,pin) for item in value.values())
    return isinstance(value,list) and any(receipt_binds(item,pin) for item in value)


def validate_case(raw, source):
    require(raw.get('schema') == 1 and raw.get('nonlinear_policy') == 'ScalarFeasibleR3' and
        raw.get('action_policy') == 'VerifiedR3Action' and raw.get('eta') == 0 and raw.get('ranks') == 1,
        'Changed native action policy/forcing/MPI scope')
    audit = Audit(raw, source)
    constraints = audit.constraints(raw['constraints'])
    state,history = np.asarray(raw['state'],dtype=CD),np.asarray(raw['history'],dtype=CD)
    oracle = ValueOracle(raw, cutoff_relative_step='1e-4')
    base = oracle.value(state,history,raw['alpha'])
    base_error = np.asarray(np.real(base),float)-np.asarray(raw['residual'],float)
    base_norms = audit.norm_record(base_error,constraints)
    # This checks construction of the same native value problem. Near-root
    # resolution is separately reported in physical norms; it is not relabelled
    # a roundoff floor or a successful 2e-12 root check.
    component_error=float(np.max(np.abs(base_error)/np.maximum(1.,np.abs(raw['residual'])),initial=0.))
    numerical_failures=[]
    if component_error>2e-10:
        numerical_failures.append(dict(gate='independent_complete_value_componentwise',error=component_error,tolerance=2e-10))
    matrix_path = Path(raw['matrix_file'])
    require(matrix_path.stat().st_size == 8*len(state)**2, 'Native matrix byte length differs')
    matrix = np.memmap(matrix_path,dtype=np.float64,mode='r',shape=(len(state),len(state)))
    directions = []
    for row in raw['directions']:
        direction = np.asarray(row['direction'],dtype=CD)
        require(direction.shape == state.shape and np.all(np.isfinite(direction)), 'Invalid native input direction')
        native = np.asarray(row['stage_action'],float)
        matrix_error = audit.norm_record(np.asarray(matrix@np.asarray(direction.real,float))-native,constraints)
        if matrix_error['maximum']>ABSOLUTE_ACTION_GATE:
            numerical_failures.append(dict(direction=row['name'],gate='native_matrix_output_consistency',error=matrix_error))
        levels,actions = [],[]
        for h_string in COMPLEX_STEPS:
            h=LD(h_string)
            action=np.asarray(np.imag(oracle.value(state+CD(1j)*h*direction,history,raw['alpha']))/h,dtype=LD)
            error=audit.norm_record(np.asarray(action,float)-native,constraints)
            levels.append(dict(complex_step=h_string,absolute_physical_error=error,
                absolute_gate=ABSOLUTE_ACTION_GATE,gate_met=bool(error['maximum']<=ABSOLUTE_ACTION_GATE)))
            actions.append(action)
        ladder=audit.norm_record(np.asarray(actions[0]-actions[1],float),constraints)
        if not all(level['gate_met'] for level in levels):
            numerical_failures.append(dict(direction=row['name'],gate='independent_absolute_full_value_action',levels=levels))
        if ladder['maximum']>ABSOLUTE_ACTION_GATE:
            numerical_failures.append(dict(direction=row['name'],gate='complex_step_ladder_absolute_action',error=ladder))
        # Independent exact mass uses the old Jacobi physical mass, external
        # column mask and actual constrained output-row mask.
        effective=np.asarray(direction,dtype=CD).copy()
        for fixed in oracle.external: effective[fixed]=0
        exact=np.zeros(len(state),dtype=LD)
        for block in range(4):
            ids=np.arange(oracle.fields[block].size)+oracle.off[block]
            exact[ids]=LD(raw['alpha'])*(np.asarray(oracle.mass[block],dtype=LD)@effective[ids].real)
        for fixed in constraints: exact[fixed]=0
        mass_error=audit.norm_record(np.asarray(exact,float)-np.asarray(row['exact_mass_action'],float),constraints)
        if mass_error['maximum']>ABSOLUTE_ACTION_GATE:
            numerical_failures.append(dict(direction=row['name'],gate='independent_exact_mass_action',error=mass_error))
        trace_rows=np.asarray(sorted(oracle.off[3]+index for index in oracle.trace_ids),int)
        directions.append(dict(name=row['name'],levels=levels,complex_step_ladder_absolute_error=ladder,
            native_matrix_output_consistency=matrix_error,independent_exact_mass_error=mass_error,
            independent_matching_trace_action_maximum=float(np.max(np.abs(actions[-1][trace_rows]),initial=0.)),
            no_large_load_subtraction=True))
    require(directions,'No actual native action directions')
    return dict(success=not numerical_failures,scope='Supplied fixed state and supplied input directions only; no nonlinear path or global constitutive qualification',
        checks=audit.checks,errors=audit.errors,action_gate=ABSOLUTE_ACTION_GATE,
        independent_value_operator='Separate clongdouble Modified_A rows, indexed tensors, moving normal quadrature/cutoff/normalization and constrained trace value projection',
        base_value_error_physical_norms=base_norms,base_componentwise_value_error=component_error,
        base_native_residual_physical_norms=audit.norm_record(raw['residual'],constraints),
        base_independent_residual_physical_norms=audit.norm_record(np.asarray(np.real(base),float),constraints),
        directions=directions,fixed_absolute_gate_met=not numerical_failures,numerical_gate_failures=numerical_failures,
        oracle_activity_and_component_ladder=dict(oracle.statistics,
            counting_scope='Repeated finite oracle evaluations; counts are not unique physical points or volume measures'),
        generalized_branch_scope=dict(outer_invariant='Bound native S/R selection, independently checked away from a binary64 tie',
            Yap='Bound native active/inactive selection; inactive derivative at exact knot',
            zero_norm='Selected zero norm subgradient for outer invariant; linear inner response at positive k',
            limiter='Inactive saturated side at threshold',cap='Uncapped side at b=delta',
            moving_panel_ties='Stable base ordering with complex movement retained',
            joint_zero_k_zero_strain='Rejected; requires a separate admissible directional value check',
            proof_of_smoothness_across_switches=False),
        numerical_precision=dict(longdouble_epsilon=str(np.finfo(LD).eps),complex_dtype=str(np.dtype(CD))),
        production_derivative_calls=0,production_residual_calls=0,trajectory_calls=0)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--expected-analysis-source',required=True)
    p.add_argument('--producer-receipt',type=Path,required=True);p.add_argument('--runtime-pin',type=Path,required=True)
    p.add_argument('--expected-producer-source',required=True)
    p.add_argument('--case',action='append',required=True,help='Unique NAME=REPORT.json fixed-state binding')
    p.add_argument('--required-case',action='append',default=[])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();require(platform.system()=='Linux','Scientific R3 action verification belongs to authoritative Linux')
    require(not a.output.exists(),'Use a fresh validation output')
    source=a.source.resolve();git=lambda *words:subprocess.check_output(['git','-C',str(source),*words],text=True).strip()
    require(git('rev-parse','HEAD')==a.expected_analysis_source and not git('status','--porcelain'),'Wrong/unclean analysis source')
    started=time.monotonic();runtime_pin=identity(a.runtime_pin);producer_pin=identity(a.producer_receipt)
    runtime=load(a.runtime_pin);producer=load(a.producer_receipt)
    require(runtime.get('success') is True and runtime.get('clean') is True and
        runtime.get('source_commit')==a.expected_producer_source,'Unbound producing native runtime')
    key='build/RANSNonpressureActionProbe';require(key in runtime['runtime_resolved_paths'],'No producing action-probe image')
    image=identity(runtime['runtime_resolved_paths'][key]);require(image['sha256']==runtime['runtime_artifacts'][key],'Producing action probe changed')
    inputs={};cases=[];names=set();failures=[]
    for relative,expected in runtime['source_inputs'].items():
        if Path(relative).suffix.lower() in ('.cpp','.h','.hpp'):
            current=identity(source/relative);require(current['sha256']==expected,'Analysis/producer C++ differs: '+relative)
    for entry in a.case:
        name,path=entry.split('=',1);require(name and name not in names,'Repeated or empty fixed-state case name');names.add(name)
        report_pin=identity(path);raw=load(path);matrix_pin=identity(raw['matrix_file'])
        require(receipt_binds(producer,report_pin) and receipt_binds(producer,matrix_pin),
            'Producer receipt does not bind this exact native report/matrix')
        inputs[report_pin['path']]=report_pin;inputs[matrix_pin['path']]=matrix_pin
        row=dict(name=name,report=report_pin,matrix=matrix_pin)
        try:
            row['independent']=validate_case(raw,source)
            if not row['independent']['success']:
                failures.append(dict(case=name,type='NumericalGateFailure',message='Unchanged independent action/value/mass gate failed',
                    gates=row['independent']['numerical_gate_failures']))
        except Exception as error:
            row['independent']=dict(success=False,failure=dict(type=type(error).__name__,message=str(error)))
            failures.append(dict(case=name,**row['independent']['failure']))
        cases.append(row)
    for path,pin in inputs.items():require(identity(path)==pin,'Native evidence changed during independent audit')
    require(identity(a.runtime_pin)==runtime_pin and identity(a.producer_receipt)==producer_pin and
        identity(image['path'])==image,'Producing receipt/runtime/program changed during audit')
    require(git('rev-parse','HEAD')==a.expected_analysis_source and not git('status','--porcelain'),'Analysis source changed')
    imports={str(Path(module.__file__).resolve()):identity(module.__file__) for module in list(sys.modules.values())
        if getattr(module,'__file__',None) and str(Path(module.__file__).resolve()).startswith(str(source))}
    required_covered=bool(a.required_case) and set(a.required_case).issubset(names)
    result=dict(schema=1,success=not failures,analysis_success=not failures,
        source_commit=a.expected_analysis_source,producing_source_commit=a.expected_producer_source,
        source_clean_start=True,source_clean_end=True,source_commit_end=a.expected_analysis_source,
        producer_receipt=producer_pin,runtime_pin=runtime_pin,producing_program=image,
        validator=identity(__file__),imported_source_files=imports,numerical_threads={name:os.environ[name] for name in THREADS},
        cases=cases,failures=failures,absolute_action_tolerance=ABSOLUTE_ACTION_GATE,
        fixed_absolute_gate_met=not failures,required_cases=a.required_case,
        required_state_action_gates_met=bool(required_covered and not failures),
        active_dependencies_covered=None,
        active_dependencies_scope='Review measured fixture activity explicitly; passing passive directions does not establish active sensor/limiter/cap coverage',
        native_calls=0,production_derivative_calls=0,production_residual_calls=0,trajectory_calls=0,
        completed_utc=datetime.now(timezone.utc).isoformat(),seconds=time.monotonic()-started)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(success=result['success'],cases=len(cases),failures=failures,output=str(a.output))))
    return 0 if result['success'] else 1


if __name__=='__main__':
    raise SystemExit(main())
