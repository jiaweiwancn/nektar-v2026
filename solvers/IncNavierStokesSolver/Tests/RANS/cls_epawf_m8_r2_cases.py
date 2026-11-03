#!/usr/bin/env python3
"""Frozen M8-R2 input preparation; never launches a solver.

Only the two preregistered outer refinements are available. The actual
archived unforced control supplies every physical input and function.
The inner mesh, degrees, native GLL rule, samplers and certificate stay fixed.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m8_r1_cases import (X, QUADRATURE, cost as tensor_cost,
    parameters, solver_info, physics, validate_original)


REVIEWED_SOURCE = '5186d63bd825aadd198488d93cff02fc9496a232'
INNER = [0., .01, .1, .25, .5, 1., 2.]
RADIAL = {
    'outer6': INNER+[2+j/3 for j in range(1, 7)],
    'outer9': INNER+[2+2*j/9 for j in range(1, 10)],
}
EXPECTED_SIZES = {
    'outer6': [592, 592, 592, 304, 42, 42],
    'outer9': [736, 736, 736, 448, 42, 60],
}
EXPECTED_UNKNOWNS = {'outer6': 2165, 'outer9': 2759}
INNER_LAWS = ('Linear', 'SaturatedR2')
TARGET_DT = .0025
ALPHA_LEVELS = [3200., 1600., 800., 400.]
POLICY = 'ScalarFeasibleR2'
INPUT_CONTROLS = {
    'RANSWallNonlinearPolicy': POLICY,
    'RANSWallPhysicalRangeGuard': 'Off',
    'RANSWallSolveDiagnostics': 'Failure',
    'RANSWallSolveDiagnosticReplay': 'Off',
    'RANSWallSolveDiagnosticFile': 'case.nonlinear.json',
}
STAGE_TESTS = ('Off', 'MatchedR2', 'FrozenScalarR2')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cost(representation):
    if representation not in RADIAL:
        raise ValueError('Only frozen outer6 and outer9 refinements are available')
    result = tensor_cost(RADIAL[representation])
    if result['field_sizes'] != EXPECTED_SIZES[representation] or result['unknowns'] != EXPECTED_UNKNOWNS[representation]:
        raise RuntimeError('Frozen tensor-space count changed')
    baseline_unknowns, baseline_cells = 1571, 45
    result.update(representation=representation,
        FD_volume_work_ratio_to_R1_resolved=result['unknowns']/baseline_unknowns*result['cells']/baseline_cells,
        dense_factor_work_ratio_to_R1_resolved=(result['unknowns']/baseline_unknowns)**3,
        one_sided_old_witness_radius=8/3,
        derivative_error_planning_factor=1/(8 if representation == 'outer6' else 27),
        measured_timing_or_error_bound=False)
    return result


def normalized_physics(root):
    """Protected data, apart from explicitly authorized R2 experiment controls."""
    result = physics(root, omit_quadrature=True)
    result['parameters'].pop('TIMESTEP')
    result['parameters'].pop('RANSWALLMAXITERATIONS')
    for name in ('RANSWallInnerStress', 'RANSWallStageTest', *INPUT_CONTROLS):
        result['solver_info'].pop(name.upper(), None)
    return result


def session(original, destination, representation, quadrature='base',
            inner_stress='SaturatedR2', stage_test='Off'):
    if quadrature not in QUADRATURE or inner_stress not in INNER_LAWS or stage_test not in STAGE_TESTS:
        raise ValueError('Unregistered quadrature or inner law')
    specification = cost(representation)
    source = ET.parse(original).getroot()
    validate_original(source)
    root = deepcopy(source)
    tensor_mesh(root, RADIAL[representation], 3, 2, 6., 1., points=[8, 8], xknots=X)
    set_parameter(root, 'TimeStep', repr(TARGET_DT))
    set_parameter(root, 'NumSteps', 0)
    set_parameter(root, 'RANSWallMaxIterations', 12)
    for name, value in zip(('RANSWallTangentialOrder', 'RANSWallNormalOrder'), QUADRATURE[quadrature]):
        set_parameter(root, name, value)
    set_info(root, 'RANSWallInnerStress', inner_stress)
    set_info(root, 'RANSWallStageTest', stage_test)
    for name, value in INPUT_CONTROLS.items():
        set_info(root, name, value)
    if normalized_physics(root) != normalized_physics(source):
        raise RuntimeError('R2 preparation changed a protected physical input or function')
    # Explicitly retain zero-step initialization. Stage probes consume this
    # same initial XML and control alpha externally without advancing history.
    if int(parameters(root)['NUMSTEPS']) != 0 or not specification['within_dense_guard']:
        raise RuntimeError('Invalid R2 zero-step input')
    write(root, destination)
    return Path(destination)


def prepare(original, output, representation, laws=('SaturatedR2',)):
    """Prepare one representation only; outer9 execution remains conditional."""
    if not laws or len(set(laws)) != len(laws) or any(law not in INNER_LAWS for law in laws):
        raise ValueError('Requested law list is not registered')
    original, output = Path(original).resolve(), Path(output)
    source = ET.parse(original).getroot()
    validate_original(source)
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for law in laws:
        for quadrature in ('base', 'high'):
            name = representation+'-'+law.lower()+'-'+quadrature
            directory = output/name
            directory.mkdir()
            xml = session(original, directory/'case.xml', representation, quadrature, law)
            rows.append(dict(name=name, representation=representation,
                inner_stress=law, quadrature=quadrature,
                source_orders=list(QUADRATURE[quadrature]), native_GLL_points=[8, 8],
                xml=str(xml), xml_sha256=sha(xml), cost=cost(representation),
                steps=0, target_dt=TARGET_DT, scalar_policy=POLICY,
                physical_guard='Off for complete initial observation; no qualification implied'))
    payload = dict(schema=2, milestone='M8-R2', role='E representation preparation',
        reviewed_source=REVIEWED_SOURCE, generator_sha256=sha(__file__),
        original=dict(path=str(original), sha256=sha(original)),
        protected_physics=normalized_physics(source), cases=rows,
        no_solver_invoked=True, allowed_representation_order=['outer6', 'outer9 only if needed'],
        maximum_representation_changes=2, maximum_unknowns=3000,
        native_quadrature=QUADRATURE, independent_same_state_orders=[[8, 16, 8], [8, 24, 8], [12, 24, 12]],
        fixed_stage=dict(dt=TARGET_DT, time=TARGET_DT, alpha_levels=ALPHA_LEVELS,
            original_target_alpha=400., history='same actual represented initial state at every level',
            auxiliary_native_history_advance=False, maximum_Newton_iterations_per_level_ceiling=16),
        authorization='Preparation only; parent pins producing source/runtime and authorizes all Linux execution',
        qualification='Complete finite-set physical and fixed-box preflight is required separately; no new box or model retuning')
    (output/'preflight-plan.json').write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n')
    return payload


def prepare_stages(original, output, representation, frozen_scalar=False):
    """Prepare the matched native stage pair or one conditional localization.

    The equation system owns the unchanged initial history and actual implicit
    stage. NumSteps=0 prevents an accidental evolution loop. No auxiliary
    parameter is represented by a changed dt or a restarted history.
    """
    original, output = Path(original).resolve(), Path(output)
    source = ET.parse(original).getroot()
    validate_original(source)
    output.mkdir(parents=True, exist_ok=False)
    stage_test = 'FrozenScalarR2' if frozen_scalar else 'MatchedR2'
    rows = []
    for law in (('SaturatedR2',) if frozen_scalar else INNER_LAWS):
        name = representation+'-'+law.lower()+'-'+stage_test.lower()
        directory = output/name
        directory.mkdir()
        xml = session(original, directory/'case.xml', representation, 'base', law, stage_test)
        rows.append(dict(name=name, xml=str(xml), xml_sha256=sha(xml),
            representation=representation, inner_stress=law,
            stage_test=stage_test, source_orders=list(QUADRATURE['base']),
            native_GLL_points=[8, 8], cost=cost(representation), steps=0,
            target_dt=TARGET_DT, target_time=TARGET_DT, scalar_policy=POLICY))
    payload = dict(schema=2, milestone='M8-R2', generator_sha256=sha(__file__),
        reviewed_source=REVIEWED_SOURCE, no_solver_invoked=True,
        original=dict(path=str(original), sha256=sha(original)),
        protected_physics=normalized_physics(source), cases=rows,
        stage_test=stage_test, alpha_levels=[400.] if frozen_scalar else ALPHA_LEVELS,
        target_alpha=400., maximum_Newton_updates_per_level=12,
        maximum_levels=1 if frozen_scalar else 4,
        conditional_flow=('Only if both coupled target roots unresolved; one diagnostic, not a coupled root'
            if frozen_scalar else 'On auxiliary failure skip remaining auxiliary levels and attempt400 once from last converged scalar-feasible iterate or initial; no failed400 retry'),
        actual_history='identical initial h0 at every level and for both laws',
        time_and_native_history_advanced=False,
        authorization='Preparation only; root owns operator/action gate, coverage decision, source/runtime pinning and every Linux solve')
    (output/'stage-plan.json').write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n')
    return payload


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--representation', choices=tuple(RADIAL), required=True)
    parser.add_argument('--laws', choices=('saturated', 'matched'), default='saturated')
    parser.add_argument('--family', choices=('preflight', 'stage', 'frozen-scalar'), default='preflight')
    args = parser.parse_args()
    result = (prepare(args.original, args.output, args.representation,
        INNER_LAWS if args.laws == 'matched' else ('SaturatedR2',)) if args.family == 'preflight' else
        prepare_stages(args.original, args.output, args.representation, args.family == 'frozen-scalar'))
    print(json.dumps(dict(cases=[dict(name=row['name'], unknowns=row['cost']['unknowns'])
        for row in result['cases']], solver_invoked=False)))
