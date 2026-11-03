"""Frozen M8-R1 representation sessions; contains no solver invocation.

Both representations copy the actual archived control XML. Only NumSteps and
the explicitly listed native geometry/expansion metadata may change. Physical
inputs, functions, model selections, time-step and acceptance policy stay exact.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_info, set_parameter, write

BASELINE = 'ab5c3c6fee19761e8d88692e8e6c331e5889234e'
X = [-4., -3., -2., -1., 1., 2.]
RADIAL = {
    'original': [0., 2., 4.],
    'resolved': [0., .01, .1, .25, .5, 1., 2., 8/3, 10/3, 4.],
}
QUADRATURE = {'base': (8, 16), 'high': (12, 24)}
DIAGNOSTIC_INFO = {
    'RANSWallPhysicalRangeGuard': 'R1DiagnosticV1',
    'RANSWallSolveDiagnostics': 'Failure',
    'RANSWallSolveDiagnosticReplay': 'Off',
    'RANSWallSolveDiagnosticFile': 'case.nonlinear.json',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parameters(root):
    return {node.text.split('=', 1)[0].strip().upper(): node.text.split('=', 1)[1].strip()
            for node in root.findall('CONDITIONS/PARAMETERS/P')}


def solver_info(root):
    return {node.get('PROPERTY').upper(): node.get('VALUE')
            for node in root.findall('CONDITIONS/SOLVERINFO/I')}


def physics(root, omit_quadrature=False):
    """Exact physical strings; explicitly selected quadrature is separate."""
    result = dict(parameters=parameters(root), solver_info=solver_info(root),
                  functions=[dict(name=node.get('NAME'), entries=[dict(entry.attrib) for entry in node])
                             for node in root.findall('CONDITIONS/FUNCTION')],
                  time_integration=[[node.tag, node.text] for node in
                                    root.find('CONDITIONS/TIMEINTEGRATIONSCHEME')])
    result['parameters'].pop('NUMSTEPS')
    if omit_quadrature:
        for key in ('RANSWALLTANGENTIALORDER', 'RANSWALLNORMALORDER'):
            result['parameters'].pop(key)
    return result


def topology(root):
    if root.find('GEOMETRY').get('DIM') != '2':
        raise ValueError('M8-R1 preserves the actual 2D archived control')
    vertices = [list(map(float, node.text.split())) for node in root.findall('GEOMETRY/VERTEX/V')]
    x = sorted({point[0] for point in vertices})
    r = sorted({point[1] for point in vertices})
    expansions = root.findall('EXPANSIONS/E')
    velocity = [node for node in expansions if 'u' in node.get('FIELDS', '').split(',')]
    pressure = [node for node in expansions if node.get('FIELDS') == 'p']
    if len(velocity) != 1 or len(pressure) != 1:
        raise ValueError('Ambiguous archived tensor expansions')
    modes = list(map(int, velocity[0].get('NUMMODES').split(',')))
    pressure_modes = list(map(int, pressure[0].get('NUMMODES').split(',')))
    points = list(map(int, velocity[0].get('NUMPOINTS').split(',')))
    if modes != [4, 4] or pressure_modes != [2, 2] or points != [8, 8]:
        raise ValueError('M8-R1 expects the actual archived p3/p1, GLL8 control')
    if x != X or r != RADIAL['original']:
        raise ValueError('The source is not the original failed flow radial mesh')
    return dict(x=x, radial=r, degree=3, pressure_degree=1, native_GLL_points=points)


def cost(radial):
    nx, ni, no = len(X)-1, sum(b <= 2 for b in radial[1:]), sum(a >= 2 for a in radial[:-1])
    p, pp = 3, 1
    nv = (nx*p+1)*((ni+no)*p+1)
    ne = (nx*p+1)*(no*p+1)
    pi, po = (nx*pp+1)*(ni*pp+1), (nx*pp+1)*(no*pp+1)
    size = 3*nv+ne+pi+po+1
    return dict(dimension=2, degree=p, pressure_degree=pp, x=X, radial=radial,
                cells=nx*(ni+no), inner_cells=nx*ni, outer_cells=nx*no,
                field_sizes=[nv, nv, nv, ne, pi, po], gauge=1, unknowns=size,
                dense_guard=3000, within_dense_guard=size <= 3000,
                central_FD_residual_calls_per_full_Jacobian=2*size,
                single_dense_matrix_bytes=8*size*size,
                source_volume_outer_Gauss_points=nx*no*8**2,
                estimate_scope='Exact tensor-space counts before native verification; assembly/factor scaling is indicative, not a measured timing')


def validate_original(root):
    mesh = topology(root)
    p, info = parameters(root), solver_info(root)
    expected = {'RANSWALLHEIGHT': 2., 'RANSWALLROUGHNESS': .1, 'RANSWALLSAMPLEHEIGHT': 1.,
                'RANSWALLTAILSTART': 1., 'RANSWALLETA': 0., 'RANSWALLCHI': 1.,
                'RANSWALLCUBICBETA': .4, 'RANSWALLREGULARIZATIONGAMMA': 0.,
                'RANSWALLREGULARIZATIONSCALE': 1., 'RANSWALLPOSITIVITYTHETA': .2,
                'RANSWALLPOSITIVITYDELTA': .1, 'TIMESTEP': .005, 'TIME': 0.,
                'RANSWALLTRANSITIONSTART': -3., 'RANSWALLTRANSITIONEND': -2.25,
                'RANSWALLBASELINESTART': -2., 'RANSWALLSENSORDISTANCEA': 2.5,
                'RANSWALLSENSORDISTANCEB': 3.5, 'RANSWALLSENSORON': .02, 'RANSWALLSENSOROFF': .1,
                'RANSWALLTANGENTIALORDER': 8, 'RANSWALLNORMALORDER': 16}
    for key, value in expected.items():
        if float(p[key]) != value:
            raise ValueError('Archived physical/control parameter differs: '+key)
    for key, value in {'RANSWALLVARIANT': 'Regularized', 'RANSWALLDESTRUCTION': 'Capped',
                       'RANSWALLDESTRUCTIONSMOOTHING': 'Linear', 'RANSWALLREFERENCE': 'Neutral',
                       'RANSWALLMAP': 'Flat', 'RANSWALLTRANSITION': 'Upstream',
                       'RANSWALLRESTART': 'None', 'RANSWALLCHECKPOINT': 'None'}.items():
        if info[key] != value:
            raise ValueError('Archived solver selection differs: '+key)
    names = {node.get('NAME') for node in root.findall('CONDITIONS/FUNCTION')}
    if names != {'InitialPerturbationInner', 'InitialPerturbationOuter'}:
        raise ValueError('Unforced control contains an unexpected function')
    for node in root.findall('CONDITIONS/FUNCTION/E'):
        expression = node.get('VALUE')
        if expression != '0' and not expression.startswith(('(0)', '(-0)')):
            raise ValueError('Archived control has a nonzero perturbation amplitude')
    return mesh


def session(original, destination, representation, quadrature='base'):
    if representation not in RADIAL or quadrature not in QUADRATURE:
        raise ValueError('Exactly two representations and base/high quadrature are preregistered')
    source = ET.parse(original).getroot()
    validate_original(source)
    root = deepcopy(source)
    set_parameter(root, 'NumSteps', 0)
    if representation == 'resolved':
        tensor_mesh(root, RADIAL[representation], 3, 2, 6., 1., points=[8, 8], xknots=X)
    for name, value in zip(('RANSWallTangentialOrder', 'RANSWallNormalOrder'), QUADRATURE[quadrature]):
        set_parameter(root, name, value)
    if physics(root, True) != physics(source, True):
        raise RuntimeError('Zero-step representation changed a physical input or solver policy')
    write(root, destination)
    return Path(destination)


def prepare(original, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    original = Path(original).resolve()
    root = ET.parse(original).getroot()
    mesh = validate_original(root)
    records = []
    base = cost(RADIAL['original'])
    for representation in ('original', 'resolved'):
        for quadrature in ('base', 'high'):
            name = representation+'-'+quadrature
            directory = output/name
            directory.mkdir()
            xml = session(original, directory/'case.xml', representation, quadrature)
            record = dict(name=name, representation=representation, quadrature=quadrature,
                          source_quadrature=list(QUADRATURE[quadrature]), xml=str(xml),
                          xml_sha256=sha(xml), cost=cost(RADIAL[representation]), steps=0)
            record['cost']['FD_volume_assembly_ratio_to_original'] = (
                record['cost']['unknowns']/base['unknowns']*record['cost']['cells']/base['cells'])
            record['cost']['dense_factor_work_ratio_to_original'] = (record['cost']['unknowns']/base['unknowns'])**3
            records.append(record)
    payload = dict(schema=1, milestone='M8-R1', phase='C representation preflight',
        baseline_reviewed=BASELINE, no_solver_invoked=True,
        original=dict(path=str(original), sha256=sha(original), mesh=mesh),
        physics=physics(root, True), cases=records, generator_sha256=sha(__file__),
        protocol=dict(native_initializations=4, trajectories=0, permitted_mesh_choices=['original', 'resolved'],
            native_quadrature_levels=QUADRATURE, native_GLL_points=[8, 8],
            quadrature_audits=[[8, 16, 8], [8, 24, 8], [12, 24, 12]],
            quadrature_interpretation='Independent loads at the identical initial coefficients; no additional CFD time map',
            qualification_stop='First declared physical-range exit; preserve state. Finite coverage alone is not a continuum or polynomial-everywhere certificate.',
            preserved='Physical domain, H, h, samplers, transition, model/source parameters, eta, dt, tolerances and solver policy',
            authorization='Solver calls await the parent phase-A archived-ledger gate'))
    (output/'preflight-plan.json').write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n')
    return payload


def trajectory_session(original, destination, dt):
    """Prepare one of two root-authorized diagnostic controls; never run it.

    The range guard records the already exited local region at t=0, permits
    explicitly diagnostic continuation, and stops its first declared finite
    physical violation before history advance. Initial conditions/history are
    rebuilt for each dt; no stored checkpoint or old evolving field is used.
    """
    if dt not in (.005, .0025):
        raise ValueError('Only the two approved independent time steps are allowed')
    session(original, destination, 'resolved', 'base')
    root = ET.parse(destination).getroot()
    for name, value in DIAGNOSTIC_INFO.items():
        set_info(root, name, value)
    set_parameter(root, 'TimeStep', repr(dt))
    set_parameter(root, 'NumSteps', 12 if dt == .005 else 24)
    before, after = physics(ET.parse(original).getroot(), True), physics(root, True)
    before['parameters'].pop('TIMESTEP')
    after['parameters'].pop('TIMESTEP')
    for name in DIAGNOSTIC_INFO:
        if name.upper() in before['solver_info']:
            raise ValueError('Historical control unexpectedly contains an opt-in guard property')
        after['solver_info'].pop(name.upper())
    if before != after:
        raise RuntimeError('Diagnostic control changed an unapproved physical/function/solver-policy input')
    return write(root, destination)


def prepare_trajectories(original, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    original = Path(original).resolve()
    validate_original(ET.parse(original).getroot())
    records = []
    for name, dt in (('resolved-dt005', .005), ('resolved-dt0025', .0025)):
        directory = output/name
        directory.mkdir()
        xml = trajectory_session(original, directory/'case.xml', dt)
        records.append(dict(name=name, xml=str(xml), xml_sha256=sha(xml), dt=dt,
            requested_steps=12 if dt == .005 else 24, maximum_time=.06,
            cost=cost(RADIAL['resolved']), fresh_initialization=True, restart='None', checkpoint='None',
            interpretation='Explicit diagnostic continuation; outside frozen local gradient region from t=0'))
    result = dict(schema=1, milestone='M8-R1', phase='C diagnostic control preparation',
        baseline_reviewed=BASELINE, original=dict(path=str(original), sha256=sha(original)),
        generator_sha256=sha(__file__), no_solver_invoked=True, cases=records,
        diagnostic_solver_info=DIAGNOSTIC_INFO,
        qualification_from_t0=False, local_range_exit_expected_at_t0=True,
        frozen_region=dict(gradient_entry_radius=.001, k_over_K0=[.99, 1.01],
            epsilon_over_epsilon0=[.98, 1.02], minimum_lambda=.5),
        rule='Exactly two sequential root-launched controls; preserve first declared physical violation/state/history; no physical-input, mesh, region, tolerance or solver-policy retuning',
        authorization='Preparation only. Root owns build pinning, native evolution allocation and execution.')
    (output/'trajectory-plan.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--family', choices=('preflight', 'trajectory'), default='preflight')
    args = parser.parse_args()
    result = (prepare if args.family == 'preflight' else prepare_trajectories)(args.original, args.output)
    print(json.dumps(dict(cases=[dict(name=row['name'], unknowns=row['cost']['unknowns'])
                                for row in result['cases']], solver_invoked=False)))
