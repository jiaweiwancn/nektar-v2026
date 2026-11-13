#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT6 immutable manufactured inputs. No CT4/CT5 physical data are changed."""
import argparse, json, os, xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_info, set_parameter, write
from ct_accepted_feasibility import parameters
from rans_request_manifest import sha

X = [-1.2, -.8, -.4, 0., .4, .8, 1.2]
RADIAL = [0., .05, .125, .25, .5, .75, 1.]
QUADRATURE = [32, 8, 12, 12, 12]


def make(folder, runtime, model='A95', degree=3, dt=.002, steps=0,
         amplitude=.005, small=False, static=False, restart=None):
    """Each changed space/timestep starts independently; never edit histories."""
    if degree not in (3, 4) or model not in ('A93', 'A95', 'E94E95'):
        raise ValueError('outside CT6 bounded space/model contract')
    path = case_xml(folder, model, 'Natural', 'Accepted', steps=steps, dt=dt,
                    descriptor_probe=runtime/'RANSNewModelsProbe', degree=degree,
                    radial_knots=[0., .5, 1.] if small else RADIAL,
                    xknots=[-1.2, 0., 1.2] if small else X,
                    quadrature=QUADRATURE)
    root = ET.parse(path).getroot()
    conditions = root.find('CONDITIONS')
    for function in list(conditions.findall('FUNCTION')):
        conditions.remove(function)
    settings = dict(RANSWallCurvedCapability='CT6StaticTransmissionV1' if static else 'CT6ManufacturedV1',
                    RANSWallBoundaryPolicy='Closed' if static else 'AcceptedManufacturedOpenV1',
                    RANSWallMap='Gaussian', RANSWallManufactured='Off' if static else 'CT6ContinuousV1',
                    RANSWallReference='Functions', RANSWallFunctionCoordinates='Chart',
                    RANSOpenTop='DirichletCompatibleV1', RANSOpenBackflow='RejectV1',
                    RANSWallPimpleSnapshots='On', RANSWallEvaluationReuse='Off')
    for k, v in settings.items():
        set_info(root, k, v)
    p = parameters(model)
    for k, v in dict(RANSWallHeight=.5, RANSWallSampleHeight=.5, RANSWallTailStart=.5,
                     RANSWallRoughness=.05, Kinvis=1.5e-5, RANSWallFrictionVelocity=.4,
                     RANSWallChi=0., RANSWallEta=0., RANSWallMapAmplitude=amplitude,
                     RANSWallMapLengthX=1., RANSWallMapHeightAmplitude=.03 if static else 0.,
                     RANSWallMapHeightWaveX=3.141592653589793/2,
                     K0=p['q']*.16, Qstar=p['q'], Kappa=p['kappa'],
                     RANSOpenVelocityScale=.2).items():
        set_parameter(root, k, v)
    set_info(root, 'RANSWallRestart', str(restart) if restart else 'None')
    set_info(root, 'RANSWallCheckpoint', str(folder/'checkpoint'))
    set_info(root, 'RANSWallPimpleSnapshotDirectory', str(folder/'pimple-states'))
    write(root, path)
    construction = dict(schema='CT6ContinuousManufacturedCaseV1', model=model,
        wall='AcceptedBudgetAWF', policy='Natural', provider=settings['RANSWallManufactured'],
        capability=settings['RANSWallCurvedCapability'], coordinates='chart (q,r), Cartesian U and gradients',
        x=[-1.2, 0., 1.2] if small else X, radial=[0., .5, 1.] if small else RADIAL,
        degree=degree, pressure_degree=degree-2, quadrature=QUADRATURE,
        dt=dt, steps=steps, H=.5, nu=1.5e-5, roughness=.05,
        map_amplitude=amplitude, map_length=1., model_parameters=p,
        pressure='absolute Pi; traction-determined level; no physical mean constraint',
        configuration_sha256=sha(path), manufactured_forcing=True,
        physical_experimental_case=False)
    write_json(folder/'construction.json', construction)
    return path


def grid():
    # Both exact interface sides and genuine interiors; no ambiguous H row.
    points, sides = [], []
    for q in (-1.1, -.9, -.7, -.3, 0., .3, .7, .9, 1.1):
        for r in (0., .025, .075, .2, .4, .65, .9, 1.):
            points.append([q, 0., r]); sides.append(0)
        for side in (-1, 1):
            points.append([q, 0., .5]); sides.append(side)
    return points, sides


def request(small=False):
    g, sides = grid()
    return dict(operation='functions', export_basis=True, export_test_basis=True,
        export_constitutive=True, export_resolution=True, export_tensor_validation=True,
        export_manufactured=True, export_open_verification=small,
        export_mixed_space=small, grid=g, grid_sides=sides)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runtime', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--model', choices=['A93', 'A95', 'E94E95'], default='A95')
    p.add_argument('--degree', type=int, choices=[3, 4], default=3)
    p.add_argument('--dt', type=float, default=.002)
    p.add_argument('--steps', type=int, default=0)
    p.add_argument('--small', action='store_true')
    p.add_argument('--static', action='store_true')
    a = p.parse_args()
    make(a.output, a.runtime, a.model, a.degree, a.dt, a.steps,
         small=a.small, static=a.static)
    write_json(a.output/'initial-request.json', request(a.small))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
