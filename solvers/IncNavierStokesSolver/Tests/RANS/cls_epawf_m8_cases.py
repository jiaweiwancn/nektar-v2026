"""M8 sessions for explicitly selected CLS research variants.

The original generators remain unchanged.  Geometry, physical wall inputs and
time integration are inherited from the qualified native wall backend; only
named model parameters below select the approved M8 extension.  Numerical
execution and scientific interpretation belong on the Linux server.
"""
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True

from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m4_cases import _evolution, tensor_mesh
from cls_epawf_m6_cases import transition_session
from cls_epawf_m7_cases import perturbation_session as original_perturbation
from validate_cls_epawf_m7_transition_flow import session as original_flow


VARIANTS = ('Original', 'CubicScaled', 'Regularized')
DEFAULT_THETA, DEFAULT_DELTA = .2, .1


def apply_variant(root, variant='CubicScaled', beta=.4, gamma=0.,
                  theta=DEFAULT_THETA, delta=DEFAULT_DELTA, smooth=False,
                  gamma_scale=0.):
    if variant not in VARIANTS:
        raise ValueError('Unknown explicit M8 variant: '+str(variant))
    if not all(math.isfinite(value) for value in (beta, gamma, theta, delta, gamma_scale)):
        raise ValueError('Nonfinite M8 model parameter')
    if not 0 <= beta <= 1 or gamma < 0 or gamma_scale < 0 or not 0 < theta < 1 or not 0 < delta < 1:
        raise ValueError('M8 requires beta in [0,1], gamma>=0, theta/delta in (0,1)')
    if variant != 'Regularized' and (gamma != 0 or gamma_scale != 0):
        raise ValueError('A persistent gamma stress requires Regularized selection')
    if variant == 'Original' and beta != 1:
        raise ValueError('The Original control must retain beta=1')
    set_info(root, 'RANSWallVariant', variant)
    set_info(root, 'RANSWallDestruction', 'Original' if variant == 'Original' else 'Capped')
    set_info(root, 'RANSWallDestructionSmoothing', 'C2' if smooth else 'Linear')
    for name, value in {'RANSWallCubicBeta': beta,
                        'RANSWallRegularizationGamma': gamma,
                        'RANSWallRegularizationScale': gamma_scale,
                        'RANSWallPositivityTheta': theta,
                        'RANSWallPositivityDelta': delta,
                        'RANSWallDiagnostics': 1}.items():
        set_parameter(root, name, repr(value))
    return root


def perturbation_session(template, destination, variant='CubicScaled', beta=.4,
                         gamma=0., theta=DEFAULT_THETA, delta=DEFAULT_DELTA,
                         smooth=False, gamma_scale=0., maximum_iterations=80, **kwargs):
    original_perturbation(template, destination, **kwargs)
    root = ET.parse(destination).getroot()
    apply_variant(root, variant, beta, gamma, theta, delta, smooth, gamma_scale)
    set_parameter(root, 'RANSWallMaxIterations', maximum_iterations)
    return write(root, destination)


def depletion_session(template, destination, depletion=.95, **kwargs):
    """Positive polynomial source-repair mechanism test, not a covariance cure.

    At x=L/2,r=0 the inner k factor is 1-depletion; at the physical
    rK=H/2 sample it is 1-depletion/4.  The polynomial is represented exactly
    by p>=3, vanishes at side boundaries/H and never depletes outer epsilon.
    A depleted sampled-viscosity inner stress can be nonrealizable.  The
    native runner records that fact separately and never promotes this case
    to a robust-flow qualification merely because its source test passes.
    """
    if not 0 <= depletion < 1:
        raise ValueError('Depletion must keep the initial polynomial positive')
    kwargs.setdefault('kind', 'none')
    perturbation_session(template, destination, **kwargs)
    root = ET.parse(destination).getroot()
    if any(e.text and 'RANSWallHeight' in e.text for e in root.findall('CONDITIONS/PARAMETERS/P')):
        # The native M7 bounded family uses H=2 and L=8 unless explicitly
        # overridden by the generator's physical matching-height study.
        height = float(kwargs.get('height', 2.))
    else:
        raise ValueError('Missing native matching-height parameter')
    node = root.find("CONDITIONS/FUNCTION[@NAME='InitialPerturbationInner']/E[@VAR='k']")
    if node is None:
        raise ValueError('Missing explicit initial perturbation function')
    z = 'z' if kwargs.get('dimension', 2) == 3 else 'y'
    expression = f'-({depletion:.17g})*Qstar*0.16*4*(x/8)*(1-x/8)*(1-({z})/({height:.17g}))^2'
    if kwargs.get('dimension', 2) == 3:
        expression += '*(sin(PI*y/2))^2'
    node.set('VALUE', expression)
    return write(root, destination)


def manufactured_session(template, destination, variant='CubicScaled', beta=.4,
                         gamma=0., theta=DEFAULT_THETA, delta=DEFAULT_DELTA,
                         smooth=False, gamma_scale=0., **kwargs):
    """The zero-velocity M6 MMS remains exact for common cubic scaling.

    A Regularized stress changes the zero-velocity manufactured equation
    whenever Sref is nonzero.  It requires a separately derived body source,
    so a nonzero gamma must not silently reuse the original MMS forcing.
    """
    if gamma or gamma_scale:
        raise ValueError('Regularized MMS needs its independently derived persistent-stress body force')
    transition_session(template, destination, **kwargs)
    root = ET.parse(destination).getroot()
    apply_variant(root, variant, beta, gamma, theta, delta, smooth, gamma_scale)
    return write(root, destination)


def flow_session(template, destination, amplitude=.001, dt=.005, steps=1,
                 variant='Regularized', beta=.4, gamma=0.,
                 theta=DEFAULT_THETA, delta=DEFAULT_DELTA, smooth=False,
                 gamma_scale=0.,
                 degree=3, refine=1, hill=False, active=True,
                 maximum_iterations=80, checkpoint=None, restart=None):
    """Paired unforced finite-flow fixture with fixed physical wall inputs.

    The optional small hill is a numerical qualification fixture, not the
    Cao--Tamura physical geometry.  The original M7 pulse, ordinary neutral
    inflow, physical sensors and one-way upstream mask are retained.
    """
    original_flow(template, destination, amplitude, dt, steps)
    root = ET.parse(destination).getroot()
    x = [-4., -3., -2., -1., 1., 2.]
    z = [0., 2., 4.]
    if refine > 1:
        x = [a+(b-a)*i/refine for a, b in zip(x, x[1:])
             for i in range(refine)]+[x[-1]]
        z = [a+(b-a)*i/refine for a, b in zip(z, z[1:])
             for i in range(refine)]+[z[-1]]
    tensor_mesh(root, z, degree, 2, 6., 1., xknots=x)
    _evolution(root, steps, 2, dt, 0., False, restart, checkpoint)
    set_info(root, 'RANSWallMap', 'CosineHill' if hill else 'Flat')
    set_parameter(root, 'RANSWallMapAmplitude', .02)
    set_parameter(root, 'RANSWallMapLengthX', 1.)
    set_parameter(root, 'RANSWallChi', int(active))
    set_parameter(root, 'RANSWallMaxIterations', maximum_iterations)
    apply_variant(root, variant, beta, gamma, theta, delta, smooth, gamma_scale)
    return write(root, destination)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('destination', type=Path)
    p.add_argument('--variant', choices=VARIANTS, default='CubicScaled')
    p.add_argument('--beta', type=float, default=.4)
    p.add_argument('--gamma', type=float, default=0.)
    p.add_argument('--gamma-scale', type=float, default=0.)
    p.add_argument('--theta', type=float, default=DEFAULT_THETA)
    p.add_argument('--delta', type=float, default=DEFAULT_DELTA)
    p.add_argument('--smooth', action='store_true')
    p.add_argument('--steps', type=int, default=1)
    p.add_argument('--family', choices=('flat', 'depletion', 'flow'), default='flat')
    a = p.parse_args()
    generator = {'flat': perturbation_session, 'depletion': depletion_session,
                 'flow': flow_session}[a.family]
    generator(Path(__file__).with_name('EquilibriumWallEvolution.xml'),
              a.destination, variant=a.variant, beta=a.beta, gamma=a.gamma,
              gamma_scale=a.gamma_scale, theta=a.theta, delta=a.delta, smooth=a.smooth,
              steps=a.steps)
