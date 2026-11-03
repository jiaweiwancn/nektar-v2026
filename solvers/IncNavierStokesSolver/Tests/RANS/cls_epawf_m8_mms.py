"""Genuine flat manufactured solution of the persistent-gamma wall equations.

This source is differentiated from the actual selected stress and scalar
equations.  No residual anchoring or reference reset is used.  Both subdomain
velocities are polynomial, with one common interface velocity and traction;
the dependent epsilon trace and zero k flux are compatible exactly.
"""
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True

from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m4_cases import _evolution, tensor_mesh


def flat_session(template, destination, variant='Regularized', beta=.4,
                 gamma=0., gamma_scale=.25, theta=.2, delta=.1, smooth=False,
                 degree=3, refine=1, steps=2, dt=.005, order=2,
                 maximum_iterations=80, checkpoint=None, restart=None,
                 quadrature=12, normal_order=16):
    from cls_epawf_m8_cases import apply_variant
    root = ET.parse(template).getroot()
    x = [float(i)/refine for i in range(refine+1)]
    z = [a+(b-a)*i/refine for a, b in zip((0., .1, 2.), (.1, 2., 4.))
         for i in range(refine)]+[4.]
    tensor_mesh(root, z, degree, 2, 1., 1., xknots=x)
    _evolution(root, steps, order, dt, 0., False, restart, checkpoint)
    set_info(root, 'RANSWallReference', 'Functions')
    set_info(root, 'RANSWallTransition', 'Fixed')
    set_info(root, 'RANSWallMap', 'Flat')
    set_info(root, 'RANSWallFunctionCoordinates', 'Chart')
    set_info(root, 'TurbulenceModel', 'CLSHighRe')
    for key, value in {'RANSWallHeight': 2., 'RANSWallRoughness': .1,
                       'RANSWallSampleHeight': 1., 'RANSWallTailStart': 1.,
                       'Kinvis': .000015, 'RANSWallFrictionVelocity': .4,
                       'RANSWallDirectionX': 1., 'RANSWallDirectionY': 0.,
                       'RANSWallEta': 0, 'RANSWallChi': 1,
                       'RANSWallTangentialOrder': quadrature,
                       'RANSWallNormalOrder': normal_order,
                       'RANSWallMaxIterations': maximum_iterations,
                       'RANSWallExpressionDimension': 2}.items():
        set_parameter(root, key, repr(value))
    apply_variant(root, variant, beta, gamma, theta, delta, smooth, gamma_scale)
    # Parameters use the independently documented dimensional formulas.  In
    # parallel shear all active cubic terms cancel for any common beta.
    additions = {
        'M8MMSUstar': '.4',
        'M8MMSTop': '4',
        'M8MMSG': 'M8MMSUstar/(Kappa*RANSWallHeight)',
        'M8MMSCmu': '0.3*(1-exp(-0.36*exp(0.75*Qstar)))/(1+0.35*Qstar^1.5)',
        'M8MMSNuT': 'M8MMSCmu*K0^2/E0',
        'M8MMSTraction': 'M8MMSNuT*M8MMSG',
        'M8MMSDH': 'Kinvis+Slope*max(0,RANSWallHeight-Zv)',
        'M8MMSInnerCoefficient': 'M8MMSTraction/(RANSWallHeight*M8MMSDH)',
        'M8MMSUH': 'M8MMSTraction*RANSWallHeight/(2*M8MMSDH)',
        'M8MMSJump': '(11/75)*M8MMSNuT*K0*M8MMSG^2/E0',
        'M8MMSPi': '-M8MMSJump*(M8MMSTop-RANSWallHeight)/M8MMSTop',
        'M8MMSPo': 'M8MMSJump*RANSWallHeight/M8MMSTop',
        'M8MMSOuterForceAtH': 'RANSWallRegularizationScale*Kappa*M8MMSUstar*M8MMSG+RANSWallRegularizationGamma*M8MMSUstar/(Kappa*RANSWallHeight^2)',
    }
    for key, value in additions.items():
        set_parameter(root, key, value)
    r = 'y'
    nt = f'(Slope*max((0*t),({r}-Zv)))'
    nt_derivative = f'(Slope*(1+sign({r}-Zv))/2)'
    D = f'(Kinvis+{nt})'
    below_roughness = f'((1-sign({r}-RANSWallRoughness))/2)'
    gradient_i = f'(M8MMSInnerCoefficient*{r})'
    ui = f'(M8MMSInnerCoefficient*{r}^2/2)'
    uo = f'(M8MMSUH+M8MMSG*({r}-RANSWallHeight))'
    gamma_r = f'(RANSWallRegularizationGamma+RANSWallRegularizationScale*Kappa*M8MMSUstar*{r})'
    # tau_o = nut*g + gamma(r)*(g-us/(kappa*r)); hence -div(tau_o)
    # is the required prescribed mean body force.  Its value at H also
    # determines the actual auxiliary roughness force, which must be
    # cancelled explicitly in the prescribed inner body.
    outer_shear_derivative = f'(RANSWallRegularizationScale*Kappa*M8MMSUstar*M8MMSG+RANSWallRegularizationGamma*M8MMSUstar/(Kappa*{r}^2))'
    body_o = f'(-{outer_shear_derivative})'
    body_i = f'(-M8MMSInnerCoefficient*({D}+{r}*{nt_derivative})-{below_roughness}*M8MMSOuterForceAtH)'
    production_i = f'({nt}*{gradient_i}^2)'
    destruction_i = f'((K0/Qstar)^2*{nt}/{D}^2)'
    production_o = f'((M8MMSNuT*M8MMSG+{gamma_r}*(M8MMSG-M8MMSUstar/(Kappa*{r})))*M8MMSG)'
    ratio = f'(K0^1.5/(2.5*E0*{r}))'
    yap = f'(0.83*E0^2/K0*max(({ratio}-1),(0*t))*{ratio}^2)'
    inner = dict(u=ui, v='0', k='K0', p='M8MMSPi')
    outer = dict(u=uo, v='0', k='K0', eps_tilde='E0', p='M8MMSPo')
    functions = {}
    for side, values in (('Inner', inner), ('Outer', outer)):
        for prefix in ('InitialConditions', 'ExactSolution'):
            functions[prefix+side] = values
        functions['BoundaryValues'+side] = {key: value for key, value in values.items() if key != 'p'}
    functions['BodyForceInner'] = dict(u=body_i, v='0')
    functions['BodyForceOuter'] = dict(u=body_o, v='0')
    functions['ScalarSourceInner'] = dict(k=f'({destruction_i}-{production_i})')
    functions['ScalarSourceOuter'] = dict(k=f'(E0-{production_o})',
        eps_tilde=f'(-1.44*E0/K0*{production_o}+1.92*E0^2/K0-{yap})')
    conditions = root.find('CONDITIONS')
    for node in list(conditions.findall('FUNCTION')):
        conditions.remove(node)
    for name, values in functions.items():
        node = ET.SubElement(conditions, 'FUNCTION', {'NAME': name})
        for variable, value in values.items():
            ET.SubElement(node, 'E', {'VAR': variable, 'VALUE': value})
    return write(root, destination)


def coefficients(raw):
    """Independent dimensional coefficients, evaluated only on Linux."""
    nu, H, rough, _, _, K, q, kappa = raw['parameters']
    us = raw['ustar']
    epsilon = us**3/(kappa*H)
    g = us/(kappa*H)
    cmu = .3*(-math.expm1(-.36*math.exp(.75*q)))/(1+.35*q**1.5)
    nt = cmu*K*K/epsilon
    friction = rough*math.sqrt(K)/(70*nu)
    exponent = max(.5-.4*friction**.7, 1-.79*friction**(-.28)) if friction else 0.
    zv = 10.7*nu*(1-friction**exponent)/math.sqrt(K) if friction else 10.7*nu/math.sqrt(K)
    slope = .2295*math.sqrt(K)
    DH = nu+slope*max(H-zv, 0.)
    traction = nt*g
    a = traction/(H*DH)
    jump = 11/75*nt*K*g*g/epsilon
    return dict(nu=nu, H=H, roughness=rough, K=K, q=q, kappa=kappa, us=us,
                epsilon=epsilon, g=g, nt=nt, zv=zv, slope=slope, inner_coefficient=a,
                UH=traction*H/(2*DH), traction=traction, pressure_jump=jump)


def verify(raw):
    """Independent all-field polynomial and interface accuracy, not a reset."""
    import numpy as np
    from cls_epawf_m7_energy_oracle import Audit
    from validate_cls_epawf_m4 import spaces
    from cls_epawf_m8_native_oracle import selected_stress
    audit = Audit()
    audit.require(raw['dimension'] == 2 and not raw['mapped'] and raw['eta'] == 0 and
                  not raw['neutral_reference'] and raw['chi'] == 1,
                  'Wrong genuine flat nonzero-gamma MMS mode')
    c = coefficients(raw)
    fields, state, offsets = spaces(raw), np.asarray(raw['state']), raw['offsets']
    top = float(fields[0].upper[2])
    pressure_i = -c['pressure_jump']*(top-c['H'])/top
    pressure_o = c['pressure_jump']*c['H']/top
    scaled_errors = []
    for block, field in enumerate(fields):
        square = 0.
        for element in field.elements:
            points, weights = field.quadrature(element, element['degree']+3)
            ids, phi, _ = field.rows(element, points)
            value = phi@state[offsets[block]+ids]
            r = points[:, 2]
            exact = (np.where(r <= c['H'], c['inner_coefficient']*r*r/2,
                              c['UH']+c['g']*(r-c['H'])) if block == 0 else
                     np.full_like(r, 0. if block == 1 else c['K'] if block == 2 else
                                  c['epsilon'] if block == 3 else pressure_i if block == 4 else pressure_o))
            square += float(weights@((value-exact)**2))
        scale = c['us'] if block < 2 else c['K'] if block == 2 else c['epsilon'] if block == 3 else c['us']**2
        error = math.sqrt(max(0., square))/(scale*math.sqrt(field.volume))
        audit.require(error <= 1e-9, f'Genuine flat MMS physical field {block} error {error}')
        scaled_errors.append(error)
    # Direct independently indexed stress at the true interface verifies the
    # analytical pressure offset and shear match for the actual gamma model.
    G = np.zeros((3, 3)); G[0, 2] = c['g']
    outer, _ = selected_stress(raw, G, c['K'], c['epsilon'], [0.5, 0., c['H']])
    audit.near(-outer[0, 2], c['traction'], 'm8_mms_exact_shared_tangential_traction', 2e-12)
    audit.near(2*c['K']/3-outer[2, 2], c['pressure_jump'], 'm8_mms_exact_normal_pressure_offset', 2e-12)
    face = raw['state_audits'][-1]['interface']
    audit.near(face[1]/face[0], -c['pressure_jump'], 'm8_mms_native_normal_pressure_offset', 1e-9)
    audit.near(face[2:6], np.zeros(4), 'm8_mms_native_total_traction_and_zero_k_flux', 1e-9)
    source = raw['state_audits'][-1]['robust_source']
    audit.near(source[4:7], [0., 0., 1.], 'm8_mms_inactive_repair_plateau', 1e-10)
    audit.require(max(abs(value) for value in raw['step_audits'][-1][4:]) <= 1e-9,
                  'Genuine MMS did not remain stationary without anchoring')
    return dict(schema=1, success=True, checks=audit.checks, errors=audit.errors,
                coefficients=c, field_errors=scaled_errors,
                persistent_gamma=raw['robust_parameters'][1] != 0 or raw['robust_parameters'][5] != 0,
                scope='Actual nonzero-gamma flat operator; polynomial complete state, independently derived forcing, shared traction and dependent epsilon trace',
                curved_mms='not established by this flat construction')


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('destination', type=Path)
    p.add_argument('--gamma', type=float, default=0.)
    p.add_argument('--gamma-scale', type=float, default=.25)
    p.add_argument('--steps', type=int, default=2)
    a = p.parse_args()
    flat_session(Path(__file__).with_name('EquilibriumWallEvolution.xml'), a.destination,
                 gamma=a.gamma, gamma_scale=a.gamma_scale, steps=a.steps)
