"""M7 bounded perturbations, preserving the frozen CLS/AWF laws and traces."""
import sys
sys.dont_write_bytecode = True
import math
import xml.etree.ElementTree as ET
from cls_epawf_m3_cases import set_info, set_parameter, write
from cls_epawf_m4_cases import neutral_session, tensor_mesh
from cls_epawf_m6_cases import transition_session, arc_offset


def perturbation_session(template, destination, kind='joint', amplitude=.001,
                         wave=1, degree=3, refine=1, dt=.0025, steps=4,
                         dimension=2, eta=1, chi=1, height=2., transition=False,
                         ranks_shift=0, order=2, checkpoint=None, restart=None):
    if transition:
        transition_session(template, destination, family='feet', degree=degree,
                           steps=steps, dt=dt, order=order, refine=refine,
                           quadrature=8, normal_order=16,checkpoint=checkpoint,restart=restart)
        root = ET.parse(destination).getroot()
        shift=arc_offset();L=6.-shift
        x = f'(x+{4.-shift})/{L}'; z = 'y'; Z = 4.
    else:
        neutral_session(template, destination, degree=degree, dimension=dimension,
                        inner_knots=[0.,height] if dimension==3 else [0., .1, 1., height],
                        outer_knots=[height,8.] if dimension==3 else [height, 4., 8.],
                        length=8., span=2., steps=steps, order=order, dt=dt,
                        eta=eta, chi=chi, tangential_order=8, normal_order=16,
                        checkpoint=checkpoint,restart=restart)
        root = ET.parse(destination).getroot()
        if refine > 1:
            knots = [0., .1, 1., height, 4., 8.]
            knots = sorted(set(a+(b-a)*j/refine for a,b in zip(knots,knots[1:]) for j in range(refine)) | {8.})
            tensor_mesh(root, knots, degree, dimension, 8., 2.,
                        xknots=[8.*i/refine for i in range(refine+1)])
        x = 'x/8'; z = 'y' if dimension==2 else 'z'; L = 8.; Z = 8.
    set_parameter(root, 'RANSWallDiagnostics', 1)
    set_parameter(root, 'RANSWallOwnershipShift', ranks_shift)
    set_parameter(root, 'RANSWallMaxIterations', 35)
    # All factors vanish on the ordinary essential boundaries. The stream
    # function pair is divergence-free before spatial projection.
    bx = f'(sin({wave}*PI*({x})))^2'
    bz = f'(sin(PI*{z}/{Z}))^2'
    osc = f'sin({2*wave}*PI*({x}))'
    scalar = f'({bx})*({bz})'
    if kind == 'interface':
        scalar += f'*exp(-(({z}-{height})/0.5)^2)'
    variables = ['u','v'] if dimension==2 else ['u','v','w']
    conditions = root.find('CONDITIONS')
    for region in ('Inner','Outer'):
        node = ET.SubElement(conditions, 'FUNCTION', {'NAME':'InitialPerturbation'+region})
        values = {v:'0' for v in variables+['k','eps_tilde','p']}
        if kind in ('velocity','joint','interface','transition'):
            values['u'] = f'{amplitude*.4}*({bx})*sin(2*PI*{z}/{Z})'
            values[variables[-1]] = f'{-amplitude*.4*wave*Z/L}*({osc})*({bz})'
        if kind in ('k','joint','interface','transition'):
            values['k'] = f'{amplitude}*Qstar*0.16*({scalar})'
        if region == 'Outer' and kind in ('epsilon','joint','interface','transition'):
            values['eps_tilde'] = f'{amplitude*.064}/(Kappa*{height})*({bx})*(sin(PI*({z}-{height})/({Z}-{height})))^2'
        if kind == 'node_alias':
            values['k'] = f'-1.001*Qstar*0.16*4*({x})*(1-({x}))*(1-({z}/{Z})^2)'
        if kind == 'invalid_sampler':
            values['k'] = f'-4*Qstar*0.16*4*({x})*(1-({x}))*(1-({z}/{Z}))^2'
        if dimension == 3:
            for var in values:
                if values[var] != '0':values[var] += '*(sin(PI*y/2))^2'
        for var,value in values.items():ET.SubElement(node,'E',{'VAR':var,'VALUE':value})
    return write(root,destination)
