"""Independent BDF2 flat AWF momentum inventory, using accepted native traces.

Pointwise flux residuals need not vanish in CG. Check that their independently
reconstructed CV balance equals the reported residual, and retain its magnitude.
"""
import csv
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def flat_momentum_budget(case):
    case=Path(case);root=ET.parse(case).getroot()
    info={n.attrib['PROPERTY']:n.attrib['VALUE'] for n in root.findall('CONDITIONS/SOLVERINFO/I')}
    params={n.text.split('=',1)[0].strip():float(n.text.split('=',1)[1]) for n in root.findall('CONDITIONS/PARAMETERS/P')}
    if (info.get('RANSWallGeometry')!='Plane' or params.get('RANSTerrainSlope',0.)!=0
            or info.get('RANSWallModel','None')=='None'):
        return dict(supported=False,reason='Only horizontal plane AWF inventories are assessed')
    if int(root.find('CONDITIONS/TIMEINTEGRATIONSCHEME/ORDER').text)!=2:
        return dict(supported=False,reason='The independent time stencil is BDF2')
    forcing=root.find('FORCING')
    if forcing is not None and len(forcing):raise ValueError('This budget oracle requires no external forcing')
    dt=params['TimeStep'];groups={}
    for path in case.parent.glob(case.stem+'.terrain-points-r*.csv'):
        with path.open(newline='') as stream:
            for r in csv.DictReader(stream):
                t=float(r['time'])
                if r['phase']=='1' and abs(t/dt-round(t/dt))<1e-7:
                    groups.setdefault(int(r['point']),{})[round(t/dt)]=r
    maximum=physical=0.;count=0
    for timeline in groups.values():
        for n,r in timeline.items():
            if n-1 not in timeline or n-2 not in timeline:continue
            t=float(r['time'])
            if t<params.get('RANSStartupDuration',0.):continue
            a,b=timeline[n-1],timeline[n-2]
            derivative=(3*float(r['uP'])-4*float(a['uP'])+float(b['uP']))/(2*dt)
            convection=(float(r['convection_t']) if info.get('RANSAdvectionTreatment','Explicit')=='Implicit' else
                        2*float(a['convection_t'])-float(b['convection_t']))
            balance=float(r['H'])*(derivative+convection+float(r['pressure_t']))+float(r['traction_x'])-float(r['resolved_shear'])
            reported=float(r['tangent_flux_residual'])
            if not math.isfinite(balance+reported):raise ValueError('Nonfinite strip inventory')
            maximum=max(maximum,abs(balance-reported));physical=max(physical,abs(balance));count+=1
    if not count:raise ValueError('At least three post-startup BDF2 trace times are required')
    if maximum>5e-8:raise ValueError('Independent BDF2 strip momentum mismatch: '+str(maximum))
    return dict(supported=True,samples=count,maximum_inventory_identity_error=maximum,
        maximum_pointwise_flux_residual=physical,
        partition=info.get('RANSAdvectionTreatment','Explicit'),
        definition='H*d_t(U_P)+H*(convection_P+mechanical_pressure_gradient_P)+tau_wall-T_H; no external forcing',
        pointwise_zero_required=False)
