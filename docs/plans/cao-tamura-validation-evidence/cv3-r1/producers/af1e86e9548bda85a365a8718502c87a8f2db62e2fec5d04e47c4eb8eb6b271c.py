#!/usr/bin/env python3
"""Cao--Tamura experiment preparation and observation (not a calibrated case).

Only the frozen flat-floor data can enter inlet fitting or calibration. All
native sampling is delegated to RANSObservationAudit, never nearest-node data.
"""
import argparse
import csv
import hashlib
import json
import math
from functools import lru_cache
import os
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
from generate_terrain import plane, terrain_mesh
from generate_m3 import element
from generate_m4 import info, parameter
from generate_m6 import original_j0
from validate_terrain import bounded_solver
from cao_tamura_observables import speedup, reattachment
from cao_tamura_routes import envelope, resolve_route, route_metadata

ROOT = Path(__file__).resolve().parents[4]
DATA = ROOT/'docs/plans/cao-tamura-validation-data'
VERSION = 'cao-tamura-workflow-v3'
THREADS = {k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
MODELS = ('standard', 'standard_original', 'standard_fixed', 'standard_switched',
          'cls_original', 'cls_fixed', 'cls_switched')
DEFAULTS = dict(surface='S', geometry='floor', model='standard_fixed',
                diffusion='Molecular', dimension=2, degree=4, refinement=1,
                layers=10, span=.16, span_layers=2, hm=.011, point_p=.0055,
                point_n=.022, ks=0., xmin=-1., xmax=1.2, ramp_end=-.2,
                ramp_width_h=2., dt=1e-5, steps=2, checkpoint=1,
                outer_length=.025, outer_intensity=.003, failure_policy='UseOriginal',
                inlet_delta=.25, k_scale=1., epsilon_scale=1., x_columns=None,
                velocity_scaling='inlet_fit', z0=None,
                inlet_reference_q=None,inlet_reference_kappa=None,outer_blend='Geometric',mean_outer_cap=False,pressure_degree=None,time_order=2,epsilon_profile='NeutralLog',startup_flow_times=0.,
                inlet_wall_collar=False,k_profile='NeutralConstant',variance_factor=1.,outer_blend_cells=10,
                advection='Explicit',quadrature_points=None,wall_route=None)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def frozen_data():
    snapshot = json.loads((DATA/'snapshot.json').read_text())
    for name, expected in snapshot['files'].items():
        if sha(DATA/name) != expected:
            raise ValueError('Frozen CV1 data changed: '+name)
    return snapshot


def records(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, columns=None):
    if not rows and columns is None:
        raise ValueError('Empty CSV requires columns')
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns or list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def solve_linear(a, b):
    a = [list(row)+[v] for row, v in zip(a, b)]
    for i in range(len(a)):
        pivot = max(range(i, len(a)), key=lambda j: abs(a[j][i]))
        a[i], a[pivot] = a[pivot], a[i]
        if abs(a[i][i]) < 1e-15:
            raise ValueError('Singular inlet fit')
        a[i] = [v/a[i][i] for v in a[i]]
        for j in range(len(a)):
            if j != i:
                v = a[j][i]; a[j] = [x-v*y for x, y in zip(a[j], a[i])]
    return [row[-1] for row in a]


def profile_fit(surface, hm):
    """Small weighted log-polynomial initialized using Fig. 3a alone."""
    points = [p for p in records(DATA/'profiles.csv')
              if p['series_id'] == 'f3a_'+surface and float(p['weight']) > 0
              and .04*float(p['y_over_H']) <= .25]
    x = [math.log(.04*float(p['y_over_H'])/hm) for p in points]
    y = [float(p['value']) for p in points]
    w = [float(p['weight']) for p in points]
    coeff = solve_linear([[sum(t*z**(i+j) for z, t in zip(x, w)) for j in range(4)]
                         for i in range(4)],
                        [sum(t*v*z**i for z, v, t in zip(x, y, w)) for i in range(4)])
    return coeff, points


def blend(r):
    r = max(0., min(1., r))
    return 1-10*r**3+15*r**4-6*r**5


@lru_cache(maxsize=32)
def rms_fit(surface,hm,point_n):
    """Allowed Fig.3c only; no reserved vertical variance enters inferred k."""
    points=[p for p in records(DATA/'profiles.csv') if p['series_id']=='f3c_'+surface
            and float(p['weight'])>0 and float(p['value'])>0 and .04*float(p['y_over_H'])>=point_n]
    x=[.04*float(p['y_over_H'])/.25 for p in points]
    y=[math.log(float(p['value'])) for p in points];w=[float(p['weight']) for p in points]
    coeff=solve_linear([[sum(t*z**(i+j) for z,t in zip(x,w)) for j in range(4)] for i in range(4)],
                       [sum(t*v*z**i for z,v,t in zip(x,y,w)) for i in range(4)])
    return coeff,points


def inlet_reference(c):
    standard=c['model'].startswith('standard')
    return (c['inlet_reference_q'] or (10/3 if standard else 3.190675081624288),
            c['inlet_reference_kappa'] or (math.sqrt(.48*1.3/(10/3)) if standard else .440553905186659))


def near_k_over_star2(c):
    if c['k_profile']=='NeutralConstant':return c['k_scale']*inlet_reference(c)[0]
    coeff,_=rms_fit(c['surface'],c['hm'],c['point_n'])
    l=c['point_n']/c['inlet_delta']
    return c['k_scale']*c['variance_factor']*math.exp(2*sum(a*l**i for i,a in enumerate(coeff)))


def collar_data(c,coeff):
    """Exact primitive of a convex blend of canonical and measured shear."""
    q,kap=inlet_reference(c);ratio=c['point_n']/c['hm']
    a=(.033 if c['surface']=='S' else .051)*math.sqrt(near_k_over_star2(c)/q)/kap
    b=[1.]+[0.]*5
    for power,value in ((3,-10.),(4,15.),(5,-6.)):
        for j in range(power+1):b[j]+=value*math.comb(power,j)*(-1.)**(power-j)/(ratio-1)**power
    l=math.log(.25/c['inlet_delta'])
    d=[a-coeff[1]-2*coeff[2]*l-3*coeff[3]*l*l,-2*coeff[2]-6*coeff[3]*l,-3*coeff[3]]
    return a,b,d


def collar_primitive(x,b,d):
    l=math.log(x)
    value=b[0]*(d[0]*l+d[1]*l*l/2+d[2]*l**3/3)
    for m in range(1,6):
        value+=b[m]*x**m/m*(d[0]+d[1]*(l-1/m)+d[2]*(l*l-2*l/m+2/m**2))
    return value


def mean_raw(eta, c, coeff):
    delta=c['inlet_delta'];l = math.log(eta*.25/(delta*c['hm']))
    w = blend((eta-.8*delta)/(.4*delta))
    value=sum(a*l**i for i,a in enumerate(coeff))
    if c['mean_outer_cap']:
        t=max(0.,min(1.,(value-.98)/.04))
        value=.98+.04*(t-t**3+.5*t**4)+min(value-.98,0.)
    return w*value+(1-w)


def mean_raw_derivative(eta,c,coeff):
    delta=c['inlet_delta'];l=math.log(eta*.25/(delta*c['hm']))
    value=sum(a*l**i for i,a in enumerate(coeff))
    derivative=sum(i*a*l**(i-1)/eta for i,a in enumerate(coeff) if i)
    if c['mean_outer_cap']:
        t=max(0.,min(1.,(value-.98)/.04))
        derivative*=1-3*t*t+2*t**3
        value=.98+.04*(t-t**3+.5*t**4)+min(value-.98,0.)
    t=max(0.,min(1.,(eta-.8*delta)/(.4*delta)))
    return blend(t)*derivative-30*t*t*(1-t)**2/(.4*delta)*(value-1)


def mean_ratio(eta,c,coeff):
    value=mean_raw(eta,c,coeff)
    if c['inlet_wall_collar'] and eta<c['point_n']:
        _,b,d=collar_data(c,coeff)
        value+=collar_primitive(max(1.,eta/c['hm']),b,d)-collar_primitive(c['point_n']/c['hm'],b,d)
    return value


def mean_derivative(eta,c,coeff):
    derivative=mean_raw_derivative(eta,c,coeff)
    if c['inlet_wall_collar'] and eta<c['point_n']:
        a,_,_=collar_data(c,coeff);w=blend((eta-c['hm'])/(c['point_n']-c['hm']))
        derivative=w*a/eta+(1-w)*derivative
    return derivative


def config(value):
    unknown = set(value)-set(DEFAULTS)-{'assumed_nu'}
    if unknown:
        raise ValueError('Unknown configuration keys: '+str(sorted(unknown)))
    c = {**DEFAULTS, **value}
    if c['surface'] not in ('S', 'R') or c['geometry'] not in ('floor', 'hill') or c['model'] not in MODELS:
        raise ValueError('Unsupported surface, geometry or model/wall combination')
    c['wall_route']=resolve_route(c['model'],c['wall_route'])
    if c['diffusion'] not in ('Molecular', 'HighRe') or c['dimension'] not in (2, 3):
        raise ValueError('Unsupported diffusion or dimension')
    if c['failure_policy'] not in ('UseOriginal', 'RequirePreservation'):
        raise ValueError('Unsupported AWF failure policy')
    if not 0<=c['startup_flow_times']<=10:
        raise ValueError('Use 0..10 flow-through times for startup continuation')
    if c['startup_flow_times'] and c['geometry']!='floor':
        raise ValueError('Startup continuation is for empty floors; set startup_flow_times=0 for a hill')
    if c['k_profile'] not in ('NeutralConstant','MeasuredRms') or not 0<c['variance_factor']<=3:
        raise ValueError('Unsupported inferred variance profile or factor')
    if type(c['inlet_wall_collar']) is not bool:
        raise ValueError('inlet_wall_collar must be a boolean')
    if c['inlet_wall_collar'] and not (1.5<=c['point_n']/c['hm']<=3 and c['point_n']<.8*c['inlet_delta']):
        raise ValueError('Wall-compatible collar requires 1.5<=N/H<=3 below the outer blend')
    if c['epsilon_profile'] not in ('NeutralLog','LocalShear','ShearLength'):
        raise ValueError('Unsupported inferred inlet dissipation profile')
    if c['advection'] not in ('Explicit','Implicit'):
        raise ValueError('Unsupported RANS advection partition')
    if type(c['mean_outer_cap']) is not bool:
        raise ValueError('mean_outer_cap must be a boolean')
    if c['outer_blend'] not in ('Geometric','Arithmetic'):
        raise ValueError('Unsupported outer turbulence interpolation')
    if c['velocity_scaling'] not in ('inlet_fit','incident_ReH'):
        raise ValueError('Unsupported velocity dimensionalization')
    for name in ('degree', 'refinement', 'layers', 'span_layers', 'steps', 'checkpoint','outer_blend_cells'):
        if type(c[name]) is not int or c[name] < 1:
            raise ValueError(name+' must be a positive integer')
    if type(c['time_order']) is not int or c['time_order'] not in (1,2):
        raise ValueError('time_order must be 1 or 2')
    if c['pressure_degree'] is not None and (type(c['pressure_degree']) is not int or not 1<=c['pressure_degree']<c['degree']):
        raise ValueError('Require 1<=pressure_degree<velocity degree')
    if c['quadrature_points'] is not None and (type(c['quadrature_points']) is not int or c['quadrature_points']<c['degree']+2):
        raise ValueError('Require at least degree+2 quadrature points')
    if c['degree'] < 3 or c['layers'] < 5:
        raise ValueError('Use degree >=3 and at least five outer layers')
    for name, v in c.items():
        if isinstance(v, (int, float)) and not math.isfinite(v):
            raise ValueError('Nonfinite '+name)
    H, P, N = (c[k] for k in ('hm', 'point_p', 'point_n'))
    if not 0 < P < H < N < .1 or not 0 <= c['ks'] < H:
        raise ValueError('Require 0<P<hm<N<0.1 and 0<=ks<hm')
    if .04*math.pi**2/(2*.1**2)*H > .25:
        raise ValueError('hm violates the hill curvature guard, including paired floors')
    if c['surface'] == 'S' and c['ks'] != 0:
        raise ValueError('Smooth AWF uses ks=0')
    if not c['xmin'] < -.52 or not c['xmax'] >= .8:
        raise ValueError('Provide approach fetch upstream of widest ramp and a wake outlet')
    if not c['xmin'] < c['ramp_end']-.04*c['ramp_width_h'] < c['ramp_end'] < -.1:
        raise ValueError('The complete transition must lie upstream of the hill foot')
    if min(c['dt'], c['outer_length'], c['outer_intensity'], c['span']) <= 0:
        raise ValueError('Time, length and intensity parameters must be positive')
    if not .08<=c['inlet_delta']<=.45 or not 0<c['k_scale']<=5 or not 0<c['epsilon_scale']<=10:
        raise ValueError('Inlet thickness/variance/length parameters are outside declared admissible bounds')
    if c['x_columns'] is not None:
        x=c['x_columns']
        if c['geometry']!='floor' or not isinstance(x,list) or len(x)<3 or not all(isinstance(a,(int,float)) and math.isfinite(a) for a in x) or x!=sorted(set(x)) or x[0]!=c['xmin'] or x[-1]!=c['xmax'] or 0. not in x:
            raise ValueError('Explicit floor x_columns must increase from xmin to xmax and contain x=0')
    if c['z0'] is not None and not 0<c['z0']<c['hm']:
        raise ValueError('Require 0<z0<hm')
    for key in ('inlet_reference_q','inlet_reference_kappa'):
        if c[key] is not None and (not isinstance(c[key],(int,float)) or not math.isfinite(c[key]) or c[key]<=0):
            raise ValueError(key+' must be positive')
    if 'assumed_nu' in c and c['assumed_nu'] <= 0:
        raise ValueError('assumed_nu must be positive')
    return c


def reference(c):
    standard = c['model'].startswith('standard')
    q = 10/3 if standard else 3.190675081624288
    kap = math.sqrt(.48*1.3/q) if standard else .440553905186659
    coeff, points = profile_fit(c['surface'], c['hm'])
    if c['inlet_wall_collar'] and c['mean_outer_cap'] and any(
            mean_raw(c['hm']+(c['point_n']-c['hm'])*i/32,
                     {**c,'mean_outer_cap':False},coeff)>=.98 for i in range(33)):
        raise ValueError('The exact shear collar requires an inactive mean cap below N')
    re_h = 12000 if c['surface'] == 'S' else 10500
    uh = re_h*c['assumed_nu']/.04 if 'assumed_nu' in c else 1.
    if c['velocity_scaling']=='incident_ReH':
        pts=sorted((.04*float(p['y_over_H']),float(p['value'])) for p in points)
        a,b=next((a,b) for a,b in zip(pts,pts[1:]) if a[0]<=.04<=b[0])
        ratio=a[1]+(.04-a[0])*(b[1]-a[1])/(b[0]-a[0])
    else:ratio=mean_ratio(.04,c,coeff)
    uinf = uh/ratio
    star = uinf*(.033 if c['surface'] == 'S' else .051)
    k_outer = 1.5*(c['outer_intensity']*uinf)**2
    inlet_q=c['inlet_reference_q'] if c['inlet_reference_q'] is not None else q
    inlet_kappa=c['inlet_reference_kappa'] if c['inlet_reference_kappa'] is not None else kap
    result = dict(q=q, kappa=kap,inlet_q=inlet_q,inlet_kappa=inlet_kappa, coefficients=coeff, U_H=uh, Uinf=uinf,
                  ustar=star, nu=.04*uh/re_h, re_h=re_h, k_near=near_k_over_star2(c)*star**2,
                  k_outer=k_outer, epsilon_outer=.09**.75*k_outer**1.5/c['outer_length'],
                  displacement=0. if c['surface'] == 'S' else .0035,
                  z0=c['z0'] if c['z0'] is not None else 4e-6 if c['surface'] == 'S' else .0002,
                  fit_series='f3a_'+c['surface'], fit_points=len(points))
    if min(mean_ratio(c['hm']*(.8/c['hm'])**(i/200), c, coeff) for i in range(201)) <= 0:
        raise ValueError('Fitted inlet mean is nonpositive')
    result['rms_coefficients'],rms_points=(rms_fit(c['surface'],c['hm'],c['point_n']) if c['k_profile']=='MeasuredRms' else ([],[]))
    result['rms_fit_points']=len(rms_points)
    result['epsilon_inlet_near']=epsilon_seed(c['hm'],c,result)
    result['epsilon_unit_near']=epsilon_seed(c['hm'],c,result,1.)
    H = c['hm']; nu = result['nu']; K = result['k_near']
    j0 = original_j0(K, H, c['ks'], nu)
    beta = .09*2.55*math.sqrt(K)
    rh = c['ks']*math.sqrt(K)/(70*nu)
    intercept = 10.7 if rh == 0 else 10.7*(1-rh**max(.5-.4*rh**.7, 1-.79*rh**(-.28)))
    zv = nu*intercept/math.sqrt(K)
    dh = nu+beta*max(0., H-zv)
    if H <= zv:
        l = 0.
    elif zv > 0:
        l = (math.log(dh/nu)-1+nu/dh)/beta
    else:
        d0 = nu-beta*zv
        l = (math.log(dh/d0)-nu/d0+nu/dh)/beta
    anchor_star=math.sqrt(K/q);A = anchor_star**3/kap
    lam = anchor_star**4*l/A
    result.update(original_anchor=anchor_star**2*j0, ep_lambda=lam,
                  ep_cutoff=H*math.exp(1-lam), wall_re=math.sqrt(K)*H/nu,
                  wall_epsilon_near=K**1.5/(2.55*H))
    if c['model'] == 'standard' or c['model'].endswith(('_fixed', '_switched')):
        result['wall_epsilon_near'] = A/H
    return result


def inlet_k(eta,c,r):
    delta=c['inlet_delta'];w=blend((eta-.8*delta)/(.4*delta))
    seed=r['k_near']
    if c['k_profile']=='MeasuredRms':
        l=max(c['point_n'],min(1.2*delta,eta))/delta
        raw=c['k_scale']*c['variance_factor']*r['ustar']**2*math.exp(2*sum(a*l**i for i,a in enumerate(r['rms_coefficients'])))
        b=blend((eta-c['point_n'])/c['point_n'])
        seed=math.exp(b*math.log(seed)+(1-b)*math.log(raw))
    return (w*seed+(1-w)*r['k_outer'] if c['outer_blend']=='Arithmetic' else
            math.exp(w*math.log(seed)+(1-w)*math.log(r['k_outer'])))


def epsilon_seed(eta,c,r,scale=None):
    scale=c['epsilon_scale'] if scale is None else scale
    if c['epsilon_profile']=='NeutralLog':return scale*r['ustar']**3/(r['inlet_kappa']*eta)
    k=inlet_k(eta,c,r)
    outer=(.09**.75*k**1.5/c['outer_length'] if c['epsilon_profile']=='ShearLength' else r['epsilon_outer'])
    return math.hypot(scale*k*r['Uinf']*mean_derivative(eta,c,r['coefficients'])/r['inlet_q'],outer)


def wall_epsilon(c,r,x=None):
    """Initial wall trace of the actual spatially blended AWF policy."""
    if resolve_route(c['model'],c.get('wall_route'))!='EPTransition':
        return r['wall_epsilon_near']
    x=c['xmin'] if x is None else x
    w=envelope(x,c['ramp_end']-.04*c['ramp_width_h'],c['ramp_end'])
    original=r['k_near']**1.5/(2.55*c['hm'])
    return original+w*(r['wall_epsilon_near']-original)


def wall_epsilon_expression(c):
    if resolve_route(c['model'],c.get('wall_route'))!='EPTransition':return 'EpsilonWall'
    a=c['ramp_end']-.04*c['ramp_width_h'];width=.04*c['ramp_width_h']
    t=f'((x-({a:.17g}))/({width:.17g}))'
    positive=f'(({t}+abs({t}))/2)'
    clip=f'(({positive}+1-abs({positive}-1))/2)'
    w=f'(1-10*{clip}^3+15*{clip}^4-6*{clip}^5)'
    return f'(EpsilonOriginalH+{w}*(EpsilonWall-EpsilonOriginalH))'


def inlet_values(eta,c,r,x=None):
    delta=c['inlet_delta'];w=blend((eta-.8*delta)/(.4*delta))
    collar=blend((eta-c['hm'])/(c['point_n']-c['hm']))
    k=inlet_k(eta,c,r);epsilon=epsilon_seed(eta,c,r)
    if c['epsilon_profile']=='ShearLength':
        inner=epsilon_seed(eta,c,r,1.)*wall_epsilon(c,r,x)/r['epsilon_unit_near']
        epsilon=math.exp(collar*math.log(inner)+(1-collar)*math.log(epsilon))
    else:
        epsilon+=collar*(wall_epsilon(c,r,x)-r['epsilon_inlet_near'])
        if epsilon<=0:raise ValueError('Inferred epsilon collar is nonpositive')
        if c['epsilon_profile']=='NeutralLog':
            epsilon=(w*epsilon+(1-w)*r['epsilon_outer'] if c['outer_blend']=='Arithmetic' else
                     math.exp(w*math.log(epsilon)+(1-w)*math.log(r['epsilon_outer'])))
    if epsilon<=0:raise ValueError('Inferred epsilon profile is nonpositive')
    return dict(u=r['Uinf']*mean_ratio(eta,c,r['coefficients']),v=0.,w=0.,k=k,epsilon=epsilon,p=2*k/3)


def expressions(c, r):
    z = 'y' if c['dimension'] == 2 else 'z'
    def smooth(a, b):
        t = f'(({z}-({a:.17g}))/({b-a:.17g}))'
        positive = f'(({t}+abs({t}))/2)'
        clip = f'(({positive}+1-abs({positive}-1))/2)'
        return f'(1-10*{clip}^3+15*{clip}^4-6*{clip}^5)'
    delta=c['inlet_delta'];w = smooth(.8*delta,1.2*delta)
    l = f'log({z}*0.25/({delta:.17g}*H))'
    poly = '('+'+'.join(f'({a:.17g})*({l})^{i}' for i, a in enumerate(r['coefficients']))+')'
    derivative='('+'+'.join(f'({i*a:.17g})*({l})^{i-1}/{z}' for i,a in enumerate(r['coefficients']) if i)+')'
    if c['mean_outer_cap']:
        t=f'(({poly}-.98)/.04)'
        positive=f'(({t}+abs({t}))/2)'
        t=f'(({positive}+1-abs({positive}-1))/2)'
        below=f'(({poly}-.98-abs({poly}-.98))/2)'
        derivative=f'({derivative}*(1-3*{t}^2+2*{t}^3))'
        poly=f'(.98+.04*({t}-{t}^3+.5*{t}^4)+{below})'
    def clip(value,lower,upper):
        positive=f'((({value})-({lower:.17g})+abs(({value})-({lower:.17g})))/2+({lower:.17g}))'
        return f'(({positive}+({upper:.17g})-abs({positive}-({upper:.17g})))/2)'
    t=clip(f'({z}-({.8*delta:.17g}))/({.4*delta:.17g})',0.,1.)
    gradient=f'({w}*{derivative}-30*{t}^2*(1-{t})^2/({.4*delta:.17g})*({poly}-1))'
    mean=f'({w}*{poly}+(1-{w}))'
    collar=smooth(c['hm'],c['point_n'])
    if c['inlet_wall_collar']:
        aa,bb,dd=collar_data(c,r['coefficients']);ratio=c['point_n']/c['hm']
        x=clip(f'{z}/H',1.,ratio);lx=f'log({x})'
        terms=[f'({bb[0]:.17g})*(({dd[0]:.17g})*{lx}+({dd[1]/2:.17g})*{lx}^2+({dd[2]/3:.17g})*{lx}^3)']
        for m in range(1,6):
            terms.append(f'({bb[m]/m:.17g})*{x}^{m}*(({dd[0]:.17g})+({dd[1]:.17g})*({lx}-({1/m:.17g}))+({dd[2]:.17g})*({lx}^2-({2/m:.17g})*{lx}+({2/m**2:.17g})))')
        mean=f'({mean}+('+ '+'.join(terms)+f')-({collar_primitive(ratio,bb,dd):.17g}))'
        gradient=f'({collar}*({aa:.17g})/{z}+(1-{collar})*{gradient})'
    seed='Knear'
    if c['k_profile']=='MeasuredRms':
        eta=clip(z,c['point_n'],1.2*delta)
        lr=f'({eta}/({delta:.17g}))'
        pk='('+'+'.join(f'({v:.17g})*({lr})^{i}' for i,v in enumerate(r['rms_coefficients']))+')'
        raw=f'(({c["k_scale"]*c["variance_factor"]*r["ustar"]**2:.17g})*exp(2*{pk}))'
        wk=smooth(c['point_n'],2*c['point_n'])
        seed=f'exp({wk}*log(Knear)+(1-{wk})*log({raw}))'
    k=(f'exp({w}*log({seed})+(1-{w})*log(Kouter))' if c['outer_blend']=='Geometric' else
       f'({w}*{seed}+(1-{w})*Kouter)')
    if c['epsilon_profile'] in ('LocalShear','ShearLength'):
        outer=f'((0.09)^0.75*{k}^1.5/OuterLength)' if c['epsilon_profile']=='ShearLength' else 'EpsilonOuter'
        def shear_seed(scale):
            return f'sqrt((({scale:.17g})*{k}*Uinf*{gradient}/InletQ)^2+{outer}^2)'
        base=shear_seed(c['epsilon_scale'])
        if c['epsilon_profile']=='ShearLength':
            inner=f'({shear_seed(1.)}*EpsilonWall/EpsilonUnitH)'
            eps=f'exp({collar}*log({inner})+(1-{collar})*log({base}))'
        else:eps=f'({base}+{collar}*(EpsilonWall-EpsilonInletH))'
    else:
        eps=f'(({c["epsilon_scale"]:.17g})*Ustar^3/(InletKappa*{z})+{collar}*(EpsilonWall-EpsilonInletH))'
        eps=(f'exp({w}*log({eps})+(1-{w})*log(EpsilonOuter))' if c['outer_blend']=='Geometric' else
             f'({w}*{eps}+(1-{w})*EpsilonOuter)')
    eps=eps.replace('EpsilonWall',wall_epsilon_expression(c))
    return dict(u=f'Uinf*{mean}',v='0',w='0',k=k,eps=eps,p=f'2/3*{k}')


def generate(value, output):
    frozen = frozen_data(); c = config(value); r = reference(c)
    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=True)
    path = output/'case.xml'
    if path.exists():
        raise ValueError('Refusing to overwrite an existing case')
    plane(path, model=c['model'], degree=c['degree'], layers=4, nx=2,
          dt=c['dt'], steps=c['steps'], dimension=c['dimension'], height=c['hm'], diffusion=c['diffusion'])
    tree = ET.parse(path); root = tree.getroot(); conditions = root.find('CONDITIONS')
    root.find('CONDITIONS/TIMEINTEGRATIONSCHEME/ORDER').text=str(c['time_order'])
    if c['pressure_degree'] is not None:
        root.find("EXPANSIONS/E[@FIELDS='p']").set('NUMMODES',','.join([str(c['pressure_degree']+1)]*c['dimension']))
    if c['quadrature_points'] is not None:
        for e in root.findall('EXPANSIONS/E'):
            e.set('NUMPOINTS',','.join([str(c['quadrature_points'])]*c['dimension']))
    # Every width uses the same complete mesh and feet. Refine the whole domain.
    xs = sorted(set([c['xmin'], c['xmax']]+[x for x in
        (-.8, -.68, -.56, -.52, -.48, -.44, -.40, -.36, -.32, -.28, -.24, -.20,
         -.16, -.12, -.1, -.075, -.05, -.025, 0., .025, .05, .075, .1, .14,
         .2, .3, .45, .6, .8, 1.) if c['xmin'] < x < c['xmax']]))
    def refine(values):
        n = c['refinement']
        return [a+(b-a)*i/n for a,b in zip(values, values[1:]) for i in range(n)]+[values[-1]]
    xs = refine(c['x_columns'] if c['x_columns'] is not None else xs)
    begin,end=.8*c['inlet_delta'],1.2*c['inlet_delta']
    etas = sorted(set([c['hm'], c['point_n'], begin, end, .8]+[begin+i*(end-begin)/c['outer_blend_cells'] for i in range(1,c['outer_blend_cells'])]+
                     [c['point_n']*(.8/c['point_n'])**(i/c['layers']) for i in range(c['layers']+1)]))
    terrain_mesh(root, xs, refine(etas), .04 if c['geometry'] == 'hill' else 0.,
                 top=.8, dimension=c['dimension'], span=c['span'],
                 normal_n=c['point_n'], span_layers=c['span_layers'])
    info(root, 'RANSWallGeometry', 'CosineHill' if c['geometry'] == 'hill' else 'Plane')
    info(root, 'RANSAWFFailurePolicy', c['failure_policy'])
    info(root, 'RANSAWFSpatialWeight', 'Ramp' if c['wall_route']=='EPTransition' else 'None')
    info(root, 'RANSExperimentalRoute',c['wall_route'])
    info(root, 'AdvectionType', 'Convective'); info(root, 'SpectralHPDealiasing', 'False')
    info(root, 'RANSAdvectionTreatment', c['advection'])
    for key, val in dict(H=c['hm'], Z=.8, L=c['xmax']-c['xmin'], Ustar=r['ustar'],
            Qstar=r['q'], Kappa=r['kappa'],InletQ=r['inlet_q'],InletKappa=r['inlet_kappa'], Uinf=r['Uinf'], Uref=r['Uinf'],
            MatchingVelocity=inlet_values(c['hm'], c, r)['u'], Knear=r['k_near'],
            Kouter=r['k_outer'], EpsilonOuter=r['epsilon_outer'], EpsilonInletH=r['epsilon_inlet_near'], EpsilonUnitH=r['epsilon_unit_near'], OuterLength=c['outer_length'], EpsilonWall=r['wall_epsilon_near'],
            EpsilonOriginalH=r['k_near']**1.5/(2.55*c['hm']),
            Kinvis=r['nu'], ResolvedNu=r['nu'] if c['diffusion']=='Molecular' else 0.,
            Span=c['span'], RANSRoughnessHeight=c['ks'], RANSHeightShift=0.,
            RANSStripHeight=c['hm'], RANSStripPointP=c['point_p'], RANSStripPointN=c['point_n'],
            RANSTerrainHeight=.04 if c['geometry'] == 'hill' else 0., RANSTerrainHalfLength=.1,
            RANSTerrainSlope=0., RANSAWFTransitionStart=c['ramp_end']-.04*c['ramp_width_h'],
            RANSAWFTransitionEnd=c['ramp_end'], RANSExperimentDisplacement=r['displacement'],
            RANSExperimentZ0=r['z0'],
            RANSExperimentBlendBegin=begin, RANSExperimentBlendEnd=end,
            RANSReferenceMonitor=0, RANSStageAudit=0, RANSAWFPointDiagnostics=1,
            RANSStartupDuration=c['startup_flow_times']*(c['xmax']-c['xmin'])/r['Uinf'],
            RANSWallGeometryTolerance=1e-7, RANSNonlinearTolerance=1e-8,
            RANSPressureTolerance=1e-8, IO_CheckSteps=c['checkpoint'],
            TimeStep=c['dt'], NumSteps=c['steps']).items():
        parameter(root, key, val)
    if c['model'] == 'standard':
        parameter(root, 'RANSZ0', r['z0'])
    expr = expressions(c, r)
    names = [e.text for e in root.findall('CONDITIONS/VARIABLES/V')]
    values = {n: expr['eps' if n == 'eps_tilde' else n] for n in names}
    for group in list(conditions.findall('FUNCTION')):
        if group.attrib['NAME'] in ('InitialConditions', 'ExactSolution', 'RANSExperimentalInlet'):
            conditions.remove(group)
    for label in ('InitialConditions', 'RANSExperimentalInlet'):
        group = element(conditions, 'FUNCTION', NAME=label)
        for name in names:
            element(group, 'E', VAR=name, VALUE=values[name])
    for region in root.findall('CONDITIONS/BOUNDARYCONDITIONS/REGION'):
        b = int(region.attrib['REF'])
        if b == 0:
            # VCS reapplies Dirichlet data before the first implicit wall
            # update. The startup trace must match the prepared inlet collar,
            # including the actual K amplitude and the selected wall policy.
            for entry in region:
                if entry.attrib.get('VAR') in ('eps','eps_tilde'):
                    entry.set('VALUE',wall_epsilon_expression(c))
            continue
        if b >= 4:
            continue
        region.clear(); region.attrib['REF'] = str(b)
        for name in names:
            if name == 'p':
                element(region, 'N', VAR=name, VALUE='0', USERDEFINEDTYPE='RANSPressure')
            elif b == 3:
                element(region, 'D', VAR=name, VALUE=values[name])
            elif b == 1:
                attrs = {'USERDEFINEDTYPE': 'RANSOutflow'} if name in names[:c['dimension']] else {}
                element(region, 'N', VAR=name, VALUE='0', **attrs)
            else:
                roof = r['Uinf'] if name == 'u' else r['k_outer'] if name == 'k' else r['epsilon_outer'] if name in ('eps', 'eps_tilde') else 0.
                element(region, 'D', VAR=name, VALUE=f'{roof:.17g}')
    fingerprint = hashlib.sha256(json.dumps(dict(version=VERSION, config={k:v for k,v in c.items() if k not in ('steps', 'checkpoint')},
                           reference=r, data_snapshot=sha(DATA/'snapshot.json')), sort_keys=True).encode()).hexdigest()
    info(root, 'RANSExperimentalSignature', fingerprint)
    tree.write(path, encoding='utf-8', xml_declaration=True)
    receipt = dict(version=VERSION, status='provisional_unfitted', config=c, reference=r,
                   wall_route=route_metadata(path),
                   physics_fingerprint=fingerprint, case_sha256=sha(path),
                   cv1_version=frozen['version'], cv1_snapshot_sha256=sha(DATA/'snapshot.json'),
                   inference=dict(epsilon=c['epsilon_profile']+' prescribed inference; transported epsilon is unmodified and unmeasured',
                                  k=('variance_factor*sigma_u_fit^2, with common k_scale and a constant N..H collar; unmeasured span variance'
                                     if c['k_profile']=='MeasuredRms' else 'Q*u_tau^2 near wall; outer isotropy assumption; no reserved vertical-rms fit'),
                                  scale='assumed physical viscosity' if 'assumed_nu' in c else 'analysis units U_H=1; Re_H fixed; absolute airspeed unknown',
                                  displacement='observational only; physical Yap distance and wall origin unchanged'),
                   roof='0.8m Dirichlet free stream with nonzero inferred turbulence; facility/domain sensitivity pending CV4',
                   outlet='zero prescribed residual traction RANSOutflow; natural scalar flux',
                   preparation_required=True, calibration_complete=False)
    write_json(output/'manifest.json', receipt)
    return receipt


def environment(build):
    build = Path(build).resolve(); os.environ.update(THREADS)
    library = build/'geometry-lib/libSpatialDomains.so'
    if not library.exists():
        raise ValueError('Private full-precision SpatialDomains library required')
    os.environ['LD_LIBRARY_PATH'] = ':'.join(dict.fromkeys(
        [str(library.parent)]+[p for p in os.environ.get('LD_LIBRARY_PATH','').split(':') if p]))
    return build


def execute(command, directory, label, timeout=600):
    directory = Path(directory); command = list(map(str, command))
    with (directory/(label+'.log')).open('w') as log:
        code, expired = bounded_solver(command, directory, log, timeout)
    result = dict(command=command, exit_code=code, timeout=expired, threads=THREADS)
    write_json(directory/(label+'.json'), result)
    if code:
        raise RuntimeError('Command failed; see '+str(directory/(label+'.log')))
    return result


def allocation(ranks, cpus):
    ids = list(map(int, cpus.split(',')))
    if not 1 <= ranks <= 16 or len(set(ids)) != len(ids) or not ranks <= len(ids) <= 16:
        raise ValueError('Explicit allocation: 1..16 ranks, at most sixteen distinct CPUs')
    return ['taskset', '-c', cpus]+(['mpiexec', '-np', str(ranks), '--bind-to', 'none'] if ranks > 1 else [])


def prepare(case, build):
    build = environment(build); case = Path(case).resolve()
    manifest = json.loads((case.parent/'manifest.json').read_text())
    if sha(case) != manifest['case_sha256']:
        raise ValueError('Case changed before preparation')
    initial = case.parent/'initial.fld'
    command = [build/'RANSWallFieldAudit', case, '--experimental-initial', initial]
    receipt = execute(command, case.parent, 'prepare', 300)
    receipt['geometry'] = json.loads(Path(str(initial)+'.geometry.json').read_text())
    # Independent quadrature and streamfunction audit; SciPy is a validation dependency only.
    from scipy.integrate import quad
    c, ref = manifest['config'], manifest['reference']
    mismatch = normal = 0.; minimum = {'k': math.inf, 'epsilon': math.inf}
    span_values = {}; expected_cache = {}
    def expected_values(x, z, eta, nx, nz, foot_x):
        amp = .04 if c['geometry'] == 'hill' else 0.
        a = math.pi*x/.1
        h, hp = (amp/2*(1+math.cos(a)), -amp*math.pi/.2*math.sin(a)) if abs(x) < .1 else (0., 0.)
        e = z-h*((.8-z)/(.8-h))**2
        b = inlet_values(e, c, ref,foot_x); v = inlet_values(eta, c, ref,foot_x)
        t = max(0., min(1., (eta-c['point_n'])/c['point_n']))
        chi = blend(t); derivative = -30*t*t*(1-t)**2/c['point_n']
        integral = quad(lambda s: inlet_values(s, c, ref)['u'], e, eta,
                        epsabs=1e-12, epsrel=2e-12, points=[p for p in (.2, .3, c['point_n']) if min(e, eta)<p<max(e, eta)])[0] if derivative else 0.
        a = chi*v['u']+derivative*integral
        gx = (1-chi)*b['u']*(-hp*(.8-z)**2*(.8+h)/(.8-h)**3)+a*nx
        gz = (1-chi)*b['u']*(1+2*h*(.8-z)/(.8-h)**2)+a*nz
        expected = dict(u=gz, v=0., w=0.)
        expected['v' if c['dimension']==2 else 'w'] = -gx
        for name in ('k', 'p', 'eps', 'eps_tilde'):
            key = 'epsilon' if name in ('eps', 'eps_tilde') else name
            expected[name] = chi*v[key]+(1-chi)*b[key]
        return expected, abs(gz*nx-gx*nz) if abs(eta-c['hm']) < 1e-9 else 0.
    # All fields and periodic span copies share the same independent 2D
    # streamfunction integral. Stream the large audit and compute it once.
    def input_rows():
        with Path(str(initial)+'.csv').open(newline='') as stream:
            yield from csv.DictReader(stream)
    for row in input_rows():
        name = row['field']; x = float(row['x']); z = float(row['y' if c['dimension'] == 2 else 'z'])
        eta = float(row['distance']); nx, nz = float(row['normal_x']), float(row['normal_z'])
        foot_x=float(row['foot_x'])
        key = tuple(round(v, 14) for v in (x, z, eta, nx, nz,foot_x))
        if key not in expected_cache:
            expected_cache[key] = expected_values(x, z, eta, nx, nz,foot_x)
        expected, defect = expected_cache[key]; expected = expected[name]
        normal = max(normal, defect)
        raw, projected = float(row['input']), float(row['projected'])
        if not math.isfinite(raw+projected):
            raise RuntimeError('Nonfinite prepared field')
        mismatch = max(mismatch, abs(raw-expected)/max(1., abs(expected)))
        if name in ('k', 'eps', 'eps_tilde'):
            key = 'k' if name == 'k' else 'epsilon'; minimum[key] = min(minimum[key], projected)
        if c['dimension'] == 3:
            key = (name, round(x, 10), round(z, 10))
            low, high = span_values.get(key, (projected, projected))
            span_values[key] = min(low, projected), max(high, projected)
    span_error = max(((b-a)/max(1., abs(a), abs(b)) for a,b in span_values.values()), default=0.)
    if mismatch > 5e-10 or normal > 1e-8 or min(minimum.values()) <= 0 or span_error > 5e-10:
        raise RuntimeError(f'Prepared-field audit failed: {mismatch}, {normal}, {minimum}, {span_error}')
    tree = ET.parse(case); group = tree.getroot().find("CONDITIONS/FUNCTION[@NAME='InitialConditions']")
    group.clear(); group.attrib['NAME'] = 'InitialConditions'
    for v in tree.getroot().findall('CONDITIONS/VARIABLES/V'):
        element(group, 'F', VAR=v.text, FILE=initial.name)
    tree.write(case, encoding='utf-8', xml_declaration=True)
    receipt.update(maximum_relative_input_error=mismatch, raw_interface_normal_max=normal,
                   minimum_projected=minimum, scaled_span_variation=span_error, initial_sha256=sha(initial))
    manifest.update(preparation_required=False, case_sha256=sha(case), preparation=receipt)
    write_json(case.parent/'manifest.json', manifest)
    return receipt


def probes(surface, scope, output, tail_heights=None):
    frozen_data(); series = {s['id']:s for s in json.loads((DATA/'series.json').read_text())}
    rows, skipped = [], []
    if scope=='dividing':
        heights=tail_heights or ([.002,.004,.006] if surface=='S' else [.007,.009,.011])
        if len(heights)<3 or len(set(heights))!=len(heights) or any(not math.isfinite(h) or h<=0 for h in heights):
            raise ValueError('Declare at least three distinct positive tail heights')
        rows=[dict(id=f'dividing:{j}:{i}',x=.11+.01*i,span=0,z=h)
              for j,h in enumerate(sorted(heights)) for i in range(30)]
    for p in records(DATA/'profiles.csv'):
        s = series[p['series_id']]
        if scope=='dividing' or s['case'] != (surface if scope == 'flat' else surface*2):
            continue
        if scope == 'flat' and s['role'] not in ('calibration', 'flat_crosscheck'):
            continue
        if scope != 'flat' and s['role'] not in ('validation', 'dimensional_holdout', 'wake_mean_crosscheck'):
            continue
        if 'delta75' in s['height_kind']:
            skipped.append(dict(id=p['series_id']+':'+p['point_id'], reason='experimental_delta75_unknown')); continue
        z = .04*float(p['y_over_H']); x = .04*float(p['x_over_H'])
        if scope in ('reference','local-reference'):
            z=.04*float(p['yprime_over_H'])
            if scope=='reference':x=0.
        rows.append(dict(id=p['series_id']+':'+p['point_id'], x=x, span=0, z=z))
    if scope in ('hill','reference','local-reference'):
        rows.append(dict(id='crest_speedup',x=0.,span=0,z=.05 if scope=='hill' else .01))
    write_csv(output, rows, ['id', 'x', 'span', 'z'])
    result = dict(version=VERSION, scope=scope, surface=surface, points=len(rows),
                  unavailable=skipped, csv_sha256=sha(output),
                  cv1_snapshot_sha256=sha(DATA/'snapshot.json'),
                  normalization='global laboratory X/span/Z; physical global heights from frozen CV1 conversion',
                  reference_station_x=0. if scope=='reference' else None,
                  tail_heights_m=sorted({r['z'] for r in rows}) if scope=='dividing' else None,
                  dividing_method='x=0.11..0.40m, step 0.01m; declared low-tail height sensitivity required' if scope=='dividing' else None)
    write_json(str(output)+'.json', result)
    return result


def run(case, build, ranks, cpus, timeout=1800):
    case = Path(case).resolve(); build = environment(build)
    manifest = json.loads((case.parent/'manifest.json').read_text())
    if manifest.get('wall_route',route_metadata(case))!=route_metadata(case):
        raise ValueError('Wall route differs from its producing manifest')
    if manifest['preparation_required'] or sha(case) != manifest['case_sha256']:
        raise ValueError('Prepare the unmodified case before running')
    solver_hash = sha(build/'IncNavierStokesSolver')
    geometry_hash = sha(build/'geometry-lib/libSpatialDomains.so')
    command = allocation(ranks, cpus)+[build/'IncNavierStokesSolver', case.name, '--use-scotch']
    provenance=dict(ranks=ranks, cpus=cpus, case_sha256=sha(case),
        solver_sha256=solver_hash, geometry_sha256=geometry_hash,
        physics_fingerprint=manifest['physics_fingerprint'], cv1_snapshot_sha256=sha(DATA/'snapshot.json'))
    write_json(case.parent/'run-start.json',dict(command=list(map(str,command)),**provenance))
    try:
        result = execute(command, case.parent, 'run', timeout)
    except RuntimeError:
        result=json.loads((case.parent/'run.json').read_text());result.update(provenance,complete=False)
        result['binaries_unchanged_during_run']=(solver_hash==sha(build/'IncNavierStokesSolver') and geometry_hash==sha(build/'geometry-lib/libSpatialDomains.so'))
        try:result['partial_stages']=stage_metrics(case)
        except (ValueError,FileNotFoundError):pass
        write_json(case.parent/'run.json',result)
        raise
    result.update(provenance,complete=True)
    result['stages'] = stage_metrics(case)
    result['binaries_unchanged_during_run'] = (
        solver_hash == sha(build/'IncNavierStokesSolver') and
        geometry_hash == sha(build/'geometry-lib/libSpatialDomains.so'))
    write_json(case.parent/'run.json', result)
    if not result['binaries_unchanged_during_run']:
        raise RuntimeError('Validation binaries changed during the run; retain the run as unaccepted')
    return result


def stage_metrics(case):
    case = Path(case)
    stages = records(case.with_suffix('.rans.csv'))
    wall = records(case.with_suffix('.terrain-wall.csv'))
    if not stages or not wall:
        raise ValueError('Missing completed RANS stages/wall diagnostics')
    return dict(completed_stages=len(stages), final_time=float(stages[-1]['time']),
                maximum_weak_divergence=max(float(r['weak_divergence_l2']) for r in stages),
                minimum_k=min(float(r['min_k']) for r in stages),
                minimum_epsilon=min(float(r['min_eps']) for r in stages),
                maximum_guard_activations=max(float(r['guard_activations']) for r in stages),
                maximum_nonpositive_linear_images=max(float(r.get('nonpositive_linear_images',0)) for r in stages),
                maximum_startup_weight=max(float(r.get('startup_weight',0)) for r in stages),
                final_startup_weight=float(stages[-1].get('startup_weight',0)),
                last_forced_stage_time=max((float(r['time']) for r in stages if float(r.get('startup_weight',0))>0),default=0.),
                maximum_projected_normal=max(float(r['normal_projected_max']) for r in wall),
                maximum_strong_normal=max(float(r['normal_pointwise_max']) for r in wall),
                maximum_normal_integral_defect=max(abs(float(r['normal_integral_defect'])) for r in wall),
                maximum_resolved_mass_defect=max(abs(float(r['resolved_mass_defect'])) for r in wall),
                maximum_N_inverse_residual=max(float(r['N_inverse_residual']) for r in wall),
                units='raw solver units; U_H=1 analysis scaling unless an assumed physical nu was declared')


def observe(case, field, points, output, build, ranks=1, cpus='1'):
    case, field, points, output = (Path(p).resolve() for p in (case, field, points, output))
    build = environment(build); output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = case.parent/'manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text())['case_sha256'] != sha(case):
        raise ValueError('Producing XML differs from its case manifest')
    work = output.parent/(output.stem+'-native')
    work.mkdir(exist_ok=False)
    shutil.copy2(case, work/'sample.xml')
    # Native postprocessing is read-only: no solver diagnostic, history or partition files are reused.
    command = allocation(ranks, cpus)+[build/'RANSObservationAudit', work/'sample.xml', field, points, output]
    result = execute(command, work, 'observe')
    native = json.loads(Path(str(output)+'.native.json').read_text())
    fingerprint = ET.parse(case).getroot().find("CONDITIONS/SOLVERINFO/I[@PROPERTY='RANSExperimentalSignature']")
    result.update(case_sha256=sha(case), field_sha256=sha(field) if field.is_file() else None,
                  wall_route=route_metadata(case),
                  points_sha256=sha(points), output_sha256=sha(output),
                  observer_sha256=sha(build/'RANSObservationAudit'),
                  geometry_sha256=sha(build/'geometry-lib/libSpatialDomains.so'),
                  observation_configuration=native['configuration'],
                  experimental_fingerprint=fingerprint.attrib['VALUE'] if fingerprint is not None else None,
                  manifest_sha256=sha(manifest_path) if manifest_path.exists() else None, ranks=ranks, cpus=cpus)
    # Bind exclusion of temporary startup forcing to the producing XML, not
    # merely an optional description in a sidecar manifest.
    duration=[n.text.split('=',1)[1].strip() for n in ET.parse(case).getroot().findall('CONDITIONS/PARAMETERS/P')
              if n.text.split('=',1)[0].strip()=='RANSStartupDuration']
    result['startup_end_time']=float(duration[0]) if duration else 0.
    write_json(str(output)+'.json', result)
    return result


def statistics(inputs, output, verification=False):
    """Trapezoidal physical-time averages, modelled stress plus resolved covariance."""
    if not verification:
        provenance = [json.loads(Path(str(p)+'.json').read_text()) for p in inputs]
        if any(r.get('output_sha256') != sha(p) for r,p in zip(provenance,inputs)):
            raise ValueError('Observation receipt/hash mismatch')
        for key in ('observation_configuration','experimental_fingerprint','points_sha256','observer_sha256','geometry_sha256'):
            if len({r[key] for r in provenance}) != 1:
                raise ValueError('Inconsistent observation provenance: '+key)
        if len({json.dumps(r.get('wall_route'),sort_keys=True) for r in provenance})!=1:
            raise ValueError('Inconsistent observation provenance: wall_route')
        if not provenance[0]['experimental_fingerprint']:
            raise ValueError('Experimental statistics require a bound experimental case')
    snapshots = [records(p) for p in inputs]
    snapshots.sort(key=lambda s: float(s[0]['time']))
    times = [float(s[0]['time']) for s in snapshots]
    if len(times) < 2 or any(b <= a for a,b in zip(times, times[1:])):
        raise ValueError('Supply at least two distinct physical times')
    if not verification and times[0] <= 0:
        raise ValueError('Experimental statistics exclude zero-time/startup fields')
    if not verification and any(t<r.get('startup_end_time',0.) for t,r in
            zip([float(records(p)[0]['time']) for p in inputs],provenance)):
        raise ValueError('Experimental statistics exclude the forced startup continuation')
    ids = [r['id'] for r in snapshots[0]]
    for s in snapshots:
        if [r['id'] for r in s] != ids or len({r['time'] for r in s}) != 1:
            raise ValueError('Inconsistent probe IDs or times')
        if any(r['origin'] != 'checkpoint' for r in s):
            raise ValueError('Statistics require evolved native checkpoints, not prepared/manufactured fields')
    weights = [(times[1]-times[0])/2]+[(times[i+1]-times[i-1])/2 for i in range(1,len(times)-1)]+[(times[-1]-times[-2])/2]
    weights = [w/(times[-1]-times[0]) for w in weights]
    def average(values):
        return sum(w*v for w,v in zip(weights, values))
    result = []
    for i in range(len(ids)):
        samples = [s[i] for s in snapshots]; first = samples[0]
        if any(any(r[k] != first[k] for k in ('x', 'span', 'z', 'kind')) for r in samples):
            raise ValueError('Probe positions or support changed during averaging')
        row = {k:first[k] for k in ('id', 'x', 'span', 'z', 'kind')}
        row.update(time_begin=times[0], time_end=times[-1], samples=len(times),
                   origin='verification_statistics' if verification else 'checkpoint_statistics')
        means = {}
        for key in ('u', 'v', 'w', 'k', 'epsilon_model'):
            row[key] = average([float(r[key]) for r in samples]) if all(r[key] != '' for r in samples) else ''
            means[key] = row[key]
        for key, a, b in (('R_xx','u','u'), ('R_yy','v','v'), ('R_zz','w','w'),
                          ('R_xy','u','v'), ('R_xz','u','w'), ('R_yz','v','w')):
            row[key] = ''
            if all(r[key] != '' for r in samples):
                row[key] = average([float(r[key])+(float(r[a])-means[a])*(float(r[b])-means[b]) for r in samples])
        result.append(row)
    write_csv(output, result)
    receipt = dict(version=VERSION, verification_only=verification, time_begin=times[0], time_end=times[-1], samples=len(times),
                   weighting='trapezoidal physical time; total R=mean(R_model)+covariance(U_resolved)',
                   stationarity='not assessed; CV4 must establish a stationary window',
                   input_sha256={str(p):sha(p) for p in inputs}, output_sha256=sha(output),
                   negative_diagonal_count=sum(float(r[k]) < 0 for r in result for k in ('R_xx','R_yy','R_zz') if r[k] != ''))
    if not verification:
        receipt.update(experimental_fingerprint=provenance[0]['experimental_fingerprint'],
                       observation_configuration=provenance[0]['observation_configuration'],
                       wall_route=provenance[0].get('wall_route'))
    write_json(str(output)+'.json', receipt)
    return receipt


def compare(stats, manifest_path, output, calibration=False, flat_stats=None, flat_manifest=None, dividing_stats=None):
    frozen_data(); manifest = json.loads(Path(manifest_path).read_text())
    if manifest.get('cv1_snapshot_sha256') != sha(DATA/'snapshot.json'):
        raise ValueError('Case does not bind the frozen CV1 snapshot')
    surface = manifest['config']['surface']; ref = manifest['reference']
    if not calibration and (not manifest.get('calibration_complete') or not manifest.get('calibration_version')):
        raise ValueError('Freeze a versioned CV3 calibration before inspecting reserved or hill comparison scores')
    statistics_receipt = json.loads(Path(str(stats)+'.json').read_text())
    if statistics_receipt.get('output_sha256') != sha(stats) or statistics_receipt.get('experimental_fingerprint') != manifest['physics_fingerprint']:
        raise ValueError('Statistics differ from the candidate manifest')
    predicted = {p['id']:p for p in records(stats)}
    flat = {p['id']:p for p in records(flat_stats)} if flat_stats else {}
    if flat_stats:
        if flat_manifest is None:
            raise ValueError('Matched-flat statistics require --flat-manifest')
        fm=json.loads(Path(flat_manifest).read_text());fr=json.loads(Path(str(flat_stats)+'.json').read_text())
        transfer=lambda c:{k:v for k,v in c.items() if k not in ('geometry','steps','checkpoint','startup_flow_times')}
        if fm['config']['geometry']!='floor' or transfer(fm['config'])!=transfer(manifest['config']) or fm['reference']!=ref:
            raise ValueError('Flat reference must transfer the identical surface/inlet/model/wall configuration')
        if fr.get('output_sha256')!=sha(flat_stats) or fr.get('experimental_fingerprint')!=fm['physics_fingerprint']:
            raise ValueError('Unbound matched-flat statistics')
    def matched(r, ident):
        f=flat.get(ident)
        if f is None:return None
        x=float(r['x']);h=.02*(1+math.cos(math.pi*x/.1)) if abs(x)<.1 else 0.
        if abs(float(f['x']))>1e-10 or abs(float(f['z'])-(float(r['z'])-h))>1e-9 or f['origin']!='checkpoint_statistics':
            raise ValueError('Flat reference requires x=0 and the same physical vertical yprime')
        if any(abs(float(r[k])-float(f[k]))>1e-12 for k in ('time_begin','time_end')):
            raise ValueError('Matched hill/floor statistical windows differ')
        return f
    series = json.loads((DATA/'series.json').read_text())
    source = records(DATA/'profiles.csv'); results = []
    floor_case=manifest['config']['geometry']=='floor'
    if calibration and manifest['config']['geometry'] != 'floor':
        raise ValueError('Calibration accepts empty-floor cases only')
    for s in series:
        if s['case'] != (surface if calibration or floor_case else surface*2):
            continue
        if calibration and s['role'] != 'calibration':
            continue
        if not calibration and s['role'] != ('flat_crosscheck' if floor_case else 'validation'):
            continue
        points = [p for p in source if p['series_id'] == s['id'] and float(p['weight']) > 0]
        used = []; missing = []; negative = []
        for p in points:
            ident = s['id']+':'+p['point_id']; r = predicted.get(ident)
            value = None; quantity = s['quantity']; unsupported_reason=None
            if r is not None:
                if r.get('origin') != 'checkpoint_statistics':
                    raise ValueError('Scoring requires native checkpoint statistics')
                if quantity == 'U_over_Uinf' and r['u'] != '':
                    value = float(r['u'])/ref['Uinf']
                elif quantity == 'minus_uv_over_Uinf2' and r['R_xz'] != '':
                    value = -float(r['R_xz'])/ref['Uinf']**2
                elif quantity in ('sigma_u_over_utau', 'sigma_u_over_Uinf', 'sigma_vertical_over_utau'):
                    key = 'R_zz' if 'vertical' in quantity else 'R_xx'
                    if r[key] != '':
                        if float(r[key]) < 0:
                            negative.append(ident)
                        else:
                            value = math.sqrt(float(r[key]))/(ref['ustar'] if 'utau' in quantity else ref['Uinf'])
                elif quantity in ('speedup','sigma_u_over_flat_same_height','uv_over_flat_same_height'):
                    f=matched(r,ident)
                    if f is not None:
                        if quantity=='speedup' and r['u']!='' and f['u']!='':
                            answer=speedup(float(r['u']),float(f['u']),ref['Uinf'])
                            if answer['available']:value=answer['value']
                            else:unsupported_reason=answer['reason']
                        elif quantity=='sigma_u_over_flat_same_height' and r['R_xx']!='' and f['R_xx']!='':
                            if float(r['R_xx'])>=0 and float(f['R_xx'])>1e-20*ref['Uinf']**2:
                                value=math.sqrt(float(r['R_xx'])/float(f['R_xx']))
                            elif min(float(r['R_xx']),float(f['R_xx']))<0:negative.append(ident)
                        elif quantity=='uv_over_flat_same_height' and r['R_xz']!='' and f['R_xz']!='':
                            if abs(float(f['R_xz']))>1e-12*ref['Uinf']**2:value=float(r['R_xz'])/float(f['R_xz'])
            if value is None:
                missing.append(dict(id=ident, reason=unsupported_reason or ('missing' if r is None else r['kind'] if r['kind'] != 'resolved_native' else 'operator_or_reference_unavailable')))
            else:
                used.append((float(p['weight']), value-float(p['value']), float(p['value'])))
        weight = sum(w for w,_,_ in used); total = sum(float(p['weight']) for p in points)
        rmse = math.sqrt(sum(w*e*e for w,e,_ in used)/weight) if weight else None
        results.append(dict(series=s['id'], role=s['role'], quantity=s['quantity'],
                            complete=not missing, supported_weight=weight, total_weight=total,
                            conditional_rmse=rmse, complete_range_rmse=rmse if not missing else None,
                            conditional_relative_rms_error=rmse/math.sqrt(sum(w*q*q for w,_,q in used)/weight)
                                if weight and rmse is not None and s['quantity']=='sigma_u_over_Uinf' and sum(w*q*q for w,_,q in used)>0 else None,
                            conditional_bias=sum(w*e for w,e,_ in used)/weight if weight else None,
                            missing=missing, negative_diagonals=negative))
    scalars={}
    if not calibration and not floor_case:
        target={r['id']:float(r['value']) for r in records(DATA/'reported_scalars.csv')
                if r['id'].startswith(('crest_S_','reattachment_'))}
        r=predicted.get('crest_speedup');f=matched(r,'crest_speedup') if r else None
        if r and (abs(float(r['x']))>1e-10 or abs(float(r['z'])-.05)>1e-10):
            raise ValueError('Crest speedup requires x=0, yprime/H=0.25')
        s=speedup(float(r['u']),float(f['u']),ref['Uinf']) if r and f and r['u']!='' and f['u']!='' else dict(available=False,reason='matching crest/flat probes required')
        if s['available']:
            s.update(reported=target['crest_S_'+surface*2],absolute_error=abs(s['value']-target['crest_S_'+surface*2]))
        scalars['crest_speedup']=s
        xr=dict(available=False,reason='dividing-line statistical probes required')
        if dividing_stats:
            dr=json.loads(Path(str(dividing_stats)+'.json').read_text())
            if dr.get('output_sha256')!=sha(dividing_stats) or dr.get('experimental_fingerprint')!=manifest['physics_fingerprint']:
                raise ValueError('Unbound dividing-line statistics')
            xr=reattachment(records(dividing_stats))
            if xr['available']:xr.update(reported=target['reattachment_'+surface*2],absolute_error_over_H=abs(xr['x_R_over_H']-target['reattachment_'+surface*2]))
        scalars['reattachment']=xr
    receipt = dict(version=VERSION, operation='flat_candidate_objective' if calibration else 'reserved_flat_crosscheck' if floor_case else 'heldout_profile_comparison',
                   calibration_complete=bool(manifest.get('calibration_complete')),calibration_version=manifest.get('calibration_version'),
                   stationarity_accepted=False,
                   note='Observable evaluation; developed CV3 calibration and CV4 stationarity/convergence are separate gates.',
                   cv1_snapshot_sha256=sha(DATA/'snapshot.json'), statistics_sha256=sha(stats),
                   manifest_sha256=sha(manifest_path), results=results,scalars=scalars,
                   deferred=['wall-shear zero requires a separately averaged native traction trace; it is not the dividing-line x_R',
                             'delta75-dependent physical heights', 'dimensional near-floor amplitudes',
                             'turbulence ratios without matched flat profile'])
    write_json(output, receipt)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='command', required=True)
    g = sub.add_parser('generate'); g.add_argument('--config', type=Path, required=True); g.add_argument('--output', type=Path, required=True)
    g = sub.add_parser('prepare'); g.add_argument('case', type=Path); g.add_argument('--build', type=Path, required=True)
    g = sub.add_parser('probes'); g.add_argument('--surface', choices=('S','R'), required=True); g.add_argument('--scope', choices=('flat','hill','reference','local-reference','dividing'), required=True); g.add_argument('--tail-heights',help='Comma-separated physical vertical tail heights in metres'); g.add_argument('--output', type=Path, required=True)
    g = sub.add_parser('run'); g.add_argument('case', type=Path); g.add_argument('--build', type=Path, required=True); g.add_argument('--ranks', type=int, default=1); g.add_argument('--cpus', required=True); g.add_argument('--timeout', type=int, default=1800)
    g = sub.add_parser('observe'); g.add_argument('case', type=Path); g.add_argument('field', type=Path); g.add_argument('points', type=Path); g.add_argument('--output', type=Path, required=True); g.add_argument('--build', type=Path, required=True); g.add_argument('--ranks', type=int, default=1); g.add_argument('--cpus', required=True)
    g = sub.add_parser('statistics'); g.add_argument('inputs', type=Path, nargs='+'); g.add_argument('--output', type=Path, required=True)
    for name in ('calibrate', 'compare'):
        g = sub.add_parser(name); g.add_argument('statistics', type=Path); g.add_argument('--manifest', type=Path, required=True); g.add_argument('--output', type=Path, required=True); g.add_argument('--flat-statistics', type=Path)
        g.add_argument('--flat-manifest',type=Path);g.add_argument('--dividing-statistics',type=Path)
    a = p.parse_args()
    if a.command == 'generate': result = generate(json.loads(a.config.read_text()), a.output)
    elif a.command == 'prepare': result = prepare(a.case, a.build)
    elif a.command == 'probes': result = probes(a.surface, a.scope, a.output, list(map(float,a.tail_heights.split(','))) if a.tail_heights else None)
    elif a.command == 'run': result = run(a.case, a.build, a.ranks, a.cpus, a.timeout)
    elif a.command == 'observe': result = observe(a.case, a.field, a.points, a.output, a.build, a.ranks, a.cpus)
    elif a.command == 'statistics': result = statistics(a.inputs, a.output)
    else: result = compare(a.statistics, a.manifest, a.output, a.command == 'calibrate', a.flat_statistics,a.flat_manifest,a.dividing_statistics)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
