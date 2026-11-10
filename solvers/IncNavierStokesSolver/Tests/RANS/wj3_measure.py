#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Physical comparison and independent replay for immutable native WJ packets.

No governing solve or projection is performed here. SVD is used only to inspect
already completed frozen systems, never to generate a trajectory correction.
"""
import math
from functools import lru_cache
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.linalg import svdvals
from scipy.optimize import minimize_scalar
from wj2_cases import normal_knots
from test_wj2_native import wall,independent,QE,KAPPA,CSTAR
from run_cls_epawf_picard import load,sha,require

FIELDS=('u','v','k','outer_epsilon','reduced_pressure','mechanical_pressure')
K0=QE*.25


def observations(order=6):
    # Identical physical rule on every mesh, split at the finest supported
    # elements AND physical breakpoints; zero-weight interface/end samples
    # augment sampled Linf without changing physical L2 integrals.
    _,zv,_,_,zc=wall(K0)
    cuts=sorted(set(normal_knots(2)+[p for p in (zc,zv,.1,.5,1.) if 0<p<2]))
    gx,wx=np.polynomial.legendre.leggauss(4);gz,wz=np.polynomial.legendre.leggauss(order)
    rows=[]
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        for x,w in zip((gx+1)/2,wx/2):
            for z,v in zip(lo+(gz+1)*(hi-lo)/2,wz*(hi-lo)/2):rows.append([x,0.,z,w*v,0.])
    for z in sorted(set(cuts+[0.,1.,2.])):
        for x in (0.,.25,.5,.75,1.):
            for side in ((-1,1) if z==1 else (0,)):rows.append([x,0.,z,0.,side])
    return rows


@lru_cache(maxsize=None)
def continuum(z,inner,policy):
    a,zv,t,A,c=wall(K0);nu=1.5e-5;d=nu-a*zv
    def J(s):return math.log1p(a*s/d)/a
    def Q(s):return A*(s/c if s<c else 1+math.log(s/c))-t*t*(J(s)-nu*s/(d*(d+a*s)))
    k=K0
    if inner and policy=='Budget':
        cuts=sorted(set([z,1.]+[p for p in (c,.1) if z<p<1]))
        k-=sum(quad(lambda s:Q(s)/(d+a*s),lo,hi,epsabs=1e-12,epsrel=3e-12)[0] for lo,hi in zip(cuts[:-1],cuts[1:]))
    u=t*J(z) if inner else t*J(1)+.5/KAPPA*math.log(z)
    jump=-K0/(81/20);pi=jump/2 if inner else -jump/2
    return np.array([u,0.,k,0. if inner else .125/(KAPPA*z),pi,pi-2*k/3])


def arrays(measurement,policy):
    rows=measurement['observations']
    return (np.array([r['value'] for r in rows]),np.array([r['represented_initial'] for r in rows]),
        np.array([continuum(r['position'][2],r['inner'],policy) for r in rows]),
        np.array([r['weight'] for r in rows]),np.array([not r['inner'] for r in rows]))


def norms(difference,weights,outer):
    d=np.array(difference,copy=True);w=np.array(weights);volume=sum(w)
    # ONE overall gauge for both pressure traces and mechanical pressure.
    shift=float(np.dot(w,d[:,4])/volume)
    d[:,4]-=shift;d[:,5]-=shift
    result={}
    for j,name in enumerate(FIELDS):
        active=outer if j==3 else np.ones(len(w),dtype=bool)
        result[name]=dict(L2=float(np.sqrt(np.dot(w[active],d[active,j]**2))),
            sampled_Linf=float(max(abs(d[active,j]))))
    return dict(fields=result,overall_pressure_gauge_shift=shift)


def errors(measurement,policy):
    value,initial,reference,w,outer=arrays(measurement,policy)
    flux=measurement['interfaces'];jump=-K0/(81/20)
    return dict(continuum_error=norms(value-reference,w,outer),initial_drift=norms(value-initial,w,outer),
        initial_representation=norms(initial-reference,w,outer),
        interface=dict(traction_max=max(float(max(abs(np.array(f['inner_outward_traction'])+f['outer_outward_traction']))) for f in flux),
            k_flux_max=max(abs(sum(f['opposing_k_diffusive_flux'])) for f in flux),
            k_value_max=max(abs(np.diff(f['inner_outer_k'])[0]) for f in flux),
            epsilon_trace_max=max(abs(np.diff(f['epsilon_trace_target'])[0]) for f in flux),
            pressure_jump_error=max(abs(f['inner_outer_mechanical_pressure'][0]-f['inner_outer_mechanical_pressure'][1]-jump) for f in flux),
            sampled_K_range=[min(f['sampled_K'] for f in flux),max(f['sampled_K'] for f in flux)],
            cutoff_range=[min(f['cutoff'] for f in flux),max(f['cutoff'] for f in flux)]),
        strong_divergence_sampled_Linf=measurement['strong_divergence_sampled_Linf'])


def compare(a,b,policy):
    require([(r['position'],r['weight'],r['side']) for r in a['observations']]==
            [(r['position'],r['weight'],r['side']) for r in b['observations']], 'unequal physical observation rules')
    av,_,_,w,outer=arrays(a,policy);bv,_,_,_,_=arrays(b,policy)
    return norms(av-bv,w,outer)


def replay_witnesses(packet):
    replay={}
    for name,w in packet['principal_witnesses'].items():
        k,e=w['k'],w['epsilon'];G=np.array(w['raw_G_row_major']).reshape(3,3)
        xi=np.array(w['xi']);a=np.array(w['polarization']);H=np.outer(a,xi)
        require(abs(a@xi)<2e-15 and abs(xi@xi-1)<2e-15 and abs(a@a-1)<2e-15,'rank-one witness lost incompressibility')
        c=independent(k,e,G);step=1e-30
        d=independent(k,e,G.astype(complex)+1j*step*H,c['scalar'][12])['R'].imag/step
        form=-np.sum(H*d)/c['scalar'][1]
        assert_allclose(form,w['minimum_over_nu_eff'],rtol=2e-10,atol=3e-12)
        assert_allclose(d.ravel(),w['dR_row_major'],rtol=2e-10,atol=3e-12)
        assert_allclose(np.linalg.eigvalsh(c['R']),w['covariance_eigenvalues'],rtol=3e-12,atol=3e-12)
        # Independent local angular minimization through complex differentiation,
        # not the production quadratic-form construction. This replays finite
        # in-plane scope; it is not a global-state stability certificate.
        def angular(theta):
            direction=np.array([math.cos(theta),0.,math.sin(theta)])
            pol=np.array([-math.sin(theta),0.,math.cos(theta)]) if name=='in_plane' else np.array([0.,1.,0.])
            h=np.outer(pol,direction)
            action=independent(k,e,G.astype(complex)+1j*step*h,c['scalar'][12])['R'].imag/step
            return -float(np.sum(h*action)/c['scalar'][1])
        theta=w['angle_radians'];opt=minimize_scalar(angular,bounds=(theta-math.pi/32,theta+math.pi/32),method='bounded',options={'xatol':1e-13})
        assert_allclose(opt.fun,w['minimum_over_nu_eff'],rtol=2e-10,atol=3e-12)
        replay[name]=dict(complex_step_form=float(form),independent_local_minimum=float(opt.fun),
            production_minimum=w['minimum_over_nu_eff'],angle=float(opt.x),state_sha256=packet['state_sha256'])
    return replay


def arithmetic(report,measurement):
    systems={};raw_floor=np.zeros(6);eps=np.finfo(float).eps
    for name in ('momentum','pressure','k','epsilon'):
        entry=report['pimple_attempts'][-1][name+'_arithmetic'];require(sha(Path(entry['path']))==entry['sha256'],'arithmetic packet changed')
        packet=load(Path(entry['path']));mt=packet['matrix'];metric=packet['mass_or_correction_map']
        require(sha(Path(mt['path']))==mt['sha256'] and sha(Path(metric['path']))==metric['sha256'],'arithmetic matrix changed')
        A=np.fromfile(mt['path'],dtype='<f8').reshape(tuple(map(int,mt['shape'])))
        M=np.fromfile(metric['path'],dtype='<f8').reshape(tuple(map(int,metric['shape'])))
        x=np.array(packet['solution']);rhs=np.array(packet['rhs']);r=A.astype(np.longdouble)@x.astype(np.longdouble)-rhs.astype(np.longdouble)
        singular=svdvals(A,check_finite=True);n=len(x);gamma=n*eps/(1-n*eps)
        residual=float(np.linalg.norm(r));formation=gamma*float(np.linalg.norm(abs(A)@abs(x)+abs(rhs)))
        coefficient=(residual+formation)/singular[-1]
        if name=='pressure':
            # Pressure coefficient uncertainty also propagates through already
            # retained BG (transpose in this file) to velocity. Euclidean-to-
            # physical conversion is bounded with the native mass spectra.
            amplification=float(svdvals(M,check_finite=True)[0]);velocity=coefficient*amplification
            raw_floor[0:2]+=velocity
            raw_floor[4:6]+=coefficient
        else:
            physical=coefficient*math.sqrt(max(np.linalg.eigvalsh(M)[-1],0))
            if name=='momentum':raw_floor[:2]+=physical
            else:raw_floor[2 if name=='k' else 3]+=physical
        systems[name]=dict(size=n,sigma_min=float(singular[-1]),sigma_max=float(singular[0]),
            condition_2=float(singular[0]/singular[-1]),long_double_residual_2=residual,
            rounding_formation_reserve_2=formation,coefficient_error_estimate=float(coefficient),
            matrix_sha256=mt['sha256'],packet_sha256=entry['sha256'])
    # A declared diagnostic estimate, not a perturbation theorem: linear growth
    # over N sequential sweeps, factor two for two compared trajectories and a
    # 100x safety reserve. Evaluation sensitivity is measured independently.
    N=int(report['steps']);floor=100*(2*N*raw_floor+np.sqrt(2)*np.array(measurement['arithmetic_evaluation_sensitivity']))
    floor[5]+=2*floor[2]/3
    return dict(systems=systems,estimated_comparison_floor=dict(zip(FIELDS,map(float,floor))),
        raw_one_step_physical_error_estimate=dict(zip(FIELDS,map(float,raw_floor))),safety_factor=100,
        propagation='2*N accumulated worst observed terminal frozen-system uncertainty; finite short-trajectory diagnostic estimate, not a rigorous nonlinear perturbation bound',
        excludes=['physical spatial residual','nonlinear one-sweep stage defect'],
        evaluation_sensitivity=measurement['arithmetic_evaluation_sensitivity'])
