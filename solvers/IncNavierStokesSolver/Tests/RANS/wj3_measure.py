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
from scipy.linalg import svd,svdvals
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
    speed=np.hypot(d[:,0],d[:,1])
    result['velocity']=dict(L2=float(np.sqrt(np.dot(w,speed**2))),sampled_Linf=float(max(speed)))
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
    systems={};raw_floor=np.zeros(6);growth=np.ones(6);velocity_metric=1.;eps=np.finfo(float).eps
    for name in ('momentum','pressure','k','epsilon'):
        entry=report['pimple_attempts'][-1][name+'_arithmetic'];require(sha(Path(entry['path']))==entry['sha256'],'arithmetic packet changed')
        packet=load(Path(entry['path']));mt=packet['matrix'];metric=packet['mass_or_correction_map']
        require(sha(Path(mt['path']))==mt['sha256'] and sha(Path(metric['path']))==metric['sha256'],'arithmetic matrix changed')
        A=np.fromfile(mt['path'],dtype='<f8').reshape(tuple(map(int,mt['shape'])))
        M=np.fromfile(metric['path'],dtype='<f8').reshape(tuple(map(int,metric['shape'])))
        x=np.array(packet['solution']);rhs=np.array(packet['rhs']);r=A.astype(np.longdouble)@x.astype(np.longdouble)-rhs.astype(np.longdouble)
        left,singular,right=svd(A,check_finite=True,full_matrices=False);n=len(x);gamma=n*eps/(1-n*eps)
        residual=float(np.linalg.norm(r));formation=gamma*float(np.linalg.norm(abs(A)@abs(x)+abs(rhs)))
        coefficient=(residual+formation)/singular[-1]
        if name=='pressure':
            # Pressure coefficient uncertainty also propagates through already
            # retained BG (transpose in this file) to velocity. Euclidean-to-
            # physical conversion is bounded with the native mass spectra.
            amplification=float(svdvals(M,check_finite=True)[0]);velocity=coefficient*amplification*velocity_metric
            raw_floor[0:2]+=velocity
            raw_floor[4:6]+=coefficient*math.sqrt(packet['physical_mass_operator_bound'])
        else:
            physical_spectrum=np.linalg.eigvalsh(M)
            physical=coefficient*math.sqrt(max(physical_spectrum[-1],0))
            if name=='momentum':
                velocity_metric=math.sqrt(packet['physical_mass_operator_bound']);raw_floor[:2]+=physical
            else:raw_floor[2 if name=='k' else 3]+=physical
        amplification=1.
        if name!='pressure':
            inverse=(right.T/singular)@left.T  # offline SVD condition analysis only
            free=np.array([i for i in range(n) if not (A[i,i]==1 and np.count_nonzero(A[i])==1)])
            Mf=M[np.ix_(free,free)];lam,V=np.linalg.eigh(Mf)
            require(min(lam)>0,'nonpositive native physical metric')
            half=(V*np.sqrt(lam))@V.T;invhalf=(V/np.sqrt(lam))@V.T
            transfer=half@inverse[np.ix_(free,free)]@(packet['alpha_mass']*Mf)@invhalf
            amplification=max(1.,float(svdvals(transfer,check_finite=True)[0]))
            indices=range(2) if name=='momentum' else [2 if name=='k' else 3]
            for i in indices:growth[i]=max(growth[i],amplification)
        systems[name]=dict(size=n,frozen_history_mass_norm_amplification=amplification,sigma_min=float(singular[-1]),sigma_max=float(singular[0]),
            condition_2=float(singular[0]/singular[-1]),long_double_residual_2=residual,
            rounding_formation_reserve_2=formation,coefficient_error_estimate=float(coefficient),
            matrix_sha256=mt['sha256'],packet_sha256=entry['sha256'],
            physical_metric_min_eigenvalue=float(physical_spectrum[0]) if name!='pressure' else None,
            physical_mass_operator_bound=packet['physical_mass_operator_bound'])
    # A declared frozen-system diagnostic estimate, not a nonlinear theorem:
    # mass-metric history-map propagation, factor two for compared trajectories,
    # 100x safety reserve. Evaluation sensitivity is measured independently.
    N=int(report['steps'])
    growth[4:6]=max(growth[:4])
    propagation=np.array([N if g<=1+1e-14 else math.expm1(N*math.log(g))/(g-1) for g in growth])
    floor=100*(2*propagation*raw_floor+np.sqrt(2)*np.array(measurement['arithmetic_evaluation_sensitivity']))
    floor[5]+=2*floor[2]/3
    floors=dict(zip(FIELDS,map(float,floor)));floors['velocity']=math.hypot(floors['u'],floors['v'])
    sampled=None;point_factors=None
    if 'evaluation_basis_row_norm_max' in measurement:
        rows=np.array(measurement['evaluation_basis_row_norm_max'])
        point_factors=np.ones(6)
        for indices,name in (([0,1],'momentum'),([2],'k'),([3],'epsilon')):
            minimum=systems[name]['physical_metric_min_eigenvalue']
            require(minimum>0,'nonpositive full frozen physical metric')
            for i in indices:point_factors[i]=rows[i]/math.sqrt(minimum)
        point_factors[4]=rows[4]/math.sqrt(systems['pressure']['physical_mass_operator_bound'])
        point_factors[5]=point_factors[4]
        values=floor*point_factors;values[5]=values[4]+2*values[2]/3
        sampled=dict(zip(FIELDS,map(float,values)));sampled['velocity']=math.hypot(values[0],values[1])
    return dict(systems=systems,estimated_comparison_floor=floors,estimated_sampled_comparison_floor=sampled,
        sampled_conversion_factors=point_factors.tolist() if point_factors is not None else None,
        raw_one_step_physical_error_estimate=dict(zip(FIELDS,map(float,raw_floor))),safety_factor=100,
        propagation='terminal frozen-system uncertainty summed over N steps using terminal frozen history-map mass-norm amplification; factor2 for a pair of trajectories; finite diagnostic estimate, not a nonlinear perturbation theorem',
        accumulated_propagation_factors=propagation.tolist(),frozen_history_map_amplifications=growth.tolist(),
        excludes=['physical spatial residual','nonlinear one-sweep stage defect'],
        evaluation_sensitivity=measurement['arithmetic_evaluation_sensitivity'])


def checkpoint_payload(directory):
    """Read immutable native binary64 values; the compiled reader independently
    validates schema/checksum/signature and all slots before any measurement.
    This function never creates or modifies a checkpoint.
    """
    import struct
    data=(Path(directory)/'State.bin').read_bytes();position=0
    def scalar(fmt):
        nonlocal position
        value=struct.unpack_from('<'+fmt,data,position)[0];position+=struct.calcsize('<'+fmt);return value
    def vector():
        nonlocal position
        n=scalar('Q');require(n<100000000,'invalid native vector size')
        result=list(struct.unpack_from('<'+str(n)+'d',data,position));position+=8*n;return result
    magic,endian,schema=[scalar('Q') for _ in range(3)];n=scalar('Q')
    signature=data[position:position+n].decode();position+=n
    order=scalar('Q');steps=scalar('Q');t,dt,eta=[scalar('d') for _ in range(3)]
    primary=vector();history=[vector() for _ in range(scalar('Q'))];times=vector();state=vector();reference=vector();defect=vector()
    require(position==len(data)-8,'checkpoint trailing data or unsupported schema')
    return dict(signature=signature,order=order,steps=steps,time=t,dt=dt,eta=eta,state=state,
        native_differential_state=primary,native_history=history,native_history_times=times,checkpoint=str(directory),
        reference_sha256=__import__('run_cls_epawf_picard').state_hash(reference))
