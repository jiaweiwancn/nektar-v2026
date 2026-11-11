#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded compiled geometry/coupled-terrain verification, independent WJ loads."""
import argparse
import functools
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time
import uuid
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.optimize import brentq
from validate_wj1 import independent, CSTAR, QE, KAPPA
from wj_h1_cases import case_xml, provisional_roughness, NU, US, HW, X, Z

THREADS={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def require(b,s):
    if not b:raise AssertionError(s)

def run(argv,folder,timeout=900,expected_failure=None):
    folder=Path(folder);log=folder/'native.log';receipt=folder/'command.json'
    require(not log.exists() and not receipt.exists(),'immutable native command output already exists')
    env=dict(os.environ,**THREADS);start=time.monotonic()
    request=dict(argv=list(map(str,argv)),cwd=str(folder.resolve()),threads=THREADS,
        timeout_seconds=timeout,rss_ceiling_bytes=8*1024**3,output_ceiling_bytes=4*1024**3)
    save(folder/'request.json',request)
    def limits():
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
        resource.setrlimit(resource.RLIMIT_FSIZE,(4*1024**3,4*1024**3))
    with log.open('wb') as f:
        try:
            proc=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=folder,env=env,
                stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits,start_new_session=True)
            code=proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=10)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            code=124
    text=log.read_text(errors='replace');rss=re.findall(r'Maximum resident set size \(kbytes\): (\d+)',text)
    result=dict(request,exit=code,seconds=time.monotonic()-start,log_sha256=sha(log),
        peak_rss_bytes=int(rss[-1])*1024 if rss else None,
        output_bytes=sum(p.stat().st_size for p in folder.rglob('*') if p.is_file()))
    save(receipt,result)
    if expected_failure is None:require(code==0,f'native command failed ({code}); retained {log}')
    else:require(code not in (0,124) and expected_failure in text,'wrong native rejection reason: '+str(log))
    require(result['output_bytes']<=4*1024**3,'native output ceiling')
    return result

def graph(s,kind,side=0):
    if kind==2:return .2*s,.2,0.
    if kind==0:return 0.,0.,0.
    a=math.pi/2.5
    if abs(s)==2.5:
        return 0.,0.,(.5*a*a if (s<0 and side>0) or (s>0 and side<0) else 0.)
    if abs(s)>2.5:return 0.,0.,0.
    return .5*(1+math.cos(a*s)),-.5*a*math.sin(a*s),-.5*a*a*math.cos(a*s)

def geometry(s,r,kind,foot=0,region=-1):
    f,fp,fpp=graph(s,kind,foot);g=math.hypot(1,fp);c=-fpp/g**3;n=np.array([-fp/g,0,1/g])
    normal=(r<HW or (r==HW and region<=0));D=.0875+(r if normal else HW)
    a=1+D*c;base=np.array([s,0,f])+D*n
    jac=np.array([[a,0,n[0]],[0,1,0],[fp*a,0,n[2]]])
    eta=0.
    if not normal:
        eta=(r-HW)/(20-HW)
        top=np.array([s,0,20.]);column=(top-base)/(20-HW)
        jac[:,0]=(1-eta)*jac[:,0]+eta*np.array([1,0,0])
        jac[:,2]=column;base=(1-eta)*base+eta*top
    return base,jac,c/(1+.0875*c),g*(1+.0875*c),n,eta

def check_geometry(packet):
    out=[]
    for m in packet['maps']:
        errors=[]
        for row in m['samples']:
            s,_,r=row['q'];expected,j,c,b,n,eta=geometry(s,r,m['kind'],row['foot'],row['region'])
            assert_allclose(row['point'],expected,rtol=0,atol=8e-14)
            assert_allclose(row['jacobian'],j,rtol=2e-13,atol=4e-14)
            assert_allclose(row['inverse'],np.linalg.inv(j),rtol=3e-13,atol=4e-14)
            assert_allclose([row['det'],row['base_area'],row['kappa_eff'],row['eta']],
                [np.linalg.det(j),b,c,eta],rtol=3e-13,atol=4e-14)
            # Independent centered/one-sided physical-gradient action, including
            # the two computational H derivatives and both C1 foot owners.
            grad=[]
            def value(a,b):
                q,*_=geometry(a,b,m['kind']);x,_,z=q
                return .7*x*x+.2*x*z+.1*z
            x,_,z=expected
            for axis,v in ((0,s),(2,r)):
                h=2e-5
                side=row['foot'] if axis==0 and row['foot'] else 0
                if axis==2 and r in (0.,HW,20.):side=1 if r==0 or (r==HW and row['region']>0) else -1
                def ev(delta):return value(s+delta if axis==0 else s,r+delta if axis==2 else r)
                if side:
                    d=side*(-3*ev(0)+4*ev(side*h)-ev(side*2*h))/(2*h)
                else:d=(ev(h)-ev(-h))/(2*h)
                grad.append(d)
            actual=np.array([grad[0],0,grad[1]])@np.array(row['inverse'])
            error=float(np.max(abs(actual-np.array([1.4*x+.2*z,0,.2*x+.1]))))
            errors.append(error)
        require(max(errors)<2e-7,'independent finite-difference physical gradient')
        require(m['inverse_error']<2e-10,'composite physical inverse')
        def v(s):
            f,fp,fpp=graph(s,m['kind']);g=math.hypot(1,fp);c=-fpp/g**3;D=.0875+HW
            return g*(HW+c*(.0875*HW+HW*HW/2))+.5*((1+D*c)*(20-f-D*g)+(20-f-D/g))
        volume=sum(quad(v,a,b,epsabs=2e-12,epsrel=2e-13)[0] for a,b in zip(X[:-1],X[1:]))
        for integral in m['integrals']:
            assert_allclose(integral['volume'],volume,rtol=2e-12,atol=2e-11)
            assert_allclose(integral['boundary_constant_flux'],0,rtol=0,atol=2e-11)
            assert_allclose(integral['boundary_x_flux'],[volume,0,0],rtol=0,atol=3e-10)
            assert_allclose(integral['boundary_z_flux'],[0,0,volume],rtol=0,atol=3e-10)
        out.append(dict(kind=m['kind'],physical_volume=volume,gradient_fd_max=max(errors),inverse_error=m['inverse_error']))
    require(packet['invalid_geometry_rejections']==5,'invalid geometry rejection coverage')
    assert_allclose(packet['frame_check']['R_error'],0,atol=2e-16)
    assert_allclose(packet['frame_check']['tangential_outer_traction'],US*US,rtol=3e-13)
    assert_allclose(packet['frame_check']['production'],packet['frame_check']['epsilon'],rtol=3e-13)
    for row in packet['signed_forcing']:
        assert_allclose(row['traction'],(row['U_H']-row['F']*row['J1'])/row['J0'],rtol=3e-13,atol=2e-15)
        require(row['mean_production']>=0,'turbulent production lost its signed-traction square')
    require(packet['signed_forcing'][0]['traction']>0 and packet['signed_forcing'][1]['traction']<0
        and packet['signed_forcing'][2]['traction']<0,'pressure-gradient/reversal witnesses absent')
    for row in packet['calibration']:
        assert_allclose(row['h_s'],provisional_roughness(row['H']),rtol=1e-12)
        require(row['z_v']<=0 and 0<row['h_s']<row['H'] and row['gate']>=1 and row['minimum_reference_k']>0,
            'conditional flat calibration/reference gates')
    K=QE*US*US;proof=[]
    for row in packet['curved_columns']:
        c=row['kappa_eff'];a,zv,A,rc,P=wall(K,c);d=NU-a*zv;t=K/QE
        assert_allclose(row['cutoff'],rc,rtol=3e-11,atol=3e-14)
        def Q(r):
            y=d+a*r
            production=t*t/a*((c/a)*(y-d)+(1-c*d/a-NU*c/a)*math.log(y/d)
                +NU*(1-c*d/a)*(1/y-1/d))
            loss=A/rc*(r+c*r*r/2) if r<=rc else A*(1+math.log(r/rc)+c*(r-rc/2))
            return loss-production
        values=[]
        for r in np.linspace(0,HW,5):
            knots=sorted(set([float(r),HW]+([rc] if r<rc<HW else [])))
            values.append(K-sum(quad(lambda z:Q(z)/((1+c*z)*(d+a*z)),l,u,epsabs=2e-14,epsrel=2e-12)[0]
                for l,u in zip(knots[:-1],knots[1:])))
        assert_allclose(row['reference_k'],values,rtol=2e-10,atol=3e-13)
        assert_allclose(row['Q_endpoints'],[0,0],atol=2e-14)
        require(row['minimum_k_101']>0,'curved auxiliary reference nonpositive')
        proof.append(dict(kappa_eff=c,cutoff=rc,gate=row['gate'],K_ground=values[0],
            initial_defect_margin=t*t*a*(-zv)/d**2-A/rc))
    require(min(v['initial_defect_margin'] for v in proof)>0,'rough-column single-crossing positivity certificate fails')
    return dict(maps=out,calibration=packet['calibration'],curved_columns=proof,invalid_geometry_rejections=5,frame=packet['frame_check'],
        positivity_scope='At fixed selected K/h/H, p_ref is unimodal, E-p_ref is positive after cutoff and negative initially; Q_J<=0 by the budget. Linear-in-curvature budgets bracket every intermediate cutoff between endpoint roots. Thus k_col>=K for the certified curvature interval; no evolving-k positivity theorem.')

@functools.lru_cache(maxsize=2048)
def wall(K,c):
    h=provisional_roughness();a=KAPPA*math.sqrt(K/QE);x=h*math.sqrt(K)/(70*NU)
    b=10.7*(1-x**max(.5-.4*x**.7,1-.79*x**(-.28)))
    zv=NU*b/math.sqrt(K);A=(K/QE)**1.5/KAPPA;t=K/QE
    def prod(r):
        nt=a*max(r-zv,0);return (1+c*r)*t*t*nt/(NU+nt)**2
    points=sorted(set([0.,h,HW]+([zv] if 0<zv<HW else [])))
    P=sum(quad(prod,l,r,epsabs=1e-14,epsrel=1e-12)[0] for l,r in zip(points[:-1],points[1:]))
    def budget(rc):return A*(1+math.log(HW/rc)+c*(HW-rc/2))
    require(P>=budget(HW),'metric budget gate failed in independent oracle')
    rc=brentq(lambda z:budget(z)-P,HW*1e-12,HW,xtol=2e-15)
    return a,zv,A,rc,P

def independent_weak(packet):
    alpha=packet['alpha_mass'];expected=np.zeros((2,3,6));parts={}
    def add(name,mode,block,value):
        parts.setdefault(name,np.zeros((3,6)))[mode,block]+=value
    for row in packet['samples']:
        inner=row['inner'];w=row['weight'];u=np.array(row['U']);G=np.array(row['G']);M=G+G.T
        k,e,pi=row['k'],row['epsilon'],row['pi'];gk=np.array(row['grad_k']);ge=np.array(row['grad_epsilon'])
        _,j,curvature,_,normal,_=geometry(row['chart'][0],row['chart'][2],packet['geometry_kind'])
        assert_allclose(w,row['canonical_weight']*np.linalg.det(j),rtol=5e-13,atol=2e-15)
        divergence=np.trace(G);force=np.zeros(3)
        if inner:
            _,_,curv,_,n,_=geometry(row['chart'][0],HW,packet['geometry_kind'])
            t=np.array([n[2],0,-n[0]]);uh=np.array(row['H_U']);gh=np.array(row['H_G'])
            curvH=curv/(1+HW*curv)
            scalar=np.dot(t,np.array(row['H_grad_pi'])-2*np.array(row['H_grad_k'])/3)
            scalar+=2*np.dot(t,uh)*(t@gh@t-curvH*np.dot(n,uh))
            expected_force=scalar*t
            assert_allclose(row['force'],expected_force,rtol=3e-11,atol=3e-12)
            if row['chart'][2]<=provisional_roughness():force=expected_force
            a,zv,A,rc,P=wall(row['K'],row['kappa_eff'])
            assert_allclose(row['cutoff'],rc,rtol=3e-11,atol=3e-14)
            nt=a*max(row['chart'][2]-zv,0);dk=NU+nt
            production=.5*nt*np.sum(M*M)-2*k*divergence/3
            destruction=A/max(row['chart'][2],rc);sk=production-destruction
            T=dk*M-pi*np.eye(3);de=se=0.
        else:
            law=independent(k,e,G);R=law['R'];production=-np.sum(R*G)
            destruction=e;sk=production-e;dk=CSTAR*k*k/e;de=dk/1.3
            se=(1.44*production-1.82*e)*e/k
            T=-R+(2*k/3-pi)*np.eye(3)
        for mode,test in enumerate(row['modes']):
            phi=test['v'];grad=np.array(test['grad_v']);ephi=test['e'];egrad=np.array(test['grad_e']);hist=row['history']
            for b,axis in enumerate((0,2)):
                storage=alpha*phi*(u[axis]-hist[b]);stress=np.dot(grad,T[axis]);adv=-np.dot(grad,u)*u[axis];rough=-phi*force[axis]
                expected[0 if inner else 1,mode,b]+=w*(storage+stress+adv+rough)
                add('time_mass',mode,b,w*storage);add('inner_momentum_stress' if inner else 'outer_momentum_stress',mode,b,w*stress)
                add('conservative_advection',mode,b,w*adv);add('roughness_force',mode,b,w*rough)
            storage=alpha*phi*(k-hist[2]);adv=-np.dot(grad,u)*k;diff=dk*np.dot(grad,gk)
            expected[0 if inner else 1,mode,2]+=w*(storage+adv+diff-phi*sk)
            add('time_mass',mode,2,w*storage);add('conservative_advection',mode,2,w*adv);add('k_diffusion',mode,2,w*diff)
            add('inner_k_production' if inner else 'outer_k_production',mode,2,-w*phi*production)
            add('inner_k_destruction' if inner else 'outer_k_destruction',mode,2,w*phi*destruction)
            if not inner:
                storage=alpha*ephi*(e-hist[3]);adv=-np.dot(egrad,u)*e;diff=de*np.dot(egrad,ge)
                expected[1,mode,3]+=w*(storage+adv+diff-ephi*se)
                add('time_mass',mode,3,w*storage);add('conservative_advection',mode,3,w*adv);add('epsilon_diffusion',mode,3,w*diff)
                add('epsilon_production',mode,3,-w*ephi*1.44*production*e/k)
                add('epsilon_destruction',mode,3,w*ephi*1.82*e*e/k)
            pb=4 if inner else 5
            expected[0 if inner else 1,mode,pb]+=w*test['p']*(divergence+row['lambda']/packet['volume'])
            add('continuity',mode,pb,w*test['p']*divergence)
            add('gauge',mode,pb,w*test['p']*row['lambda']/packet['volume'])
    errors={}
    component_names=('inner_momentum_stress','outer_momentum_stress','conservative_advection',
        'k_diffusion','epsilon_diffusion','inner_k_production','inner_k_destruction','outer_k_production',
        'outer_k_destruction','epsilon_production','epsilon_destruction','epsilon_yap','prescribed_scalar_source',
        'roughness_force','prescribed_body','continuity','gauge','time_mass')
    for name,native in packet['native_terms'].items():
        if name not in component_names:continue
        oracle=parts.get(name,np.zeros((3,6)));native=np.array(native)
        # Different independent quadrature; absolute floors cover nearly zero
        # mode loads. Production/source and pressure errors remain separate.
        assert_allclose(native,oracle,rtol=2e-8,atol=2e-9,err_msg=name)
        errors[name]=float(np.max(abs(native-oracle)))
    actual=np.array([packet['native_inner'],packet['native_outer']])
    assert_allclose(actual,expected,rtol=2e-8,atol=4e-9,err_msg='regional physical weak loads')
    return dict(samples=len(packet['samples']),component_max_errors=errors,
        regional_max_error=float(np.max(abs(actual-expected))),regional_native=actual.tolist(),regional_independent=expected.tolist(),
        scope='three declared test modes, actual represented fields; raw regional/exterior reaction loads, not a nonlinear Jacobian')

def check_native(report,folder):
    attempts=report['pimple_attempts'];require(report['steps'] in (1,3),'actual accepted terrain steps missing')
    require(len(attempts)==report['steps'] and len(report['callbacks'])==report['steps'],'native callback/attempt counts')
    for attempt in attempts:
        for key in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):
            require(attempt[key]==1,'one-sweep solve/bound count: '+key)
        require(attempt['bound_activations']==0 and attempt['additional_governing_equation_solves']==0,'correction/hidden solve')
        snap=json.loads(Path(attempt['k_snapshot']).read_text())
        require(np.asarray(snap['raw'],dtype='<f8').tobytes()==np.asarray(snap['corrected'],dtype='<f8').tobytes(),'raw/corrected k state changed')
    audits=[report['startup_audit']]+[v['audit'] for v in report['callbacks']]
    for audit in audits:require(audit['wj_hard_admissible']==[1],'native hard model-domain predicate')
    interfaces=[report['startup_interface']]+[e['interface'] for e in report['events']]
    K=[];cut=[];ep=[]
    for packet in interfaces:
        require(packet['shared_value_test_row_defect']<2e-12,'common interface basis mismatch')
        for row in packet['samples']:
            require(row['inner_geometry']!=row['outer_geometry'],'one-sided element ownership lost')
            assert_allclose(row['mass_inner']+row['mass_outer'],0,atol=2e-12)
            assert_allclose(row['epsilon'],row['epsilon_target'],rtol=0,atol=1e-11)
            require(row['cutoff_gate']>=1 and row['metric_stretch']>.1,'cutoff/geometry gate')
        K.append([v['K'] for v in packet['samples']]);cut.append([v['cutoff'] for v in packet['samples']]);ep.append([v['epsilon_target'] for v in packet['samples']])
    require(np.max(abs(np.array(K[-1])-K[0]))>1e-13,'endpoint K did not evolve')
    require(np.max(abs(np.array(cut[-1])-cut[0]))>1e-15,'metric cutoff did not respond to current K')
    require(np.max(abs(np.array(ep[-1])-ep[0]))>1e-14,'dependent epsilon trace did not evolve')
    initial=independent_weak(json.loads((folder/'report.json.startup-weak.json').read_text()))
    final=independent_weak(json.loads(Path(report['final_weak']).read_text()))
    return dict(status='PASS',steps=report['steps'],state_size=report['state_size'],time=report['time'],dt=report['dt'],
        maximum_K_change=float(np.max(abs(np.array(K[-1])-K[0]))),maximum_cutoff_change=float(np.max(abs(np.array(cut[-1])-cut[0]))),
        maximum_epsilon_trace_change=float(np.max(abs(np.array(ep[-1])-ep[0]))),
        startup_independent_loads=initial,final_independent_loads=final,
        signed_wall_traction_range=[min(v for event in report['events'] for v in event['physical_samples']['wall_stress']),
            max(v for event in report['events'] for v in event['physical_samples']['wall_stress'])],
        corrections=0,governing_solves_each=report['steps'],logical_bounds=report['steps'])

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--geometry-probe',type=Path);parser.add_argument('--native-probe',type=Path)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--kind',choices=('hill','zero','plane'),default='hill')
    parser.add_argument('--mode',choices=('geometry','native','guards'),required=True);args=parser.parse_args()
    folder=args.output/('run-'+uuid.uuid4().hex[:12]);folder.mkdir(parents=True,exist_ok=False)
    save(folder/'design.json',dict(mode=args.mode,kind=args.kind,xknots=X,zknots=Z,
        nu=NU,ustar=US,matching_height=HW,roughness=provisional_roughness(),dt=1e-6,
        steps=1 if args.kind=='plane' else 3,driver_sha256=sha(__file__)))
    try:
        if args.mode=='geometry':
            receipt=run([args.geometry_probe.resolve(),folder.resolve()/'geometry.json'],folder)
            result=dict(status='PASS',geometry=check_geometry(json.loads((folder/'geometry.json').read_text())),command=receipt)
        elif args.mode=='native':
            xml=case_xml(folder/'input',args.kind,1 if args.kind=='plane' else 3)
            receipt=run([args.native_probe.resolve(),xml.resolve(),folder.resolve()/'report.json'],folder)
            result=check_native(json.loads((folder/'report.json').read_text()),folder);result['command']=receipt
            if args.kind=='plane':require(result['signed_wall_traction_range'][0]<0,'reversed inclined native traction was lost')
        else:
            import xml.etree.ElementTree as ET
            from cls_epawf_m3_cases import set_info,set_parameter,write
            cases=[('missing_route','RANSWallNativeGeometry','FlatV1','composite native evolution requires explicit BudgetTerrainV1'),
                ('EP','RANSWallFormulation','WJCalibratedEPAWF','BudgetTerrainV1 requires WJCalibratedBudgetAWF'),
                ('capacity','RANSWallDenseEvolutionCapacity','4331','dense_capacity_route'),
                ('verify','RANSWallPimpleVerify','On','additional governing verification solves'),
                ('fold','RANSWallHeight',2.,'composite inner normal reach/minimum stretch')]
            results=[]
            for name,key,value,reason in cases:
                sub=folder/name;xml=case_xml(sub/'input','hill',0);root=ET.parse(xml).getroot()
                (set_parameter if key=='RANSWallHeight' else set_info)(root,key,value)
                if name=='fold':
                    set_parameter(root,'RANSWallSampleHeight',value);set_parameter(root,'RANSWallTailStart',value)
                write(root,xml)
                results.append(run([args.native_probe.resolve(),xml.resolve(),sub.resolve()/'report.json'],sub,expected_failure=reason))
            result=dict(status='PASS',reason_specific_rejections=len(results),commands=results)
        save(folder/'result.json',result);print(json.dumps(dict(result=str(folder/'result.json'),status='PASS')))
    except BaseException as error:
        save(folder/'failure.json',dict(status='OBSTRUCTED',error=repr(error)))
        raise
if __name__=='__main__':main()
