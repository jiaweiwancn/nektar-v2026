# SPDX-License-Identifier: MIT
"""Independent actual-state conservative loads on the H2 flat composite map."""
import numpy as np
from numpy.testing import assert_allclose
from wj_h2_profiles import wall,constitutive,CD,NU,D0,TOP
from validate_wj_h1 import require
def independent_weak(packet,report,enforce=True):
    H=report["parameters"][1];h=report["parameters"][2]
    alpha=packet['alpha_mass'];expected=np.zeros((2,3,6));parts={}
    def add(name,mode,block,value):
        parts.setdefault(name,np.zeros((3,6)))[mode,block]+=value
    for row in packet['samples']:
        inner=row['inner'];w=row['weight'];u=np.array(row['U']);G=np.array(row['G']);M=G+G.T
        k,e,pi=row['k'],row['epsilon'],row['pi'];gk=np.array(row['grad_k']);ge=np.array(row['grad_epsilon'])
        determinant=1. if inner else (TOP-D0-H)/(TOP-H)
        assert_allclose(w,row['canonical_weight']*determinant,rtol=5e-13,atol=2e-15)
        divergence=np.trace(G);force=np.zeros(3)
        if inner:
            uh=np.array(row['H_U']);gh=np.array(row['H_G'])
            expected_force=np.array([row['H_grad_pi'][0]-2*row['H_grad_k'][0]/3+2*uh[0]*gh[0,0],0,0])
            assert_allclose(row['force'],expected_force,rtol=3e-11,atol=3e-12)
            if row['chart'][2]<=h:force=expected_force
            selected=wall(row['K'],H,h);a,zv,A,rc=selected['a'],selected['zv'],selected['A'],selected['rc']
            assert_allclose(row['cutoff'],rc,rtol=3e-11,atol=3e-14)
            nt=a*max(row['chart'][2]-zv,0);dk=NU+nt
            production=.5*nt*np.sum(M*M)-2*k*divergence/3
            destruction=A/max(row['chart'][2],rc);sk=production-destruction
            T=dk*M-pi*np.eye(3);de=se=0.
        else:
            R,_,_=constitutive(k,e,G);production=-np.sum(R*G)
            destruction=e;sk=production-e;dk=CD*k*k/e;de=dk/1.3
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
    exterior=np.zeros_like(expected)
    support=min(float(v['upper'][2]) for v in report['fields'][0]['mesh'] if float(v['lower'][2])==H)
    for row in report['open_boundary']['boundary_rows']:
        x,_,rho=row['chart'];w=row['weight'];inner=row['inner'];side=0 if inner else 1
        hat=0 if rho<H/2 or rho>support else (rho-H/2)/(H/2) if rho<H else (support-rho)/(support-H)
        # Recompute each physical flux, without using the returned native flux.
        n=np.array(row['normal']);u=np.array(row['U']);un=u@n;k=row['k'];e=row['epsilon'];G=np.array(row['G'])
        if inner:
            v=wall(row['K'],H,h);nt=v['a']*max(rho-v['zv'],0);dk=NU+nt;de=0
        else:dk=CD*k*k/e;de=dk/1.3
        jk=-dk*np.array(row['grad_k'])@n if row['role']=='inlet' else 0
        je=-de*np.array(row['grad_epsilon'])@n if row['role']=='inlet' and not inner else 0
        if row['role']=='top' and 'ReferenceTractionV1' in report['signature']:
            K=report['parameters'][5];je=CD*K*K/1.3/(TOP-D0)
        traction=np.array(row['applied_traction'])
        flux=[u[0]*un-traction[0],u[2]*un-traction[2],k*un+jk,e*un+je]
        for mode,phi in enumerate((1.,x/5,hat)):
            for block in range(3 if inner else 4):exterior[side,mode,block]+=w*phi*flux[block]
    expected+=exterior
    parts['exterior_boundary']=exterior.sum(axis=0)
    errors={};conditions={}
    component_names=('inner_momentum_stress','outer_momentum_stress','conservative_advection',
        'k_diffusion','epsilon_diffusion','inner_k_production','inner_k_destruction','outer_k_production',
        'outer_k_destruction','epsilon_production','epsilon_destruction','epsilon_yap','prescribed_scalar_source',
        'roughness_force','prescribed_body','continuity','gauge','time_mass','exterior_boundary')
    for name,native in packet['native_terms'].items():
        if name not in component_names:continue
        oracle=parts.get(name,np.zeros((3,6)));native=np.array(native)
        # Different independent quadrature; absolute floors cover nearly zero
        # mode loads. Production/source and pressure errors remain separate.
        conditions[name]=bool(np.all(np.isfinite(native)) and np.all(np.isfinite(oracle)) and
            np.all(abs(native-oracle)<=2e-9+2e-8*abs(oracle)))
        if enforce:assert_allclose(native,oracle,rtol=2e-8,atol=2e-9,err_msg=name)
        errors[name]=float(np.max(abs(native-oracle)))
    actual=np.array([packet['native_inner'],packet['native_outer']])
    conditions['regional_sum']=bool(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)) and
        np.all(abs(actual-expected)<=4e-9+2e-8*abs(expected)))
    if enforce:assert_allclose(actual,expected,rtol=2e-8,atol=4e-9,err_msg='regional physical weak loads')
    return dict(status='PASS' if all(conditions.values()) else 'UNRESOLVED',conditions=conditions,
        outer_quadrature_order=packet.get('outer_quadrature_order',24),samples=len(packet['samples']),component_max_errors=errors,
        independent_components={k:v.tolist() for k,v in parts.items()},native_components=packet['native_terms'],
        regional_max_error=float(np.max(abs(actual-expected))),regional_native=actual.tolist(),regional_independent=expected.tolist(),
        scope='three declared test modes, actual represented fields; raw regional/exterior reaction loads, not a nonlinear Jacobian')
