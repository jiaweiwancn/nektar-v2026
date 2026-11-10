# SPDX-License-Identifier: MIT
"""Independent forward tensors/production compared with actual native samples."""
import numpy as np
from validate_a93_a95_e94e95 import oracle
from rans_preservation_oracle import Reference
from rans_duration_contract import require


def analyze(packet, config, weights):
    observed=packet['observation'];native=packet['native_constitutive']
    require(native is not None and native['descriptor']==config['descriptor'],'native_constitutive_unavailable_or_foreign')
    ref=Reference(config['model'],config['policy'],config['wall']);defects={};values=[];w=[];margins=[]
    for sample in native['samples']:
        j=sample['sample_index'];A=np.array(sample['raw_gradient']);k=sample['k'];e=sample['epsilon'];expect=oracle(config['model'],k,e,A);R=np.array(sample['R']);Q=np.array(sample['Q']);point=observed['grid'][j]
        require(sample['coordinate']==point and sample['side']==observed['grid_sides'][j],'constitutive_sample_identity_mismatch')
        entries=dict(R=np.max(np.abs(R-expect['R'])),Q=np.max(np.abs(Q-expect['Q'])),P=abs(sample['P']-expect['production']),epsilon_production=abs(sample['epsilon_production']-expect['eps_production']),epsilon_destruction=abs(sample['epsilon_destruction']-expect['eps_destruction']),
            normal_traction=abs(observed['grid_traction'][j][2]-(-R[2,2]+2*k/3-sample['Pi'])),shear_traction=abs(observed['grid_traction'][j][0]+R[0,2]))
        sd=(A+A.T)/2-np.trace(A)*np.eye(3)/3
        isotropic=-2*k*np.trace(A)/3;nonlinear=-np.sum(Q*A);linear=2*sample['nu_t']*np.sum(sd*sd)
        entries['signed_full_production']=abs(sample['P']-(isotropic+nonlinear+linear))
        if sample['a95_radicand'] is not None:
            rad=1-9*expect['cmu']**2*np.sum(sd*sd)*(k/e)**2;entries['A95_radicand']=abs(rad-sample['a95_radicand']);margins.append(dict(value=sample['a95_radicand'],coordinate=point,side=sample['side']))
        for key,value in entries.items():defects[key]=max(defects.get(key,0.),float(value))
        eig=np.linalg.eigvalsh(R);defects['covariance_eigenvalues']=max(defects.get('covariance_eigenvalues',0.),float(np.max(np.abs(eig-np.array(sample['covariance_eigenvalues'])))));analytic=ref.at(point[2],1);R0=np.array(analytic['R'])
        values.append([R[0,0],R[1,1],R[2,2],R[0,2],Q[0,2],np.linalg.norm(Q),sample['P'],nonlinear,isotropic,np.trace(A),eig[0],R[2,2]-R0[2,2]])
        w.append(weights[j])
    require(values and all(np.isfinite(v) and v<=5e-11 for v in defects.values()),'independent_constitutive_mismatch')
    a=np.array(values);w=np.array(w);labels=['Rxx','Ryy','Rzz','Rxz','Qxz','Q_frobenius','P','minus_Q_contract_raw_A','isotropic_production','raw_divergence','covariance_minimum','Rzz_minus_reference']
    stats={name:dict(minimum=float(a[:,j].min()),maximum=float(a[:,j].max()),weighted_RMS=float(np.sqrt(np.sum(w*a[:,j]**2)/w.sum()))) for j,name in enumerate(labels)}
    return dict(status='pass',configuration=config,time=observed['time'],samples=len(a),defects=defects,absolute_comparison_limit=5e-11,statistics=stats,minimum_A95_radicand=min(margins,key=lambda x:x['value']) if margins else None,
        neutral_identity='analytical simple-shear Qxz=0; Rxz=-u_star^2. Actual represented/evolved Qxz is separately observed and not forced to zero.',
        scope='independent A-polynomial/eigenvalue oracle; actual native R,Q,P, signed production including raw-divergence isotropic work, full normal/shear traction')
