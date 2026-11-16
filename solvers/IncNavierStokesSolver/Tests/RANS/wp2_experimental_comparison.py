# SPDX-License-Identifier: MIT
"""Retained-state flat-profile discrepancy, with all frozen points and weights.

Only actual affine polynomial primitive fields/gradients are interpolated.
Selected stresses are recomputed; no stress-profile interpolation is used.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss,legvander
from cao_tamura_h2g_constitutive import parts
from run_cls_epawf_picard import require,sha,state_hash

US=.07426793792419734


def primitive_at(report,x,z):
    a=report['wp2_common_samples'];require(a['order']==8,'native GL8 primitive data required')
    zk=np.asarray(a['physical_z_partition']);xk=np.array([-.2,-.1,-.05,0,.05,.1,.6]);nz=len(zk)-1
    require(xk[0]<=x<=xk[-1] and zk[0]<=z<=zk[-1],'physical observation outside the flat domain')
    ix=min(len(xk)-2,max(0,int(np.searchsorted(xk,x,side='right')-1)))
    iz=min(nz-1,max(0,int(np.searchsorted(zk,z,side='right')-1)))
    old,w=leggauss(8);V=legvander(old,3);projection=(np.arange(4)*2+1)[:,None]/2*V.T*w[None,:]
    px=legvander(2*(x-xk[ix])/(xk[ix+1]-xk[ix])-1,3).ravel()@projection
    pz=legvander(2*(z-zk[iz])/(zk[iz+1]-zk[iz])-1,3).ravel()@projection
    v=np.asarray(a['values']).reshape(6,8,nz,8,14)[ix,:,iz,:,:]
    return np.einsum('i,ijf,j->f',px,v[:,:,[0,1,2,3,4,10,11,12,13]],pz)


def compare(report,reference_file,stations=(-.15,0.,.6)):
    from cao_tamura_h2d_state import load
    reference=load(reference_file);require(reference['dataset']=='cao-tamura-2006-v1.1','frozen data version changed')
    law='calibrated_standard' if report['model_specification']['name']=='StandardKEpsilon' else 'regularized_CLS'
    kap=report['parameters'][7];Uinf=US/.051;outputs={}
    for x in stations:
        series={}
        for name,s in reference['series'].items():
            rows=[]
            for record in s['records']:
                z=record['z'];require(z>report['parameters'][1],'this upstream dataset must use the actual outer law')
                v=primitive_at(report,x,z);g=np.zeros((1,3,3));g[0,0,0],g[0,0,2],g[0,2,0],g[0,2,2]=v[5:9]
                c=parts(g,np.array([v[2]]),np.array([v[3]]),np.array([False]),np.array([0.]),np.array([[x,0,z]]),US,kap,law)
                R=c['tensor'][0]
                prediction=v[0]/Uinf if name=='f3a_R' else -R[0,2]/Uinf**2 if name=='f3b_R' else np.sqrt(R[0,0])/US
                rows.append(dict(z=z,weight=record['weight'],observed=record['observed'],extraction_bound=record['extraction_bound'],
                                 prediction=float(prediction),residual=float(prediction-record['observed'])))
            weight=sum(q['weight'] for q in rows)
            series[name]=dict(count=len(rows),records=rows,weighted_RMSE=float(np.sqrt(sum(q['weight']*q['residual']**2 for q in rows)/weight)))
            if name=='f3a_R':
                band=[q for q in rows if .01<=q['z']<=.051]
                series[name]['unchanged_seven_point_band_RMSE']=float(np.sqrt(sum(q['weight']*q['residual']**2 for q in band)/sum(q['weight'] for q in band)))
        outputs[str(x)]=series
    return dict(version='WP2PreservedExperimentalDiscrepancyV1',model=report['model_specification'],actual_time=report['time'],
        state_sha256=state_hash(report['state']),frozen_reference_file=str(reference_file),frozen_reference_file_sha256=sha(reference_file),
        snapshot_sha256=reference['snapshot_sha256'],stations=outputs,Uinf=Uinf,
        point_lookup='physical flat coordinates; exact local cubic reconstruction from native primitive samples; right native limit at x=0, interior limit at outlet; stresses recomputed',
        interpretation='numerical flat stations are declared comparison locations, not an asserted identification of an unresolved experimental streamwise coordinate; no new fit',
        held_out='Fig3d and all hill observations remain unscored; full-height surface-layer numerical preservation is separate from finite-boundary-layer experimental agreement')
