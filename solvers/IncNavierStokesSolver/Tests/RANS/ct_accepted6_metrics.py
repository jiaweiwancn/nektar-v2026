# SPDX-License-Identifier: MIT
"""Physical CT6 metrics on unchanged native rows; absolute augmented pressure."""
import math,os
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import legfit,legval
from ct_accepted_feasibility import parameters,require
from ct_accepted_open import row_value
from validate_a93_a95_e94e95 import oracle


def endpoint_trace(packet,state):
    o=packet['observation'];off=o['offsets'];degree=o['elements'][0]['degree']
    bounds=sorted({x for e in o['elements'] for x in (e['lower'][0],e['upper'][0])})
    samples={}
    for point in o['volume']:
        q=point['coordinate'][0];v=row_value(point['CT6_prescribed_data']['H_velocity_k_row'],state,off[2])[0]
        if q in samples:require(abs(v-samples[q])<1e-12,'inconsistent_actual_endpoint_trace')
        samples[q]=v
    polynomials=[];defect=0.
    for a,b in zip(bounds[:-1],bounds[1:],strict=True):
        data=sorted((q,v) for q,v in samples.items() if a<q<b)
        require(len(data)>=degree+1,'insufficient_native_endpoint_basis_samples')
        x=np.array([2*(q-a)/(b-a)-1 for q,v in data]);y=np.array([v for q,v in data])
        coefficients=legfit(x,y,degree);defect=max(defect,float(max(abs(legval(x,coefficients)-y))))
        polynomials.append((a,b,coefficients))
    require(defect<=2e-12,'endpoint_polynomial_restriction_not_exact')
    def evaluate(x):
        a,b,c=next(p for p in polynomials if p[0]<=x<=p[1])
        return float(legval(2*(x-a)/(b-a)-1,c))
    return evaluate,defect


def wall_nut(model,K,r):
    p=parameters(model);nu=1.5e-5;hr=.05;x=hr*math.sqrt(K)/(70*nu)
    B=10.7*(1-x**max(.5-.4*x**.7,1-.79*x**(-.28)))
    return p['alpha']*math.sqrt(K)*max(0.,r-nu*B/math.sqrt(K))


def metrics(packet,model,initial_state=None):
    o=packet['observation'];s=np.asarray(o['state']);off=o['offsets'];p=parameters(model)
    K0=p['q']*.16;scales=np.array([.2,.2,K0,.4**3/(p['kappa']*.5),.1])
    squared=np.zeros(5);gradient=np.zeros(5);volumes=np.zeros(5);drift=np.zeros(5)
    for point in o['volume']:
        region=point['region'];actual=[];derivatives=[]
        for b,row in ((0,point['velocity_k']),(1,point['velocity_k']),(2,point['velocity_k']),
                      (3,point['epsilon']),(4+region,point['pressure'])):
            if row is None:actual.append(0.);derivatives.append(np.zeros(3))
            else:
                value,d=row_value(row,s,off[b]);actual.append(value);derivatives.append(d)
        exact=point['CT6_prescribed_data'];a=np.array(actual)-np.asarray(exact['continuous_values'])
        dg=np.array(derivatives)-np.asarray(exact['continuous_gradients']);w=point['weight']
        use=np.ones(5);use[3]=region
        squared+=w*use*a*a;gradient+=w*use*np.sum(dg*dg,axis=1);volumes+=w*use
        if initial_state is not None:
            base=[]
            for b,row in ((0,point['velocity_k']),(1,point['velocity_k']),(2,point['velocity_k']),
                          (3,point['epsilon']),(4+region,point['pressure'])):
                base.append(0. if row is None else row_value(row,initial_state,off[b])[0])
            drift+=w*use*(np.array(actual)-base)**2
    rms=np.sqrt(squared/volumes)/scales;grms=np.sqrt(gradient/volumes)/scales
    fields={name:dict(scaled_weighted_RMS=float(rms[j]),true_L2=float(math.sqrt(squared[j])),
            scaled_gradient_RMS=float(grms[j]),volume=float(volumes[j]),scale=float(scales[j]),
            scaled_drift_from_actual_initial=float(math.sqrt(drift[j]/volumes[j])/scales[j]) if initial_state is not None else None)
            for j,name in enumerate(('U','w','k','outer_epsilon','Pi'))}
    trace,trace_defect=endpoint_trace(packet,s);wall=0.;wall_measure=0.;wall_mean=np.zeros(3);wall_exact=np.zeros(3);interface={};outlet=0.;outlet_measure=0.
    interface_error=[0.,0.];interface_measure=[0.,0.]
    for face in o['faces']:
        region=face['region'];q=face['coordinate'];v=face['velocity_k'];U=np.zeros(3);G=np.zeros((3,3))
        for b,axis in enumerate((0,2)):U[axis],G[axis]=row_value(v,s,off[b])
        k,gk=row_value(v,s,off[2]);pi,_=row_value(face['pressure'],s,off[4+region])
        if region:
            e,_=row_value(face['epsilon'],s,off[3]);law=oracle(model,k,e,G);D=law['nut'];T=-law['R']+np.eye(3)*(2*k/3-pi)
        else:
            D=1.5e-5+wall_nut(model,trace(q[0]),q[2]);T=D*(G+G.T)-np.eye(3)*pi
        area=np.asarray(face['area']);norm=np.linalg.norm(area);normal=area/norm;weight=face['weight']*norm
        traction=T@normal;truth=np.asarray(face['CT6_continuous']['T'])@normal
        if q[2]==0.:
            wall+=weight*float((traction-truth)@(traction-truth));wall_measure+=weight
            wall_mean+=weight*traction;wall_exact+=weight*truth
        if face['interface']:
            key=tuple(q);interface.setdefault(key,[]).append(dict(traction=traction,flux=float((U*k-D*gk)@normal),weight=weight))
            interface_error[region]+=weight*float((traction-truth)@(traction-truth));interface_measure[region]+=weight
        elif q[0]==max(e['upper'][0] for e in o['elements']):
            outlet+=weight*float((traction-truth)@(traction-truth));outlet_measure+=weight
    jump=0.;flux=0.;iw=0.
    for key,rows in interface.items():
        require(len(rows)==2,'two_sided_interface_pair_missing')
        a,b=rows;require(abs(a['weight']-b['weight'])<1e-13,'interface_measures_differ')
        d=a['traction']+b['traction'];jump+=a['weight']*float(d@d);flux+=a['weight']*(a['flux']+b['flux'])**2;iw+=a['weight']
    traction=dict(ground_scaled_RMS=math.sqrt(wall/wall_measure)/.16,
        ground_signed_mean=(wall_mean/wall_measure).tolist(),ground_reference_mean=(wall_exact/wall_measure).tolist(),
        full_interface_jump_scaled_RMS=math.sqrt(jump/iw)/.16,
        inner_interface_reference_error_scaled_RMS=math.sqrt(interface_error[0]/interface_measure[0])/.16,
        outer_interface_reference_error_scaled_RMS=math.sqrt(interface_error[1]/interface_measure[1])/.16,
        conservative_k_interface_jump_scaled_RMS=math.sqrt(flux/iw)/(K0*.2),
        outlet_reconstruction_scaled_RMS=math.sqrt(outlet/outlet_measure)/.16,
        current_endpoint_trace_restriction_defect=trace_defect)
    limits=np.array([.005,.005,.005,.005,.01]);passed=bool(np.all(rms<=limits)) and max(traction[k] for k in
        ('ground_scaled_RMS','full_interface_jump_scaled_RMS','inner_interface_reference_error_scaled_RMS','outer_interface_reference_error_scaled_RMS'))<=.005
    return dict(schema='CT6ManufacturedPhysicalMetricsV1',status='pass' if passed else 'fail',
        fields=fields,traction=traction,pressure='absolute augmented Pi, no fitted gauge subtraction',
        field_limits=dict(zip(('U','w','k','outer_epsilon','Pi'),limits.tolist(),strict=True)),
        traction_limit=.005,scope='current FE coefficients and actual endpoint K; projection/error, drift and reconstructed physical traction remain separate')
