#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Initial target inner-shear Stokes/physical-face pressure-image diagnostic."""
import argparse,os,resource,time
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted12_pressure import Frozen,maximum
from ct_accepted9_pressure_run import sealed
from ct_accepted_feasibility import parameters,write_json,require
from ct_accepted_open import row_value
from rans_request_manifest import RequestManifest,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct10','ct11','pressure','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));resource.setrlimit(resource.RLIMIT_CORE,(0,0));start=time.process_time()
    primary=sealed(a.pressure);target=sealed(a.ct11/'target-v1');c=target['case'];ini=n.read(a.ct10/'target-preflight-v1/target/readiness.json');off=ini['observation']['offsets'];nm=off[2];state=np.asarray(ini['observation']['state'])
    row=primary['members']['0.001']['steps'][0];require(sha(Path(row['snapshot']['path']))==row['snapshot']['sha256'],'changed initial captured operator')
    s=n.read(row['snapshot']['path']);require(s['input']==state.tolist(),'foreign initial frozen state')
    manifest=RequestManifest(a.output,dict(schema='CT12InitialShearFaceV1',source=sha(Path(__file__)),contract=sha(a.contract),pressure=sha(a.pressure/'results.json'),target=sha(a.ct11/'target-v1/results.json'),initial_sha256=exact_vector_hash(state),snapshot=row['snapshot'],scope='read-only initial inner-shear face/volume and pressure-image identity; no native process, trajectory or governing correction'))
    mom=ini['ct10_momentum'];delta=sum((np.asarray(mom['actual'][k])-np.asarray(mom['analytic'][k]) for k in ('inner/volume/linear_shear','inner/volume/molecular_shear')),np.zeros(nm))
    ends={}
    for point,side,q in zip(ini['observation']['grid'],ini['observation']['grid_sides'],ini['test_basis'],strict=True):
        if point[2]==c['H'] and side==1:ends[round(point[0],14)]=q['velocity_k']
    alpha=parameters('A95')['alpha'];face=np.zeros(nm);outlet=np.zeros(nm);samples=[]
    for f in ini['observation']['faces']:
        if f['interface'] or f['region'] or not f['area'][0]:continue
        x,_,z=f['coordinate'];K=row_value(ends[round(x,14)],state,off[2])[0];require(K>0,'invalid actual initial endpoint K')
        rough=c['roughness_input']*np.sqrt(K)/(70*c['nu']);B=10.7*(1-rough**max(.5-.4*rough**.7,1-.79*rough**(-.28)))
        D=c['nu']+alpha*np.sqrt(K)*max(0.,z-c['nu']*B/np.sqrt(K))
        G=np.zeros((3,3))
        for block,axis in ((0,0),(1,2)):G[axis]=row_value(f['velocity_k'],state,off[block])[1]
        tau=D*(G[0,2]+G[2,0]);q=f['velocity_k'];ids=np.asarray(q['ids'],int)+off[1]
        load=f['weight']*f['area'][0]*np.asarray(q['phi'])*(tau-c['ustar']**2)
        np.add.at(face,ids,load)
        if f['area'][0]>0:np.add.at(outlet,ids,load)
        samples.append(dict(coordinate=f['coordinate'],geometry=q['geometry'],owner_element=q['element'],area=f['area'],weight=f['weight'],actual_tau=tau,reference_tau=c['ustar']**2,actual_K_H=K,raw_gradient=G.tolist()))
    normal=delta.copy();normal[:off[1]]=0.;streamwise=delta-normal
    free=np.asarray(mom['free_rows'],int);fixed=np.asarray(mom['fixed_rows'],int)
    defects=dict(all_row=maximum(face-normal),free_row=maximum((face-normal)[free]),prescribed_row=maximum((face-normal)[fixed]))
    write_json(a.output/'partial-stokes.json',dict(defects=defects,maximum_face_load=maximum(face),maximum_native_normal_shear_difference=maximum(normal),scope='partial read-only witness,not qualified'))
    require(max(defects.values())<=2e-10,'independent shear physical-face/native-volume Stokes identity')
    fr=Frozen(s,ini,c,n.inventory(c)[2],len(n.inventory(c)[2]));pn=fr.image(normal,'initial-inner-shear/normal');pu=fr.image(streamwise,'initial-inner-shear/streamwise');pf=fr.image(face,'initial-inner-shear/physical-faces')
    repeat=fr.measure(fr.image(normal,'initial-inner-shear/normal-tighter',2e-14)-pn);ferror=fr.measure(pf-pn)
    require(max(repeat['maximum_scaled_value'],ferror['maximum_scaled_value'])<1e-8,'shear face pressure-image uncertainty')
    native=np.asarray(s['pressure_increment']);den=fr.inner(native,native)
    wrong=fr.image(-2*outlet,'negative/reversed-outlet-orientation');signal=fr.measure(wrong)
    require(signal['maximum_scaled_value']>1e-8,'wrong outlet orientation has no pressure signal')
    out=dict(schema='CT12InitialShearFaceV1',status='pass',native_stokes_defects=defects,
        interpretation='initial radial reference is tangentially constant; normal-momentum inner-shear weak load equals actual-minus-analytical inlet/outlet shear flux. Prescribed inlet reactions retained separately; free inlet tests vanish. This does not add a governing boundary load.',
        normal=dict(image=fr.measure(pn),signed_alignment_with_native_increment=fr.inner(pn,native)/den),streamwise=dict(image=fr.measure(pu),signed_alignment_with_native_increment=fr.inner(pu,native)/den),
        physical_face_image=fr.measure(pf),face_image_error=ferror,tighter_repeat=repeat,
        negative_reversed_outlet=dict(status='pass',pressure_defect=signal,raw_load_defect=maximum(2*outlet),intended='wrong orientation gives resolved pressure-image error; raw absolute row bound alone need not resolve a small load'),
        solves=fr.receipts,CPU=time.process_time()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        sample_count=len(samples),sample_witness=max(samples,key=lambda q:abs(q['actual_tau']-q['reference_tau'])),physical_horizon=0.,governing_change=False)
    write_json(a.output/'physical-face-samples.json',samples);write_json(a.output/'results.json',out);manifest.seal(a.output,['results.json','physical-face-samples.json']);return 0
if __name__=='__main__':raise SystemExit(main())
