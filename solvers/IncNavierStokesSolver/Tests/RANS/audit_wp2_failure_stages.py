#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Zero-solve forensic replay of saved physical stages at the WP2 rejection."""
import argparse,time,resource
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss,legvander
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require,sha,state_hash
from continue_cls_epawf_pimple import atomic
from wp2_all_state import Rows,wall
from wp2_common_quadrature import resample


def main():
    parser=argparse.ArgumentParser()
    for name in ['root','output']:parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();r=args.root;started=time.monotonic()
    path=r/'preservation-recovery-standard/step390/case.evolution.json';report=load(path)
    packet_path=Path(report['attempt_packet_status']['last_written_packet']);packet=load(packet_path)
    require(sha(packet_path)==report['attempt_packet_status']['sha256'],'packet changed')
    require(sha(Path(packet['matrix_file']))==packet['matrix_sha256'],'frozen k matrix changed')
    attempt=report['pimple_attempts'][0];rows=Rows(r/'geometry-rows-selected/result.json')
    audit=report['physical_range_events'][0]['audit'];point=audit['accepted_awf_inner_minimum_point'][:3]
    geometry=int(audit['native_inner_minimum_geometry_element_side_block'][0])
    bounds=next(e for e in report['fields'][0]['mesh'] if e['geometry']==geometry)
    knots=np.asarray(rows.meta['physical_z_partition']);j=int(np.searchsorted(knots,point[2],side='right')-1)
    mask=(rows.a['geometry'][:,0]==geometry)&(rows.physical[:,2]>knots[j])&(rows.physical[:,2]<knots[j+1])
    require(sum(mask)==64,'native one-sided panel identity')
    z,w=leggauss(8);V=legvander(z,3);P=(2*np.arange(4)+1)[:,None]/2*V.T*w
    xhat=2*(point[0]-bounds['lower'][0])/(bounds['upper'][0]-bounds['lower'][0])-1
    zhat=2*(point[2]-knots[j])/(knots[j+1]-knots[j])-1
    ix=legvander([xhat],3)[0]@P;iz=legvander([zhat],3)[0]@P
    def at(data):return np.einsum('i,ijf,j->f',ix,np.asarray(data)[mask].reshape(8,8,-1),iz)
    names={'committed':'incoming_sha256','post_flow':'post_flow_sha256','pre_k':'pre_k_sha256','raw':'raw_k_sha256','corrected':'corrected_k_sha256'}
    result=[];nu,H,hs=report['parameters'][:3];q,kappa=report['parameters'][6:8];K0=report['parameters'][5]
    for name,identity in names.items():
        state=np.asarray(packet[name]);require(state_hash(state)==attempt[identity],'stage hash differs: '+name)
        v=rows.tuple(state,report);values=at(v['values']);G=at(v['G'].reshape(-1,9)).reshape(3,3)
        endpoint=np.repeat(v['endpoint_K'],len(v['values'])//len(v['endpoint_K']));K=float(at(endpoint)[0])
        profile=wall(np.asarray([K]),np.asarray([point[2]]),nu,H,hs,q,kappa);nt=float(profile['nt'][0]);k=float(values[2])
        tensor=2*k/3*np.eye(3)-nt*(G+G.T);ev,vec=np.linalg.eigh(tensor)
        dense=resample(dict(report,wp2_common_samples=dict(report['wp2_common_samples'],values=v['values']),physical_samples=dict(wall_K=v['endpoint_K'])),28)
        R=dense['R'];minimum=np.minimum(R[:,1,1],.5*(R[:,0,0]+R[:,2,2]-np.hypot(R[:,0,0]-R[:,2,2],2*R[:,0,2])))
        witnesses=[]
        for region in [0,1]:
            ids=np.flatnonzero(dense['regions']==region);lim=int(ids[np.argmin(minimum[ids])]);witnesses.append(dict(region=region,minimum=float(minimum[lim]),physical=dense['physical'][lim].tolist()))
        result.append(dict(stage=name,state_sha256=state_hash(state),native_witness=dict(physical=point,geometry=geometry,tangential_side=1,radial_side=0,k=k,endpoint_K=K,dependent_epsilon=float(profile['E'][0]),eddy_viscosity=nt,gradient=G.tolist(),tensor=tensor.tolist(),eigenvalues=ev.tolist(),critical_eigenvector=vec[:,0].tolist(),actual_production=float(-np.sum(tensor*G))),GL28_volume_witnesses=witnesses))
    require(state_hash(packet['raw'])==state_hash(packet['corrected']),'unexpected bound change')
    require(state_hash(packet['corrected'])==state_hash(report['rejected_state']),'terminal candidate differs')
    require(state_hash(packet['fixed_reference'])==state_hash(report['fixed_reference']),'fixed reference differs')
    first=next((x['stage'] for x in result if x['native_witness']['eigenvalues'][0]<0 or any(w['minimum']<0 for w in x['GL28_volume_witnesses'])),None)
    residual=np.asarray(packet['transport_residual']);increment=np.asarray(packet['raw_increment']);n=len(residual)
    matrix=np.fromfile(packet['matrix_file'],dtype='<f8').reshape(n,n)
    action=matrix@increment+residual
    require(np.max(abs(action))<2e-12,'independent frozen k linear solve replay differs')
    atomic(args.output,dict(version='WP2SavedStageSequenceV1',stages=result,first_negative_on_declared_diagnostic_sets=first,
        project_margin=1e-4*K0,independent_k_linear_residual_max=float(np.max(abs(action))),native_k_linear_residual=attempt['k_frozen_linear_residual'],
        raw_corrected_bitwise_identical=True,report=dict(path=str(path),sha256=sha(path)),packet=dict(path=str(packet_path),sha256=sha(packet_path)),
        source_sha256=sha(Path(__file__)),seconds=time.monotonic()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Read-only saved stages, exact native-sided limiting point and seen GL28 volume grid. No new PDE solve or complete production-inventory certification; no state, history or source modification.'))


if __name__=='__main__':main()
