#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Decompose sealed A95 response samples without changing fields or equations."""
import argparse, gzip, json, os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest, sha
from rans_duration_contract import require
from rans_preservation_oracle import Reference
from rans_qualification import write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparison',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();comparison=load(a.comparison/'results.json');cfg=comparison['configuration']
    require(cfg['model']=='A95' and cfg['wall']=='Accepted' and comparison['time']==.016,'stress_decomposition_foreign_experiment')
    ref=Reference(cfg['model'],cfg['policy'],cfg['wall']);inputs={};responses=[];max_defect=0.
    vector_file=a.comparison/'response-vectors.npz'
    require(sha(vector_file)==comparison['response_vectors_sha256'],'stress_decomposition_vectors_changed')
    with np.load(vector_file,allow_pickle=False) as vectors:weights=vectors['Ryy_weights'].copy()
    for member in range(4):
        arrays=[]
        for kind in ('perturbed','baseline'):
            path=a.comparison/(str(member)+'-'+kind)/'packet.json.gz';inputs[str(path)]=sha(path)
            with gzip.open(path,'rt') as stream:packet=json.load(stream)
            rows=[];identities=[]
            for sample in packet['native_constitutive']['samples']:
                index=sample['sample_index'];g=np.array(sample['raw_gradient']);R=np.array(sample['R']);Q=np.array(sample['Q'])
                require(g.shape==(3,3) and np.max(np.abs(g[1]))==0 and np.max(np.abs(g[:,1]))==0,
                        'stress_decomposition_requires_declared_flat_2D')
                div=float(np.trace(g));k=sample['k'];nt=sample['nu_t'];sd=(g+g.T)/2-div*np.eye(3)/3
                scalar=2*k/3;divergence=2*nt*div/3
                iso=-2*k*div/3;nonlinear=-float(np.sum(Q*g));linear=2*nt*float(np.sum(sd*sd))
                max_defect=max(max_defect,abs(R[1,1]-(scalar+divergence)),abs(Q[1,1]),
                    abs(div-sample['raw_divergence']),abs(sample['P']-(iso+nonlinear+linear)))
                rows.append([R[1,1]/ref.us**2,scalar/ref.us**2,divergence/ref.us**2,
                    div/(ref.us/ref.H),sample['P']/(ref.A/ref.H),iso/(ref.A/ref.H),
                    nonlinear/(ref.A/ref.H),linear/(ref.A/ref.H)])
                identities.append((index,sample['coordinate'],sample['side'],sample['geometry']))
            if kind=='perturbed':previous=identities
            else:require(previous==identities,'stress_decomposition_baseline_sample_mismatch')
            require(len(rows)==len(weights),'stress_decomposition_weights_differ')
            arrays.append(np.array(rows))
        responses.append(arrays[0]-arrays[1])
    require(max_defect<=5e-14,'stress_decomposition_native_identity_defect')
    intent=dict(schema='DirectionalStressResponseDecompositionV1',comparison_sha256=sha(a.comparison/'results.json'),
        packet_sha256=inputs,response_vectors_sha256=sha(vector_file),configuration=cfg,time=.016,
        scope='read-only sealed native samples; not a new evolution or changed response gate')
    seal=RequestManifest(a.output,intent);labels=('Ryy','Ryy_scalar_part','Ryy_deviatoric_divergence_part',
        'raw_divergence','P','P_isotropic','P_nonlinear','P_linear')
    def rms(v):return float(np.sqrt(np.dot(weights,v*v)/weights.sum()))
    rows={name:dict(response_RMS=[rms(v[:,j]) for v in responses],
        spatial_difference=rms(responses[0][:,j]-responses[1][:,j]),
        timestep_differences=[rms(responses[0][:,j]-responses[2][:,j]),rms(responses[1][:,j]-responses[3][:,j])])
        for j,name in enumerate(labels)}
    original=comparison['fields']['Ryy']
    require(np.allclose(rows['Ryy']['response_RMS'],original['response_RMS'],rtol=1e-13,atol=1e-18),
            'stress_decomposition_original_response_changed')
    result=dict(status='pass',maximum_dimensional_identity_defect=max_defect,statistics=rows,
        normalization='Ryy by u_star^2, raw divergence by u_star/H, production by epsilon(H); weighted physical RMS on original outer samples',
        identity='A95 flat-2D: Qyy=0 and Ryy=2k/3+2*nu_t*trace(raw A)/3 under DeviatoricStrainV1; signed P includes all three terms',
        interpretation='weak discrete continuity does not require pointwise trace(raw A)=0; this decomposition diagnoses its measured stress contribution, without changing the extension or imposing an extra constraint')
    write_json(a.output/'results.json',result);seal.seal(a.output,['results.json']);print('stress/production decomposition reproduced sealed native samples');return 0


if __name__=='__main__':raise SystemExit(main())
