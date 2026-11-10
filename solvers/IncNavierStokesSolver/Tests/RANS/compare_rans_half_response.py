#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Pointwise half-amplitude response check, with fixed physical RMS floors."""
import argparse,os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from qualify_rans_preservation import load
from rans_qualification import write_json
from rans_request_manifest import sha
from rans_duration_contract import require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--full',type=Path,required=True);p.add_argument('--half',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();full,half=[load(x/'results.json') for x in (a.full,a.half)];require(full['configuration']==half['configuration'] and full['timesteps']==half['timesteps'],'half_response_configuration_mismatch');require(abs(full['time']-half['time'])<1e-13 and full['amplitude']==2*half['amplitude'],'half_response_amplitude_or_time_mismatch')
    for path,d in ((a.full,full),(a.half,half)):require(sha(path/'response-vectors.npz')==d['response_vectors_sha256'],'half_response_vectors_changed')
    f,h=[np.load(x/'response-vectors.npz') for x in (a.full,a.half)];rows=[]
    for name in full['fields']:
        w=f[name+'_weights'];require(np.array_equal(w,h[name+'_weights']),'half_response_weights_mismatch')
        for member in range(4):
            x,y=f[str(member)+'_'+name],h[str(member)+'_'+name];require(x.shape==y.shape==w.shape,'half_response_array_shape');signal=float(np.sqrt(np.sum(w*y*y)/w.sum()));difference=float(np.sqrt(np.sum(w*(x/2-y)**2)/w.sum()));limit=max(5e-7,.1*signal);rows.append(dict(field=name,member=member,half_signal_RMS=signal,departure_from_linear_scaling_RMS=difference,limit=limit,status='pass' if difference<=limit else 'unresolved',signal_status='resolved' if signal>5e-7 else 'unresolved_floor'))
    result=dict(status='pass' if all(x['status']=='pass' for x in rows) else 'unresolved',configuration=full['configuration'],time=full['time'],full_amplitude=full['amplitude'],half_amplitude=half['amplitude'],checks=rows,scope='bounded amplitude response; no exact linearity or decay theorem');require(not a.output.exists(),'half_response_output_exists');write_json(a.output,result);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
