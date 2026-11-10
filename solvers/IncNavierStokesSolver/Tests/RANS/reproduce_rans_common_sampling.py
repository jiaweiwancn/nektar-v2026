#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Portable one-row witness for the corrected raw/effective native test mismatch."""
import argparse,gzip,json,os,hashlib
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from rans_spatial_resolution import basis_matrix
from rans_qualification import write_json

def read(p):
 with gzip.open(p,'rt') as f:return json.load(f)
def verify(w):
 column=np.array(w['prolongation_column']);fine=w['fine_test'];raw=w['raw_common_test'];fixed=w['effective_common_test'];j=w['common_mode']
 mapped=float(np.array(fine['phi'])@column[fine['ids']]);value=lambda r:float(sum(v for i,v in zip(r['ids'],r['phi'],strict=True) if i==j))
 bad=abs(mapped-value(raw));good=abs(mapped-value(fixed));shift=max(abs(a-b) for a,b in zip(w['raw_coordinate'],w['effective_coordinate'],strict=True))
 assert bad>w['value_tolerance'] and good<=w['value_tolerance'] and shift<=1.01e-11
 return dict(status='pass',original_sampling_status='fails_identity',corrected_sampling_status='pass',raw_identity_defect=bad,effective_identity_defect=good,raw_effective_coordinate_shift=shift,value_tolerance=w['value_tolerance'])
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--witness',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.witness:d=verify(json.loads(a.witness.read_text()));assert not a.output.exists();write_json(a.output,d);return 0
 assert a.root is not None;r=a.root;old=r/'evolved-audits/previous-v2-E-N-EP-h0p5';new=r/'evolved-audits/previous-v3-E-N-EP-h0p5';paths=[old/'T-0.0016/packet.json.gz',old/'T-0.0016-common-cross/packet.json.gz',new/'T-0.0016-common-cross/packet.json.gz'];f,raw,fixed=[read(p) for p in paths];rows=f['observation']['volume'];nf=f['observation']['sizes'][0];nc=read(r/'diagnosis/initial-EP/common/packet.json.gz')['observation']['sizes'][0];vf=basis_matrix([q['velocity_k'] for q in rows],nf);vc=basis_matrix([q['velocity_k'] for q in fixed['test_basis']],nc);vr=basis_matrix([q['velocity_k'] for q in raw['test_basis']],nc);weights=np.array([q['weight'] for q in rows]);P=cho_solve(cho_factor(np.array(f['resolution']['mass'][0]),lower=True),(vf.T@vc.multiply(weights[:,None])).toarray());bad=vf@P-vr.toarray();i,j=np.unravel_index(np.argmax(np.abs(bad)),bad.shape)
 w=dict(schema='NativeCommonSamplingWitnessV1',model=f['observation']['resolved_descriptor'],source_files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],raw_coordinate=raw['observation']['grid'][i],effective_coordinate=fixed['observation']['grid'][i],fine_test=rows[i]['velocity_k'],raw_common_test=raw['test_basis'][i]['velocity_k'],effective_common_test=fixed['test_basis'][i]['velocity_k'],common_mode=int(j),prolongation_column=P[:,j].tolist(),value_tolerance=2e-10,source_operator_modified=False)
 assert not a.output.exists();a.output.mkdir();write_json(a.output/'witness.json',w);write_json(a.output/'results.json',verify(w));return 0
if __name__=='__main__':raise SystemExit(main())
