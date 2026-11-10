#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Isolated one/two-column p4 sensitivity on a partition-resolving common grid."""
import argparse,os,shutil,sys
from pathlib import Path
from rans_qualification import THREADS,compare_samples,write_json
from rans_request_manifest import RequestManifest
from rans_preservation_oracle import Reference
from qualify_rans_preservation import load,probe,replay_request
from run_rans_neutral_pilot import scaled

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path);a=p.parse_args();r=a.root.resolve();os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'));out=a.output or r/'streamwise-sensitivity';out.mkdir(exist_ok=False)
 labels=['E-N-Accepted-h0p4','E-N-Accepted-h0p4-x2'];folders=[r/'campaigns'/x for x in labels];requests=[load(f/'request-manifest.json')['request'] for f in folders];T=.0016
 for key in ('descriptor','knots','degree','pressure_degree','dt','runtime','boundary_policy','velocity_amplitude'):assert requests[0][key]==requests[1][key],key
 for key in ('wall_normal','wall_tangential','outer_volume','mass_override','mass_default','interface','reference_projection'):assert requests[0]['quadrature'][key]==requests[1]['quadrature'][key],key
 assert [x['columns'] for x in requests]==[1,2]
 for f,req in zip(folders,requests,strict=True):assert RequestManifest(f,req,True).completed(f/'T-0.0016')
 two=load(folders[1]/'T-0.0016/endpoint-samples.json');weights=two['grid_weights'];shutil.copy2(folders[0]/'T-0.0016/case.xml',out/'case.xml')
 seed=load(folders[0]/'initial-packet.json');report=load(folders[0]/'T-0.0016/case.evolution.json');request=replay_request(seed,two['grid'],state=report['state'],time=T);request['fixed_boundary']=report['initialization_receipt']['live_fixed_boundary'];one=probe(r/'runtime/RANSPreservationProbe',out,request)['observation'];one['grid_weights']=weights
 difference=scaled(compare_samples(one,two,weights),Reference('E94E95','Natural','Accepted'));assert all(v['weighted_RMS']<=(1e-4 if k=='v' else 1e-3) for k,v in difference.items())
 result=dict(status='pass',time=T,theta=T*.2,columns=[1,2],degree=4,common_x_partitions=[-1.,0.,1.],physical_differences_scaled=difference,source_requests=[load(f/'request-manifest.json')['sha256'] for f in folders],scope='short sensitivity only; p4 fails later radial/interface accuracy; no multi-column .1 qualification')
 write_json(out/'results.json',result);print('pass');return 0
if __name__=='__main__':sys.exit(main())
