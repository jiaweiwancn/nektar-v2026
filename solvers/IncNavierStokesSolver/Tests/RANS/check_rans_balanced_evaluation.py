#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Same-coefficient GL10/8 versus GL14/10 physical evaluation uncertainty."""
import argparse,os,shutil
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load,replay_request
from rans_preservation_oracle import Reference,physical_grid,representation
from diagnose_rans_balanced_spatial import run_probe
from run_rans_neutral_pilot import scaled,FIELDS

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'evaluation-uncertainty';out.mkdir(exist_ok=False);os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'));rows=[]
 for wall in ('EP','Accepted'):
  ref=Reference('E94E95','Natural',wall);coarse=r/'campaigns'/('E-N-'+wall+'-B5');fine=r/'campaigns'/('E-N-'+wall+'-B11');baseline=[load(f/'initial-metrics.json')['accuracy']['reference_error'] for f in (coarse,fine)]
  for index,folder in enumerate((coarse,fine)):
   seed=load(folder/'initial-packet.json');request=load(folder/'request-manifest.json')['request'];panels=sorted(set(request['knots']+ref.knots()));grid,weights=physical_grid(panels,14,10);target=out/(wall+'-'+str(index));target.mkdir();shutil.copy2(folder/'T-0.0016/case.xml',target/'case.xml')
   packet,_=run_probe(r/'runtime',target,replay_request(seed,grid,state=seed['reference']));receipt=packet['initialization_receipt'];old=seed['initialization_receipt']
   for key in ('raw_initial_sha256','effective_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256'):assert receipt[key]==old[key]
   high=scaled(representation(ref,packet['observation'],weights),ref);checks={}
   for f in FIELDS:
    delta=abs(high[f]['weighted_RMS']-baseline[index][f]['weighted_RMS']);spatial=abs(baseline[0][f]['weighted_RMS']-baseline[1][f]['weighted_RMS']);limit=max(1e-10,.01*spatial);checks[f]=dict(low_RMS=baseline[index][f]['weighted_RMS'],high_RMS=high[f]['weighted_RMS'],evaluation_change=delta,limit=limit,status='pass' if delta<=limit else 'fail')
   rows.append(dict(wall=wall,space='Balanced5V1' if index==0 else 'Balanced11V1',status='pass' if all(x['status']=='pass' for x in checks.values()) else 'fail',checks=checks,effective_inputs_identical=True,source_grid=request['evaluation']['grid'],target_grid_points=len(grid)))
 status='pass' if all(x['status']=='pass' for x in rows) else 'fail';write_json(out/'results.json',dict(status=status,checks=rows,metric='scaled weighted physical RMS, common global pressure gauge'));return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
