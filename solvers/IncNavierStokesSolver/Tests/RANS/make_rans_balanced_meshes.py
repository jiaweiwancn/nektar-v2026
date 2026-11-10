#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed localized hierarchy selected from signed support-energy localization."""
import argparse
from pathlib import Path
from rans_preservation_oracle import Reference
from rans_qualification import write_json

def specifications(ref):
 base=[x for x in ref.knots() if x!=ref.kernel.ze];pairs=list(zip(base[:-1],base[1:]));inner=[v for v in pairs if v[1]<=ref.H];outer=[v for v in pairs if v[0]>=ref.H]
 assert len(inner)==11 and len(outer)==5
 # Leave the largest remaining bulk panels for the finer member; refine the
 # ground/source/endpoint neighbors in the candidate. Both are fixed meshes.
 coarse_inner=[0,4,5,10];fine_inner=[0,1,4,5,6,7,8,9,10]
 rows={}
 for name,ii,oi in [('Balanced5V1',coarse_inner,[0]),('Balanced11V1',fine_inner,[0,2])]:
  chosen=[inner[i] for i in ii]+[outer[i] for i in oi]
  rows[name]=dict(schema='BalancedFixedMeshV1',degree=5,pressure_degree=3,knots=sorted(set(base+[(a+b)/2 for a,b in chosen])),split_intervals=chosen,
   omitted_inactive_mesh_break=ref.kernel.ze,retained_active_breaks=sorted(set([0.,ref.H,ref.Z,.1,ref.zc]+([ref.kernel.zv] if 0<ref.kernel.zv<ref.H else []))),
   policy='mesh only; original epsilon breakpoint is inactive for EP plateau and Accepted; native physical profile panels remain unchanged',
   rationale='full signed energy: remove inactive tiny-panel test artifact; resolve ground/cutoff, bulk 0.02..1, interface and outer contributions under 3000 unknowns')
 assert set(rows['Balanced5V1']['knots'])<set(rows['Balanced11V1']['knots'])
 return rows

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--model',default='E94E95',choices=('E94E95','A95','A93'));p.add_argument('--policy',default='Natural',choices=('Natural','Prescribed'));p.add_argument('--wall',required=True,choices=('EP','Accepted'));a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 rows=specifications(Reference(a.model,a.policy,a.wall));write_json(a.output/'spaces.json',rows)
 for name,s in rows.items():write_json(a.output/(name+'.json'),s)
 return 0
if __name__=='__main__':raise SystemExit(main())
