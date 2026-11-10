#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Partition native dual-norm energy by exact basis-support mass fractions.

Signed element contributions sum to the full squared norm. This is localization,
not an additive norm of separately integrated element loads (which cancel).
"""
import argparse,gzip,json,os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from rans_qualification import write_json
from rans_preservation_oracle import Reference
from qualify_rans_preservation import load

def energy(packet,ref):
 o=packet['observation'];r=packet['resolution'];bc={int(i) for i,v in r['constraints']};out={}
 for block,name,scale in ((0,'u',ref.us**2/ref.H),(2,'k',ref.K*ref.us/ref.H)):
  n=o['sizes'][block];off=o['offsets'][block];mass=np.array(r['mass'][block]);free=np.array([i for i in range(n) if i+off not in bc]);load=np.array(r['native']['spatial_unmasked'])[off+free];inverse=cho_solve(cho_factor(mass[np.ix_(free,free)],lower=True),load);density=np.zeros(n);density[free]=load*inverse/(scale**2*8.)
  diagonal={}
  for q in o['volume']:
   row=q['velocity_k'];g=row['geometry'];a=diagonal.setdefault(g,np.zeros(n));a[np.array(row['ids'],int)]+=q['weight']*np.square(row['phi'])
  total=sum(diagonal.values());relative=float(np.max(np.abs(total-np.diag(mass)))/np.max(np.diag(mass)));assert relative<1e-11
  elements=[]
  for g,part in diagonal.items():
   e=next(x for x in o['elements'] if x['geometry']==g);value=float(np.sum(density*part/total));elements.append(dict(geometry=g,interval=[e['lower'][2],e['upper'][2]],signed_squared_norm=value,absolute_row_energy=float(np.sum(np.abs(density)*part/total))))
  squared=float(load@inverse)/(scale**2*8.);assert abs(sum(x['signed_squared_norm'] for x in elements)-squared)<1e-12
  out[name]=dict(full_squared_norm=squared,elements=sorted(elements,key=lambda x:x['absolute_row_energy'],reverse=True),mass_diagonal_relative_check=relative)
 return out

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--wall',choices=('EP','Accepted'),required=True);a=p.parse_args();ref=Reference('E94E95','Natural',a.wall);rows=[]
 for x in load(a.directory/'results.json')['checks']:
  if x['status']!='pass':continue
  with gzip.open(Path(x['source'])/'packet.json.gz','rt') as f:packet=json.load(f)
  rows.append(dict(name=x['name'],energy=energy(packet,ref)))
 if a.output.exists():raise FileExistsError(a.output)
 write_json(a.output,dict(status='pass',checks=rows,scope=__doc__));return 0
if __name__=='__main__':raise SystemExit(main())
