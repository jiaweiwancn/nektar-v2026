#!/usr/bin/env python3
"""Physical all-slot history recheck of archived native M4 comparisons."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source',type=Path);p.add_argument('qualification',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();sys.path.insert(0,str(a.source/'solvers/IncNavierStokesSolver/Tests/RANS'))
import numpy as np
from validate_cls_epawf_m4 import spaces
pairs=[('evolving-condensed','evolving-o2-dt0.1'),('evolving-mpi2','evolving-o2-dt0.1'),('evolving-mpi4','evolving-o2-dt0.1'),
       ('stationary-restart-mpi4','stationary-o2-dt0.05')]
for direction in ('3D-straight','3D-rotated'):
 for ranks in (2,4):pairs.append((direction+'-mpi'+str(ranks),direction))
for order,split in ((1,1),(2,0),(2,1),(2,2),(2,3)):
 pairs.append((f'evolving-o{order}-restart{split}',f'evolving-o{order}-dt0.1'))
checks=0;records=[]
for first,second in pairs:
 paths=[a.qualification/name/'case.evolution.json' for name in (first,second)]
 A,B=(json.loads(path.read_text()) for path in paths)
 assert A['dimension']==B['dimension'] and A['order']==B['order'] and A['offsets']==B['offsets'] and A['dt']==B['dt']
 assert np.max(np.abs(np.asarray(A['native_history_times'])-np.asarray(B['native_history_times'])))<=1e-8
 d,us,H=A['dimension'],A['ustar'],A['parameters'][1];K,kappa=A['parameters'][5],A['parameters'][7]
 scales=[us]*d+[K,us**3/(kappa*H)];errors=[]
 for block,(field,scale,offset) in enumerate(zip(spaces(A)[:d+2],scales,A['offsets'])):
  mass=field.mass()
  for slot,(x,y) in enumerate(zip(A['native_history'],B['native_history'])):
   delta=np.asarray(x)[offset:offset+field.size]-np.asarray(y)[offset:offset+field.size]
   stored_scale=scale if slot<A['order'] else scale*A['dt']/(H/us)
   error=math.sqrt(max(0.,float(delta@mass@delta)))/(stored_scale*math.sqrt(field.volume))
   checks+=1
   if error>1e-8:raise RuntimeError(f'{first}/{second} block {block}, native slot {slot}: {error}')
   errors.append(dict(block=block,slot=slot,scaled_physical_mass_error=error))
 records.append(dict(first=first,second=second,maximum_scaled_error=max(item['scaled_physical_mass_error'] for item in errors),
                     native_slot_errors=errors,sha256=[hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]))
result=dict(schema=1,status='passed',source_commit=subprocess.check_output(['git','-C',str(a.source),'rev-parse','HEAD'],text=True).strip(),
            comparison_count=len(records),checks=checks,tolerance=1e-8,
            norm='Physical mass by field; value slots use field scale; native dt*RHS slots use dt*field scale/(H/u_star)',
            maximum_scaled_history_error=max(row['maximum_scaled_error'] for row in records),comparisons=records,
            runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='comparisons'},indent=2))
