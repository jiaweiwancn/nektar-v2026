from pathlib import Path
import json,shutil,copy
import numpy as np
from numpy.polynomial.chebyshev import chebder
from resources import identity,atomic
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/reference-inverse-v19';o.mkdir();t=r/'reference-inverse-project-v19';(t/'reference_inverse').mkdir(parents=True)
for name,target in [('reference-inverse-probe-v19.cpp','reference_inverse_probe.cpp'),('reference-inverse-CMakeLists-v19.txt','reference_inverse/CMakeLists.txt'),('reference-inverse-check-v19.py','reference_inverse_check.py'),('reference-inverse-admission-v19.md','admission.md')]:shutil.copy2(r/'inputs'/name,t/target)
table=r/'science/reference-tables-v1.json';data=json.loads(table.read_text());lo,hi=data['K_interval'];mutations=[]
for name in('constant_cutoff','increasing_cutoff','two_inverse_roots','stale_derivative'):
 d=copy.deepcopy(data);node=d['cutoff_table'];value=node['coefficients'][0]
 if name=='stale_derivative':node['logcutoff_dK_coefficients'][0]+=.01
 else:
  node['coefficients']=[value,0.]if name=='constant_cutoff'else[value,.02]if name=='increasing_cutoff'else[value,0.,.02]
  node['degree']=len(node['coefficients'])-1;node['logcutoff_dK_coefficients']=(chebder(node['coefficients'])*2/(hi-lo)).tolist()
 p=o/(name+'.json');atomic(p,d);mutations.append(dict(name=name,file=identity(p),expected_reason='analytic derivative coefficients differ'if name=='stale_derivative'else'cutoff table has no certified strictly decreasing inverse'))
atomic(o/'mutations.json',dict(schema='WJH27ReferenceInverseMutationsV1',table=identity(table),mutations=mutations,synthetic_unqualified_controls=True,TASK3_READY=False))
libs=[r/'build-v19/libWJNativeMixedFlow.a',r/'build-v19/h26-kernel/libWJNativeFlow.a'];atomic(t/'manifest.json',dict(schema='WJH27ReferenceInverseProjectV1',parent='d205335c80458009387430d2fa54f3c3e39ef560',source=identity(r/'source-v19.json'),runtime=identity(r/'science/runtime-v19.json'),libraries=[identity(p)for p in libs],files=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()],mutations=identity(o/'mutations.json')))
print('REFERENCE_INVERSE_FIXTURES_READY')
