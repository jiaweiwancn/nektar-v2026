from pathlib import Path
import sys,json
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');sys.path.insert(0,str(r/'tools-terminal-validation-v14/h27'))
import terminal_output_check as check
from resources import atomic,identity
import numpy as np
o=r/'science/coupled-terminal-v14';prior=check.read(o/'terminal-independent.json');result={}
for label,a,b,limits in [('scheduled_vs_dense','scheduled','dense',check.LIMITS),('cold_vs_uninterrupted','stitched','scheduled',{k:1e-8 for k in check.LIMITS})]:
 pa=o/('terminal-independent-'+a+'.npz');pb=o/('terminal-independent-'+b+'.npz')
 check.require(identity(pa)['sha256']==prior['reductions'][a]['artifact']['sha256']and identity(pb)['sha256']==prior['reductions'][b]['artifact']['sha256'],'reduction identity changed')
 norm=check.po.Norms()
 with np.load(pa,allow_pickle=False)as aa,np.load(pb,allow_pickle=False)as bb:
  groups=sorted({key.split('__')[0]for key in aa.files})
  for group in groups:
   region,field=group.rsplit('_',1)
   for moment,key in [('mean','mean_scaled'),('sd','standard_deviation_scaled')]:
    for meta in ['scale','x','r','weight','owner']:check.require(np.array_equal(aa[group+'__'+meta],bb[group+'__'+meta]),'observation inventory changed')
    norm.add(moment+'/'+region,field,aa[group+'__'+key]-bb[group+'__'+key],aa[group+'__weight'])
 result[label]=norm.report(limits)
packet=dict(schema='WJH27TerminalStatisticFailureDetailV1',status='DIAGNOSTIC_ONLY',prior=identity(o/'terminal-independent.json'),sources=dict(checker=identity(check.__file__)),comparisons=result)
atomic(o/'terminal-failure-detail.json',packet)
print(json.dumps({label:{k:v for k,v in rows.items()if v.get('status')!='MET'}for label,rows in result.items()}))
