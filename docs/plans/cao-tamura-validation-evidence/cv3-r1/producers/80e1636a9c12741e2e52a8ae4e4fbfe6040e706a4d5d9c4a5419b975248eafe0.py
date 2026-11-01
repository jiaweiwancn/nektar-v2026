"""Read-only decomposition of accepted wall-normal errors at trace nodes."""
from pathlib import Path
import csv,hashlib,json,math,shutil

B=Path('/media/jay/新加卷/nektar-v2026-runs/cao-tamura-calibrated-validation')
R=B/'calibration/cv3/ep-transition-v2/r1'
out=R/'normal-trace-diagnostic-v1';out.mkdir(exist_ok=False)
shutil.copy2(__file__,out/'producer.py')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
result=dict(status='COMPLETE',scope='Accepted final-frame trace diagnostics only; no instability attribution or development acceptance',
 producer_sha256=sha(out/'producer.py'),cases={})
cases={'R-centre-fieldmaximum':R/'development-residual-map/R-cls_fixed/stage-00',
 'R-lower-corner-pointwise':R/'development-lower-corner-pointwise/R-cls_fixed/stage-00',
 'S-centre-pointwise':R/'development-pointwise-scalars/S-cls_fixed/stage-01'}
for name,folder in cases.items():
 with (folder/'case.rans.csv').open(newline='') as f:stages=list(csv.DictReader(f))
 stages=[r for r in stages if r.get('time') and None not in r and all(v is not None for v in r.values())]
 time=float(stages[-1]['time']);manifest=json.loads((folder/'manifest.json').read_text())
 nodes=manifest['config']['x_columns'];rows=[];sources={}
 for p in folder.glob('case.terrain-points-r*.csv'):
  sources[str(p)]=sha(p)
  with p.open(newline='') as f:
   for row in csv.DictReader(f):
    if row.get('phase')=='1' and abs(float(row['time'])-time)<1e-10:
     r={k:float(row[k]) for k in ('x','normal_target','normal_error','normal_projected_error','mask','chi','kP')}
     if not all(math.isfinite(v) for v in r.values()):raise ValueError('Nonfinite accepted trace')
     r.update(normal_velocity=r['normal_target']+r['normal_error'],
       at_element_node=any(abs(r['x']-x)<1e-10 for x in nodes))
     rows.append(r)
 if not rows:raise ValueError('No final accepted trace')
 groups={}
 for row in rows:
  if row['at_element_node']:groups.setdefault(round(row['x'],10),[]).append(row)
 jumps=[]
 for x,values in sorted(groups.items()):
  jumps.append(dict(x=x,copies=len(values),
   normal_target_jump=max(r['normal_target'] for r in values)-min(r['normal_target'] for r in values),
   normal_velocity_jump=max(r['normal_velocity'] for r in values)-min(r['normal_velocity'] for r in values),
   maximum_pointwise_error=max(abs(r['normal_error']) for r in values)))
 result['cases'][name]=dict(folder=str(folder),time=time,source_sha256=sources,
  stage_csv_sha256=sha(folder/'case.rans.csv'),manifest_sha256=sha(folder/'manifest.json'),
  maximum_error=max(rows,key=lambda r:abs(r['normal_error'])),
  maximum_interior_error=max((abs(r['normal_error']) for r in rows if not r['at_element_node']),default=None),
  maximum_projected_error=max(abs(r['normal_projected_error']) for r in rows),
  one_sided_node_values=jumps)
 print(name,'time',time,'maximum',result['cases'][name]['maximum_error'],flush=True)
(out/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
