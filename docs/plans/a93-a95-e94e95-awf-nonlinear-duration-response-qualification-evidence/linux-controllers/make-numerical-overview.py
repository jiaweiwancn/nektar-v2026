from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
result=dict(status='pass',cases=[])
for model in ('A93','A95'):
 for member in ('B5','B11','B5-dt2','B11-dt2'):
  label=model+'-N-Accepted-duration-0.002-'+member;folder=r/'campaigns'/label;request=json.loads((folder/'request-manifest.json').read_text());q=request['request'];row=dict(label=label,configuration={k:q[k] for k in ('model','policy','wall','descriptor')},dt=q['dt'],request_sha256=request['sha256'],mesh=q['declared_mesh'],milestones=[])
  for T in (.016,.1,.5):
   p=folder/('T-'+format(T,'.4f'))/'metrics.json';d=json.loads(p.read_text());accuracy=d['all_audited_accuracy'];maxima=[max(accuracy,key=lambda v:abs(v['values'][j])) for j in range(12)];initial={k:v for k,v in d['initialization_receipt'].items() if k not in ('live_fixed_boundary','effective_initial')}
   row['milestones'].append(dict(time=d['time'],theta=d['theta'],status=d['status'],reference_error=d['reference_error'],change_from_actual_initial=d['change_from_initial'],global_pressure_gauge=d['global_pressure_gauge'],accuracy_maxima=[dict(component=j,value=v['values'][j],time=v['time'],stage='accepted') for j,v in enumerate(maxima)],domain_extrema=d['extrema'],source_and_constraint_maxima=d['source_and_constraint_maxima'],final_diagnostic=d['final_audit'],initialization=initial,history=d['history_identity'],state_sha256=d['final_state_sha256'],configuration_sha256=d['configuration_sha256'],metrics_sha256=sha(p),native_sha256=sha(p.parent/'case.evolution.json')))
  result['cases'].append(row)
(r/'numerical-overview.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
