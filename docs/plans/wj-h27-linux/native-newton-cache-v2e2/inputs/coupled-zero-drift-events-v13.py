from pathlib import Path
import json,hashlib,time,sys
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/coupled-fixtures-v13';began=time.monotonic();rows=[json.loads(s,parse_int=lambda t:-0.0 if t=='-0'else int(t))for s in(o/'poly-driftzero.steps.jsonl').read_text().splitlines()];assert [x['step']for x in rows]==[1,2,3,4]
previous=None;drift_events=0;events=0
for row in rows:
 s=row['coupled_newton'];assert s['preconditioner_policy']=='OperatorDriftV1'and s['relative_drift_limit']==0 and s['maximum_uses']==8 and s['counts_complete']and row['retries']==0
 for event in s['preconditioner_events']:
  alpha=(1 if row['step']==1 else 1.5)/row['dt'];assert event['alpha']==alpha
  age=0 if previous is None else previous['age_used'];assert event['age_before']==age and event['attempt']==0
  reason='cold'if previous is None else'alpha'if previous['alpha']!=alpha else'drift'if event['drift']>0 else'age'if age>=8 else'reuse'
  assert event['reason']==reason and event['refreshed']==(reason!='reuse')and event['age_used']==(age+1 if reason=='reuse'else 1)
  drift_events+=reason=='drift';events+=1;previous=event
 assert row['rebuilds']==3*s['refreshes']and max(row['residuals'])<=3e-12
assert drift_events>0
out=o/'driftzero-event-independent.json';assert not out.exists();out.write_text(json.dumps(dict(schema='WJH27ZeroDriftEventControlV1',status='PASS_NUMERICAL_REFRESH_EVENT_CONTROL_ONLY',full_H27_qualification=False,events=events,nonzero_drift_refreshes=drift_events,native_log=dict(path=str(o/'poly-driftzero.steps.jsonl'),sha256=hashlib.sha256((o/'poly-driftzero.steps.jsonl').read_bytes()).hexdigest()),elapsed_seconds=time.monotonic()-began),indent=2)+'\n');print('ZERO_DRIFT_EVENTS_PASS',drift_events)
