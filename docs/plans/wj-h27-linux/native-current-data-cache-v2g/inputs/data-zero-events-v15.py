from pathlib import Path
import json,hashlib,time,sys,math
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/coupled-fixtures-v15';sys.path.insert(0,str(r/'tools-data-validation-v15/h27'))
import implicit_algorithm_check as check
from resources import atomic,identity
began=time.monotonic();rows=[check.native_json(s)for s in(o/'poly-driftzero.steps.jsonl').read_text().splitlines()]
assert[x['step']for x in rows]==[1,2,3,4]
previous=None;drift_events=0;events=0
for row in rows:
    stats=check.cache_policy_statistics(row,'CurrentDataDriftV1')
    assert stats['relative_drift_limit']==0 and stats['maximum_uses']==8 and stats['counts_complete']and row['retries']==0
    for event in stats['preconditioner_events']:
        alpha=(1 if row['step']==1 else 1.5)/row['dt'];assert event['alpha']==alpha and math.isfinite(event['drift'])and event['drift']>=0
        age=0 if previous is None else previous['age_used'];assert event['age_before']==age and event['attempt']==0
        reason='cold'if previous is None else'alpha'if previous['alpha']!=alpha else'drift'if event['drift']>0 else'age'if age>=8 else'reuse'
        assert event['reason']==reason and event['refreshed']==(reason!='reuse')and event['age_used']==(age+1 if reason=='reuse'else 1)
        drift_events+=reason=='drift';events+=1;previous=event
    assert row['rebuilds']==3*stats['refreshes']and max(row['residuals'])<=3e-12
assert drift_events>0
atomic(o/'driftzero-event-independent.json',dict(schema='WJH27ZeroDriftEventControlV1',status='PASS_NUMERICAL_REFRESH_EVENT_CONTROL_ONLY',
    full_H27_qualification=False,source=identity(__file__),checker=identity(check.__file__),native_log=identity(o/'poly-driftzero.steps.jsonl'),
    preconditioner_policy='CurrentDataDriftV1',drift_measure='ScaledCurrentQuadratureDataL2V1',events=events,
    nonzero_drift_refreshes=drift_events,source_identity=rows[0]['source_identity'],elapsed_seconds=time.monotonic()-began))
print('ZERO_CURRENT_DATA_DRIFT_EVENTS_PASS',drift_events)
