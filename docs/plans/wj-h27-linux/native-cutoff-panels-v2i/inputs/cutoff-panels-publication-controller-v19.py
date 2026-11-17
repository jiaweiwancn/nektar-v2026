from pathlib import Path
import json,subprocess
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
assert json.loads((r/'operations/cutoff-panels-radial-zero-v19/receipt.json').read_text())['status']=='PASS'
assert json.loads((r/'science/cutoff-panels-v19/radial-zero-independent.json').read_text())['four_binary_histories_exact']
for label,phase,name in [('cutoff-panels-assessment-v19','analysis','cutoff-panels-assessment-v19.py'),('cutoff-panels-publication-v19','evidence','cutoff-panels-publication-v19.py')]:
 subprocess.run(base+['--label',label,'--phase',phase,'--seconds','60','--','python3','-B',str(r/'inputs'/name)],check=True)
print('SOURCE19_EVIDENCE_PUBLISHED',flush=True)
