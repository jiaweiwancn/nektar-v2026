from pathlib import Path
import subprocess
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
subprocess.run(base+['--label','data-history-overlay-v15b','--phase','analysis','--seconds','20','--','python3','-B',str(r/'inputs/data-history-overlay-v15b.py')],check=True)
o=r/'science/coupled-fixtures-v15';t=r/'tools-data-validation-v15b/h27'
subprocess.run(base+['--label','physical-data-cache-history-independent-v15b','--phase','analysis','--seconds','80','--','python3','-B',str(t/'physical_cache_history_check.py'),'--manifest',str(o/'physical-history-manifest.json'),'--table',str(r/'science/reference-tables-v1.json'),'--out',str(o/'physical-history-independent-v2.json'),'--seconds','75'],check=True)
print('SOURCE15_REPAIRED_HISTORY_PHASE_FINISHED',flush=True)
