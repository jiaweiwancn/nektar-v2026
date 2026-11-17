from pathlib import Path
import shutil
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');t=r/'tools-data-validation-v15b'
shutil.copytree(r/'tools-data-validation-v15',t,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(r/'inputs/physical_cache_history_check-v15b.py',t/'h27/physical_cache_history_check.py')
for name in ['physical-history-v15b.md','admission-v15b.md']:shutil.copy2(r/'inputs'/name,t/name)
atomic(t/'history-overlay.json',dict(scope='Read-only descriptor/numbered-frame classification and binding repair; no native or physical data changes',parent=identity(r/'tools-data-validation-v15/data-policy-overlay.json'),sources=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()and p.suffix in ['.py','.csv','.md']]))
