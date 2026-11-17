from pathlib import Path
import shutil,json
from resources import identity,atomic
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');old=r/'tools-cutoff-recovery-v19';new=r/'tools-cutoff-recovery-v19b';assert not new.exists();shutil.copytree(old,new)
for p in sorted(old.rglob('*')):
 if p.is_file():assert identity(p)['sha256']==identity(new/p.relative_to(old))['sha256']
source=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/cao-tamura-fig3/rough-fig3.csv';target=new/'cao-tamura-fig3/rough-fig3.csv';target.parent.mkdir();shutil.copy2(source,target)
manifest=json.loads((r/'source-v19.json').read_text());row=next(v for v in manifest['files']if v['path']==str(source));assert identity(target)['sha256']==row['sha256']
atomic(new/'overlay-v19b.json',dict(parent=identity(old/'overlay.json'),source=identity(r/'source-v19.json'),fixture=identity(target),source_fixture=row,all_original_helper_bytes_identical=True,scope='Missing report-only roughness CSV packaging repair; no equation/model/native input/helper code change'))
print('FROZEN_ROUGHNESS_FIXTURE_RESTORED')
