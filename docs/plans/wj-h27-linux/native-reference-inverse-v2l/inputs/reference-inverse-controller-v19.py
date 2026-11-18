from pathlib import Path
import json,subprocess
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
assert not any(e['status']=='RESERVED'for e in json.loads((r/'accounting.json').read_text())['entries'].values())
run('reference-inverse-fixtures-v19','analysis',30,['python3','-B',r/'inputs/reference-inverse-fixtures-v19.py'])
t=r/'reference-inverse-project-v19';b=r/'build-reference-inverse-v19';o=r/'science/reference-inverse-v19'
run('reference-inverse-configure-v19','build',40,['cmake','-S',t/'reference_inverse','-B',b,'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON','-DH27_PINNED_SOURCE='+str(r/'source-v19'),'-DH27_PINNED_BUILD='+str(r/'build-v19')])
run('reference-inverse-build-v19','build',180,['cmake','--build',b,'--parallel','4'])
old=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime');libs=str(old/'glibc')+':'+str(old/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64'
run('reference-inverse-native-v19','science',20,[old/'glibc/ld-linux-x86-64.so.2','--library-path',libs,b/'WJH27ReferenceInverseProbe',r/'science/reference-tables-v1.json',o/'mutations.json',o/'native.json'])
run('reference-inverse-independent-v19','analysis',40,['python3','-B',t/'reference_inverse_check.py','--table',r/'science/reference-tables-v1.json','--native',o/'native.json','--manifest',o/'mutations.json','--out',o/'independent.json'])
assert json.loads((o/'independent.json').read_text())['inverse_guard_gate_pass']is True
run('reference-inverse-reuse-audit-v19','analysis',30,['python3','-B',r/'inputs/verification-reuse-audit-v19.py'])
print('SOURCE19_REFERENCE_INVERSE_QUALIFICATION_FINISHED',flush=True)
