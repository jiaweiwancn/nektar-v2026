from pathlib import Path
import json,subprocess,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');s=r/'trace-geometry-v17';base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
run('trace-geometry-configure-v17','build',30,['cmake','-S',s,'-B',s/'build','-DCMAKE_BUILD_TYPE=Release','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake'])
run('trace-geometry-build-v17','build',100,['cmake','--build',s/'build','--target','WJH27TraceGeometryProbe','--parallel','1'])
rt=json.loads((r/'science/runtime-v15.json').read_text());exe=s/'build/WJH27TraceGeometryProbe'
manifest=dict(scope='Read-only H point/native-face diagnostic linked against exact source15 static archives; no production code changed',binary=dict(path=str(exe),sha256=hashlib.sha256(exe.read_bytes()).hexdigest()),
 static_libraries=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in [r/'build-v15/libWJNativeMixedFlow.a',r/'build-v15/h26-kernel/libWJNativeFlow.a']],runtime_source=str(r/'science/runtime-v15.json'))
with(s/'manifest.json').open('x')as f:json.dump(manifest,f,indent=2);f.write('\n')
for name,case in [('r2',r/'science/physical-reference-r2-v17/reference-r2.xml'),('h2',r/'science/physical-reference-h2-v1/reference-h2.xml')]:
 run('trace-geometry-'+name+'-native-v17','science',40,[rt['loader']['path'],'--library-path',rt['libraries'],exe,case,r/'science'/('trace-geometry-'+name+'-v17.json')])
print('TRACE_GEOMETRY_DIAGNOSTICS_FINISHED',flush=True)
