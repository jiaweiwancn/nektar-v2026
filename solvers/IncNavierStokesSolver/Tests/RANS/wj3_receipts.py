#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Command receipts and immutable source/runtime/artifact index for WJ-3."""
import argparse
from datetime import datetime,timezone
import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
from run_cls_epawf_picard import THREADS,sha,write_json


def command(args):
    root=args.root.resolve();root.mkdir(exist_ok=True,parents=True);file=root/'command-receipts.json'
    receipts=json.loads(file.read_text()) if file.exists() else []
    if any(r['label']==args.label for r in receipts):raise ValueError('immutable receipt label already exists')
    env=os.environ.copy();env.update(THREADS);start=time.time();log=root/(args.label+'.log')
    argv=args.command[1:] if args.command[:1]==['--'] else args.command
    with log.open('w') as stream:result=subprocess.run(argv,stdout=stream,stderr=subprocess.STDOUT,env=env)
    receipts.append(dict(label=args.label,command=argv,started_utc=datetime.fromtimestamp(start,timezone.utc).isoformat(),
        seconds=time.time()-start,exit_code=result.returncode,log=str(log),log_sha256=sha(log),threads=THREADS))
    current=receipts[-1]
    import fcntl
    with (root/'.receipt.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        merged=json.loads(file.read_text()) if file.exists() else []
        if any(r['label']==args.label for r in merged):raise ValueError('immutable receipt label already exists')
        merged.append(current);write_json(file,merged)
    print(json.dumps(current));return result.returncode


def collect(args):
    root=args.root.resolve();source=root/'source';build=root/'build'
    git=lambda *v:subprocess.check_output(['git','-c','core.filemode=false','-C',str(source),*v],text=True).strip()
    files=subprocess.check_output(['git','-C',str(source),'ls-files','-z']).split(b'\0')
    content={os.fsdecode(f):sha(source/os.fsdecode(f)) for f in files if f and (source/os.fsdecode(f)).is_file()}
    write_json(root/'source-content.json',content)
    runtime={}
    for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','RANSReturnedFieldsProbe'):
        binary=build/name
        ldd=subprocess.check_output(['ldd',str(binary)],text=True);libs={}
        for path in re.findall(r'(?:=>\s+)?(/[^\s()]+)',ldd):
            p=Path(path)
            if p.is_file():libs[str(p.resolve())]=sha(p.resolve())
        runtime[name]=dict(path=str(binary),sha256=sha(binary),ldd=ldd,resolved_libraries=libs)
    artifacts={}
    for scope in [root/'study',root/'study-final',build]:
        if not scope.exists():continue
        for f in scope.rglob('*'):
            if f.is_file() and (f.suffix in ('.json','.xml') or f.name=='State.bin'):
                artifacts[str(f)]=dict(bytes=f.stat().st_size,sha256=sha(f))
    write_json(root/'artifact-index.json',artifacts)
    specs={name:content['docs/literature/'+name] for name in ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')}
    flags={str(p.relative_to(build)):p.read_text() for p in [build/'CMakeFiles/IncNavierStokesSolver.dir/flags.make',build/'CMakeFiles/RANSWallResolutionProbe.dir/flags.make',build/'CMakeFiles/RANSKernels.dir/flags.make']}
    studyroot=root/'study-final' if (root/'study-final').exists() else root/'study'
    manifest=dict(schema='WJ3LinuxManifestV1',root=str(root),source_revision=git('rev-parse','HEAD'),source_tree=git('rev-parse','HEAD^{tree}'),
        source_status=git('status','--short'),source_content_sha256=sha(root/'source-content.json'),specification_sha256=specs,
        runtime=runtime,effective_flags=flags,cmake_cache_sha256=sha(build/'CMakeCache.txt'),
        dependencies=dict(python=sys.version,numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__,sympy=__import__('sympy').__version__),
        compiler=subprocess.check_output(['c++','--version'],text=True),platform=subprocess.check_output(['uname','-a'],text=True),
        allocation=dict(initial_build_workers=1,subsequent_build_workers=4,main_policy_concurrency=2,test_parallel=1,ranks=1,threads=THREADS,server_total_worker_limit=72),
        command_receipts=json.loads((root/'command-receipts.json').read_text()),command_receipts_sha256=sha(root/'command-receipts.json'),
        inventory=dict(path=str(root/'inventory.log'),sha256=sha(root/'inventory.log'),count=len(json.loads((root/'inventory.log').read_text())['tests'])),
        artifacts=dict(path=str(root/'artifact-index.json'),sha256=sha(root/'artifact-index.json'),count=len(artifacts)),
        summaries={p.name:dict(path=str(studyroot/p.name/'summary.json'),sha256=sha(studyroot/p.name/'summary.json')) for p in studyroot.iterdir() if (p/'summary.json').exists()})
    write_json(root/'manifest.json',manifest);print(root/'manifest.json')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--collect',action='store_true');p.add_argument('--label');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    if a.collect:collect(a)
    else:sys.exit(command(a))

if __name__=='__main__':main()
