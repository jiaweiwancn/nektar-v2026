# SPDX-License-Identifier: MIT
"""Same native work/limits/memory sampling; avoid statting known mappings repeatedly."""
import os,resource,subprocess,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from rans_request_manifest import sha

def native_run(command,folder,log,cpu,wall):
    """Bound the native process itself, observe its actual maps, always reap it."""
    observer_cpu=time.process_time();started=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);maps={};rss=vm=0;stop=None
    def limits():
        resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu+5))
    with (folder/log).open('x') as f:
        child=subprocess.Popen([str(x) for x in command],cwd=folder,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,**THREADS),preexec_fn=limits)
        try:
            while child.poll() is None:
                proc=Path('/proc')/str(child.pid)
                try:
                    for line in (proc/'maps').read_text().splitlines():
                        word=line.split()[-1]
                        if word.startswith('/') and word not in maps:
                            path=Path(word)
                            if path.is_file():maps[word]=sha(path)
                    for line in (proc/'status').read_text().splitlines():
                        if line.startswith('VmHWM:'):rss=max(rss,int(line.split()[1]))
                        if line.startswith('VmPeak:'):vm=max(vm,int(line.split()[1]))
                except (FileNotFoundError,ProcessLookupError):pass
                if time.monotonic()-started>wall:stop='wall_work_limit';child.terminate();break
                time.sleep(.05)
        finally:
            try:child.wait(timeout=10)
            except subprocess.TimeoutExpired:child.kill();child.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(observer_policy='OncePerNewMappedPathV1',observer_self_CPU_seconds=time.process_time()-observer_cpu,command=list(map(str,command)),exit_code=child.returncode,wall_seconds=time.monotonic()-started,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,peak_RSS_KiB=rss,peak_address_space_KiB=vm,actual_mapped_files=maps,CPU_soft=cpu,CPU_hard=cpu+5,wall_limit=wall,stop_reason=stop)
