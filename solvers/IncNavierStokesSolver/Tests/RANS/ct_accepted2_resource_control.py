#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Corrective private-allocation control; does not erase any historical peak."""
import argparse,os,subprocess,time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_open import run,require
from ct_accepted_feasibility import write_json
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();a.output.mkdir(exist_ok=False);source=r/'source-after-reboot';peak=0
    with (a.output/'git-status.log').open('x') as f:
        process=subprocess.Popen(['git','-c','index.threads=1','-c','pack.threads=1','status','--porcelain'],cwd=source,stdout=f,stderr=subprocess.STDOUT)
        while process.poll() is None:
            try:
                for line in (Path('/proc')/str(process.pid)/'status').read_text().splitlines():
                    if line.startswith('Threads:'):peak=max(peak,int(line.split()[1]))
            except FileNotFoundError:pass
            time.sleep(.005)
        require(process.returncode==0 and peak<=1,'bounded_git_bookkeeping')
    # Force this one small existing native test target to compile with real
    # compiler children, rather than accepting a no-op incremental build.
    e=run(['cmake','--build',r/'reporting/build','--parallel','2','--target','RANSNativeDiagnosticChecks','--','-B'],a.output,'two-job-build.log')
    require(e['exit_code']==0,'two_job_corrective_build_failed')
    result=dict(status='pass',git_peak_threads=peak,git_options=['index.threads=1','pack.threads=1'],compiler_jobs=2,build=e,history='original 4-job/Git helper allocation overruns remain failed; hard aggregate72 gate evaluated separately')
    write_json(a.output/'results.json',result);print('corrective build control pass');return 0
if __name__=='__main__':raise SystemExit(main())
