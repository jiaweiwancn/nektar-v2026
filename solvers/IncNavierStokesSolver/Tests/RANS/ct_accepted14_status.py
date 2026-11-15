#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only CT14 progress plus append-only owned-native memory observations."""
import argparse,hashlib,json,os,sys,time
from pathlib import Path

def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--pair',default='pair-v2');p.add_argument('--control',default='resources-pair-v2');a=p.parse_args()
    tests=a.root/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(tests))
    from ct_accepted14_budget import processes,family
    control=a.root/a.control;pair=a.root/a.pair;out=dict(members={},whole_elapsed_wall=time.time()-1794751358);rows=processes();owned=set();live=0.
    if not (control/'execution.json').exists():
        leases=read(control/'worker-leases.json');outer=rows[leases[0]]['parent'];owned=family(rows,outer);live=sum(rows[pid]['ticks'] for pid in owned if pid in rows)/os.sysconf('SC_CLK_TCK')
    out['sampled_package_CPU']=live+sum(read(path)['process_cpu_seconds'] for path in a.root.glob('resources-*/execution.json'))
    for dt in ('0.00025','0.000125'):
        folder=pair/dt;row={};path=folder/'sparse-receipts.jsonl'
        if path.exists():
            receipts=[]
            for line in path.read_text().splitlines():
                try:receipts.append(json.loads(line))
                except json.JSONDecodeError:pass # a concurrently appended last line supplies no claim
            fills=[r['factor_nnz'] for r in receipts if r.get('event')=='factor' and r.get('ordering')=='VectorCOLAMDV1'];row.update(predictors_produced=len(fills),maximum_predictor_fill=max(fills,default=0))
        path=folder/'case.evolution.json'
        if path.exists():
            q=read(path);row['native']=dict(time=q['time'],accepted_steps=q['steps'],failure=q['failure'])
        path=folder/'budget-progress.json'
        if path.exists():row['budgets']={k:{n:v.get(n) for n in ('reconstruction_status','accuracy_status','reason')} for k,v in read(path)['budgets'].items()}
        if row:out['members'][dt]=row
    path=pair/'results.json'
    if path.exists():
        q=read(path);out.update(qualification=q['qualification'],pressure=q.get('temporal',{}).get('Pi'))
    observations=control/'live-memory-observations.jsonl';prior=[json.loads(s) for s in observations.read_text().splitlines()] if observations.exists() else [];maximum={(r['pid'],r['process_start_ticks']):r['VmPeak_KiB'] for r in prior};added=[]
    for pid in sorted(owned):
        proc=Path('/proc')/str(pid)
        try:
            args=[s.decode() for s in (proc/'cmdline').read_bytes().split(b'\0') if s]
            if not args or Path(args[0]).name not in ('RANSPreservationProbe','IncNavierStokesSolver'):continue
            status={s.split(':',1)[0]:s.split(':',1)[1].strip() for s in (proc/'status').read_text().splitlines() if ':' in s};vm=int(status['VmPeak'].split()[0]);key=(pid,rows[pid]['start'])
            if vm<=maximum.get(key,-1):continue
            added.append(dict(time_unix=time.time(),pid=pid,process_start_ticks=key[1],command=args,executable_sha256=sha(Path(args[0])),VmPeak_KiB=vm,VmRSS_KiB=int(status['VmRSS'].split()[0]),OS_threads=int(status['Threads']),inspector_sha256=sha(Path(__file__))))
        except (FileNotFoundError,ProcessLookupError,PermissionError,KeyError,ValueError):pass
    if added:
        with observations.open('a') as stream:
            for row in added:stream.write(json.dumps(row)+'\n')
        all_rows=prior+added;winner=max(all_rows,key=lambda r:r['VmPeak_KiB']);summary=dict(schema='CT14LiveNativeMemoryV1',maximum_observed_address_space_KiB=winner['VmPeak_KiB'],witness=winner,samples=len(all_rows),observations_sha256=sha(observations),scope='actual /proc VmPeak observations of owned native processes; not a claim to sample their final peak; hard16GiB RLIMIT_AS remains enforced by the numerical controller')
        (control/'live-memory-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if (control/'live-memory-summary.json').exists():out['memory']=read(control/'live-memory-summary.json')
    print(json.dumps(out),flush=True)

if __name__=='__main__':main()
