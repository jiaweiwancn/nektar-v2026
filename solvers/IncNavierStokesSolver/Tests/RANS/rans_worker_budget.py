# SPDX-License-Identifier: MIT
"""Reserve complete serial native process families, including startup and MPI helpers."""
import fcntl,json,subprocess


def spawn(root, command, **kwargs):
    with (root/'worker-launch.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        processes={}
        for line in subprocess.check_output(['ps','-eo','pid,ppid,nlwp,comm,args'],text=True).splitlines()[1:]:
            p=line.split(None,4)
            if len(p)==5:processes[int(p[0])]=(int(p[1]),int(p[2]),p[3],p[4])
        path=root/'worker-leases.json';leases=json.loads(path.read_text()) if path.exists() else []
        leases=[p for p in leases if p in processes]
        covered=set(leases)
        while True:
            more={pid for pid,(parent,*_) in processes.items() if parent in covered}-covered
            if not more:break
            covered.update(more)
        other=sum(threads for pid,(parent,threads,name,args) in processes.items() if pid not in covered and (name.startswith(('IncNavier','RANS','Nek','cc1plus','f951','orted','mpiexec')) or (name.startswith('python') and 'nektar' in args)))
        # A runner/probe family needs <=8 threads (Python + solver + ORTE).
        # Keep sixteen threads unallocated for unrelated campaign startup.
        booked=other+8*len(leases)
        if booked+8+16>72:
            path.write_text(json.dumps(leases)+'\n');return None
        process=subprocess.Popen(command,**kwargs);leases.append(process.pid);path.write_text(json.dumps(leases)+'\n')
        with (root/'worker-bookings.jsonl').open('a') as out:out.write(json.dumps(dict(pid=process.pid,external_threads=other,reserved_families=len(leases),booked_threads=other+8*len(leases),unallocated_reserve=16,limit=72))+'\n')
        return process
