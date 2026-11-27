#!/usr/bin/env bash
set -euo pipefail
snapshot=$1
limit=$2
# Include one thread per rank plus all potential computational workers.
# Exclude only observed MPI I/O/progress threads blocked in poll/epoll.
ps -eLo tgid,tid,stat,comm,wchan:36 > "$snapshot"
awk -v cap="$limit" '
$3 !~ /^[ITZ]/ && $4 ~ /^(cc1plus|simpleFoam|mpirun|mpiexec|ninja|cmake|python|WJ|IncNavier|orted|f951|gcc|g\+\+|time)/ {
    if ($1 != $2 && $4 ~ /^(simpleFoam|mpirun|mpiexec|IncNavier|orted)/ &&
        $5 ~ /^(ep_poll|do_poll|do_select)(\.|$)/) { mpi_io++; next }
    by_command[$4]++
    workers++
}
END {
    for (name in by_command) print name, by_command[name]
    print "counted workers",workers+0,"MPI poll services",mpi_io+0,"limit",cap
    if (workers > cap) exit 86
}' "$snapshot"
free -h
