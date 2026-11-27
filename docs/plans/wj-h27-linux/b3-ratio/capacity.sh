#!/usr/bin/env bash
set -euo pipefail
snapshot=$1
requested=$2
ps -eLo tgid,tid,stat,comm,wchan:36 > "$snapshot"
status=0
awk -v extra="$requested" '
$3 !~ /^[ITZ]/ && $4 ~ /^(cc1plus|simpleFoam|mpirun|mpiexec|ninja|cmake|python|WJ|IncNavier|orted|f951|gcc|g\+\+|time)/ {
    if ($3 ~ /^S/ && $4 ~ /^(time|timer|timeout)$/) { wrappers++; next }
    if ($1 != $2 && $4 ~ /^(simpleFoam|mpirun|mpiexec|IncNavier|orted)/ &&
        $5 ~ /^(ep_poll|do_poll|do_select)(\.|$)/) { mpi_io++; next }
    by_command[$4]++
    workers++
}
END {
    for (name in by_command) print name, by_command[name]
    print "counted workers",workers+0,"requested",extra,"MPI poll services",mpi_io+0,"wait wrappers",wrappers+0
    if (workers+extra > 72) exit 86
}' "$snapshot" || status=$?
free -h
exit "$status"
