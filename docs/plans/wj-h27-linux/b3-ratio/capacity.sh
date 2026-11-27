#!/usr/bin/env bash
set -euo pipefail
snapshot=$1
requested=$2
ps -eLo tgid,tid,stat,comm,wchan:36,ppid > "$snapshot"
status=0
awk -v extra="$requested" '
$1 == $2 && $3 !~ /^[ITZ]/ && $4 ~ /^(bash|timeout|time|mpirun|IncNavier|simpleFoam|python)/ {
    has_child[$6]++
}
$3 !~ /^[ITZ]/ && $4 ~ /^(cc1plus|simpleFoam|mpirun|mpiexec|ninja|cmake|python|WJ|IncNavier|orted|f951|gcc|g\+\+|time)/ {
    if ($3 ~ /^S/ && $4 ~ /^(time|timer|timeout)$/) { wrappers++; next }
    # Count rank leaders, even when they wait. A launcher or daemon stopped
    # in MPI progress I/O is a service, as are the corresponding rank threads.
    if ($3 ~ /^S/ && $5 ~ /^(ep_poll|do_poll|do_select)(\.|$)/ &&
        ($4 ~ /^(mpirun|mpiexec|orted)$/ || ($1 != $2 && $4 ~ /^(simpleFoam|IncNavier)/))) {
        mpi_io++
        next
    }
    if ($1 == $2 && $3 ~ /^S/ && $4 ~ /^python/ && $5 ~ /^do_(wait|select)(\.|$)/) {
        possible_wrapper[$1] = $4
        next
    }
    by_command[$4]++
    workers++
}
END {
    # A single caller blocked on a separately counted child is orchestration.
    # Any library compute threads of that process are still counted above.
    for (pid in possible_wrapper) {
        if (has_child[pid]) { controllers++; continue }
        by_command[possible_wrapper[pid]]++
        workers++
    }
    for (name in by_command) print name, by_command[name]
    print "counted workers",workers+0,"requested",extra,"MPI poll services",mpi_io+0,"wait wrappers",wrappers+0,
          "child controllers",controllers+0
    if (workers+extra > 72) exit 86
}' "$snapshot" || status=$?
free -h
exit "$status"
