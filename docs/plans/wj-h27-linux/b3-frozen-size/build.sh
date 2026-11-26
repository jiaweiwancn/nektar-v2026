#!/usr/bin/env bash
set -euo pipefail
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source threads.sh
git -C source fetch "$r/source.bundle" HEAD
git -C source checkout --detach FETCH_HEAD
git -C source rev-parse HEAD > tested-source.txt
test -z "$(git -C source status --porcelain)"
ps -eo stat,ppid,nlwp,comm | awk '$1 !~ /^[ITZ]/ && ($4~/^(mpirun|cc1plus|IncNavier|WJH27|python|orted)/ || $4=="time") {n+=$3} END{print "tracked active threads before build",n+0; if(n>62)exit 86}'
free -h
/usr/bin/time -v -o configure.cost cmake -S "$r/source/solvers/IncNavierStokesSolver/Tests/RANS/h27" -B "$r/build" \
  -DCMAKE_BUILD_TYPE=Release "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
/usr/bin/time -v -o build.cost timeout --kill-after=15s 620s cmake --build build --parallel 8 \
  --target IncNavierStokesSolver WJH27SteadyKrylovProbe
sha256sum build/IncNavierStokesSolver build/WJH27SteadyKrylovProbe > executable.sha256
sed -n '/^CXX_DEFINES/p' build/CMakeFiles/WJNativeMixedFlow.dir/flags.make > source-fingerprint.txt
/usr/bin/time -v -o probe.cost timeout 55s build/WJH27SteadyKrylovProbe controls.json
test -z "$(git -C source status --porcelain)"
