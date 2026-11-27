#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
stage=$1
parallel=$2
[[ "$parallel" =~ ^([1-9]|1[0-6])$ ]]
target="$r/bin/$stage"
mkdir -p "$target"
latest=4300
if test "$stage" = momentum; then latest=7050; fi
while ! bash capacity.sh "$target/workers.ps" "$((parallel+2))"; do
    test "$(($(date +%s)-$(cat start.epoch)))" -lt "$latest"
    sleep 20
done
date +%s > "$target/start.epoch"
sha256sum "$stage.bundle" > "$target/bundle.sha256"
git -C source fetch "$r/$stage.bundle" HEAD
git -C source checkout --detach FETCH_HEAD
git -C source rev-parse HEAD > "$target/source.txt"
test -z "$(git -C source status --porcelain)"
/usr/bin/time -v -o "$target/configure.cost" cmake -S source/solvers/IncNavierStokesSolver/Tests/RANS/h27 -B build \
    -DCMAKE_BUILD_TYPE=Release "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
/usr/bin/time -v -o "$target/build.cost" timeout --kill-after=15s 980s cmake --build build --parallel "$parallel" \
    --target IncNavierStokesSolver WJH27SteadyKrylovProbe
cp build/IncNavierStokesSolver build/WJH27SteadyKrylovProbe "$target/"
(cd "$target" && sha256sum IncNavierStokesSolver WJH27SteadyKrylovProbe > executable.sha256)
sed -n '/^CXX_DEFINES/p' build/CMakeFiles/WJNativeMixedFlow.dir/flags.make > "$target/source-fingerprint.txt"
/usr/bin/time -v -o "$target/probe.cost" timeout 60s "$target/WJH27SteadyKrylovProbe" "$target/controls.json"
date +%s > "$target/finish.epoch"
test -z "$(git -C source status --porcelain)"
