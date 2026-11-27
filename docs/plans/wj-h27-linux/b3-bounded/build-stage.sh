#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")" && pwd)"
cd "$r"
source "/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace/env.sh"
source threads.sh
stage=$1
bundle="$r/$stage.bundle"
target="$r/bin/$stage"
mkdir -p "$target"
bash workers.sh "$target/workers.ps" 54
date +%s > "$target/start.epoch"
sha256sum "$bundle" > "$target/bundle.sha256"
git -C source fetch "$bundle" HEAD
git -C source checkout --detach FETCH_HEAD
git -C source rev-parse HEAD > "$target/source.txt"
test -z "$(git -C source status --porcelain)"
/usr/bin/time -v -o "$target/configure.cost" cmake -S source/solvers/IncNavierStokesSolver/Tests/RANS/h27 -B build \
  -DCMAKE_BUILD_TYPE=Release "-DNektar++_DIR=$NEKTAR_INSTALL/lib64/nektar++/cmake"
/usr/bin/time -v -o "$target/build.cost" timeout --kill-after=15s 760s cmake --build build --parallel 16 \
  --target IncNavierStokesSolver WJH27SteadyKrylovProbe
cp build/IncNavierStokesSolver build/WJH27SteadyKrylovProbe "$target/"
(cd "$target" && sha256sum IncNavierStokesSolver WJH27SteadyKrylovProbe > executable.sha256)
sed -n '/^CXX_DEFINES/p' build/CMakeFiles/WJNativeMixedFlow.dir/flags.make > "$target/source-fingerprint.txt"
/usr/bin/time -v -o "$target/probe.cost" timeout 60s "$target/WJH27SteadyKrylovProbe" "$target/controls.json"
date +%s > "$target/finish.epoch"
test -z "$(git -C source status --porcelain)"
