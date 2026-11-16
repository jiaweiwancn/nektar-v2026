#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only actual-row contraction prototype, never a native speedup claim."""
import argparse, ctypes, json, os, resource, subprocess, time
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest, sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ct14', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_AS, (16*1024**3, 16*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (100, 105))
    source = Path(__file__).with_name('CT16PointContractionBench.cpp')
    header = Path(__file__).parents[2]/'EquationSystems/RANSWallSharedVector.h'
    initial = a.ct14/'static-v1/target/frozen.json'
    members = [a.ct14/'pair-v2'/dt/'case.evolution.json' for dt in ('0.00025', '0.000125')]
    m = RequestManifest(a.output, dict(schema='CT16ActualRowContractionPrototypeV1',
        source_sha256=sha(source), controller_sha256=sha(Path(__file__)), shared_vector_header_sha256=sha(header),
        initial_sha256=sha(initial), member_sha256=[sha(f) for f in members],
        maximum_rows=65536, deterministic_sample_per_inventory=4096, maximum_coefficients=128,
        repetitions=16, timed_pairs_per_state=3, physical_advance=False,
        qualification_scope='Sampled read-only primitive exactness/throughput; not governing or whole-workflow admission'))
    libfile = a.output/'contraction.so'
    command = ['c++', '-std=c++17', '-O3', '-DNDEBUG', '-shared', '-fPIC', '-I'+str(Path(__file__).parents[3]), str(source), '-o', str(libfile)]
    tick = time.monotonic(); child = resource.getrusage(resource.RUSAGE_CHILDREN)
    compiled = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (a.output/'compiler.log').write_text(compiled.stdout)
    require(compiled.returncode == 0, 'prototype compilation failed')
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    build = dict(command=command, CPU=after.ru_utime+after.ru_stime-child.ru_utime-child.ru_stime, wall=time.monotonic()-tick,
                 compiler=subprocess.check_output(['c++', '--version'], text=True), binary_sha256=sha(libfile))
    packet = json.loads(initial.read_text()); o = packet['observation']; offsets = o['offsets']
    queries = []
    def append(rows, regions):
        indices = np.linspace(0, len(rows)-1, min(4096, len(rows)), dtype=int) if rows else []
        for i in indices:
            q = rows[i]; region = regions[i]
            for name, block in [('velocity_k', 0), ('velocity_k', 1), ('velocity_k', 2), ('epsilon', 3), ('pressure', 4+region)]:
                row = q.get(name)
                if row is not None:
                    require(0 < len(row['ids']) <= 128 and 'constant' in row, 'actual row lacks bounded native constant representation')
                    queries.append((row, offsets[block]))
    append(packet['test_basis'], o['grid_regions'])
    append(o['volume'], [q['region'] for q in o['volume']])
    append(o['faces'], [q['region'] for q in o['faces']])
    require(0 < len(queries) <= 65536, 'actual query count exceeds prospective bound')
    ptr = np.zeros(len(queries)+1, dtype=np.uint64)
    ptr[1:] = np.cumsum([len(r['ids']) for r, _ in queries])
    ids = np.asarray([j for r, _ in queries for j in r['ids']], np.int32)
    fields = [np.asarray([x for r, _ in queries for x in r[key]], np.float64) for key in ('phi', 'constant')]
    fields += [np.asarray([x for r, _ in queries for x in r['derivative'][d]], np.float64) for d in range(3)]
    off = np.asarray([b for _, b in queries], np.int32)
    lib = ctypes.CDLL(str(libfile)); lib.ct16_contraction_error.restype = ctypes.c_char_p
    lib.ct16_contraction_make.restype = ctypes.c_void_p
    lib.ct16_contraction_make.argtypes = [ctypes.c_size_t]+[ctypes.c_void_p]*8
    lib.ct16_contraction_run.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
    lib.ct16_contraction_free.argtypes = [ctypes.c_void_p]
    arrays = [ptr, ids, *fields, off]
    handle = lib.ct16_contraction_make(len(queries), *[v.ctypes.data for v in arrays])
    require(handle, 'prototype native row construction:'+str(lib.ct16_contraction_error()))
    states = [o['state']]+[json.loads(f.read_text())['state'] for f in members]
    out = dict(schema='CT16ActualRowContractionPrototypeV1', status='incomplete', build=build,
               rows=len(queries), coefficients=int(ptr[-1]), states=[], physical_advance=False)
    try:
        for index, values in enumerate(states):
            state = np.asarray(values, np.float64); results = [np.zeros((len(queries), 4)) for _ in range(2)]
            def run(policy, repeats):
                code = lib.ct16_contraction_run(handle, len(state), state.ctypes.data, repeats, policy, results[policy].ctypes.data)
                require(code == 0, 'prototype contraction:'+str(lib.ct16_contraction_error()))
            run(0, 1); run(1, 1)
            require(np.array_equal(results[0].view(np.uint64), results[1].view(np.uint64)), 'packed contractions differ in binary64 bits')
            times = {str(i): [] for i in (0, 1)}
            for repeat in range(3):
                for policy in ((0, 1) if repeat % 2 == 0 else (1, 0)):
                    cpu = time.process_time(); wall = time.monotonic(); run(policy, 16)
                    times[str(policy)].append(dict(CPU=time.process_time()-cpu, wall=time.monotonic()-wall))
            require(np.array_equal(results[0].view(np.uint64), results[1].view(np.uint64)), 'timed contractions changed output')
            out['states'].append(dict(index=index, binary64_equal=True, maximum_difference=0., timing=times))
            write_json(a.output/'progress.json', out)
    finally:
        lib.ct16_contraction_free(handle)
    out['status'] = 'pass'; out['interpretation'] = 'Read-only prototype exactness and sampled primitive timings only. No new governing runtime,accepted state,full action/history qualification or production admission.'
    write_json(a.output/'results.json', out); m.seal(a.output, ['results.json'])
    print('CT16 primitive contractions pass', len(queries), flush=True)


if __name__ == '__main__':
    main()
