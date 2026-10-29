# Project purpose

This project implements, tests, and applies new features and numerical methods in
Nektar++. These instructions apply throughout this repository.

## Development and validation environments

- Use this checkout, `/Users/jay/Documents/nektar-v2026`, for source changes,
  code review, and Git-based work.
- Treat the Linux server checkout at `/media/jay/新加卷/nektar-v2026` as the
  authoritative environment for compiling Nektar++, running solvers, numerical
  validation, and heavier tests.
- Connect from the Mac with `ssh linux-usb` (account `jay` on
  `jay-Super-Server`). If `ssh linux-usb` is unavailable, connect over Tailscale
  with `ssh jay@100.105.247.110`.
- Both the Mac and Linux server support GitHub SSH authentication.
- Synchronize the intended source changes to Linux before validation. Check the
  revision and working-tree state in both checkouts, and preserve any existing
  uncommitted work when synchronizing.

## Linux build and installation

The existing Release installation has MPI enabled with `NEKTAR_USE_MPI=ON`:

- Source: `/media/jay/新加卷/nektar-v2026`
- Build directory: `/home/jay/build/nektar-v2026-mpi`
- Installation: `/home/jay/.local/opt/nektar-v2026-mpi`
- Shell environment: `/home/jay/.local/opt/nektar-v2026-mpi/env.sh`

Activate the environment in an existing Linux shell with:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
```

Rebuild and reinstall on Linux after source changes as needed. Select tests that
exercise the changed behavior, and use Linux results for build, solver, numerical,
and performance conclusions. Report the tested revision, commands, and relevant
MPI/thread configuration.

## Heavy test storage

- Store heavier or storage-consuming tests and their generated data under
  `/media/jay/新加卷` on the Linux server.
- Use a separate working directory such as
  `/media/jay/新加卷/nektar-v2026-runs/<experiment>` for each experiment.
- Keep large meshes, field files, checkpoints, restart data, logs, and results in
  these run directories, outside the source checkout. Direct large temporary
  files there as well, including through `TMPDIR` when applicable.
- Keep small, reproducible regression-test definitions and necessary fixtures in
  the repository when appropriate.

## Linux worker limit

- **Never employ more than 72 workers in total on the Linux server.** This is a
  hard limit even though the server exposes more CPUs.
- Apply the limit to compiler jobs, MPI ranks, test-runner concurrency, and other
  parallel or background work. Account for concurrent jobs and threads within
  each MPI rank so their combined worker count stays at or below 72.
- Set parallelism explicitly. Avoid unbounded build parallelism or automatically
  using all CPUs reported by `nproc` or similar tools.
- For MPI jobs, default to `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.
  Use a hybrid MPI/thread configuration only with an explicit allocation that
  respects the same total limit.
- Choose parallel test-runner settings according to the MPI ranks and threads
  used by all tests that can run concurrently.

For example, this rebuild uses 32 compiler workers and must share the 72-worker
budget with any other concurrent work:

```bash
cmake --build /home/jay/build/nektar-v2026-mpi --parallel 32
cmake --install /home/jay/build/nektar-v2026-mpi
```
