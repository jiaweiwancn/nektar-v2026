Nektar++
========
Nektar++ is an open-source software framework designed to support the
development of high-performance scalable solvers for partial differential
equations (PDEs) using the spectral/hp element method.

This package consists of a set of libraries (the framework) and a number of
pre-written PDE solvers for a selection of application domains.

The software and User Guide is available for download from
<http://www.nektar.info/>.


User Guide
----------
Detailed information on compiling, installing and using the software is
available in the User Guide. This document is available as a pre-compiled PDF
from the downloads section of the project website.

The [neutral-ABL URANS example](solvers/IncNavierStokesSolver/Examples/NeutralABL/README.md)
documents the standard and cubic CLS turbulence models, original/corrected
analytical wall functions, runnable cases, full-history restart and measured
validation limits on this development branch.
The [terrain-wall example](solvers/IncNavierStokesSolver/Examples/TerrainHill/README.md)
adds physical normal-coordinate wall strips, integrated strip continuity and
an upstream transition to original AWF, with paired hill/empty-floor pilots.
Its Cao–Tamura entry point supplies versioned experimental candidates, native
physical probes and flat-only calibration objectives; developed calibration
and experimental accuracy remain separate milestones in the
[calibrated-validation plan](docs/plans/cao-tamura-calibrated-validation.md).
Its [active CV3 protocol](docs/plans/cao-tamura-cv3-calibration-protocol.md)
now prioritizes the matched EP-AWF/transition/original-AWF floor route.
[CV3-R1 execution](docs/plans/cao-tamura-cv3-r1-reference-feasibility.md) records
reference recovery, developed standard precursors and unresolved CLS feasibility.

The [revised CLS and EP-AWF implementation plan](docs/plans/cls-epawf-revised-implementation-plan.md)
defines the new spectral/hp wall-element project.
[M0 contract and baseline](docs/plans/cls-epawf-revised-m0.md) freezes its
model/discretization contract, equation mapping and exact-base provenance.
[M1 kernel verification](docs/plans/cls-epawf-revised-m1.md) reproduces its
reference identities and constitutive failures on Linux.
[M2 spectral profiles and loads](docs/plans/cls-epawf-revised-m2.md) adds
metric-balanced cutoffs, all-mode source/flux checks and the specified inner-k
positivity counterexample, with 19 passing Linux tests.
[M3 native flat wall stage](docs/plans/cls-epawf-revised-m3.md) adds actual
subdomain fields, coupled pressure/transport, analytical AWF refinement and
full/condensed nonlinear stage verification, with 22 passing Linux tests.
The [revised wall-stage guide](solvers/IncNavierStokesSolver/Examples/RevisedEquilibriumWall/README.md)
documents its runnable sessions and bounded native CG backend.
[M4 constrained reference and native evolution](docs/plans/cls-epawf-revised-m4.md)
adds native IMEX1/2 startup/history, fixed free-load anchoring, full 3D flat
reference/rotation, MPI 1/2/4 and full-history restart. The complete M4
qualification and 25 other registered Linux regressions pass. The original
stability and positivity failures remain recorded physical gates.
[M5 conservative curved patch](docs/plans/cls-epawf-revised-m5.md)
adds certified analytic normal maps, variable-height interfaces, full Cartesian
momentum/scalar coupling, conservative flux and energy audits, native curved
manufactured cases, MPI ownership and restart checks. All 28 registered tests
are covered by passing Linux qualification, note checks and regression/replay
receipts. The original physical failure gates remain unchanged.
[M6 upstream transition and adapted baseline](docs/plans/cls-epawf-revised-m6.md)
adds a physical C2 one-way mask, independent state/admissibility gates, coupled
trace derivatives and constrained baseline-face support. Linux qualification,
orientation/operator replay and earlier regressions cover 31 registered tests,
including MPI 1/2/4 and full-history restart.
[M7 admissibility and perturbation diagnosis](docs/plans/cls-epawf-revised-m7.md)
adds accepted-state/transfer audits, independent energy/covariance accounting,
native perturbations, dt/p/mesh and matching-height studies, MPI and perturbed
restart. All 33 registered gates are covered by passing source-bound Linux
verification. The original three-component CLS transverse mode grows under
refinement despite positive sampled covariance, so robust-model qualification
remains failed. M7 is closed as a faithful-model diagnosis; no M8 remedy is
adopted and M9 remains open.


Tutorials
---------
A number of tutorials are available, designed to walk the user through the
basics of spectral/hp element methods, through the use of individual solvers and
performing specific types of calculations.

The tutorials are available from <http://doc.nektar.info/tutorials/latest>.


Pre-requisites
--------------
Nektar++ requires the following software to be installed on the users system:

- CMake
- BLAS/LAPACK

Additional software is also required. This can either be installed system-wide
or it can be downloaded and compiled automatically during the build process.

For more detailed information, please see the User Guide.


Compilation
-----------
On most UNIX-based systems a default compilation can be performed using the
following commands from the top-level of the source tree:

    mkdir build
    cd build
    cmake ..
    make

To alter the build configuration (for example, to enable parallel execution
support) we recommend using the `ccmake` command instead of `cmake`. 

For more detailed operating-system specific instructions, please see the
User Guide.


Installation
------------
The default installation location is in a `dist` subdirectory of the `build`
directory. This can be changed by setting the `CMAKE_INSTALL_PREFIX` option
using `ccmake`. To install the compiled libraries, solvers and header files, on
UNIX-based systems run:

    make install
