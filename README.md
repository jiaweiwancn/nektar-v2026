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

HIT inflow extension
--------------------
This branch implements three-component homogeneous isotropic turbulence inflow
for full 3D and Homogeneous1D velocity-correction calculations. Start with the
[HIT guide and test/data index](solvers/IncNavierStokesSolver/Tests/HITInflow/README.md),
the [full 3D example](solvers/IncNavierStokesSolver/Examples/HITInflow3D/README.md),
or the [Quasi-3D example](solvers/IncNavierStokesSolver/Examples/HITInflowQuasi3D/README.md).
The [comparison tools](solvers/IncNavierStokesSolver/Tests/HITInflow/pilots/README.md)
regenerate the retained experimental comparison. Its
[report](docs/plans/hit-inflow-milestone-7.md) documents a limited-band SVV result;
full-band physical validation and 2% precision are not established.
