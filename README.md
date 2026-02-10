# OSlo — A Unified Fortran ODE Solver Toolbox

![License](https://img.shields.io/github/license/MarcoGrossi92/OSlo)
![Fortran](https://img.shields.io/badge/Fortran-90%2B-blue)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)

**OSlo** is a high-performance, Fortran-based toolbox for solving systems of ordinary differential equations (ODEs).
It provides a **unified and modular interface** to multiple state-of-the-art ODE solvers, allowing users to switch between stiff and non-stiff integrators with minimal code changes.

OSlo is designed for **scientific computing and HPC environments**, bundling well-established solvers such as **DVODE**, **FATODE**, **Hairer–Wanner methods**, **Intel ODE**, and **SUNDIALS** under a single, consistent API.

---

## Key Features

* **Unified Solver Interface**
  Switch between different ODE solvers using a common Fortran API.

* **Broad Solver Coverage**

  * Hairer & Wanner solvers for stiff and non-stiff problems
  * **DVODE** (classic and reliable)
  * **FATODE** with automatic differentiation
  * **Intel ODE** (Linux only, Intel compiler)
  * **SUNDIALS** (CVODE, IDA, etc.)

* **Parallel Execution**

  * Shared-memory parallelism via **OpenMP**
  * Distributed-memory parallelism via **MPI**

* **Compiler Flexibility**

  * GNU Fortran (`gfortran`)
  * Intel Fortran (`ifort` / `ifx`)

* **Robust Linear Algebra**

  * LAPACK integration for dense and sparse systems

* **Modular Architecture**

  * Easy to extend or embed in existing Fortran codes

---

## Requirements

### Mandatory

* **CMake** ≥ 3.23
* **Fortran compiler**: GNU or Intel
* **LAPACK**
* **C compiler** (required for SUNDIALS)

### Optional

* **OpenMP** — shared-memory parallelism
* **MPI** — distributed-memory parallelism

---

## Platform Support

OSlo has been tested on the following configurations:

| OS             | Compiler       | Status |
| -------------- | -------------- | ------ |
| macOS          | GNU (gfortran) | ✅      |
| Ubuntu Linux   | GNU            | ✅      |
| Ubuntu Linux   | Intel          | ✅      |
| OpenSUSE Linux | GNU            | ✅      |
| OpenSUSE Linux | Intel          | ✅      |

The toolbox is expected to work on most Unix-like systems with a compatible Fortran compiler.

---

## Installation

### Quick Start

```bash
git clone https://github.com/MarcoGrossi92/OSlo.git
cd OSlo
./install.sh build --compiler=gnu
```

### Build Variants

```bash
./install.sh [GLOBAL_OPTIONS] build [OPTIONS]
```

**Global options**

* `-v, --verbose` — verbose output

**Build options**

* `--compiler=gnu|intel`
* `--use-openmp`
* `--use-mpi`
* `--use-sundials` (enabled by default)

**Examples**

```bash
# GNU compiler + OpenMP
./install.sh build --compiler=gnu --use-openmp

# Intel compiler + MPI
./install.sh build --compiler=intel --use-mpi
```

---

## Manual CMake Build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DUSE_MPI=OFF \
      -DUSE_SUNDIALS=ON \
      ..
make
make install
```

---

## Usage

OSlo builds a Fortran library that can be linked into your own applications.

### Minimal Example

```fortran
use oslo
implicit none

integer, parameter :: neq = 3
real(8) :: y(neq), t, tout
real(8) :: rtol(neq), atol(neq)
integer :: err

call setup_odesolver(N=neq, solver='dvodef90', RT=rtol, AT=atol)
call run_odesolver(neq, t, tout, y, Fgeneral, err)
```

### Linking Against OSlo

```bash
gfortran -c your_program.f90 -Ibuild/modules
gfortran your_program.o -Lbuild/lib -loslo -o your_program
```

---

## Third-Party Software

### Bundled

* **Hairer & Wanner Solvers**
  © University of Geneva — BSD-style license
  [https://www.unige.ch/~hairer/software.html](https://www.unige.ch/~hairer/software.html)

* **FATODE**
  Fast and Accurate ODE solver with AD
  [https://people.cs.vt.edu/~asandu/Software/FATODE/](https://people.cs.vt.edu/~asandu/Software/FATODE/)

### Optional

* **SUNDIALS**
  © Lawrence Livermore National Laboratory — BSD license
  [https://computing.llnl.gov/projects/sundials](https://computing.llnl.gov/projects/sundials)

* **Intel ODE**
  Intel-optimized ODE solvers (Linux only)

---

## Project Structure

```text
OSlo/
├── src/            # Core library and solvers
├── lib/            # Third-party solver sources
├── bin/            # Executables and tests
├── cmake/          # CMake modules
├── scripts/        # Helper scripts
├── install.sh      # Build helper
└── CMakeLists.txt
```

---

## Contributing

Contributions are welcome — bug fixes, solver additions, documentation, or performance improvements.

Typical workflow:

1. Fork the repo
2. Create a feature branch
3. Commit changes
4. Open a pull request

---

## Support

Please use the GitHub issue tracker for:

* bug reports
* feature requests
* solver questions

👉 [https://github.com/MarcoGrossi92/OSlo/issues](https://github.com/MarcoGrossi92/OSlo/issues)

---
