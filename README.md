# OSlo - ODE Solver Toolbox

![License](https://img.shields.io/github/license/MarcoGrossi92/OSlo)
![Fortran](https://img.shields.io/badge/Fortran-90%2B-blue)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)

OSlo is a comprehensive Fortran-based ODE (Ordinary Differential Equations) solver toolbox that provides integrated access to multiple state-of-the-art numerical solvers. It bundles FATODE, DVODE, Intel ODE solvers, and SUNDIALS, offering a unified interface for solving stiff and non-stiff ODEs in high-performance computing environments.

## Features

- **Multiple ODE Solvers**:
  - Hairer and Wanner solvers for stiff ODE problems
  - FATODE: Fast and Accurate ODE solver with automatic differentiation
  - Intel ODE: Intel's optimized ODE solvers (Linux only)
  - SUNDIALS: Suite of nonlinear and differential/algebraic equation solvers
- **Flexible Parallelization**: Support for OpenMP and MPI
- **Compiler Support**: Intel and GNU Fortran compilers
- **LAPACK Integration**: Linear algebra support for robust ODE solving
- **Modular Design**: Easy to integrate into existing scientific computing workflows

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Platform Support](#platform-support)
- [License](#license)
- [Third-Party Dependencies](#third-party-dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [Support](#support)

## Requirements

### System Requirements

- **CMake** >= 3.23
- **Fortran compiler**: Intel or GNU (gfortran)
- **C compiler**: For SUNDIALS support
- **LAPACK**: Linear algebra library

### Optional Requirements

- **MPI**: For distributed-memory parallelization
- **OpenMP**: For shared-memory parallelization

## Platform Support

OSlo has been successfully tested and verified to work on:

| Operating System | Fortran Compiler | Status |
|------------------|------------------|--------|
| macOS | GNU (gfortran) | ✅ Tested |
| Ubuntu (Linux) | GNU | ✅ Tested |
| Ubuntu (Linux) | Intel | ✅ Tested |
| OpenSUSE (Linux) | GNU | ✅ Tested |
| OpenSUSE (Linux) | Intel | ✅ Tested |

The project is designed to be portable across Unix-like systems with compatible Fortran compilers.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Third-Party Dependencies

This project includes adapted versions of the following third-party
libraries:

- **Hairer and Wanner ODE Solvers**  
  Used for high-order ODE and DAE integration  
  © 2004, University of Geneva (UNIGE)  
  Licensed under a BSD-style license  
  https://www.unige.ch/~hairer/software.html

- **FATODE**  
  Fast and Accurate ODE solver with automatic differentiation  
  https://people.cs.vt.edu/~asandu/Software/FATODE/index.html

This project optionally uses the following third-party libraries:

- **SUNDIALS**  
  Suite of nonlinear and differential/algebraic equation solvers  
  © Lawrence Livermore National Laboratory  
  Licensed under a BSD-style license  
  https://computing.llnl.gov/projects/sundials

- **Intel ODE**  
  Intel's optimized ODE solvers for Linux systems

## Installation

### Quick Start

Clone the repository and run the installation script:

```bash
git clone https://github.com/MarcoGrossi92/OSlo.git
cd OSlo
./install.sh build --compiler=gnu
```

### Build Options

The `install.sh` script provides several build options:

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global Options:**
- `-v, --verbose`: Enable verbose output

**Commands:**

- **build**: Perform a full build
  - `--compiler=<name>`: Set compiler suite (intel, gnu)
  - `--use-openmp`: Enable OpenMP parallelization
  - `--use-mpi`: Enable MPI parallelization
  - `--use-sundials`: Enable SUNDIALS solver (default: on)

**Examples:**

```bash
# Build with GNU compiler and OpenMP
./install.sh build --compiler=gnu --use-openmp

# Build with Intel compiler and MPI
./install.sh build --compiler=intel --use-mpi
```

### Manual CMake Build

If you prefer to use CMake directly:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DUSE_OPENMP=ON \
      -DUSE_SUNDIALS=ON \
      ..
make
make install
```

**CMake Options:**
- `USE_MPI` (ON/OFF): Enable MPI parallelization
- `USE_OPENMP` (ON/OFF): Enable OpenMP parallelization
- `USE_SUNDIALS` (ON/OFF): Enable SUNDIALS solver suite

### Running Tests

After a successful build, run the test suite:

```bash
./bin/test/robertson
```

## Usage

OSlo provides a Fortran library that can be linked to your applications. The library is built during installation in the `build/lib/` directory.

### Basic Example

```fortran
use oslo
implicit none
integer, parameter :: neq = 3
real(8) :: Y(neq), t, tout
real(8) :: RT(neq), AT(neq)

call setup_odesolver(N=neq, solver='dvodef90', RT=RT, AT=AT)
call run_odesolver(neq, t, tout, Y, Fgeneral, err)
```

### Integrating OSlo in Your Project

Link the OSlo library when compiling your Fortran code:

```bash
gfortran -c your_program.f90 -Ibuild/modules
gfortran your_program.o -Lbuild/lib -loslo -o your_program
```

## Project Structure

```
OSlo/
├── bin/                    # Compiled executables
├── build/                  # Build directory (created after build)
│   ├── lib/                # Compiled libraries
│   └── modules/            # Fortran module files
├── src/                    # Source code
│   ├── lib/                # Fortran library implementations
│   └── test/               # Test programs
├── lib/                    # External ODE solver packages
│   ├── FATODE/             # FATODE solver
│   ├── Intel-ODE/          # Intel ODE solvers (Linux only)
│   └── sundials/           # SUNDIALS suite
├── scripts/                # Utility scripts (testing, etc.)
├── cmake/                  # CMake modules and configuration
├── CMakeLists.txt          # Main CMake configuration
├── install.sh              # Installation script
└── LICENSE                 # GNU General Public License v3.0
```

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Support

For issues, feature requests, or questions, please open an issue on the [GitHub repository](https://github.com/MarcoGrossi92/OSlo/issues). 
