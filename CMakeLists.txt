cmake_minimum_required(VERSION 3.13)

project (OSlo LANGUAGES Fortran)

set(CMAKE_VERBOSE_MAKEFILE FALSE)
# seach path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING. You should review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries. There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

# Set directories
include(GNUInstallDirs)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Define third-party libs
# Intel-ODE
if(UNIX AND NOT APPLE)
  add_compile_definitions(INTEL)
  file(GLOB INTELODE "${CMAKE_CURRENT_SOURCE_DIR}/lib/Intel-ODE/lib/intel64/*")
endif()
# ODEPACK
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set ( odepackD ${CMAKE_CURRENT_SOURCE_DIR}/lib/odepack )
  add_subdirectory(${odepackD})
endif()
# FATODE
add_subdirectory(lib/FATODE)

# SUNDIALS
set(SUNDIALS_DIR "lib/sundials/build")
find_package(SUNDIALS REQUIRED
             PATHS ${SUNDIALS_DIR}
             NO_DEFAULT_PATH)

# Build library
add_subdirectory(src/lib)

# Build tests and applications (only if NerO build has not been called by another program)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message("-- Full build allowed")
    add_subdirectory(src/test)
endif()

