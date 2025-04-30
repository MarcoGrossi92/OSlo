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

# FATODE
add_subdirectory(lib/FATODE)

# SUNDIALS
set(SUNDIALS_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sundials-install")
# Check if cloned
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials/CMakeLists.txt")
    message(STATUS "Sundials not found. Attempting to initialize...")
    execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to initialize git submodules.")
    endif()
endif()
# Check if already installed
if(NOT EXISTS "${SUNDIALS_INSTALL_DIR}/lib/cmake/sundials/SUNDIALSConfig.cmake")
    message(STATUS "SUNDIALS not installed, building and installing it...")

    # Build SUNDIALS in a temporary build directory
    execute_process(COMMAND ${CMAKE_COMMAND} -S ${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials
                                           -B ${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials/build
                                           -DCMAKE_INSTALL_PREFIX=${SUNDIALS_INSTALL_DIR}
                                           -DCMAKE_BUILD_TYPE=Release
                                           -DBUILD_SHARED_LIBS=OFF
                                           -DSUNDIALS_BUILD_WITH_MONITORING=OFF
                                           -DBUILD_FORTRAN_MODULE_INTERFACE=ON 
                                           -DEXAMPLES_ENABLE_C=OFF
                                           -DEXAMPLES_ENABLE_CXX=OFF
                                           -DEXAMPLES_ENABLE_F2003=ON
                    RESULT_VARIABLE CONFIGURE_RESULT)
    if(NOT CONFIGURE_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to configure SUNDIALS")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials/build --target install
                    RESULT_VARIABLE BUILD_RESULT)
    if(NOT BUILD_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to build and install SUNDIALS")
    endif()
endif()
find_package(SUNDIALS REQUIRED
             PATHS ${SUNDIALS_INSTALL_DIR}
             NO_DEFAULT_PATH)

# Build library
add_subdirectory(src/lib)

# Build tests and applications (only if OSlo build has not been called by another program)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message("-- Full build allowed")
    add_subdirectory(src/test)
endif()

