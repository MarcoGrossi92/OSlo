# Turns on either OpenMP and/or MPI

IF (USE_OPENMP)
    # Find OpenMP
    IF (NOT OpenMP_Fortran_FLAGS)
        FIND_PACKAGE (OpenMP_Fortran)
        IF (NOT OpenMP_Fortran_FLAGS)
            MESSAGE (FATAL_ERROR "Fortran compiler does not support OpenMP")
        ENDIF (NOT OpenMP_Fortran_FLAGS)
    ENDIF (NOT OpenMP_Fortran_FLAGS)
ENDIF()

IF (USE_MPI)
    # Find MPI
    IF (NOT MPI_Fortran_FOUND)
        FIND_PACKAGE (MPI REQUIRED)
        message(STATUS "MPI enabled, using compilers: ${MPI_Fortran_COMPILER}")
    ENDIF (NOT MPI_Fortran_FOUND)
ENDIF()

IF (NOT USE_OPENMP AND NOT USE_MPI)
    # Turn off both OpenMP and MPI
    SET (OMP_NUM_PROCS 0 CACHE
         STRING "Number of processors OpenMP may use" FORCE)
    UNSET (OpenMP_Fortran_FLAGS CACHE)
    UNSET (GOMP_Fortran_LINK_FLAGS CACHE)
    UNSET (MPI_FOUND CACHE)
    UNSET (MPI_COMPILER CACHE)
    UNSET (MPI_LIBRARY CACHE)
ENDIF ()
