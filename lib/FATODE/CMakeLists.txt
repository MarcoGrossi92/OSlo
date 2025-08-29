# FATODE

file(GLOB FATODE_FILES *.F90)
add_compile_definitions(FULL_ALGEBRA)

if (${LAPACK_FOUND})
  add_library(fatode OBJECT ${FATODE_FILES})
else()
  # We compile lapack files downloaded from netlib
  add_library(fatode OBJECT ${FATODE_FILES} lapack/*.f)
endif()
target_include_directories(fatode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
IF(USE_OPENMP)
  SET_TARGET_PROPERTIES(fatode PROPERTIES
                        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                        LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF()
