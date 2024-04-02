# Hairer solvers
file(GLOB HSRC hairer*)
add_library(hairer ${HSRC})

# OSlo
add_library(oslo OSlo.f90 dvodeOMP.f90)
target_link_libraries(oslo hairer fatode)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(oslo odepack_mod)
endif()
if (UNIX AND NOT APPLE)
  target_link_libraries(oslo ${INTELODE})
endif()
IF(USE_OPENMP)
  SET_TARGET_PROPERTIES(oslo PROPERTIES
                        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                        LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF()
