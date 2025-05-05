# Compile OSlo lib

file(GLOB HSRC hairer*)
add_library(OSlo OSlo.f90 dvodeOMP.f90 ${HSRC})
target_link_libraries(OSlo PRIVATE fatode)
if (UNIX AND NOT APPLE)
  target_link_libraries(OSlo PRIVATE ${INTELODE})
endif()
target_link_libraries(OSlo
  PRIVATE
    SUNDIALS::fcvode_mod
    SUNDIALS::fcore_mod
    SUNDIALS::farkode_mod
    SUNDIALS::fsunnonlinsolnewton_mod
)
IF(USE_OPENMP)
   SET_TARGET_PROPERTIES(OSlo PROPERTIES
                         COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                         LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF()

