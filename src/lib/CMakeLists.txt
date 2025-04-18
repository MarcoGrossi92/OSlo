# Compile OSlo lib

file(GLOB HSRC hairer*)
add_library(OSlo OSlo.f90 dvodeOMP.f90 ${HSRC})
target_link_libraries(OSlo fatode)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(OSlo odepack_mod)
endif()
if (UNIX AND NOT APPLE)
  target_link_libraries(OSlo ${INTELODE})
endif()

