# Hairer solvers
file(GLOB HSRC hairer*)
add_library(hairer ${HSRC})

# OSlo
add_library(oslo OSlo.f90 dvode.f dvodeOMP.f90)
target_link_libraries(oslo odepack_mod hairer fatode)
if(UNIX AND NOT APPLE)
  target_link_libraries(oslo ${INTELODE})
endif()
