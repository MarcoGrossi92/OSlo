# Hairer solvers
file(GLOB HSRC hairerToolbox.f radau5.f radau.f rodas.f sdirk4.f)
add_library(hairer ${HSRC})

# OSlo
add_library(oslo OSlo.f90 dvode.f dvodeOMP.f90)
target_link_libraries(oslo odepack_mod hairer)
if(UNIX AND NOT APPLE)
  target_link_libraries(oslo ${INTELODE})
endif()
