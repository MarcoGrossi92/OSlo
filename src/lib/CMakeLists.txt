file(GLOB DVODESRC dvodef90.f90)
add_library(DVODE ${DVODESRC})

file(GLOB HSRC *.f)
add_library(hairer ${HSRC})
