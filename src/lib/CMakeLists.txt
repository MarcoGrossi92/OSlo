# DVODE
add_library(dvode dvode.f)

# DVODE F90 OpenMP
add_library(dvodef90OMP dvode_OpenMP.f90)

# Hairer solvers
file(GLOB HSRC libradau.f radau5.f radau.f rodas.f sdirk4.f)
add_library(hairer ${HSRC})
