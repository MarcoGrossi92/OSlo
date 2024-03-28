option(REQUIRE_FAST_LAPACK 
"If ON, then cmake configuration will 
error out if e.g. OpenBLAS is not found" OFF)

add_compile_definitions(FULL_ALGEBRA)

file(GLOB SRC_FILES *.F90)
add_library(fatode ${SRC_FILES})

if (${REQUIRE_FAST_LAPACK})
  find_package(LAPACK REQUIRED)
else()
  find_package(LAPACK)
endif()
if (${LAPACK_FOUND})
  target_link_libraries(fatode ${LAPACK_LIBRARIES})
else()
  # We compile lapack files downloaded from netlib
  file(GLOB LAPACK_SRC_FILES ../LAPACK/*.f)
  add_library(lapack ${LAPACK_SRC_FILES})
  target_link_libraries(fatode lapack)
endif()
