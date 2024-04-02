option(REQUIRE_FAST_LAPACK 
"If ON, then cmake configuration will 
error out if e.g. OpenBLAS is not found" OFF)

if (${REQUIRE_FAST_LAPACK})
  find_package(LAPACK REQUIRED)
endif()

file(GLOB FATODE_FILES *.F90)
add_compile_definitions(FULL_ALGEBRA)

if (${LAPACK_FOUND})
  add_library(fatode ${FATODE_FILES})
  target_link_libraries(fatode ${LAPACK_LIBRARIES})
else()
  # We compile lapack files downloaded from netlib
  file(GLOB FATODE_LAPACK_FILES lapack/*.f)
  add_library(fatode ${FATODE_FILES} ${FATODE_LAPACK_FILES})
endif()
