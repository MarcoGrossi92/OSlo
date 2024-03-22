file(GLOB_RECURSE TESTS *.f90 *.F90)
set(TESTS ${TESTS} PARENT_SCOPE)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")


foreach(TEST ${TESTS})
  get_filename_component(EXE_NAME ${TEST} NAME_WE)
  add_executable(${EXE_NAME} ${TEST})
  add_test(${EXE_NAME} ${EXECUTABLE_OUTPUT_PATH}/${EXE_NAME} COMMAND ${EXE_NAME} -stress)
  target_link_libraries(${EXE_NAME} oslo)
  target_include_directories(${EXE_NAME}
      INTERFACE
          $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  )
  IF(USE_OPENMP)
      SET_TARGET_PROPERTIES(${EXE_NAME} PROPERTIES
                            COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                            LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
  ENDIF(USE_OPENMP)
endforeach()