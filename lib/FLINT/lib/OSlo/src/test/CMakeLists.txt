file(GLOB_RECURSE TESTS *.f90 *.F90)
set(TESTS ${TESTS} PARENT_SCOPE)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")


foreach(TEST ${TESTS})
  get_filename_component(EXE_NAME ${TEST} NAME_WE)
  add_executable(${EXE_NAME} ${TEST})
  add_test(${EXE_NAME} ${EXECUTABLE_OUTPUT_PATH}/${EXE_NAME} COMMAND ${EXE_NAME} -stress)
  target_link_libraries(${EXE_NAME} PRIVATE OSlo)
  target_include_directories(${EXE_NAME}
      INTERFACE
          $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  )
  if(USE_MPI)
    target_compile_options(${EXE_NAME} PRIVATE ${MPI_Fortran_COMPILE_FLAGS})
    target_link_options(${EXE_NAME} PRIVATE ${MPI_Fortran_LINK_FLAGS})
    target_include_directories(${EXE_NAME} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(${EXE_NAME} PRIVATE ${MPI_Fortran_LIBRARIES})
  endif()

  if(USE_OPENMP)
    target_compile_options(${EXE_NAME} PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(${EXE_NAME} PRIVATE ${OpenMP_Fortran_FLAGS})
  endif()

endforeach()