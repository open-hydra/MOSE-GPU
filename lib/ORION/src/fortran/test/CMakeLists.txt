file(GLOB_RECURSE TESTS *.f90 *.F90)
set(TESTS ${TESTS} PARENT_SCOPE)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")

foreach(TEST ${TESTS})
  get_filename_component(EXE_NAME ${TEST} NAME_WE)
  add_executable(${EXE_NAME} ${TEST})
  #set_target_properties(${EXE_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(${EXE_NAME} ORION)
  if (USE_TECIO)
    target_link_libraries(${EXE_NAME} ${LINKLIBS})
  endif()
  target_include_directories(${EXE_NAME}
      INTERFACE
          $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  )
endforeach()
