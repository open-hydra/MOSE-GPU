# FATODE

file(GLOB FATODE_FILES *.F90)
add_compile_definitions(FULL_ALGEBRA)

if (${LAPACK_FOUND})
  add_library(fatode STATIC ${FATODE_FILES})
else()
  # We compile lapack files downloaded from netlib
  file(GLOB LAPACK_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lapack/*.f")
  list(APPEND FATODE_FILES ${LAPACK_FILES})
  add_library(fatode STATIC ${FATODE_FILES})
endif()
IF(USE_OPENMP)
  SET_TARGET_PROPERTIES(fatode PROPERTIES
                        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                        LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF()
