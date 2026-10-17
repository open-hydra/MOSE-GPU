# Compile OSlo lib

file(GLOB HSRC hairer*)
add_library(OSlo OSlo.f90 dvodeOMP.f90 ${HSRC})
if(BLAS_FOUND AND LAPACK_FOUND)
  target_link_libraries(OSlo PRIVATE fatode ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
  target_link_libraries(OSlo PRIVATE fatode)
endif()
if (UNIX AND NOT APPLE)
  target_link_libraries(OSlo PRIVATE ${INTELODE})
endif()
target_include_directories(OSlo INTERFACE
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
if (USE_SUNDIALS)
  target_link_libraries(OSlo
    PRIVATE
      SUNDIALS::fcvodes_mod
      SUNDIALS::fidas_mod
      SUNDIALS::fcore_mod
      SUNDIALS::fsunnonlinsolnewton_mod
  )
endif()
IF(USE_OPENMP)
   SET_TARGET_PROPERTIES(OSlo PROPERTIES
                         COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                         LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF()

