cmake_minimum_required(VERSION 3.0)
project(YourProjectName LANGUAGES Fortran)

# Set appropriate compiler for C++ and Fortran
if(UNIX AND APPLE)
  set(CMAKE_Fortran_COMPILER gfortran)
endif()

add_definitions(-fcray-pointer)

# Set link libraries
set(LINKLIBS -lstdc++)
file(GLOB TECIOLIB ../../teciosrc/build/libtecio.a)
include_directories(../../teciosrc)

# Set executable name and sources
file(GLOB F90FILES *.f90)

# Add executable targets
add_executable(test-f90 ${F90FILES})
target_link_libraries(test-f90 ${TECIOLIB} ${LINKLIBS})
