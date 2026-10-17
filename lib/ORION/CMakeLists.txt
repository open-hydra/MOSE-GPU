# CMake project file for ORION

cmake_minimum_required(VERSION 3.13)

project (ORION LANGUAGES Fortran CXX)

set(CMAKE_VERBOSE_MAKEFILE FALSE)
# seach path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING. You should review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries. There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Set the suffix for the compiler ID based on the compilers used.
INCLUDE(${CMAKE_MODULE_PATH}/SetCompilerID.cmake)

# TecIO
OPTION(USE_TECIO "Use TecIO" ON)
if (USE_TECIO)
  if (USE_MPI)
    set(basename "teciompi")
  else()
    set(basename "tecio")
  endif()
  set(TECIO_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/TecIO/${basename}-install-${COMPILER_SUFFIX}")
  if(NOT EXISTS "${TECIO_INSTALL_DIR}")
    message(STATUS "TecIO not installed, building and installing it...")

    # Build TecIO in a temporary build directory
    execute_process(COMMAND ${CMAKE_COMMAND} -S ${CMAKE_CURRENT_SOURCE_DIR}/lib/TecIO/${basename}src
                                           -B ${CMAKE_CURRENT_SOURCE_DIR}/lib/TecIO/${basename}-build-${COMPILER_SUFFIX}
                                           -DCMAKE_INSTALL_PREFIX=${TECIO_INSTALL_DIR}
                                           -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                           -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                    RESULT_VARIABLE CONFIGURE_RESULT)
    if(NOT CONFIGURE_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to configure TecIO")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_SOURCE_DIR}/lib/TecIO/${basename}-build-${COMPILER_SUFFIX}
                    RESULT_VARIABLE BUILD_RESULT)
    if(NOT BUILD_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to build TecIO")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --install ${CMAKE_CURRENT_SOURCE_DIR}/lib/TecIO/${basename}-build-${COMPILER_SUFFIX}
                    RESULT_VARIABLE INSTALL_RESULT)
    if(NOT INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to install TecIO")
    endif()
  endif()
  find_package(${basename} REQUIRED
             PATHS ${TECIO_INSTALL_DIR}
             NO_DEFAULT_PATH)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/TecIO/${basename}src)
  if (NOT USE_MPI)
    # With CXX linker lstdc++ is not needed, but with Fortran it is.
    # set(LINKLIBS -lpthread -lm -lstdc++)
    set(LINKLIBS -lpthread)
  endif()
endif()

include(GNUInstallDirs)

# Set directories
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Build library
add_subdirectory(src/fortran/lib)

# Build tests and applications (only if ORION build has not been called by another program)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message("-- Full build allowed")
    add_subdirectory(src/fortran/test)
    add_subdirectory(src/fortran/app)
endif()

