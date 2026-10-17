# CMake project file for 'MOSE'

cmake_minimum_required(VERSION 3.23)

# Name the code
set(solver MOSE)
set(lib MOSEL)
set(batch MOSE-Batch)

project (${solver} LANGUAGES Fortran CXX C)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # seach path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Set some options the user may choose
  # Uncomment the below if you want the user to choose a parallelization library
  OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)
  OPTION(USE_OPENACC "OpenACC GPU offload (NVHPC/nvfortran only)" OFF)

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING. You should  review this file and make sure the flags 
  # are to your liking.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

  # OpenACC GPU offload (chemistry-first port, see GPU_PORTING_PLAN.md).
  # NVHPC only; the flags must also reach the link line for the device link.
  if(USE_OPENACC)
    if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
      message(FATAL_ERROR "USE_OPENACC requires nvfortran (NVHPC), got: ${CMAKE_Fortran_COMPILER_ID}")
    endif()
    set(OPENACC_GPU_OPTS "cc80" CACHE STRING "Comma-separated -gpu= options: compute capability (A30 = cc80), plus a cudaX.Y toolkit pin when building on a driverless host (e.g. cc80,cuda12.6)")
    set(OPENACC_FLAGS "-acc=gpu -gpu=${OPENACC_GPU_OPTS} -Minfo=accel")
    string(APPEND CMAKE_Fortran_FLAGS " ${OPENACC_FLAGS}")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OPENACC_FLAGS}")
  endif()
  # Locate and set parallelization libraries. There are some CMake peculiarities
  # taken care of here, such as the fact that the FindOpenMP routine doesn't know
  # about Fortran.
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # ORION
  if(MASTER STREQUAL None)
    set ( orionD ${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/ )
  elseif(MASTER STREQUAL hydra)
    set ( orionD $ENV{HYDRADIR}/lib/ORION/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/")
  endif()
  add_subdirectory(${orionD})

  # FLINT
  if(MASTER STREQUAL None)
    set ( flintD ${CMAKE_CURRENT_SOURCE_DIR}/lib/FLINT/ )
  elseif(MASTER STREQUAL hydra)
    set ( flintD $ENV{HYDRADIR}/lib/FLINT/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/FLINT/")
  endif()
  add_subdirectory(${flintD})

  # FiNeR
  if(MASTER STREQUAL None)
    set ( finerD ${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/ )
  elseif(MASTER STREQUAL hydra)
    set ( finerD $ENV{HYDRADIR}/lib/third_party/FiNeR/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/third_party/FiNeR/")
  endif()
  add_subdirectory(${finerD})

endif()

# Build library
add_subdirectory(src/lib)

OPTION(BASIC_TEST "Build basic tests" OFF)
if(BASIC_TEST)
  add_subdirectory(src/test)
endif()

# Build application (only if this file has not been called by another CMakeLists)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)
endif()
