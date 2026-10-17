cmake_minimum_required(VERSION 3.23)

# -------- Compiler overrides (for macOS) --------
# if (UNIX AND APPLE)
#     set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-14")
#     set(CMAKE_C_COMPILER   "/opt/homebrew/bin/gcc-14")
# endif ()

# -------- Say my name! --------
set(lib "FLINT")

project (FLINT LANGUAGES C CXX Fortran)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  set(CMAKE_VERBOSE_MAKEFILE FALSE)

  # Search path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 

endif()

# -------- Dependencies linking --------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(LinkDependencies)

# -------- Git submodules --------

# ORION
# If FLINT is built as part of another project (assumed that project already has ORION), skip the FLINT's ORION 
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  if ("${MASTER}" STREQUAL "hydra")
    set ( orionD $ENV{HYDRADIR}/lib/ORION/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/")
  else()
    set ( orionD ${CMAKE_SOURCE_DIR}/lib/ORION/ )
  endif()
  add_subdirectory(${orionD})
endif()

# OSlo
# if FLINT is built as part of hydra, use the hydra's OSlo, otherwise use the FLINT's one
if ("${MASTER}" STREQUAL "hydra")
  set ( oslosrcD $ENV{HYDRADIR}/lib/OSlo )
  set ( osloD $ENV{HYDRADIR}/lib/OSlo "${CMAKE_SOURCE_DIR}/build/lib/OSlo")
else()
  set ( osloD ${CMAKE_CURRENT_SOURCE_DIR}/lib/OSlo )
  set ( oslosrcD ${CMAKE_CURRENT_SOURCE_DIR}/lib/OSlo )
  # Check if cloned
  if(NOT EXISTS "${oslosrcD}/README.md")
    message(STATUS "OSlo not found. Attempting to initialize...")
    execute_process(COMMAND git submodule update --init lib/OSlo
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to initialize git submodules.")
    endif()
  endif()
endif()
add_subdirectory(${osloD})

# Cantera
# Cantera is optional and only used by FLINT
OPTION(USE_CANTERA "Use Cantera" ON)
if (USE_CANTERA)
  set(canteraD "${CMAKE_CURRENT_SOURCE_DIR}/lib/cantera")
  # Check if cloned
  if(NOT EXISTS "${canteraD}/README.rst")
    message(STATUS "Cantera not found. Attempting to initialize...")
    execute_process(COMMAND git submodule update --init lib/cantera
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to initialize git submodules.")
    endif()
  endif()
  # Check if already installed
  set(CANTERA_CONF "${canteraD}/config.log")
  if(NOT EXISTS "${CANTERA_CONF}")
    message(STATUS "Cantera not installed, building and installing it...")
    find_package(Boost)
    if (NOT Boost_FOUND)
      message(FATAL_ERROR "Boost is required to build Cantera. Please install Boost and try again.")
    endif()
    if(APPLE)
      set (sundials_libdir ${oslosrcD}/sundials-install/lib)
    elseif(UNIX)
      set (sundials_libdir ${oslosrcD}/sundials-install/lib64)
    endif()
    include(ExternalProject)
    ExternalProject_Add(
        cantera_world
        PREFIX cantera_build
        SOURCE_DIR ${canteraD}
        BINARY_DIR ${CMAKE_BINARY_DIR}/cantera_build
        CONFIGURE_COMMAND ""  # Cantera does not use configure step
        BUILD_COMMAND scons -C ${canteraD} build
            prefix=${CMAKE_BINARY_DIR}/cantera_install
            f90_interface=y
            FORTRAN=gfortran
            system_sundials=y
            sundials_include=${oslosrcD}/sundials-install/include
            sundials_libdir=${sundials_libdir}
            boost_inc_dir=${Boost_INCLUDE_DIR}
            system_eigen=n
        INSTALL_COMMAND scons -C ${canteraD} install
    )
  else()
    message(STATUS "Cantera already installed.")
  endif()
  link_directories("${canteraD}/build/lib")
  include_directories("${canteraD}/build/src/fortran/")
  add_compile_definitions(CANTERA)
endif(USE_CANTERA)

# -------- Fortran directories --------
include(GNUInstallDirs)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# -------- Define FLINT sources --------
file(GLOB LIBSRC src/lib/*.f90 src/lib/Lib_ChemMech/*.f90)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  set(FLINT_SRC ${LIBSRC})
else()
  set(FLINT_SRC ${LIBSRC} PARENT_SCOPE)
endif()

# -------- Build FLINT test suite --------
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  add_subdirectory("src/test/Fortran")
  if (USE_CANTERA)
    add_subdirectory("src/test/CXX")
  endif()
endif()