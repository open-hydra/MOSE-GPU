cmake_minimum_required(VERSION 2.8.12) # requiring 2.8.12 sets the CMP0022 policy to new (non-transitive linking, i.e. just the interfaces) by default

# -------- Project metadata and language setup --------
project(teciompi LANGUAGES C CXX Fortran)

# Default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type (Debug|Release)" FORCE)
endif ()

# -------- Boost detection and fallback --------
find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
else()
    add_custom_target(ExtractBoost
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/.. tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/../boost.tar.gz
        COMMENT "Extracting Boost library"
        VERBATIM)

    add_custom_target(CreateBoostLink ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../boost ${CMAKE_CURRENT_SOURCE_DIR}/boost
        COMMENT "Creating symbolic link to Boost directory"
        VERBATIM)

    add_dependencies(CreateBoostLink ExtractBoost)

    set(Boost_INCLUDE_DIR "boost")
    message(WARNING "Cannot find Boost. Hardcoding include directory to ${Boost_INCLUDE_DIR}")
    include_directories(SYSTEM ${Boost_INCLUDE_DIR})
endif()

FIND_PACKAGE(MPI)
IF (NOT MPI_CXX_FOUND)
    message(FATAL_ERROR "Unable to find MPI installation.")
ENDIF ()

include_directories(.)
if(NOT "${MPI_CXX_INCLUDE_PATH}" STREQUAL "")
    include_directories(SYSTEM "${MPI_CXX_INCLUDE_PATH}")
endif()

# -------- Compilation flags --------
IF (WIN32)
    set(BaseFlags "/EHsc /MP /wd\"4996\" /D MSWIN /D TP_PROJECT_USES_BOOST /D BOOST_ALL_NO_LIB /D MAKEARCHIVE /D NO_THIRD_PARTY_LIBS /D TECIOMPI /D NO_ASSERTS")
ELSE ()
    set(BaseFlags "-DLINUX -DLINUX64 -DTP_PROJECT_USES_BOOST -DBOOST_ALL_NO_LIB -DMAKEARCHIVE -DNO_THIRD_PARTY_LIBS -DTECIOMPI -DOMPI_SKIP_MPICXX -DNO_ASSERTS -fPIC -fvisibility=hidden -w")
ENDIF ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BaseFlags} ${MPI_CXX_COMPILE_FLAGS}")

# Portable AVX2 baseline (NVHPC nvc++ only). nvc++ defaults -tp to the build host arch;
# on a Cascade Lake build host that is AVX-512, which SIGILLs on AMD EPYC (fatgpu) nodes.
# Pin x86-64-v3 so the static TecIO archive runs on both Intel and AMD GPU nodes.
if (CMAKE_CXX_COMPILER_ID MATCHES "NVHPC" OR CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -tp=x86-64-v3")
endif()

# -------- Library source files --------
FILE(GLOB mainFiles "*.cpp" "*.h")
LIST(REMOVE_ITEM mainFiles "${CMAKE_CURRENT_SOURCE_DIR}/szcombine.cpp")

add_library(teciompi STATIC ${mainFiles})
add_library(teciompi::teciompi ALIAS teciompi)

if (NOT Boost_FOUND)
    add_dependencies(teciompi CreateBoostLink)
endif()

# -------- Executable --------

add_executable(szcombine "szcombine.cpp")
target_link_libraries(szcombine teciompi ${MPI_LIBRARIES})

# -------- Examples --------
set(example_dirs
    ijkpartitioned-MPI
)

foreach(example_dir ${example_dirs})
    set ( testD ../examples/${example_dir}/ "${CMAKE_CURRENT_SOURCE_DIR}/../examples/${example_dir}/build")
    add_subdirectory(${testD})
endforeach()


# -------- Install rules --------
include(GNUInstallDirs)

# Headers (if any)
install(DIRECTORY ./
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/teciompi
    FILES_MATCHING PATTERN "*.h"
)

# Library install/export
install(TARGETS teciompi
    EXPORT teciompiTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/teciompi
)

install(EXPORT teciompiTargets
    FILE teciompiTargets.cmake
    NAMESPACE teciompi::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/teciompi
)

# -------- CMake package config --------
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/teciompiConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/teciompiConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/teciompi
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/teciompiConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/teciompi
)

