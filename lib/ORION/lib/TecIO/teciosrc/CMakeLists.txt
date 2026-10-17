cmake_minimum_required(VERSION 2.8.12)

# -------- Project metadata and language setup --------
project(tecio LANGUAGES CXX Fortran)

# Default build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type (Debug|Release)" FORCE)
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

include_directories(.)

# -------- Compilation flags --------
if (WIN32)
    set(BaseFlags "/EHsc /MP /wd\"4996\" /D MSWIN /D TP_PROJECT_USES_BOOST /D BOOST_ALL_NO_LIB /D MAKEARCHIVE /D NO_THIRD_PARTY_LIBS /D NO_ASSERTS")
else()
    set(BaseFlags "-DLINUX -DLINUX64 -DTP_PROJECT_USES_BOOST -DBOOST_ALL_NO_LIB -DMAKEARCHIVE -DNO_THIRD_PARTY_LIBS -DNO_ASSERTS -fmessage-length=0 -fPIC -fvisibility=hidden -w")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BaseFlags}")

# -------- Library source files --------
file(GLOB mainFiles "*.cpp" "*.h")
list(FILTER mainFiles EXCLUDE REGEX "szcombine\\.cpp$")

add_library(tecio STATIC ${mainFiles})
add_library(tecio::tecio ALIAS tecio)

set_target_properties(tecio PROPERTIES
    LINKER_LANGUAGE Fortran
    PUBLIC_HEADER ""  # Add headers here if desired
)

if (NOT Boost_FOUND)
    add_dependencies(tecio CreateBoostLink)
endif()

# -------- Executable --------
# add_executable(szcombine szcombine.cpp)
# if (WIN32)
#     target_link_libraries(szcombine tecio)
# else()
#     target_link_libraries(szcombine tecio pthread)
# endif()

# -------- Examples --------
set(example_dirs
    rewriteszl
    ijkpartitioned
    ij_ordered
)

foreach(example_dir ${example_dirs})
    set ( testD ../examples/${example_dir}/ "${CMAKE_CURRENT_SOURCE_DIR}/../examples/${example_dir}/build")
    add_subdirectory(${testD})
endforeach()

# -------- Install rules --------
include(GNUInstallDirs)

# Headers (if any)
install(DIRECTORY ./
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tecio
    FILES_MATCHING PATTERN "*.h"
)

# Library install/export
install(TARGETS tecio
    EXPORT tecioTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tecio
)

install(EXPORT tecioTargets
    FILE tecioTargets.cmake
    NAMESPACE tecio::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tecio
)

# -------- CMake package config --------
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tecioConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/tecioConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tecio
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/tecioConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tecio
)
