# Collect all compiler IDs
set(CXX_COMPILER_ID "${CMAKE_CXX_COMPILER_ID}")
set(FORTRAN_COMPILER_ID "${CMAKE_Fortran_COMPILER_ID}")

# Define a unified suffix
if(CXX_COMPILER_ID STREQUAL "GNU" AND FORTRAN_COMPILER_ID STREQUAL "GNU")
    set(COMPILER_SUFFIX "GNU")
elseif(CXX_COMPILER_ID STREQUAL "Intel" AND FORTRAN_COMPILER_ID STREQUAL "Intel")
    set(COMPILER_SUFFIX "Intel")
elseif(CXX_COMPILER_ID STREQUAL "Clang" AND FORTRAN_COMPILER_ID STREQUAL "Clang")
    # For example, using Clang + Flang or Clang + GFortran
    set(COMPILER_SUFFIX "ClangMixed")
else()
    # Fallback: concatenate the three IDs
    set(COMPILER_SUFFIX "${CXX_COMPILER_ID}-${FORTRAN_COMPILER_ID}")
endif()