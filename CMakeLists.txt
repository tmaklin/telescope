cmake_minimum_required (VERSION 3.11)
project (telescope)
include(FetchContent)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Determine build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -g -Wall -Wextra -Wpedantic")
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xHost -fp-model fast=2 -no-prec-div -fimf-use-svml=true -ipo -m64")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -xHost -fp-model fast=2 -no-prec-div -fimf-use-svml=true -ipo -m64")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -O3 -xHost -fp-model fast=2 -no-prec-div -fimf-use-svml=true -ipo -m64")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -funroll-loops")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ffast-math -funroll-loops")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ffast-math -funroll-loops")
  endif()
endif()

## FLTO
if(CMAKE_BUILD_WITH_FLTO)
  cmake_policy(SET CMP0069 NEW)
  set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION 1)
endif()

## Native instructions
if (CMAKE_WITH_NATIVE_INSTRUCTIONS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
else()
  set(CMAKE_WITH_NATIVE_INSTRUCTIONS 0)
endif()

## PGO
if(CMAKE_PGO_GENERATE AND NOT CMAKE_PGO_USE)
  if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate")
  endif()
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
  endif()
endif()

## Set C++11 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 17)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++17)
endif()

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE TELESCOPE_BUILD_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "telescope version: ${TELESCOPE_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(TELESCOPE_BUILD_VERSION "unknown")
endif()
string(TIMESTAMP TELESCOPE_BUILD_TIMESTAMP)
### Generate a telescope_version.h file containing build version and timestamp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/telescope_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/telescope_version.h @ONLY)

## telescope library
add_library(libtelescope
${CMAKE_CURRENT_SOURCE_DIR}/src/write_alignments.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/read_themisto_alignments.cpp)

set_target_properties(libtelescope PROPERTIES OUTPUT_NAME telescope)

## telescope executable
add_executable(telescope ${CMAKE_CURRENT_SOURCE_DIR}/src/telescope.cpp)

## Dependencies
### Check supported compression types
#### zlib
if ((DEFINED ZLIB_LIBRARY AND DEFINED ZLIB_INCLUDE_DIR) AND (NOT DEFINED ZLIB_FOUND))
  message(STATUS "zlib library provided in: " ${ZLIB_LIBRARY})
  message(STATUS "zlib headers provided in: " ${ZLIB_INCLUDE_DIR})
  include_directories(${ZLIB_INCLUDE_DIR})
  target_link_libraries(telescope libtelescope ${ZLIB_LIBRARY})
else()
  find_package(ZLIB)
  if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIR})
    target_link_libraries(telescope libtelescope ${ZLIB_LIBRARY})
  endif()
endif()

#### bzip2
if (DEFINED BZIP2_LIBRARIES AND DEFINED BZIP2_INCLUDE_DIR AND (NOT DEFINED BZIP2_FOUND))
  message(STATUS "bzip2 library provided in: " ${BZIP2_LIBRARIES})
  message(STATUS "bzip2 headers provided in: " ${BZIP2_INCLUDE_DIR})
  include_directories(${BZIP2_INCLUDE_DIR})
  target_link_libraries(telescope libtelescope ${BZIP2_LIBRARIES})
else()
  find_package(BZip2)
  if (BZIP2_FOUND)
    include_directories(${BZIP2_INCLUDE_DIR})
    target_link_libraries(telescope libtelescope ${BZIP2_LIBRARIES})
  endif()
endif()

#### lzma
if (DEFINED LIBLZMA_LIBRARY AND DEFINED LIBLZMA_INCLUDE_DIR AND (NOT DEFINED LIBLZMA_FOUND))
  message(STATUS "liblzma library provided in: " ${LIBLZMA_LIBRARY})
  message(STATUS "liblzma headers provided in: " ${LIBLZMA_INCLUDE_DIR})
  include_directories(${LIBLZMA_INCLUDE_DIR})
  target_link_libraries(telescope libtelescope ${LIBLZMA_LIBRARY})
else()
  find_package(LibLZMA)
  if (LIBLZMA_FOUND)
    include_directories(${LIBLZMA_INCLUDE_DIR})
    target_link_libraries(telescope libtelescope ${LIBLZMA_LIBRARY})
  endif()
endif()

### Check OpenMP support
find_package(OpenMP)
if (OPENMP_FOUND)
  set(TELESCOPE_OPENMP_SUPPORT 1)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
  target_link_libraries(libtelescope OpenMP::OpenMP_CXX)
else()
  set(TELESCOPE_OPENMP_SUPPORT 0)
endif()
#### Configure OpenMP if it is supported on the system.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/telescope_openmp_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/include/telescope_openmp_config.hpp @ONLY)

### Check required dependenciesand and download them if not found
#### BitMagic - supplied with the project
#### BitMagic
if (DEFINED CMAKE_BITMAGIC_HEADERS)
  message(STATUS "BitMagic headers provided in: ${CMAKE_BITMAGIC_HEADERS}")
else()
  set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/BitMagic-7.12.3/src)
endif()
include_directories(${CMAKE_BITMAGIC_HEADERS})

#### bxzstr
if (DEFINED CMAKE_BXZSTR_HEADERS)
  message(STATUS "bxzstr headers provided in: ${CMAKE_BXZSTR_HEADERS}")
else()
  FetchContent_Declare(bxzstr
    GIT_REPOSITORY    https://github.com/tmaklin/bxzstr.git
    GIT_TAG           v1.1.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS        -D ZSTD_FOUND=0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(bxzstr)
  set(CMAKE_BXZSTR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr/include)
endif()
include_directories(${CMAKE_BXZSTR_HEADERS})

#### cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  message(STATUS "cxxargs headers provided in: ${CMAKE_CXXARGS_HEADERS}")
else()
  FetchContent_Declare(cxxargs
    GIT_REPOSITORY    https://github.com/tmaklin/cxxargs.git
    GIT_TAG           v1.1.4
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs"
    BUILD_IN_SOURCE   0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(cxxargs)
  set(CMAKE_CXXARGS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs/include)
endif()
include_directories("${CMAKE_CXXARGS_HEADERS}")

#### cxxio
if (DEFINED CMAKE_CXXIO_HEADERS)
  message(STATUS "cxxio headers provided in: ${CMAKE_CXXIO_HEADERS}")
else()
  FetchContent_Declare(cxxio
    GIT_REPOSITORY    https://github.com/tmaklin/cxxio.git
    GIT_TAG           v0.1.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/cxxio"
    BUILD_IN_SOURCE   0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(cxxio)
  set(CMAKE_CXXIO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxio/include)
endif()
include_directories("${CMAKE_CXXIO_HEADERS}")

#### alignment-writer
if (DEFINED CMAKE_ALIGNMENT_WRITER_HEADERS AND DEFINED CMAKE_ALIGNMENT_WRITER_LIBRARY)
  message(STATUS "alignment-writer headers provided in: ${CMAKE_ALIGNMENT_WRITER_HEADERS}")
  message(STATUS "alignment-writer library provided in: ${CMAKE_ALIGNMENT_WRITER_LIBRARY}")
else()
  FetchContent_Declare(alignment-writer
    GIT_REPOSITORY    https://github.com/tmaklin/alignment-writer.git
    GIT_TAG           v0.5.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_BXZSTR_HEADERS=${CMAKE_BXZSTR_HEADERS}
		    -D CMAKE_CXXARGS_HEADERS=${CMAKE_CXXARGS_HEADERS}
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(alignment-writer)
  add_dependencies(libtelescope libalignmentwriter)
  set_target_properties(alignment-writer PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer/include)
  set(CMAKE_ALIGNMENT_WRITER_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/lib/libalignment-writer.a)
endif()
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})
target_link_libraries(libtelescope ${CMAKE_ALIGNMENT_WRITER_LIBRARY})

#### googletest if building tests
if (DEFINED CMAKE_ALIGNMENT_WRITER_HEADERS AND DEFINED CMAKE_ALIGNMENT_WRITER_LIBRARY)
  message(STATUS "alignment-writer headers provided in: ${CMAKE_ALIGNMENT_WRITER_HEADERS}")
  message(STATUS "alignment-writer library provided in: ${CMAKE_ALIGNMENT_WRITER_LIBRARY}")
else()
  FetchContent_Declare(alignment-writer
    GIT_REPOSITORY    https://github.com/tmaklin/alignment-writer.git
    GIT_TAG           v0.5.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_BXZSTR_HEADERS=${CMAKE_BXZSTR_HEADERS}
		    -D CMAKE_CXXARGS_HEADERS=${CMAKE_CXXARGS_HEADERS}
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(alignment-writer)
  add_dependencies(libtelescope libalignmentwriter)
  set_target_properties(alignment-writer PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer/include)
  set(CMAKE_ALIGNMENT_WRITER_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/lib/libalignment-writer.a)
endif()
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})
target_link_libraries(libtelescope ${CMAKE_ALIGNMENT_WRITER_LIBRARY})

# #### Download googletest if building tests
# if(CMAKE_BUILD_TESTS)
#   if (DEFINED CMAKE_GOOGLETEST_HEADERS)
#     message(STATUS "googletest headers provided in: ${CMAKE_GOOGLETEST_HEADERS}")
#   else()
#     FetchContent_Declare(googletest
#       GIT_REPOSITORY    https://github.com/google/googletest.git
#       GIT_TAG           release-1.11.0
#       PREFIX            "external"
#       SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest"
#       BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/googletest"
#       BUILD_IN_SOURCE   0
#       CMAKE_ARGS      -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
# 		      -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
# 		      -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
# 		      -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
# 		      -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
#     INSTALL_COMMAND   ""
#     )
#   FetchContent_MakeAvailable(googletest)
#   add_dependencies(libtelescope googletest)
#   enable_testing()
#   endif()
# endif()

## Project headers
set(CMAKE_TELESCOPE_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_TELESCOPE_HEADERS})
