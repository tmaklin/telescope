cmake_minimum_required (VERSION 3.0)
project (telescope)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(ZLIB REQUIRED)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wno-sign-compare -g")

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

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
## Generate a version.h file containing build version and timestamp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/version.h @ONLY)

include_directories(include ${CMAKE_CURRENT_SOURCE_DIR}/external ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(libtelescope
${CMAKE_CURRENT_SOURCE_DIR}/src/write_alignments.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/read_themisto_alignments.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/read_kallisto_alignments.cpp)

set_target_properties(libtelescope PROPERTIES OUTPUT_NAME telescope)

add_executable(telescope ${CMAKE_CURRENT_SOURCE_DIR}/src/telescope.cpp)
target_link_libraries(telescope libtelescope z)
