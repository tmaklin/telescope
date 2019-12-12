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

include_directories(include ${CMAKE_SOURCE_DIR}/external)

add_library(libtelescope ${CMAKE_SOURCE_DIR}/src/read_alignments.cpp ${CMAKE_SOURCE_DIR}/src/write_alignments.cpp)
set_target_properties(libtelescope PROPERTIES OUTPUT_NAME telescope)

add_executable(telescope ${CMAKE_SOURCE_DIR}/src/telescope.cpp)
target_link_libraries(telescope libtelescope z)
