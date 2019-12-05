cmake_minimum_required (VERSION 2.8.12)
project (telescope)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wno-sign-compare -g")

## Set C++11 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 11)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++11)
endif()

include_directories(include ${CMAKE_SOURCE_DIR}/external)

add_executable(convert_alignments_intersection ${CMAKE_SOURCE_DIR}/src/convert_alignments_intersection.cpp)
add_executable(convert_alignments_union ${CMAKE_SOURCE_DIR}/src/convert_alignments_union.cpp)
add_executable(convert_alignments_unpaired ${CMAKE_SOURCE_DIR}/src/convert_alignments_unpaired.cpp)
add_executable(telescope ${CMAKE_SOURCE_DIR}/src/telescope.cpp
${CMAKE_SOURCE_DIR}/src/read_alignments.cpp ${CMAKE_SOURCE_DIR}/src/write_alignments.cpp)

target_link_libraries(convert_alignments_intersection z)
target_link_libraries(convert_alignments_union z)
target_link_libraries(convert_alignments_unpaired z)
target_link_libraries(telescope z)
