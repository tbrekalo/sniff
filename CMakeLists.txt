cmake_minimum_required(VERSION 3.21)
project(
  sniff
  VERSION 1.0.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# output dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)

FetchContent_Declare(
  bioparser
  GIT_REPOSITORY https://github.com/rvaser/bioparser
  GIT_TAG 3.0.15)

FetchContent_Declare(
  biosoup
  GIT_REPOSITORY https://github.com/rvaser/biosoup
  GIT_TAG 0.10.0)

FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts
  GIT_TAG v3.1.1)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt
  GIT_TAG 9.1.0)

FetchContent_Declare(
  unordered_dense
  GIT_REPOSITORY https://github.com/martinus/unordered_dense
  GIT_TAG v4.0.0
)

FetchContent_MakeAvailable(bioparser biosoup cxxopts fmt unordered_dense)

find_package(TBB)
add_library(
  sniff_lib
  src/algo.cc
  src/config.cc
  src/io.cc
  src/kmer.cc
  src/map.cc
  src/map_config.cc
  src/match.cc
  src/minimize.cc
  src/minimize_config.cc
  src/overlap.cc
  src/sketch.cc)
target_include_directories(
  sniff_lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                   $<INSTALL_INTERFACE:include>)
target_link_libraries(
  sniff_lib
  PUBLIC biosoup TBB::tbb
  PRIVATE bioparser fmt unordered_dense)

add_executable(sniff src/main.cc)
target_include_directories(sniff
                           PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
target_link_libraries(sniff PRIVATE cxxopts fmt sniff_lib)
target_compile_options(
  sniff PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
                $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
target_link_options(sniff PRIVATE
                    $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)

set_target_properties(sniff PROPERTIES VERSION ${sniff_VERSION}
                                       SOVERSION ${sniff_VERSION_MAJOR})
configure_file(${PROJECT_SOURCE_DIR}/src/version.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/include/sniff/version.h)
target_include_directories(
  sniff PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/sniff>)

if(${PROJECT_IS_TOP_LEVEL})
  include(${CMAKE_CURRENT_LIST_DIR}/test/SniffTest.cmake)
endif()
