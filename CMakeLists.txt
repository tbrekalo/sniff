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

option(SNIFF_DEV_MODE "Build sniff with development" OFF)
option(SNIFF_BUILD_ASAN "Build Debug and RelWithDebInfo with ASAN"
       ${SNIFF_DEV_MODE})
option(SNIFF_BUILD_TESTS "Build sniff unit tests" ${SNIFF_DEV_MODE})

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
  GIT_TAG 10.0.0)

FetchContent_Declare(
  quill
  GIT_REPOSITORY https://github.com/odygrd/quill
  GIT_TAG v3.2.0)

FetchContent_Declare(
  unordered_dense
  GIT_REPOSITORY https://github.com/martinus/unordered_dense
  GIT_TAG v4.0.0)

FetchContent_MakeAvailable(bioparser biosoup cxxopts fmt quill unordered_dense)

if(SNIFF_DEV_MODE)
  add_definitions(-DSNIFF_DEV_MODE=1)
  add_definitions(-DQUILL_ACTIVE_LOG_LEVEL=QUILL_LOG_LEVEL_TRACE_L3)
else()
  add_definitions(-DSNIFF_DEV_MODE=0)
  add_definitions(-DQUILL_ACTIVE_LOG_LEVEL=QUILL_LOG_LEVEL_INFO)
endif()

find_package(TBB)
add_library(
  sniff_lib
  src/algo.cc
  src/config.cc
  src/io.cc
  src/kmer.cc
  src/map.cc
  src/match.cc
  src/minimize.cc
  src/overlap.cc
  src/sketch.cc)
target_include_directories(
  sniff_lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                   $<INSTALL_INTERFACE:include>)
target_link_libraries(
  sniff_lib
  PUBLIC biosoup TBB::tbb
  PRIVATE bioparser fmt quill::quill unordered_dense)

add_executable(sniff src/main.cc)
target_include_directories(sniff
                           PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
target_link_libraries(sniff PRIVATE cxxopts fmt quill::quill sniff_lib)

if(SNIFF_BUILD_ASAN)
  target_compile_options(
    sniff PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>
                  $<$<CONFIG:Debug,RelWithDebInfo>:-fno-omit-frame-pointer>)
  target_link_options(sniff PRIVATE
                      $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address>)
endif()

set_target_properties(sniff PROPERTIES VERSION ${sniff_VERSION}
                                       SOVERSION ${sniff_VERSION_MAJOR})
configure_file(${PROJECT_SOURCE_DIR}/src/version.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/include/sniff/version.h)
target_include_directories(
  sniff PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/sniff>)

if(SNIFF_BUILD_TESTS)
  include(${CMAKE_CURRENT_LIST_DIR}/test/SniffTest.cmake)
endif()

if(SNIFF_DEV_MODE)
  include(${CMAKE_CURRENT_LIST_DIR}/tools/SniffTools.cmake)
endif()
