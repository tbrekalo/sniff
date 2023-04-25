cmake_minimum_required(VERSION 3.17)
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

find_package(edlib 1.2.7 QUIET)
if(NOT edlib_FOUND)
  FetchContent_Declare(
    edlib
    GIT_REPOSITORY https://github.com/martinsos/edlib
    GIT_TAG v1.2.7)

  FetchContent_GetProperties(edlib)
  if(NOT edlib_POPULATED)
    FetchContent_Populate(edlib)
    add_subdirectory(${edlib_SOURCE_DIR} ${edlib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

FetchContent_MakeAvailable(bioparser biosoup cxxopts)

find_package(TBB)
add_executable(sniff src/main.cc src/io.cc)
target_include_directories(sniff
                           PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
target_link_libraries(sniff PRIVATE bioparser biosoup cxxopts edlib)
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