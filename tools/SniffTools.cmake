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

add_executable(
  pairs_edit_dist
  ${CMAKE_CURRENT_LIST_DIR}/src/pairs_edit_dist.cc)
target_link_libraries(pairs_edit_dist PRIVATE cxxopts edlib fmt sniff_lib)
