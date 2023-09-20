FetchContent_Declare(
  wfa2
  GIT_REPOSITORY https://github.com/smarco/WFA2-lib
  GIT_TAG v2.3.3)

FetchContent_GetProperties(wfa2)
if(NOT wfa2_POPULATED)
  FetchContent_Populate(wfa2)
  add_subdirectory(${wfa2_SOURCE_DIR} ${wfa2_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

add_executable(pairs_edit_dist ${CMAKE_CURRENT_LIST_DIR}/src/pairs_edit_dist.cc)
target_link_libraries(
  pairs_edit_dist PRIVATE cxxopts::cxxopts fmt::fmt sniff_lib
                          unordered_dense::unordered_dense wfa2cpp_static)
