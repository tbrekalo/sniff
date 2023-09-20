find_package(Catch2 3.3 REQUIRED)

add_executable(
  sniff_test
  ${CMAKE_CURRENT_LIST_DIR}/src/kmer.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/map.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/match.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/minimize.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/overlap.cc)
target_link_libraries(sniff_test PRIVATE sniff_lib Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(sniff_test)
