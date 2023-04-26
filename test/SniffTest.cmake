find_package(Catch2 3.3 QUIET)
if(NOT Catch2_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.3.2)

  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()

add_executable(sniff_test 
  ${CMAKE_CURRENT_LIST_DIR}/src/kmer.cc
  ${CMAKE_CURRENT_LIST_DIR}/src/minimize.cc)
target_link_libraries(sniff_test PRIVATE sniff_lib Catch2::Catch2WithMain)
