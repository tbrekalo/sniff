#include "sniff/overlap.h"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("overlap-length", "[overlap]") {
  SECTION("query-longer") {
    CHECK(sniff::OverlapLength(sniff::Overlap{.query_start = 0,
                                              .query_end = 10,
                                              .target_start = 0,
                                              .target_end = 5}) == 10);
  }

  SECTION("target-longer") {
    CHECK(sniff::OverlapLength(sniff::Overlap{.query_start = 0,
                                              .query_end = 5,
                                              .target_start = 0,
                                              .target_end = 10}) == 10);
  }

  SECTION("query-target-equal") {
    CHECK(sniff::OverlapLength(sniff::Overlap{.query_start = 0,
                                              .query_end = 5,
                                              .target_start = 0,
                                              .target_end = 5}) == 5);
  }
}

TEST_CASE("overlap-error", "[overlap]") {
  SECTION("query-longer") {
    CHECK(sniff::OverlapError(sniff::Overlap{.query_start = 0,
                                             .query_end = 10,
                                             .target_start = 0,
                                             .target_end = 5}) == 0.5);
  }

  SECTION("target-longer") {
    CHECK(sniff::OverlapError(sniff::Overlap{.query_start = 0,
                                             .query_end = 5,
                                             .target_start = 0,
                                             .target_end = 10}) == 0.5);
  }

  SECTION("query-target-equal") {
    CHECK(sniff::OverlapError(sniff::Overlap{.query_start = 0,
                                             .query_end = 5,
                                             .target_start = 0,
                                             .target_end = 5}) == 0.0);
  }
}
