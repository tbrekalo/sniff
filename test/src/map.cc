#include "sniff/map.h"

#include <array>
#include <random>

#include "catch2/catch_test_macros.hpp"

static constexpr auto kMapCfg = sniff::MapConfig{
    .min_chain_length = 2, .max_chain_gap_length = 100, .kmer_len = 5};

TEST_CASE("map-one-overlap", "[map][overlap]") {
  auto rng_engine = std::mt19937{42};
  auto matches = std::vector<sniff::Match>{
      sniff::Match{.query_pos = 13, .target_pos = 1},
      sniff::Match{.query_pos = 20, .target_pos = 4},
      sniff::Match{.query_pos = 4, .target_pos = 7},
      sniff::Match{.query_pos = 9, .target_pos = 10},
      sniff::Match{.query_pos = 11, .target_pos = 13},
  };

  auto const assertions = [](std::vector<sniff::Overlap> overlaps) -> void {
    REQUIRE(overlaps.size() == 1);
    CHECK(overlaps[0].query_start == 4);
    CHECK(overlaps[0].query_end == 16);

    CHECK(overlaps[0].target_start == 7);
    CHECK(overlaps[0].target_end == 18);
  };

  SECTION("2nd-map-is-dominant") {
    std::shuffle(matches.begin(), matches.end(), rng_engine);
    assertions(sniff::Map(kMapCfg, matches));
  }

  SECTION("two-equal-dominant-maps-one-after-another") {
    matches.push_back(sniff::Match{.query_pos = 21, .target_pos = 6});
    std::shuffle(matches.begin(), matches.end(), rng_engine);
    assertions(sniff::Map(kMapCfg, matches));
  }
}

TEST_CASE("map-two-overlaps", "[map][overlap]") {
  constexpr auto kExpectedOverlaps =
      std::array<sniff::Overlap, 2>{sniff::Overlap{.query_start = 0,
                                                   .query_end = 14,
                                                   .target_start = 1,
                                                   .target_end = 12},
                                    sniff::Overlap{.query_start = 113,
                                                   .query_end = 127,
                                                   .target_start = 108,
                                                   .target_end = 127}};

  auto rng_engine = std::mt19937{42};
  auto matches = std::vector<sniff::Match>{
      sniff::Match{.query_pos = 0, .target_pos = 1},
      sniff::Match{.query_pos = 4, .target_pos = 5},
      sniff::Match{.query_pos = 9, .target_pos = 7},

      sniff::Match{.query_pos = 113, .target_pos = 108},
      sniff::Match{.query_pos = 115, .target_pos = 118},
      sniff::Match{.query_pos = 122, .target_pos = 122},
  };

  std::shuffle(matches.begin(), matches.end(), rng_engine);
  auto const overlaps = sniff::Map(kMapCfg, matches);

  REQUIRE(overlaps.size() == kExpectedOverlaps.size());
}
