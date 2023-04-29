#include "sniff/chain.h"

#include <array>

#include "catch2/catch_test_macros.hpp"

static auto kArgMatches = std::vector<sniff::Match>{
    sniff::Match{.query_pos = 13, .target_pos = 1},
    sniff::Match{.query_pos = 20, .target_pos = 4},
    sniff::Match{.query_pos = 4, .target_pos = 7},
    sniff::Match{.query_pos = 9, .target_pos = 10},
    sniff::Match{.query_pos = 11, .target_pos = 13},
};

TEST_CASE("chain", "[chain][overlap]") {
  auto chain_cfg = sniff::ChainConfig{.min_target_chain_matches = 2,
                                      .max_target_allowed_gap = 100,
                                      .kmer_len = 5};

  auto overlaps = sniff::Chain(chain_cfg, kArgMatches);
  REQUIRE(overlaps.size() == 1);

  CHECK(overlaps[0].query_start == 4);
  CHECK(overlaps[0].query_end == 16);

  CHECK(overlaps[0].target_start == 7);
  CHECK(overlaps[0].target_end == 18);
}
