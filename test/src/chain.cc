#include "sniff/chain.h"

#include "catch2/catch_test_macros.hpp"

static auto const kQueryKMers = std::vector<sniff::KMer>{
    sniff::KMer{.position = 0, .value = 0},
    sniff::KMer{.position = 5, .value = 1},
    sniff::KMer{.position = 7, .value = 2},
    sniff::KMer{.position = 10, .value = 2},
    sniff::KMer{.position = 15, .value = 7},
    sniff::KMer{.position = 19, .value = 5},
};

static auto const kTargetKMers =
    std::vector<sniff::KMer>{sniff::KMer{.position = 0, .value = 5},
                             sniff::KMer{.position = 3, .value = 1},
                             sniff::KMer{.position = 13, .value = 2},
                             sniff::KMer{.position = 17, .value = 2},
                             sniff::KMer{.position = 21, .value = 5}};

TEST_CASE("chain", "[chain][overlap]") {
  auto chain_cfg = sniff::ChainConfig{.min_target_chain_matches = 2,
                                      .max_target_allowed_gap = 100,
                                      .kmer_len = 5};

  auto overlaps = sniff::Chain(chain_cfg, kQueryKMers, kTargetKMers);
  REQUIRE(overlaps.size() == 1);

  CHECK(overlaps[0].query_start == 5);
  CHECK(overlaps[0].query_end == 15);

  CHECK(overlaps[0].target_start == 3);
  CHECK(overlaps[0].target_end == 22);
}
