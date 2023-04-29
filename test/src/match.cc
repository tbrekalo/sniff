#include "sniff/match.h"

#include <array>

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

// sorted by target position
static constexpr auto kExpectedMatches = std::array<sniff::Match, 4>{
    sniff::Match{.query_pos = 19, .target_pos = 0},
    sniff::Match{.query_pos = 5, .target_pos = 3},
    sniff::Match{.query_pos = 7, .target_pos = 13},
    sniff::Match{.query_pos = 10, .target_pos = 17},
};

TEST_CASE("match-equality", "[match]") {
  REQUIRE(sniff::Match{} == sniff::Match{});
  REQUIRE(sniff::Match{0, 1} == sniff::Match{0, 1});
  REQUIRE_FALSE(sniff::Match{1, 0} == sniff::Match{0, 1});
  REQUIRE_FALSE(sniff::Match{0, 1} == sniff::Match{1, 1});
  REQUIRE_FALSE(sniff::Match{42, 314} == sniff::Match{101, 404});
}

TEST_CASE("match-inequality", "[match]") {
  REQUIRE_FALSE(sniff::Match{} != sniff::Match{});
  REQUIRE_FALSE(sniff::Match{0, 1} == sniff::Match{1, 0});
  REQUIRE_FALSE(sniff::Match{1, 0} == sniff::Match{0, 1});
  REQUIRE_FALSE(sniff::Match{0, 1} == sniff::Match{1, 1});
  REQUIRE_FALSE(sniff::Match{42, 314} == sniff::Match{101, 404});
}

TEST_CASE("make-matches", "[match][kmer]") {
  auto const matches = sniff::MakeMatches(kQueryKMers, kTargetKMers);
  REQUIRE(matches.size() == kExpectedMatches.size());
  for (std::size_t i = 0; i < matches.size(); ++i) {
    CHECK(matches[i] == kExpectedMatches[i]);
  }
}
