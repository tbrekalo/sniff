#include "sniff/kmer.h"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("kmer equality", "[kmer][euqlaity]") {
  REQUIRE(sniff::KMer{} == sniff::KMer{});
  REQUIRE(sniff::KMer{.position = 1, .value = 5} ==
          sniff::KMer{.position = 1, .value = 5});
  REQUIRE_FALSE(sniff::KMer{.position = 1, .value = 5} ==
                sniff::KMer{.position = 1, .value = 6});
  REQUIRE_FALSE(sniff::KMer{.position = 2, .value = 5} ==
                sniff::KMer{.position = 1, .value = 5});
  REQUIRE_FALSE(sniff::KMer{.position = 2, .value = 5} ==
                sniff::KMer{.position = 1, .value = 6});
}

TEST_CASE("kmer inequality", "[kmer][inequality]") {
  REQUIRE_FALSE(sniff::KMer{} != sniff::KMer{});
  REQUIRE_FALSE(sniff::KMer{.position = 1, .value = 5} !=
                sniff::KMer{.position = 1, .value = 5});
  REQUIRE(sniff::KMer{.position = 1, .value = 5} !=
          sniff::KMer{.position = 1, .value = 6});
  REQUIRE(sniff::KMer{.position = 2, .value = 5} !=
          sniff::KMer{.position = 1, .value = 5});
  REQUIRE(sniff::KMer{.position = 2, .value = 5} !=
          sniff::KMer{.position = 1, .value = 6});
}
