#include "sniff/minimize.h"

#include <array>
#include <utility>

#include "catch2/catch_test_macros.hpp"

static constexpr auto kTestSequence =
    std::string_view{"GCGTGCCATAACCACCATATTCGACGATTCAAC"};

static constexpr auto kTestExpectedMinimizersK15W5 = std::array<sniff::KMer, 8>{
    sniff::KMer{1, 462733637},  sniff::KMer{5, 348210483},
    sniff::KMer{6, 319100111},  sniff::KMer{7, 202658621},
    sniff::KMer{8, 810634486},  sniff::KMer{12, 290256408},
    sniff::KMer{14, 349135247}, sniff::KMer{15, 322799165}};

static constexpr auto kTestExpectedMinimizersK7W7 =
    std::array<sniff::KMer, 4>{sniff::KMer{5, 5313}, sniff::KMer{10, 1300},
                               sniff::KMer{14, 5327}, sniff::KMer{19, 15750}};

static constexpr auto kTestExpectedMinimizersK5W7 = std::array<sniff::KMer, 6>{
    sniff::KMer{0, 622},  sniff::KMer{5, 332},  sniff::KMer{11, 325},
    sniff::KMer{14, 332}, sniff::KMer{15, 307}, sniff::KMer{21, 390}};

// kmer_length greate than the window_length
// kmer_length=15; window_length=5
TEST_CASE("minimizeK15W5", "[minimize]") {
  auto const minimizers =
      sniff::Minimize({.kmer_len = 15, .window_len = 5}, kTestSequence);
  REQUIRE(minimizers.size() == kTestExpectedMinimizersK15W5.size());
  for (std::size_t i = 0; i < minimizers.size(); ++i) {
    CHECK(kTestExpectedMinimizersK15W5[i] == minimizers[i]);
  }
}

// kmer length is equal to the window length
// kmer_length=7; window_length=7
TEST_CASE("minimizeK7W7", "[minimize]") {
  auto const minimizers =
      sniff::Minimize({.kmer_len = 7, .window_len = 7}, kTestSequence);
  REQUIRE(minimizers.size() == kTestExpectedMinimizersK7W7.size());
  for (std::size_t i = 0; i < minimizers.size(); ++i) {
    CHECK(kTestExpectedMinimizersK7W7[i] == minimizers[i]);
  }
}

// kmer length is less than the window length
// kmer_length=7; window_length=7
TEST_CASE("minimizeK5W7", "[minimize]") {
  auto const minimizers = sniff::Minimize({.kmer_len=5, .window_len=7}, kTestSequence);
  REQUIRE(minimizers.size() == kTestExpectedMinimizersK5W7.size());
  for (std::size_t i = 0; i < minimizers.size(); ++i) {
    CHECK(kTestExpectedMinimizersK5W7[i] == minimizers[i]);
  }
}
