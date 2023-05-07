#pragma once

#include <cstdint>

namespace sniff {

struct MapConfig {
  std::uint32_t min_chain_length = 4;
  std::uint32_t max_chain_gap_length = 100;
  std::uint32_t kmer_len;
};

}  // namespace sniff
