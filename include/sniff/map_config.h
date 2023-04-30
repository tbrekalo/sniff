#pragma once

#include <cstdint>

namespace sniff {

struct MapConfig {
  std::uint32_t min_chain_length;
  std::uint32_t max_chain_gap_length;
  std::uint32_t kmer_len;
};

}  // namespace sniff
