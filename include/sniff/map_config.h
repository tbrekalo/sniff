#pragma once

#include <cstdint>

namespace sniff {

struct MapConfig {
  std::uint32_t min_target_chain_matches;
  std::uint32_t max_target_allowed_gap;
  std::uint32_t kmer_len;
};

}  // namespace sniff
