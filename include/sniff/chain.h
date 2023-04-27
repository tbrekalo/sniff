#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "sniff/kmer.h"
#include "sniff/overlap.h"

namespace sniff {

struct ChainConfig {
  std::uint32_t min_target_chain_matches;
  std::uint32_t max_target_allowed_gap;
  std::uint32_t kmer_len;
};

auto Chain(ChainConfig cfg, std::vector<KMer> query_sketch,
           std::vector<KMer> target_sketch) -> std::vector<Overlap>;

}  // namespace sniff
