#pragma once

#include "sniff/match.h"
#include "sniff/overlap.h"

namespace sniff {

struct ChainConfig {
  std::uint32_t min_target_chain_matches;
  std::uint32_t max_target_allowed_gap;
  std::uint32_t kmer_len;
};

auto Chain(ChainConfig cfg, std::vector<Match> matches) -> std::vector<Overlap>;

}  // namespace sniff
