#pragma once

#include <span>

#include "sniff/config.h"
#include "sniff/match.h"
#include "sniff/overlap.h"

namespace sniff {

struct MapConfig {
  std::uint32_t min_chain_length = 4;
  std::uint32_t max_chain_gap_length = 100;
  std::uint32_t kmer_len;
};

auto Map(MapConfig cfg, std::span<Match const> matches) -> std::vector<Overlap>;

}  // namespace sniff
