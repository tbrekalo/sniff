#pragma once

#include <vector>

#include "sniff/kmer.h"

namespace sniff {

struct Match {
  std::uint32_t query_pos;
  std::uint32_t target_pos;

  friend constexpr auto operator==(Match const& lhs, Match const& rhs) -> bool {
    return lhs.query_pos == rhs.query_pos && lhs.target_pos == rhs.target_pos;
  }

  friend constexpr auto operator!=(Match const& lhs, Match const& rhs) -> bool {
    return !(lhs == rhs);
  }
};

// returns vector of matches sorted by target position
auto MakeMatches(std::vector<KMer> query_sketch,
                 std::vector<KMer> target_sketch) -> std::vector<Match>;

}  // namespace sniff
