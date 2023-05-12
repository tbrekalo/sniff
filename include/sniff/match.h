#pragma once

#include <compare>
#include <vector>

#include "sniff/kmer.h"

namespace sniff {

struct Match {
  std::uint32_t query_id;
  std::uint32_t query_pos;

  std::uint32_t target_id;
  std::uint32_t target_pos;

  friend constexpr auto operator<=>(Match const& lhs,
                                    Match const& rhs) = default;
};

auto MakeMatches(std::vector<KMer> query_sketch,
                 std::vector<KMer> target_sketch) -> std::vector<Match>;

}  // namespace sniff
