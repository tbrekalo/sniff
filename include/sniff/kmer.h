#pragma once

#include <compare>
#include <cstdint>

namespace sniff {

struct KMer {
  std::uint32_t position;
  std::uint64_t value;

  friend constexpr auto operator<=>(KMer const& lhs, KMer const& rhs) = default;
};

}  // namespace sniff
