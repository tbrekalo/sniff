#pragma once

#include <cstdint>

namespace sniff {

struct KMer {
  std::uint64_t position;
  std::uint64_t value;

  constexpr auto operator==(KMer const& that) const noexcept -> bool {
    return position == that.position && value == that.value;
  }

  constexpr auto operator!=(KMer const& that) const noexcept -> bool {
    return !(*this == that);
  }
};

}  // namespace sniff
