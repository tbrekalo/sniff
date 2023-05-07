#pragma once

#include <cstdint>

namespace sniff {

struct KMer {
  std::uint32_t read_id;
  std::uint32_t position;
  std::uint64_t value;

  friend constexpr auto operator==(KMer const& lhs, KMer const& rhs) noexcept
      -> bool {
    return lhs.position == rhs.position && lhs.value == rhs.value;
  }

  friend constexpr auto operator!=(KMer const& lhs, KMer const& rhs) noexcept
      -> bool {
    return !(lhs == rhs);
  }
};

}  // namespace sniff
