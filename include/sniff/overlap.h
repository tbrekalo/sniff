#pragma once

#include <cstdint>

namespace sniff {

struct Overlap {
  std::uint32_t query_start;
  std::uint32_t query_end;

  std::uint32_t target_start;
  std::uint32_t target_end;
};

auto OverlapLength(Overlap const& ovlp) -> std::uint32_t;

}  // namespace sniff
