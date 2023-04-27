#pragma once

#include <cstdint>

struct Overlap {
  std::uint32_t query_start;
  std::uint32_t query_end;

  std::uint32_t target_start;
  std::uint32_t target_end;
};
