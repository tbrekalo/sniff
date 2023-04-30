#pragma once

#include <cstdint>

namespace sniff {

struct MinimizeConfig {
  std::uint32_t kmer_len = 15;
  std::uint32_t window_len = 5;
};

}  // namespace sniff
