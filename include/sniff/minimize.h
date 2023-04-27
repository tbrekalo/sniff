#pragma once

#include <string_view>
#include <vector>

#include "sniff/kmer.h"

namespace sniff {

struct MinimizeConfig {
  std::uint32_t kmer_len = 15;
  std::uint32_t window_len = 5;
};

auto Minimize(MinimizeConfig cfg, std::string_view sequence)
    -> std::vector<KMer>;

}  // namespace sniff
