#pragma once

#include <string_view>
#include <vector>

#include "sniff/kmer.h"

namespace sniff {

auto Minimize(std::string_view sequence, std::uint32_t kmer_len,
              std::uint32_t window_len) -> std::vector<KMer>;

}  // namespace sniff
