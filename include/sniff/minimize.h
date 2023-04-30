#pragma once

#include <string_view>
#include <vector>

#include "sniff/kmer.h"
#include "sniff/minimize_config.h"

namespace sniff {

auto Minimize(MinimizeConfig cfg, std::string_view sequence)
    -> std::vector<KMer>;

}  // namespace sniff
