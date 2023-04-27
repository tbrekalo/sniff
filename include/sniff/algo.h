#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace biosoup {
class NucleicAcid;
}

namespace sniff {

struct AlgoConfig {
  double p;
  std::uint32_t length;
  std::uint32_t max_edit_distance;
};

auto FindReverseComplementPairs(
    AlgoConfig cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::pair<std::string, std::string>>;

}  // namespace sniff
