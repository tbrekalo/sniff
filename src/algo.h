#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace biosoup {
class NucleicAcid;
}

namespace sniff {

struct Config {
  double p;
  std::uint32_t length;
  std::uint32_t max_edit_distance;
};

auto FindReverseComplementPairs(
    Config cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::pair<std::string, std::string>>;

}  // namespace sniff
