#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "sniff/map_config.h"
#include "sniff/minimize_config.h"

namespace biosoup {
class NucleicAcid;
}

namespace sniff {

struct AlgoConfig {
  double p;
  std::uint32_t length;
  std::uint32_t n_neighbors;

  MapConfig map_cfg;
  MinimizeConfig minimize_cfg;
};

auto FindReverseComplementPairs(
    AlgoConfig cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::pair<std::string, std::string>>;

}  // namespace sniff
