#pragma once

#include <memory>
#include <string>
#include <vector>

#include "sniff/config.h"

namespace biosoup {
class NucleicAcid;
}

namespace sniff {

struct RcPair {
  std::string lhs;
  std::string rhs;

  auto operator<=>(RcPair const&) const noexcept = default;
};

auto FindReverseComplementPairs(
    Config const& cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<RcPair>;

}  // namespace sniff
