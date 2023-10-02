#pragma once

#include <memory>
#include <string>
#include <vector>

#include "sniff/config.h"
#include "sniff/overlap.h"

namespace biosoup {
class NucleicAcid;
}

namespace sniff {

auto FindReverseComplementPairs(
    Config const& cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<OverlapNamed>;

}  // namespace sniff
