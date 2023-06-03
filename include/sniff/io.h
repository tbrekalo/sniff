#pragma once

#include <filesystem>
#include <memory>
#include <vector>

#include "sniff/config.h"

namespace biosoup {
class NucleicAcid;
}

namespace sniff {

auto LoadReads(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>>;

}
