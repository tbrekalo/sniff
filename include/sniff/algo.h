#pragma once

#include "sniff/config.h"
#include "sniff/sketch.h"

namespace sniff {

auto FindReverseComplementPairs(Config cfg, std::vector<Sketch> sketches)
    -> std::vector<std::pair<std::string, std::string>>;

}  // namespace sniff
