#pragma once

#include <span>

#include "sniff/map_config.h"
#include "sniff/match.h"
#include "sniff/overlap.h"

namespace sniff {

auto Map(MapConfig cfg, std::span<Match const> matches) -> std::vector<Overlap>;

}  // namespace sniff
