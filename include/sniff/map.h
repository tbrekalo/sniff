#pragma once

#include "sniff/map_config.h"
#include "sniff/match.h"
#include "sniff/overlap.h"

namespace sniff {


auto Map(MapConfig cfg, std::vector<Match> matches) -> std::vector<Overlap>;

}  // namespace sniff
