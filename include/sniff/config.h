#pragma once

#include "sniff/map_config.h"
#include "sniff/minimize_config.h"

namespace sniff {

struct Config {
  double p;
  MapConfig map_cfg;
  MinimizeConfig minimize_cfg;
};

}  // namespace sniff
