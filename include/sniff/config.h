#pragma once

#include "sniff/map_config.h"
#include "sniff/minimize_config.h"

namespace sniff {

struct Config {
  double alpha_p;
  double beta_p;
  MapConfig map_cfg;
  MinimizeConfig minimize_cfg;
};

}  // namespace sniff
