#pragma once

#include <filesystem>

#include "sniff/config.h"
#include "sniff/sketch.h"

namespace sniff {

auto LoadSketches(Config cfg, std::filesystem::path const& path)
    -> std::vector<Sketch>;

}
