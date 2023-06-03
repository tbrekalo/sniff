#pragma once

#include <vector>

#include "sniff/kmer.h"

namespace sniff {

struct Sketch {
  std::uint32_t read_id;
  std::vector<KMer> minimizers;
};

}  // namespace sniff
