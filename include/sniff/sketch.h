#pragma once

#include <string>
#include <vector>

#include "sniff/kmer.h"

namespace sniff {

struct ReadIdentifier {
  std::uint32_t read_id;
  std::string read_name;
};

struct Sketch {
  ReadIdentifier read_identifier;
  std::uint32_t read_len;

  std::vector<KMer> minimizers;
  std::vector<KMer> rc_minimizers;  // reverse complement
};

}  // namespace sniff
