#pragma once

#include <cstdint>
#include <filesystem>

namespace sniff {

struct Config {
  double alpha_p;
  double beta_p;
  double filter_freq;
  std::uint32_t kmer_len;
  std::uint32_t window_len;

  std::filesystem::path log_dir;
};

}  // namespace sniff
