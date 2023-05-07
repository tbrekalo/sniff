#include "sniff/algo.h"

#include <optional>

// 3rd party
#include "ankerl/unordered_dense.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "fmt/core.h"
#include "tbb/parallel_for.h"

// sniff
#include "sniff/map.h"
#include "sniff/match.h"
#include "sniff/minimize.h"

namespace sniff {

static constexpr auto kChunkSize = 1U << 28U;  // 256 MiB

struct KMerInfo {
  std::uint32_t read_id;
  std::uint32_t read_len;
  std::uint32_t in_read_position;
};

auto FindReverseComplementPairs(Config cfg, std::vector<Sketch> read_sketches)
    -> std::vector<std::pair<std::string, std::string>> {
  auto matches = std::vector<std::optional<std::size_t>>(read_sketches.size());
  auto dst = std::vector<std::pair<std::string, std::string>>();
  auto timer = biosoup::Timer{};
  timer.Start();

  auto index =
      ankerl::unordered_dense::map<std::uint64_t, std::vector<KMerInfo>>();
  auto batch_sz = 0;
  for (std::size_t i = 0, j = i; j < read_sketches.size(); ++j) {
    batch_sz += read_sketches[j].rc_minimizers.size() * sizeof(KMer);
    if (batch_sz < kChunkSize && j + 1 < read_sketches.size()) {
      continue;
    }

    for (std::uint32_t k = i; k < j + 1; ++k) {
      for (auto const& kmer : read_sketches[k].rc_minimizers) {
        index[kmer.value].push_back(
            KMerInfo{.read_id = k,
                     .read_len = read_sketches[k].read_len,
                     .in_read_position = kmer.position});
      }
    }

    tbb::parallel_for(std::size_t(i), j, [&index](std::size_t k) -> void {

    });

    i = j;
    batch_sz = 0;
    index.clear();

    fmt::print(
        stderr,
        "\r[sniff::FindReverseComplementPairs]({:12.3f}) maped {:2.3f}% reads",
        timer.Lap(), 100. * j / read_sketches.size());
  }

  fmt::print(stderr, "\r[sniff::FindReverseComplementPairs]({:12.3f})\n",
             timer.Stop());

  return dst;
}

}  // namespace sniff
