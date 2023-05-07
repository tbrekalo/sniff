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

auto FindReverseComplementPairs(Config cfg, std::vector<Sketch> reads)
    -> std::vector<std::pair<std::string, std::string>> {
  auto matches = std::vector<std::optional<std::size_t>>(reads.size());
  auto dst = std::vector<std::pair<std::string, std::string>>();
  auto timer = biosoup::Timer{};
  timer.Start();

  auto n_mapped = std::atomic_uint32_t(0);
  auto report_ticket = std::atomic_uint32_t(0);

  auto report_status = [n_targets = reads.size(), &timer, &n_mapped,
                        &report_ticket]() -> void {
    auto const to_percent = [n_targets](std::uint32_t n) -> double {
      return (100. * n) / n_targets;
    };

    if (auto ticket = ++report_ticket; ticket == report_ticket) {
      fmt::print(stderr,
                 "\r[sniff::FindReverseComplementPairs]({:12.3f}) mapped "
                 "{:3.3f}% reads",
                 timer.Lap(), to_percent(n_mapped));
    }
  };

  timer.Start();

  fmt::print(stderr, "\r[sniff::FindReverseComplementPairs]({:12.3f})\n",
             timer.Stop());

  return dst;
}

}  // namespace sniff
