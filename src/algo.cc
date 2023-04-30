#include "sniff/algo.h"

#include <optional>

// 3rd party
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "fmt/core.h"
#include "tbb/parallel_for.h"

// sniff
#include "sniff/map.h"
#include "sniff/match.h"
#include "sniff/minimize.h"

namespace sniff {

static auto FindIndicesForLengthRelatedReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::uint32_t short_over_long_ratio_threshold) {
  auto related = std::vector<std::vector<std::size_t>>(reads.size());
  for (std::size_t i = 0; i < reads.size(); ++i) {
    for (std::size_t j = i + 1; j < reads.size(); ++j) {
      auto const len_short =
          std::min(reads[i]->inflated_len, reads[j]->inflated_len);
      auto const len_long =
          std::max(reads[i]->inflated_len, reads[j]->inflated_len);

      if (auto const ratio = static_cast<double>(len_short) / len_long;
          ratio > short_over_long_ratio_threshold) {
        related[i].push_back(j);
        related[j].push_back(i);
      }
    }
  }

  return related;
}

static auto ExtractReverseComplementQueries(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::size_t> const& indices, std::uint32_t max_sample_len)
    -> std::vector<std::string> {
  auto dst = std::vector<std::string>(indices.size());
  for (std::size_t i = 0; i < indices.size(); ++i) {
    auto tmp = biosoup::NucleicAcid("", reads[indices[i]]->InflateData());
    tmp.ReverseAndComplement();

    dst.push_back(tmp.InflateData(std::min(tmp.inflated_len, max_sample_len)));
  }

  return dst;
}

static auto FindLongestOverlap(AlgoConfig const& cfg,
                               std::vector<std::string> const& queries,
                               std::string_view target)
    -> std::optional<std::size_t> {
  std::optional<std::size_t> dst = std::nullopt;
  auto dst_ovlp = Overlap{};

  auto const calc_ovlp_target_len = [](Overlap const& ovlp) -> std::uint32_t {
    return ovlp.target_end - ovlp.target_start;
  };

  auto const target_sketch = Minimize(cfg.minimize_cfg, target);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    auto overlaps = sniff::Map(
        cfg.map_cfg,
        sniff::MakeMatches(sniff::Minimize(cfg.minimize_cfg, queries[i]),
                           target_sketch));

    if (overlaps.size() != 1) {
      continue;
    }

    if (auto const ovlp_len = calc_ovlp_target_len(overlaps.front());
        static_cast<double>(ovlp_len) / target.length() > 0.98) {
      if (ovlp_len > calc_ovlp_target_len(dst_ovlp)) {
        dst = i;
      }
    }
  }

  return dst;
}

auto FindReverseComplementPairs(
    AlgoConfig cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::pair<std::string, std::string>> {
  /* clang-format off */
  fmt::print(stderr,
             "[sniff::FindReverseComplementPairs]\n"
             "\tp: {}; len: {};\n"
             "\tk: {}; w: {}; chain: {}; gap: {};\n",
             cfg.p, cfg.length,
             cfg.minimize_cfg.kmer_len, cfg.minimize_cfg.window_len,
             cfg.map_cfg.min_chain_length, cfg.map_cfg.max_chain_gap_length);
  /* clang-format on */

  auto matches = std::vector<std::optional<std::size_t>>(reads.size());
  auto dst = std::vector<std::pair<std::string, std::string>>();
  auto timer = biosoup::Timer{};
  timer.Start();

  auto length_related_reads =
      FindIndicesForLengthRelatedReads(reads, 1.0 - cfg.p);

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

  tbb::parallel_for(std::size_t(0), length_related_reads.size(),
                    [cfg, &reads, &matches, &length_related_reads, &n_mapped,
                     &report_status](std::size_t idx) -> void {
                      auto queries = ExtractReverseComplementQueries(
                          reads, length_related_reads[idx], cfg.length);
                      auto target = reads[idx]->InflateData(0, cfg.length);
                      matches[idx] = FindLongestOverlap(cfg, queries, target);

                      ++n_mapped;
                      report_status();
                    });

  for (std::size_t i = 0; i < matches.size(); ++i) {
    if (matches[i]) {
      auto& pss = dst.emplace_back(reads[i]->name, reads[*matches[i]]->name);
      if (pss.first > pss.second) {
        std::swap(pss.first, pss.second);
      }
    }
  }

  std::sort(dst.begin(), dst.end());
  dst.erase(std::unique(dst.begin(), dst.end()), dst.end());

  fmt::print(stderr, "[sniff::FindReverseComplementPairs]({:12.3f})\n",
             timer.Stop());

  return dst;
}

}  // namespace sniff
