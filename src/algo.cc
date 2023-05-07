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

struct ReadLenIndex {
  std::uint32_t read_len;
  std::uint32_t read_index;
};

static auto FindIndicesForLengthRelatedReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::uint32_t short_over_long_ratio_threshold, std::uint32_t n_neighbors) {
  auto related = std::vector<std::vector<std::size_t>>(reads.size());
  auto len_indices = std::vector<ReadLenIndex>(reads.size());

  for (std::uint32_t i = 0; i < reads.size(); ++i) {
    len_indices[i] =
        ReadLenIndex{.read_len = reads[i]->inflated_len, .read_index = i};
  }

  std::sort(len_indices.begin(), len_indices.end(),
            [](ReadLenIndex const& lhs, ReadLenIndex const& rhs) -> bool {
              return lhs.read_len < rhs.read_len;
            });

  for (std::size_t i = 0; i < len_indices.size(); ++i) {
    for (std::size_t j = i + 1; j - i < n_neighbors && j < len_indices.size();
         ++j) {
      auto const [i_len, i_idx] = len_indices[i];
      auto const [j_len, j_idx] = len_indices[j];

      if (auto const len_ratio = static_cast<double>(j_len) / i_len;
          len_ratio < short_over_long_ratio_threshold) {
        break;
      }

      related[i_idx].push_back(j_idx);
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

    dst[i] = tmp.InflateData(0, std::min(tmp.inflated_len, max_sample_len));
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
    auto query_sketch = sniff::Minimize(cfg.minimize_cfg, queries[i]);
    auto overlaps =
        sniff::Map(cfg.map_cfg,
                   sniff::MakeMatches(std::move(query_sketch), target_sketch));

    if (overlaps.size() != 1) {
      continue;
    }

    if (auto const ovlp_len = calc_ovlp_target_len(overlaps.front());
        static_cast<double>(ovlp_len) / target.length() > 1.0 - cfg.p) {
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
      FindIndicesForLengthRelatedReads(reads, 1.0 - cfg.p, cfg.n_neighbors);

  fmt::print(
      stderr,
      "[sniff::FindReverseComplementPairs]({:12.3f}) grouped candidates\n",
      timer.Stop());

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
  tbb::parallel_for(
      std::size_t(0), length_related_reads.size(),
      [cfg, &reads, &matches, &length_related_reads, &n_mapped,
       &report_status](std::size_t idx) -> void {
        auto queries = ExtractReverseComplementQueries(
            reads, length_related_reads[idx], cfg.length);
        auto target = reads[idx]->InflateData(0, cfg.length);
        matches[idx] = [&index_mapping = length_related_reads[idx]](
                           auto opt) -> std::optional<std::size_t> {
          if (opt) {
            return index_mapping[*opt];
          }

          return std::nullopt;
        }(FindLongestOverlap(cfg, queries, target));

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

  fmt::print(stderr, "\r[sniff::FindReverseComplementPairs]({:12.3f})\n",
             timer.Stop());

  return dst;
}

}  // namespace sniff
