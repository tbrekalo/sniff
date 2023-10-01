#include "sniff/algo.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <optional>
#include <type_traits>
#include <variant>

// 3rd party
#include "ankerl/unordered_dense.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "fmt/core.h"
#include "tbb/tbb.h"

// sniff
#include "sniff/map.h"
#include "sniff/match.h"
#include "sniff/minimize.h"
#include "sniff/sketch.h"

static constexpr auto kIndexSize = 1U << 30U;

static constexpr auto kIntercept = -23.47084474;

static constexpr auto kCoefs = std::tuple{
    9.42746909, -6.64572836, -2.78147289, 16.18407094, -7.31525403, -8.86853227,
};

template <class... Args>
requires((std::is_integral_v<Args> || std::is_floating_point_v<Args>) ||
         ...) static constexpr auto PredictIsValidOvlp(std::tuple<Args...> args)
    -> bool {
  auto x = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    return (kIntercept + ... + (std::get<Is>(kCoefs) * std::get<Is>(args)));
  }
  (std::make_index_sequence<sizeof...(Args)>{});
  return (1. / (1. + std::exp(-x))) >= 0.50;
}

template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};

template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

struct Target {
  std::uint32_t read_id;
  sniff::KMer kmer;
};

struct KMerLocator {
  std::uint32_t count;
  std::variant<Target, Target const*> value;
};

using KMerLocIndex = ankerl::unordered_dense::map<std::uint64_t, KMerLocator>;

struct Index {
  KMerLocIndex locations;
  std::vector<Target> kmers;
};

static auto FlattenOverlapVec(std::vector<std::vector<sniff::Overlap>> overlaps)
    -> std::vector<sniff::Overlap> {
  auto dst = std::vector<sniff::Overlap>();
  dst.reserve(std::transform_reduce(
      overlaps.begin(), overlaps.end(), std::size_t(0), std::plus<>{},
      [](std::vector<sniff::Overlap> const& vec) -> std::size_t {
        return vec.size();
      }));

  for (auto& it : overlaps) {
    dst.insert(dst.end(), it.begin(), it.end());
    std::vector<sniff::Overlap>{}.swap(it);
  }

  return dst;
}

static auto CreateRcString(std::unique_ptr<biosoup::NucleicAcid> const& read)
    -> std::string {
  auto dst = std::string(read->inflated_len, '\0');
  for (std::uint32_t i = 0; i < dst.size(); ++i) {
    dst[i] = biosoup::kNucleotideDecoder[3 ^ read->Code(dst.size() - 1 - i)];
  }

  return dst;
}

static auto GetReadRefFromSpan(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::uint32_t read_id) -> std::unique_ptr<biosoup::NucleicAcid> const& {
  return *std::lower_bound(
      reads.begin(), reads.end(), read_id,
      [](std::unique_ptr<biosoup::NucleicAcid> const& read,
         std::uint32_t read_id) -> bool { return read->id < read_id; });
}

static auto GetFrequencyThreshold(KMerLocIndex const& index, double freq)
    -> std::uint32_t {
  if (index.size() <= 2) {
    return 0U - 1;
  }
  auto counts = std::vector<std::uint32_t>();
  counts.reserve(index.size());

  for (auto it : index) {
    counts.push_back(it.second.count);
  }

  auto const idx = static_cast<std::size_t>(counts.size() * (1. - freq));
  std::nth_element(counts.begin(), counts.begin() + idx, counts.end());
  return counts[idx];
}

// RcMinimizers -> reverse complement minimizers
static auto ExtractRcMinimizersSortedByVal(
    sniff::Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads)
    -> std::vector<Target> {
  auto dst = std::vector<Target>();
  auto const minimize_cfg = sniff::MinimizeConfig{
      .kmer_len = cfg.kmer_len, .window_len = cfg.window_len, .minhash = false};

  auto cnt = std::atomic_size_t(0);
  auto sketches = std::vector<std::vector<Target>>(reads.size());

  tbb::parallel_for(
      std::size_t(0), reads.size(),
      [&reads, &minimize_cfg, &cnt, &sketches](std::size_t const idx) -> void {
        auto const rc_kmers =
            Minimize(minimize_cfg, CreateRcString(reads[idx]));
        cnt += rc_kmers.size();

        sketches[idx].reserve(rc_kmers.size());
        for (auto const kmer : rc_kmers) {
          sketches[idx].push_back(
              Target{.read_id = reads[idx]->id, .kmer = kmer});
        }
      });

  dst.reserve(cnt);
  for (auto& sketch : sketches) {
    dst.insert(dst.end(), sketch.begin(), sketch.end());
    std::vector<Target>{}.swap(sketch);
  }

  tbb::parallel_sort(dst.begin(), dst.end(),
                     [](Target const& lhs, Target const& rhs) -> bool {
                       return lhs.kmer.value < rhs.kmer.value;
                     });

  return dst;
}

static auto IndexKMers(std::span<Target const> target_kmers) -> KMerLocIndex {
  auto dst = KMerLocIndex();
  for (std::uint32_t i = 0, j = i; i < target_kmers.size(); ++j) {
    if (j < target_kmers.size() &&
        target_kmers[i].kmer.value == target_kmers[j].kmer.value) {
      continue;
    }

    auto& locator = dst[target_kmers[i].kmer.value];
    locator.count = j - i;

    if (locator.count == 1) {
      locator.value = target_kmers[i];
    } else {
      locator.value = std::addressof(target_kmers[i]);
    }

    i = j;
  }

  return dst;
};

// Rc stands for "reverse complement"
static auto CreateRcKMerIndex(
    sniff::Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> target_reads)
    -> Index {
  auto target_kmers = ExtractRcMinimizersSortedByVal(cfg, target_reads);
  return {.locations = IndexKMers(target_kmers),
          .kmers = std::move(target_kmers)};
}

// Assumes that the input is grouped by (query_id, target_id) pairs and overlaps
// in a group are non-overlapping and sorted by ascending (query, target)
// positions.
static auto MergeOverlaps(
    sniff::Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    std::span<sniff::Overlap const> overlaps) -> std::vector<sniff::Overlap> {
  auto dst = std::vector<sniff::Overlap>();
  if (overlaps.empty()) {
    return dst;
  }

  auto buff = std::vector<sniff::Overlap>{overlaps.front()};
  auto merge_overlaps =
      [&cfg, query_reads](
          std::vector<sniff::Overlap> ovlps) -> std::optional<sniff::Overlap> {
    if (ovlps.empty()) {
      return std::nullopt;
    }

    std::sort(ovlps.begin(), ovlps.end(),
              [](sniff::Overlap const& lhs, sniff::Overlap const& rhs) -> bool {
                return lhs.query_start < rhs.query_start;
              });

    if (not std::is_sorted(
            ovlps.begin(), ovlps.end(),
            [](sniff::Overlap const& lhs, sniff::Overlap const& rhs) -> bool {
              return lhs.target_start < rhs.target_start;
            })) {
      return std::nullopt;
    }

    auto const query_len = ovlps.front().query_length;
    auto const query_score = static_cast<double>(ovlps.back().query_end -
                                                 ovlps.front().query_start) /
                             query_len;

    auto const target_len = ovlps.front().target_length;
    auto const target_score = static_cast<double>(ovlps.back().target_end -
                                                  ovlps.front().target_start) /
                              ovlps.front().target_length;

    auto const query_lhs_overhang =
        static_cast<double>(ovlps.front().query_start) / query_len;
    auto const query_rhs_overhang =
        static_cast<double>(query_len - ovlps.back().query_end) / query_len;

    auto const target_lhs_overhang =
        static_cast<double>(ovlps.front().target_start) / target_len;
    auto const target_rhs_overhang =
        static_cast<double>(target_len - ovlps.front().target_end) / target_len;

    if (query_score > cfg.beta_p && target_score > cfg.beta_p &&
        PredictIsValidOvlp(std::tuple(
            query_score, query_lhs_overhang, query_rhs_overhang, target_score,
            target_lhs_overhang, target_rhs_overhang))) {
      return sniff::Overlap{
          .query_id = ovlps.front().query_id,
          .query_length = ovlps.front().query_length,
          .query_start = ovlps.front().query_start,
          .query_end = ovlps.back().query_end,

          .target_id = ovlps.front().target_id,
          .target_length = ovlps.front().target_length,
          .target_start = ovlps.front().target_start,
          .target_end = ovlps.back().target_end,
      };
    }

    return std::nullopt;
  };

  for (std::uint32_t i = 0, j = i = 1;
       i < overlaps.size() && j < overlaps.size(); ++j) {
    if ((overlaps[i].query_id == overlaps[j].query_id &&
         overlaps[i].target_id == overlaps[j].target_id)) {
      buff.push_back(overlaps[j]);
    } else {
      if (auto const opt_ovlp = merge_overlaps(buff); opt_ovlp) {
        dst.push_back(*opt_ovlp);
      }

      buff.clear();
      buff.push_back(overlaps[j]);

      i = j;
    }
  }

  if (auto const opt_ovlp = merge_overlaps(buff); opt_ovlp) {
    dst.push_back(*opt_ovlp);
  }

  return dst;
};

static auto MapMatches(
    sniff::Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    std::vector<sniff::Match> matches) -> std::vector<sniff::Overlap> {
  std::sort(matches.begin(), matches.end(),
            [](sniff::Match const& lhs, sniff::Match const& rhs) -> bool {
              return lhs.target_id < rhs.target_id;
            });

  auto target_intervals = std::vector<std::uint32_t>{0};
  for (std::uint32_t i = 0; i < matches.size(); ++i) {
    if (i + 1 == matches.size() ||
        matches[i].target_id != matches[i + 1].target_id) {
      target_intervals.push_back(i + 1);
    }
  }

  auto ovlps_buff = std::vector<std::vector<sniff::Overlap>>(matches.size());
  tbb::parallel_for(
      std::size_t(0), target_intervals.size() - 1,
      [&cfg, query_reads, &matches, &target_intervals,
       &ovlps_buff](std::size_t read_idx) -> void {
        auto local_matches =
            std::span(matches.begin() + target_intervals[read_idx],
                      matches.begin() + target_intervals[read_idx + 1]);

        ovlps_buff[read_idx] =
            [&cfg, query_reads](std::vector<sniff::Overlap> overlaps)
            -> std::vector<sniff::Overlap> {
          using namespace std::placeholders;
          auto const get_read = std::bind(GetReadRefFromSpan, query_reads, _1);
          for (auto& ovlp : overlaps) {
            ovlp.query_length = get_read(ovlp.query_id)->inflated_len;
            ovlp.target_length = get_read(ovlp.target_id)->inflated_len;
          }

          return MergeOverlaps(cfg, query_reads, overlaps);
        }(Map({.min_chain_length = 4,
               .max_chain_gap_length = 800,
               .kmer_len = cfg.kmer_len},
              local_matches));
      });

  return FlattenOverlapVec(std::move(ovlps_buff));
}

static auto MapSketchToIndex(
    sniff::Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    sniff::Sketch const& sketch, KMerLocIndex const& index, double threshold)
    -> std::vector<sniff::Overlap> {
  auto const min_short_long_ratio = 1.0 - cfg.alpha_p;
  auto read_matches = std::vector<sniff::Match>();

  using namespace std::placeholders;
  auto const get_read = std::bind(GetReadRefFromSpan, query_reads, _1);

  auto const try_match =
      [get_read, min_short_long_ratio, &query_sketch = sketch, &read_matches,
       threshold](sniff::KMer const& query_kmer, Target const& target) -> void {
    if (query_sketch.read_id >= target.read_id) {
      return;
    }
    auto const len_ratio =
        1. *
        std::min(get_read(query_sketch.read_id)->inflated_len,
                 get_read(target.read_id)->inflated_len) /
        std::max(get_read(query_sketch.read_id)->inflated_len,
                 get_read(target.read_id)->inflated_len);
    if (len_ratio < min_short_long_ratio) {
      return;
    }

    read_matches.push_back(sniff::Match{.query_id = query_sketch.read_id,
                                        .query_pos = query_kmer.position,
                                        .target_id = target.read_id,
                                        .target_pos = target.kmer.position});
  };

  for (auto const& query_kmer : sketch.minimizers) {
    auto const cl = index.find(query_kmer.value);
    if (cl == index.end() || cl->second.count >= threshold) {
      continue;
    }

    std::visit(
        overloaded{[&try_match, &query_kmer](Target const& targe_kmer) -> void {
                     try_match(query_kmer, targe_kmer);
                   },
                   [count = cl->second.count, &try_match,
                    &query_kmer](Target const* target_kmers_ptr) -> void {
                     for (std::uint32_t i = 0; i < count; ++i) {
                       try_match(query_kmer, *(target_kmers_ptr + i));
                     }
                   }},
        cl->second.value);
  }

  return MapMatches(cfg, query_reads, std::move(read_matches));
}

static auto MapSpanToIndex(
    sniff::Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    KMerLocIndex const& target_index, double threshold)
    -> std::vector<sniff::Overlap> {
  auto const minimize_cfg = sniff::MinimizeConfig{
      .kmer_len = cfg.kmer_len, .window_len = cfg.window_len, .minhash = false};

  auto ovlps_buff =
      std::vector<std::vector<sniff::Overlap>>(query_reads.size());
  tbb::parallel_for(std::size_t(0), query_reads.size(),
                    [&cfg, query_reads, &target_index, threshold, &minimize_cfg,
                     &ovlps_buff](std::size_t idx) {
                      auto sketch = sniff::Sketch{
                          .read_id = query_reads[idx]->id,
                          .minimizers = Minimize(
                              minimize_cfg, query_reads[idx]->InflateData())};

                      ovlps_buff[idx] = MapSketchToIndex(
                          cfg, query_reads, sketch, target_index, threshold);
                    });

  return FlattenOverlapVec(std::move(ovlps_buff));
}

auto SortReadsAndReindex(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  std::sort(reads.begin(), reads.end(),
            [](std::unique_ptr<biosoup::NucleicAcid> const& lhs,
               std::unique_ptr<biosoup::NucleicAcid> const& rhs) -> bool {
              return lhs->inflated_len < rhs->inflated_len;
            });

  for (std::uint32_t idx = 0; idx < reads.size(); ++idx) {
    reads[idx]->id = idx;
  }

  return reads;
};

namespace sniff {

auto FindReverseComplementPairs(
    Config const& cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<OverlapNamed> {
  reads = SortReadsAndReindex(std::move(reads));

  auto n_ovlps = std::size_t(0);
  auto ovlps_buff = std::vector<std::vector<sniff::Overlap>>();

  auto timer = biosoup::Timer{};
  timer.Start();

  auto const scale_len =
      [p = 1.0 - cfg.alpha_p](std::uint32_t read_len) -> std::uint32_t {
    return read_len * p;
  };

  auto prev_i = std::size_t(0);
  auto batch_size = std::size_t(0);
  auto const max_batch_size = kIndexSize;
  for (std::uint32_t i = 0, j = i; j < reads.size(); ++j) {
    batch_size += reads[j]->inflated_len;
    if (batch_size < max_batch_size && j + 1U < reads.size() &&
        scale_len(reads[j]->inflated_len) < reads[i]->inflated_len) {
      continue;
    }

    auto index = CreateRcKMerIndex(
        cfg, std::span(reads.cbegin() + i, reads.cbegin() + j));

    ovlps_buff.push_back(MapSpanToIndex(
        cfg, std::span(reads.cbegin() + prev_i, reads.cbegin() + j),
        index.locations,
        GetFrequencyThreshold(index.locations, cfg.filter_freq)));
    n_ovlps += ovlps_buff.back().size();

    fmt::print(stderr, "\r[FindReverseComplementPairs]({:12.3f}) {:2.3f}%",
               timer.Lap(), 100. * j / reads.size());
    batch_size = std::size_t(0);
    prev_i = i;
    i = j + 1;
  }

  fmt::print(stderr, "\n[FindReverseComplementPairs]({:12.3f}) n pairs: {}\n",
             timer.Stop(), n_ovlps);

  auto dst = std::vector<OverlapNamed>();
  dst.reserve(n_ovlps);

  for (auto& ovlp_vec : ovlps_buff) {
    for (auto ovlp : ovlp_vec) {
      dst.push_back(OverlapNamed{
          .query_name = reads[ovlp.query_id]->name,
          .query_length = ovlp.query_length,
          .query_start = ovlp.query_start,
          .query_end = ovlp.query_end,

          .target_name = reads[ovlp.target_id]->name,
          .target_length = ovlp.target_length,
          .target_start = ovlp.target_start,
          .target_end = ovlp.target_end,
      });
    }
  }

  return dst;
}

}  // namespace sniff
