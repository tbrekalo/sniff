#include "sniff/algo.h"

#include <optional>
#include <variant>

// 3rd party
#include "ankerl/unordered_dense.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "fmt/core.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_sort.h"

// sniff
#include "sniff/map.h"
#include "sniff/match.h"
#include "sniff/minimize.h"
#include "sniff/sketch.h"

namespace sniff {

static constexpr auto kChunkSize = 1ULL << 32LLU;  // 1 GiB

template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};

template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

struct Target {
  std::uint32_t read_id;
  std::uint32_t read_len;

  KMer kmer;
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

static auto ExtractRcMinimizersSortedByVal(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads)
    -> std::vector<Target> {
  auto dst = std::vector<Target>();
  auto const minimize_cfg =
      MinimizeConfig{.kmer_len = cfg.minimize_cfg.kmer_len,
                     .window_len = cfg.minimize_cfg.window_len,
                     .minhash = false};

  auto cnt = std::atomic_size_t(0);
  auto sketches = std::vector<std::vector<Target>>(reads.size());
  tbb::parallel_for(
      std::size_t(0), reads.size(),
      [&reads, &minimize_cfg, &cnt, &sketches](std::size_t idx) -> void {
        auto rc_string = std::string(reads[idx]->InflateData());
        for (auto i = 0U; i < rc_string.size(); ++i) {
          rc_string[rc_string.size() - 1 - i] =
              biosoup::kNucleotideDecoder[3 ^ reads[idx]->Code(i)];
        }

        auto const rc_kmers = Minimize(minimize_cfg, rc_string);
        cnt += rc_kmers.size();
        for (auto const kmer : rc_kmers) {
          sketches[idx].push_back(Target{.read_id = reads[idx]->id,
                                         .read_len = reads[idx]->inflated_len,
                                         .kmer = kmer});
        }
      });

  dst.reserve(cnt);
  for (auto& sketch : sketches) {
    dst.insert(dst.end(), sketch.begin(), sketch.end());
    std::vector<Target>{}.swap(sketch);
  }

  std::sort(dst.begin(), dst.end(),
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

// Rc stands for "reverse complement"
static auto CreateRcKMerIndex(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> target_reads)
    -> Index {
  auto target_kmers = ExtractRcMinimizersSortedByVal(cfg, target_reads);
  auto index = IndexKMers(target_kmers);

  return {.locations = std::move(index), .kmers = std::move(target_kmers)};
}

static auto MapMatches(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    std::vector<Match> matches) -> std::optional<Overlap> {
  auto target_intervals = std::vector<std::uint32_t>{0};
  for (std::uint32_t i = 0; i < matches.size(); ++i) {
    if (i + 1 == matches.size() ||
        matches[i].target_id != matches[i + 1].target_id) {
      target_intervals.push_back(i + 1);
    }
  }

  auto dst_overlaps = std::vector<std::optional<Overlap>>(matches.size());
  tbb::parallel_for(std::size_t(0), target_intervals.size() - 1,
                    [&cfg, &matches, &target_intervals,
                     &dst_overlaps](std::size_t read_idx) -> void {
                      auto local_matches = std::span(
                          matches.begin() + target_intervals[read_idx],
                          matches.begin() + target_intervals[read_idx + 1]);
                      auto local_overlaps = Map(cfg.map_cfg, local_matches);
                      if (local_overlaps.size() == 1) {
                        dst_overlaps[read_idx] = local_overlaps.front();
                      }
                    });

  auto max_len = 0;
  auto dst = std::optional<Overlap>();
  for (auto const& ovlp : dst_overlaps) {
    if (ovlp &&
        // query internal overlap
        !((ovlp->query_start >
               0.125 * query_reads[ovlp->query_id]->inflated_len &&
           ovlp->query_end <
               0.875 * query_reads[ovlp->query_id]->inflated_len) ||
          // target internal overlap
          (ovlp->target_start >
               0.125 * query_reads[ovlp->target_id]->inflated_len &&
           ovlp->target_end <
               0.875 * query_reads[ovlp->target_id]->inflated_len)) &&
        OverlapLength(*ovlp) > max_len) {
      max_len = OverlapLength(*ovlp);
      dst = ovlp;
    }
  }

  return dst;
}

static auto MapSketchToIndex(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    Sketch const& sketch, KMerLocIndex const& index, double threshold)
    -> std::optional<Overlap> {
  auto const min_short_long_ratio = 1.0 - cfg.p;
  auto read_matches = std::vector<Match>();
  auto const try_match =
      [min_short_long_ratio, &query_sketch = sketch, &read_matches, threshold](
          KMer const& query_kmer, Target const& target) -> void {
    if (query_sketch.read_identifier.read_id >= target.read_id) {
      return;
    }

    if (auto const len_ratio =
            1. * std::min(query_sketch.read_len, target.read_len) /
            std::max(query_sketch.read_len, target.read_len);
        len_ratio < min_short_long_ratio) {
      return;
    }

    read_matches.push_back(
        Match{.query_id = query_sketch.read_identifier.read_id,
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
                     for (std::size_t i = 0; i < count; ++i) {
                       try_match(query_kmer, *(target_kmers_ptr + i));
                     }
                   }},
        cl->second.value);
  }

  std::sort(read_matches.begin(), read_matches.end(),
            [](Match const& lhs, Match const& rhs) -> bool {
              return lhs.target_id < rhs.target_id;
            });

  return MapMatches(cfg, query_reads, std::move(read_matches));
}

static auto MapSpanToIndex(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    KMerLocIndex const& target_index, double threshold)
    -> std::vector<std::optional<Overlap>> {
  auto const minimize_cfg =
      MinimizeConfig{.kmer_len = cfg.minimize_cfg.kmer_len,
                     .window_len = cfg.minimize_cfg.window_len,
                     .minhash = false};

  auto dst = std::vector<std::optional<Overlap>>(query_reads.size());
  tbb::parallel_for(
      std::size_t(0), query_reads.size(),
      [&cfg, query_reads, &target_index, threshold, &minimize_cfg,
       &dst](std::size_t idx) {
        auto sketch =
            Sketch{.read_identifier =
                       ReadIdentifier{.read_id = query_reads[idx]->id,
                                      .read_name = query_reads[idx]->name},
                   .read_len = query_reads[idx]->inflated_len,
                   .minimizers =
                       Minimize(minimize_cfg, query_reads[idx]->InflateData())};

        dst[idx] =
            MapSketchToIndex(cfg, query_reads, sketch, target_index, threshold);
      });

  return dst;
}

auto FindReverseComplementPairs(
    Config const& cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::pair<std::string, std::string>> {
  auto overlaps = std::vector<std::optional<Overlap>>(reads.size());
  auto dst = std::vector<std::pair<std::string, std::string>>();

  auto timer = biosoup::Timer{};
  timer.Start();

  auto batch_sz = 0ULL;
  for (std::uint32_t i = 0, j = i; j < reads.size(); ++j) {
    batch_sz += reads[j]->inflated_len;
    if (batch_sz < kChunkSize && j + 1 < reads.size()) {
      continue;
    }

    auto target_index = CreateRcKMerIndex(
        cfg, std::span(reads.cbegin() + i, reads.cbegin() + j + 1));
    auto query_reads = std::span(reads.cbegin(), reads.cbegin() + j + 1);

    auto batch_overlaps =
        MapSpanToIndex(cfg, query_reads, target_index.locations,
                       GetFrequencyThreshold(target_index.locations, 0.0002));

    for (std::uint32_t k = 0; k <= j; ++k) {
      if (!overlaps[k] ||
          (batch_overlaps[k] &&
           OverlapLength(*batch_overlaps[k]) > OverlapLength(*overlaps[k])))
        overlaps[k] = batch_overlaps[k];
    }

    i = j + 1;
    batch_sz = 0;

    fmt::print(
        stderr,
        "\r[sniff::FindReverseComplementPairs]({:12.3f}) maped {:2.3f}% reads",
        timer.Lap(), 100. * (j + 1) / reads.size());
  }

  for (auto opt_ovlp : overlaps) {
    if (opt_ovlp) {
      if (opt_ovlp->query_id > opt_ovlp->target_id) {
        opt_ovlp = ReverseOverlap(*opt_ovlp);
      }

      dst.emplace_back(reads[opt_ovlp->query_id]->name,
                       reads[opt_ovlp->target_id]->name);

      if (dst.back().first > dst.back().second) {
        std::swap(dst.back().first, dst.back().second);
      }
    }
  }

  std::sort(dst.begin(), dst.end());
  dst.erase(std::unique(dst.begin(), dst.end()), dst.end());

  fmt::print(stderr,
             "\r[sniff::FindReverseComplementPairs]({:12.3f}) found {} reverse "
             "complements\n",
             timer.Stop(), dst.size());

  return dst;
}

}  // namespace sniff
