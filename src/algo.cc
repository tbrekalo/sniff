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

namespace sniff {

static constexpr auto kChunkSize = 1ULL << 30LLU;  // 1 GiB

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

static auto ExtractRcMinimizersSortedByVal(std::span<Sketch const> sketches)
    -> std::vector<Target> {
  auto dst = std::vector<Target>();
  for (auto const& sketch : sketches) {
    for (auto const& kmer : sketch.rc_minimizers) {
      dst.push_back({.read_id = sketch.read_identifier.read_id,
                     .read_len = sketch.read_len,
                     .kmer = kmer});
    }
  }

  tbb::parallel_sort(dst.begin(), dst.end(),
                     [](Target const& lhs, Target const& rhs) -> bool {
                       return lhs.kmer.value < rhs.kmer.value;
                     });
  return dst;
}

static auto IndexKMers(std::span<Target const> kmers) -> KMerLocIndex {
  auto dst = KMerLocIndex();
  for (std::uint32_t i = 0, j = i; i < kmers.size(); ++j) {
    if (j < kmers.size() && kmers[i].kmer.value == kmers[j].kmer.value) {
      continue;
    }

    auto& locator = dst[kmers[i].kmer.value];
    locator.count = j - i;

    if (locator.count == 1) {
      locator.value = kmers[i];
    } else {
      locator.value = std::addressof(kmers[i]);
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

  std::nth_element(counts.begin(), counts.begin() + counts.size() * (1. - freq),
                   counts.end());

  auto idx = std::min(static_cast<std::size_t>(counts.size() * (1. - freq)) + 1,
                      counts.size() - 1);

  return counts[idx];
}

// Rc stands for "reverse complement"
static auto CreateRcKMerIndex(std::span<Sketch const> sketches) -> Index {
  auto kmers = ExtractRcMinimizersSortedByVal(sketches);
  auto index = IndexKMers(kmers);

  return {.locations = std::move(index), .kmers = std::move(kmers)};
}

static auto MapMatches(Config const& cfg, std::vector<Match> matches)
    -> std::optional<Overlap> {
  auto read_intervals = std::vector<std::uint32_t>{0};
  for (std::uint32_t i = 0; i < matches.size(); ++i) {
    if (i + 1 == matches.size() ||
        matches[i].target_id != matches[i + 1].target_id) {
      read_intervals.push_back(i + 1);
    }
  }

  auto dst_overlaps = std::vector<std::optional<Overlap>>(matches.size());
  tbb::parallel_for(std::size_t(0), read_intervals.size() - 1,
                    [&cfg, &matches, &read_intervals,
                     &dst_overlaps](std::size_t read_idx) -> void {
                      auto const n = read_intervals[read_idx + 1] -
                                     read_intervals[read_idx];

                      auto local_overlaps = Map(cfg.map_cfg, matches);
                      if (local_overlaps.size() == 1) {
                        dst_overlaps[read_idx] = local_overlaps.front();
                      }
                    });

  auto max_len = 0;
  auto dst = std::optional<Overlap>();
  for (auto const& ovlp : dst_overlaps) {
    if (ovlp &&
        // query internal overlap
        !((ovlp->query_start > 0.125 * cfg.sample_length &&
           ovlp->query_end < 0.875 * cfg.sample_length) ||
          // target internal overlap
          (ovlp->target_start > 0.125 * cfg.sample_length &&
           ovlp->target_end < 0.875 * cfg.sample_length)) &&
        OverlapLength(*ovlp) > max_len) {
      max_len = OverlapLength(*ovlp);
      dst = ovlp;
    }
  }

  return dst;
}

static auto MapSketchToIndex(Config const& cfg, Sketch const& sketch,
                             KMerLocIndex const& index, double threshold)
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

  return MapMatches(cfg, std::move(read_matches));
}

static auto MapSpanToIndex(Config const& cfg, std::span<Sketch const> sketches,
                           KMerLocIndex const& index, double threshold)
    -> std::vector<std::optional<Overlap>> {
  auto dst = std::vector<std::optional<Overlap>>(sketches.size());
  tbb::parallel_for(std::size_t(0), sketches.size(),
                    [&cfg, sketches, &index, threshold, &dst](std::size_t idx) {
                      dst[idx] = MapSketchToIndex(cfg, sketches[idx], index,
                                                  threshold);
                    });
  return dst;
}

auto FindReverseComplementPairs(Config cfg, std::vector<Sketch> sketches)
    -> std::vector<std::pair<std::string, std::string>> {
  auto overlaps = std::vector<std::optional<Overlap>>(sketches.size());
  auto dst = std::vector<std::pair<std::string, std::string>>();

  auto timer = biosoup::Timer{};
  timer.Start();

  auto batch_sz = 0ULL;
  for (std::uint32_t i = 0, j = i; j < sketches.size(); ++j) {
    batch_sz += sketches[j].rc_minimizers.size() * sizeof(KMer);
    if (batch_sz < kChunkSize && j + 1 < sketches.size()) {
      continue;
    }

    auto index = CreateRcKMerIndex(
        std::span(sketches.begin() + i, sketches.begin() + j + 1));
    auto batch_overlaps = MapSpanToIndex(
        cfg,
        std::span<Sketch const>(sketches.cbegin(), sketches.cbegin() + j + 1),
        index.locations, GetFrequencyThreshold(index.locations, 0.001));

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
        timer.Lap(), 100. * (j + 1) / sketches.size());
  }

  for (auto opt_ovlp : overlaps) {
    if (opt_ovlp) {
      if (opt_ovlp->query_id > opt_ovlp->target_id) {
        opt_ovlp = ReverseOverlap(*opt_ovlp);
      }

      dst.emplace_back(sketches[opt_ovlp->query_id].read_identifier.read_name,
                       sketches[opt_ovlp->target_id].read_identifier.read_name);

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
