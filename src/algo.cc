#include "sniff/algo.h"

#include <functional>
#include <optional>
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

namespace sniff {

static constexpr auto kIndexSize = 1U << 30U;

template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};

template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

struct Target {
  std::uint32_t read_id;
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

static auto CreateRcString(std::unique_ptr<biosoup::NucleicAcid> const& read)
    -> std::string {
  auto dst = std::string(read->inflated_len, '\0');
  for (auto i = 0U; i < dst.size(); ++i) {
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
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads)
    -> std::vector<Target> {
  auto dst = std::vector<Target>();
  auto const minimize_cfg = MinimizeConfig{
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
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> target_reads)
    -> Index {
  auto target_kmers = ExtractRcMinimizersSortedByVal(cfg, target_reads);
  return {.locations = IndexKMers(target_kmers),
          .kmers = std::move(target_kmers)};
}

static auto IsStrongOverlap(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    Overlap src_ovlp) -> bool {
  using namespace std::placeholders;
  auto const get_read = std::bind(GetReadRefFromSpan, query_reads, _1);

  return 1. * (src_ovlp.query_end - src_ovlp.query_start) >
             cfg.beta_p * get_read(src_ovlp.query_id)->inflated_len &&
         1. * (src_ovlp.target_end - src_ovlp.target_start) >
             cfg.beta_p * get_read(src_ovlp.target_id)->inflated_len;
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
  tbb::parallel_for(
      std::size_t(0), target_intervals.size() - 1,
      [&cfg, query_reads, &matches, &target_intervals,
       &dst_overlaps](std::size_t read_idx) -> void {
        auto local_matches =
            std::span(matches.begin() + target_intervals[read_idx],
                      matches.begin() + target_intervals[read_idx + 1]);
        auto local_overlaps = Map({.min_chain_length = 4,
                                   .max_chain_gap_length = 800,
                                   .kmer_len = cfg.kmer_len},
                                  local_matches);

        if (local_overlaps.size() == 1 &&
            IsStrongOverlap(cfg, query_reads, local_overlaps.front())) {
          dst_overlaps[read_idx] = local_overlaps.front();
        }
      });

  auto max_len = 0;
  auto dst = std::optional<Overlap>();
  for (auto const& ovlp : dst_overlaps) {
    if (ovlp && OverlapLength(*ovlp) > max_len) {
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
  auto const min_short_long_ratio = 1.0 - cfg.alpha_p;
  auto read_matches = std::vector<Match>();

  using namespace std::placeholders;
  auto const get_read = std::bind(GetReadRefFromSpan, query_reads, _1);

  auto const try_match =
      [get_read, min_short_long_ratio, &query_sketch = sketch, &read_matches,
       threshold](KMer const& query_kmer, Target const& target) -> void {
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

    read_matches.push_back(Match{.query_id = query_sketch.read_id,
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
    -> std::vector<Overlap> {
  auto const minimize_cfg = MinimizeConfig{
      .kmer_len = cfg.kmer_len, .window_len = cfg.window_len, .minhash = false};

  auto opt_ovlps = std::vector<std::optional<Overlap>>(query_reads.size());
  tbb::parallel_for(
      std::size_t(0), query_reads.size(),
      [&cfg, query_reads, &target_index, threshold, &minimize_cfg,
       &opt_ovlps](std::size_t idx) {
        auto sketch =
            Sketch{.read_id = query_reads[idx]->id,
                   .minimizers =
                       Minimize(minimize_cfg, query_reads[idx]->InflateData())};

        opt_ovlps[idx] =
            MapSketchToIndex(cfg, query_reads, sketch, target_index, threshold);
      });

  auto dst = std::vector<Overlap>();
  for (auto const& opt_ovlp : opt_ovlps) {
    if (opt_ovlp) {
      dst.push_back(*opt_ovlp);
    }
  }

  return dst;
}

static auto MakeOverlapPairs(
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads,
    std::span<std::optional<Overlap> const> opt_ovlps) -> std::vector<RcPair> {
  auto dst = std::vector<RcPair>();
  auto match_ids =
      std::vector<std::pair<std::uint32_t, std::uint32_t>>(reads.size());
  for (auto i = 0U; i < match_ids.size(); ++i) {
    match_ids[i].first = i;
  }

  for (auto const& opt_ovlp : opt_ovlps) {
    if (opt_ovlp) {
      auto const ovlp_len = OverlapLength(*opt_ovlp);
      if (match_ids[opt_ovlp->query_id].first == opt_ovlp->query_id ||
          match_ids[opt_ovlp->query_id].second < ovlp_len) {
        match_ids[opt_ovlp->query_id] = {opt_ovlp->target_id, ovlp_len};
      }

      if (match_ids[opt_ovlp->target_id].first == opt_ovlp->target_id ||
          match_ids[opt_ovlp->target_id].second < ovlp_len) {
        match_ids[opt_ovlp->target_id] = {opt_ovlp->query_id, ovlp_len};
      }
    }
  }

  for (auto lhs = 0U; lhs < match_ids.size(); ++lhs) {
    auto rhs = match_ids[lhs].first;
    if (lhs == rhs || match_ids[rhs].first != lhs || lhs > rhs) {
      continue;
    }

    auto const lhs_str = reads[lhs]->InflateData();
    auto const rhs_rc_str = CreateRcString(reads[rhs]);

    dst.push_back(RcPair{.lhs = reads[lhs]->name, .rhs = reads[rhs]->name});
    if (dst.back().lhs > dst.back().rhs) {
      std::swap(dst.back().lhs, dst.back().rhs);
    }
  }

  return dst;
}

auto ReindexAndSortReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  std::sort(reads.begin(), reads.end(),
            [](std::unique_ptr<biosoup::NucleicAcid> const& lhs,
               std::unique_ptr<biosoup::NucleicAcid> const& rhs) -> bool {
              return lhs->inflated_len < rhs->inflated_len;
            });

  for (auto idx = 0U; idx < reads.size(); ++idx) {
    reads[idx]->id = idx;
  }

  return reads;
};

auto FindReverseComplementPairs(
    Config const& cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<RcPair> {
  auto opt_ovlps = std::vector<std::optional<Overlap>>(reads.size());
  reads = ReindexAndSortReads(std::move(reads));
  auto timer = biosoup::Timer{};
  timer.Start();

  auto const scale_len =
      [p = 1.0 - cfg.alpha_p](std::uint32_t read_len) -> std::uint32_t {
    return read_len * p;
  };

  auto prev_i = std::size_t(0);
  auto batch_size = std::size_t(0);
  auto const max_batch_size = kIndexSize;
  for (auto i = 0U, j = i; j < reads.size(); ++j) {
    batch_size += reads[j]->inflated_len;
    if (batch_size < max_batch_size && j + 1U < reads.size() &&
        scale_len(reads[j]->inflated_len) < reads[i]->inflated_len) {
      continue;
    }

    auto index = CreateRcKMerIndex(
        cfg, std::span(reads.cbegin() + i, reads.cbegin() + j));

    auto batch_ovlps = MapSpanToIndex(
        cfg, std::span(reads.cbegin() + prev_i, reads.cbegin() + j),
        index.locations,
        GetFrequencyThreshold(index.locations, cfg.filter_freq));

    for (auto& batch_ovlp : batch_ovlps) {
      if (!opt_ovlps[batch_ovlp.query_id] ||
          OverlapError(batch_ovlp) <
              OverlapError(*opt_ovlps[batch_ovlp.query_id])) {
        opt_ovlps[batch_ovlp.query_id] = batch_ovlp;
      }
    }

    fmt::print(stderr, "\r[FindReverseComplementPairs]({:12.3f}) {:2.3f}%",
               timer.Lap(), 100. * j / reads.size());
    batch_size = std::size_t(0);
    prev_i = i;
    i = j + 1;
  }

  auto dst = MakeOverlapPairs(reads, opt_ovlps);
  fmt::print(stderr, "\n[FindReverseComplementPairs]({:12.3f}) n pairs: {}\n",
             timer.Stop(), dst.size());

  return dst;
}

}  // namespace sniff
