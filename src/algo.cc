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

static constexpr auto kChunkSize = 1ULL << 30LLU;  // 1 GiB

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
  auto dst = std::string(read->InflateData());
  for (auto i = 0U; i < dst.size(); ++i) {
    dst[i] = biosoup::kNucleotideDecoder[3 ^ read->Code(dst.size() - 1 - i)];
  }

  return dst;
}

// RcMinimizers -> reverse complement minimizers
static auto ExtractRcMinimizersSortedByVal(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> reads)
    -> std::vector<Target> {
  auto dst = std::vector<Target>();
  auto const minimize_cfg =
      MinimizeConfig{.kmer_len = cfg.minimize_cfg.kmer_len,
                     .window_len = cfg.minimize_cfg.window_len,
                     .minhash = true};

  auto cnt = std::atomic_size_t(0);
  auto sketches = std::vector<std::vector<Target>>(reads.size());
  tbb::parallel_for(
      std::size_t(0), reads.size(),
      [&reads, &minimize_cfg, &cnt, &sketches](std::size_t idx) -> void {
        auto const rc_kmers =
            Minimize(minimize_cfg, CreateRcString(reads[idx]));
        cnt += rc_kmers.size();
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

static auto IsSpanningOverlap(
    Config const& cfg,
    std::span<std::unique_ptr<biosoup::NucleicAcid> const> query_reads,
    Overlap src_ovlp) -> bool {
  auto const minimize_cfg =
      MinimizeConfig{.kmer_len = cfg.minimize_cfg.kmer_len,
                     .window_len = cfg.minimize_cfg.window_len,
                     .minhash = false};
  auto matches = MakeMatches(
      Minimize(minimize_cfg, query_reads[src_ovlp.query_id]->InflateData()),
      Minimize(minimize_cfg, CreateRcString(query_reads[src_ovlp.target_id])));

  auto ovlps = Map(cfg.map_cfg, matches);
  if (ovlps.size() != 1) {
    return false;
  }

  if (auto const& detail_ovlp = ovlps.front();
      !(1. * (detail_ovlp.query_end - detail_ovlp.query_start) >
            cfg.beta_p * query_reads[src_ovlp.query_id]->inflated_len &&
        1. * (detail_ovlp.target_end - detail_ovlp.target_start) >
            cfg.beta_p * query_reads[src_ovlp.target_id]->inflated_len)) {
    return false;
  }

  return true;
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
        auto local_overlaps = Map(cfg.map_cfg, local_matches);
        local_overlaps.erase(
            std::remove_if(local_overlaps.begin(), local_overlaps.end(),
                           [cfg, query_reads](Overlap const& ovlp) -> bool {
                             return !IsSpanningOverlap(cfg, query_reads, ovlp);
                           }),
            local_overlaps.end());

        if (!local_overlaps.empty()) {
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
  auto const try_match =
      [query_reads, min_short_long_ratio, &query_sketch = sketch, &read_matches,
       threshold](KMer const& query_kmer, Target const& target) -> void {
    if (query_sketch.read_id >= target.read_id) {
      return;
    }

    if (auto const len_ratio =
            1. *
            std::min(query_reads[query_sketch.read_id]->inflated_len,
                     query_reads[target.read_id]->inflated_len) /
            std::max(query_reads[query_sketch.read_id]->inflated_len,
                     query_reads[target.read_id]->inflated_len);
        len_ratio < min_short_long_ratio) {
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
            Sketch{.read_id = query_reads[idx]->id,
                   .minimizers =
                       Minimize(minimize_cfg, query_reads[idx]->InflateData())};

        dst[idx] =
            MapSketchToIndex(cfg, query_reads, sketch, target_index, threshold);
      });

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

auto FindReverseComplementPairs(
    Config const& cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<RcPair> {
  auto opt_ovlps = std::vector<std::optional<Overlap>>(reads.size());
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
    auto batch_ovlps = MapSpanToIndex(
        cfg, query_reads, target_index.locations,
        GetFrequencyThreshold(target_index.locations, cfg.map_cfg.f));

    for (std::uint32_t k = 0; k <= j; ++k) {
      if (!opt_ovlps[k] || (batch_ovlps[k] && OverlapLength(*batch_ovlps[k]) >
                                                  OverlapLength(*opt_ovlps[k])))
        opt_ovlps[k] = batch_ovlps[k];
    }

    i = j + 1;
    batch_sz = 0;

    fmt::print(
        stderr,
        "\r[sniff::FindReverseComplementPairs]({:12.3f}) maped {:2.3f}% reads",
        timer.Lap(), 100. * (j + 1) / reads.size());
  }

  auto dst = MakeOverlapPairs(reads, opt_ovlps);
  fmt::print(stderr,
             "\r[sniff::FindReverseComplementPairs]({:12.3f}) found {} reverse "
             "complements\n",
             timer.Stop(), dst.size());

  return dst;
}

}  // namespace sniff
