#include "sniff/chain.h"

#include <algorithm>
#include <iterator>
#include <limits>

#include "sniff/minimize.h"

namespace sniff {

struct Match {
  std::uint32_t query_pos;
  std::uint32_t target_pos;

  friend constexpr auto operator==(Match const& lhs, Match const& rhs) -> bool {
    return lhs.query_pos == rhs.query_pos && lhs.target_pos == rhs.target_pos;
  }

  friend constexpr auto operator!=(Match const& lhs, Match const& rhs) -> bool {
    return !(lhs == rhs);
  }
};

static auto CmpKMerByValPos(KMer const& lhs, KMer const& rhs) -> bool {
  // NOTE: can be made faster with SIMD..
  return lhs.value != rhs.value ? lhs.value < rhs.value
                                : lhs.position < rhs.position;
}

static auto CmpMatchByTargetPos(Match const& lhs, Match const& rhs) -> bool {
  return lhs.target_pos < rhs.target_pos;
}

static auto FindLongestQueryChain(std::vector<Match>::const_iterator first,
                                  std::vector<Match>::const_iterator last) {
  auto const n = static_cast<std::uint32_t>(last - first);

  auto chain = std::vector<std::uint32_t>{0};
  auto prev = std::vector<std::uint32_t>(n, n);
  for (std::uint32_t match_idx = 1; match_idx < n; ++match_idx) {
    auto chain_idx =
        std::lower_bound(
            chain.begin(), chain.end(), (first + match_idx)->query_pos,
            [first](std::uint32_t chain_head, std::uint32_t query_pos) -> bool {
              return (first + chain_head)->query_pos < query_pos;
            }) -
        chain.begin();

    if (chain_idx != chain.size()) {
      prev[match_idx] = prev[chain[chain_idx]];
      chain[chain_idx] = match_idx;
    } else {
      prev[match_idx] = chain.back();
      chain.push_back(match_idx);
    }
  }

  auto dst = std::vector<Match>(chain.size());
  for (std::uint32_t dst_idx = dst.size() - 1, match_idx = chain.back();
       match_idx != n; --dst_idx) {
    dst[dst_idx] = *(first + match_idx);
    match_idx = prev[match_idx];
  }

  return dst;
};

auto Chain(ChainConfig cfg, std::vector<KMer> query_sketch,
           std::vector<KMer> target_sketch) -> std::vector<Overlap> {
  auto matches = std::vector<Match>();

  std::sort(query_sketch.begin(), query_sketch.end(), CmpKMerByValPos);
  std::sort(target_sketch.begin(), target_sketch.end(), CmpKMerByValPos);

  for (std::size_t i = 0, j = 0;
       i < query_sketch.size() && j < target_sketch.size(); ++i) {
    for (; i < query_sketch.size() && j < target_sketch.size(); ++j) {
      if (target_sketch[j].value < query_sketch[i].value) {
        continue;
      }

      if (target_sketch[j].value > query_sketch[i].value) {
        break;
      }

      matches.push_back(Match{.query_pos = query_sketch[i].position,
                              .target_pos = target_sketch[j].position});
      ++i;
    }
  }

  std::sort(matches.begin(), matches.end(), CmpMatchByTargetPos);
  matches.push_back(
      Match{.query_pos = std::numeric_limits<std::uint32_t>::max(),
            .target_pos = std::numeric_limits<std::uint32_t>::max()});

  auto dst = std::vector<Overlap>();
  for (std::size_t i = 1, j = 0; i < matches.size(); ++i) {
    if (matches[i].target_pos - matches[j].target_pos >
        cfg.max_target_allowed_gap) {
      if (i - j >= cfg.min_target_chain_matches) {
        auto chain =
            FindLongestQueryChain(matches.begin() + j, matches.begin() + i);
        dst.push_back(
            Overlap{.query_start = chain.front().query_pos,
                    .query_end = chain.back().query_pos + cfg.kmer_len,

                    .target_start = chain.front().target_pos,
                    .target_end = chain.back().target_pos + cfg.kmer_len});
      }
      j = i;
    }
  }

  return dst;
}

}  // namespace sniff
