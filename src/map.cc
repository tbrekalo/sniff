#include "sniff/map.h"

#include <algorithm>
#include <iterator>
#include <limits>

static auto CmpMatchByTargetPos(sniff::Match const& lhs,
                                sniff::Match const& rhs) -> bool {
  return lhs.target_pos < rhs.target_pos;
}

namespace sniff {

static auto FindLongestQueryChain(std::vector<Match>::const_iterator first,
                                  std::vector<Match>::const_iterator last)
    -> std::vector<Match> {
  auto n = static_cast<std::uint32_t>(last - first);

  auto prev = std::vector<std::uint32_t>(n + 1, n);
  auto chain = std::vector<std::uint32_t>{n, 0};

  chain.reserve(n + 1);
  for (std::uint32_t match_idx = 1; match_idx < n; ++match_idx) {
    auto idx = std::lower_bound(
                   chain.begin() + 1, chain.end(), *(first + match_idx),
                   [first](std::uint32_t index, Match const& match) -> bool {
                     return (first + index)->query_pos < match.query_pos;
                   }) -
               chain.begin();

    if (idx == chain.size()) {
      chain.emplace_back();
    }

    chain[idx] = match_idx;
    prev[match_idx] = chain[idx - 1];
  }

  auto dst = std::vector<Match>(chain.size() - 1);
  for (std::uint32_t dst_idx = chain.size() - 2, curr_idx = chain.back();
       curr_idx != n; --dst_idx) {
    dst[dst_idx] = *(first + curr_idx);
    curr_idx = prev[curr_idx];
  }

  return dst;
};

auto Map(MapConfig cfg, std::vector<Match> matches) -> std::vector<Overlap> {
  std::sort(matches.begin(), matches.end(), CmpMatchByTargetPos);
  matches.push_back(
      Match{.query_pos = std::numeric_limits<std::uint32_t>::max(),
            .target_pos = std::numeric_limits<std::uint32_t>::max()});

  auto dst = std::vector<Overlap>();
  for (std::size_t i = 1, j = 0; i < matches.size(); ++i) {
    if (matches[i].target_pos - matches[j].target_pos >
        cfg.max_chain_gap_length) {
      if (i - j >= cfg.min_chain_length) {
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
