#include "sniff/map.h"

#include <algorithm>
#include <limits>
#include <span>

static auto CmpMatchByTargetPos(sniff::Match const& lhs,
                                sniff::Match const& rhs) -> bool {
  return lhs.target_pos < rhs.target_pos;
}

namespace sniff {

static auto FindLongestQueryChain(std::span<Match const> matches)
    -> std::vector<Match> {
  auto n = static_cast<std::uint32_t>(matches.size());

  auto prev = std::vector<std::uint32_t>(n + 1, n);
  auto chain = std::vector<std::uint32_t>{n, 0};

  chain.reserve(n + 1);
  for (std::uint32_t match_idx = 1; match_idx < n; ++match_idx) {
    auto idx = std::lower_bound(
                   chain.begin() + 1, chain.end(), matches[match_idx],
                   [matches](std::uint32_t index, Match const& match) -> bool {
                     return matches[index].query_pos < match.query_pos;
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
    dst[dst_idx] = matches[curr_idx];
    curr_idx = prev[curr_idx];
  }

  return dst;
};

auto Map(MapConfig cfg, std::span<Match const> src_matches)
    -> std::vector<Overlap> {
  auto matches = std::vector<Match>(src_matches.begin(), src_matches.end());
  std::sort(matches.begin(), matches.end(), CmpMatchByTargetPos);
  matches.push_back(
      Match{.query_pos = std::numeric_limits<std::uint32_t>::max(),
            .target_pos = std::numeric_limits<std::uint32_t>::max()});

  auto dst = std::vector<Overlap>();
  for (std::size_t i = 1, j = 0; i < matches.size(); ++i) {
    if (matches[i].target_pos - matches[i - 1].target_pos >
        cfg.max_chain_gap_length) {
      if (i - j >= cfg.min_chain_length) {
        auto chain = FindLongestQueryChain(
            std::span(matches.cbegin() + j, matches.cbegin() + i));
        dst.push_back(
            Overlap{.query_id = chain.front().query_id,
                    .query_start = chain.front().query_pos,
                    .query_end = chain.back().query_pos + cfg.kmer_len,

                    .target_id = chain.front().target_id,
                    .target_start = chain.front().target_pos,
                    .target_end = chain.back().target_pos + cfg.kmer_len,
                    .n_matches = static_cast<std::uint32_t>(chain.size())});
      }
      j = i;
    }
  }

  return dst;
}

}  // namespace sniff
