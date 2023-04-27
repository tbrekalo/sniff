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

static auto FindLongestChain(std::vector<Match>::const_iterator first,
                             std::vector<Match>::const_iterator last) {
  auto longest_idx = 0;
  auto chains = std::vector<std::vector<Match>>{};
  for (; first != last; ++first) {
    auto it = std::lower_bound(chains.begin(), chains.end(), first->target_pos,
                               [](std::vector<Match> const& match,
                                  std::uint32_t const target_pos) -> bool {
                                 return match.back().target_pos < target_pos;
                               });
    if (it == chains.end()) {
      chains.push_back({*first});
      if (it->size() > chains[longest_idx].size()) {
        longest_idx = it - chains.begin();
      }
    }
  }

  return chains[longest_idx];
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
    if (matches[j].target_pos - matches[i].target_pos >
        cfg.max_target_allowed_gap) {
      if (j - i > cfg.min_target_chain_matches) {
        auto chain = FindLongestChain(matches.begin() + i, matches.begin() + j);
        dst.push_back(
            Overlap{.query_start = chain.front().query_pos,
                    .query_end = chain.back().query_pos + cfg.kmer_len,

                    .target_start = chain.front().target_pos,
                    .target_end = chain.back().target_pos + cfg.kmer_len});
      }
    }
  }

  return dst;
}

}  // namespace sniff
