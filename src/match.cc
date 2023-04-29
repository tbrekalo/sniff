#include "sniff/match.h"

namespace sniff {

static auto CmpKMerByValPos(KMer const& lhs, KMer const& rhs) -> bool {
  return lhs.value != rhs.value ? lhs.value < rhs.value
                                : lhs.position < rhs.position;
}

auto MakeMatches(std::vector<KMer> query_sketch,
                 std::vector<KMer> target_sketch) -> std::vector<Match> {
  auto dst = std::vector<Match>();
  std::sort(query_sketch.begin(), query_sketch.end(), CmpKMerByValPos);
  std::sort(target_sketch.begin(), target_sketch.end(), CmpKMerByValPos);

  /* clang-format off */
  for (std::size_t query_idx = 0, target_idx = 0;
       query_idx < query_sketch.size() &&
       target_idx < target_sketch.size(); ++query_idx) {
    for (; query_idx < query_sketch.size() &&
            target_idx < target_sketch.size(); ++target_idx) {

      if (target_sketch[target_idx].value < query_sketch[query_idx].value) {
        continue;
      }

      if (target_sketch[target_idx].value > query_sketch[query_idx].value) {
        break;
      }

      dst.push_back(Match{.query_pos = query_sketch[query_idx].position,
                          .target_pos = target_sketch[target_idx].position});
      ++query_idx;
    }
  }
  /* clang-format on */

  return dst;
}

}  // namespace sniff
