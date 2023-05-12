#include "sniff/overlap.h"

#include <algorithm>

namespace sniff {

auto ReverseOverlap(Overlap const& ovlp) -> Overlap {
  return {
      .query_id = ovlp.target_id,
      .query_start = ovlp.target_start,
      .query_end = ovlp.target_end,

      .target_id = ovlp.query_id,
      .target_start = ovlp.query_start,
      .target_end = ovlp.query_end,
  };
}

auto OverlapLength(Overlap const& ovlp) -> std::uint32_t {
  return std::max(ovlp.query_end - ovlp.query_start,
                  ovlp.target_end - ovlp.target_start);
}

auto OverlapError(Overlap const& ovlp) -> double {
  return 1.0 -
         static_cast<double>(std::min(ovlp.query_end - ovlp.query_start,
                                      ovlp.target_end - ovlp.target_start)) /
             std::max(ovlp.query_end - ovlp.query_start,
                      ovlp.target_end - ovlp.target_start);
}

}  // namespace sniff
