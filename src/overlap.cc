#include "sniff/overlap.h"

#include <algorithm>

namespace sniff {

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
