#pragma once

#include <compare>
#include <cstdint>
#include <string>

namespace sniff {

struct Overlap {
  std::uint32_t query_id;
  std::uint32_t query_length;
  std::uint32_t query_start;
  std::uint32_t query_end;

  std::uint32_t target_id;
  std::uint32_t target_length;
  std::uint32_t target_start;
  std::uint32_t target_end;

  friend constexpr auto operator<=>(const Overlap& lhs,
                                    const Overlap& rhs) = default;
};

struct OverlapNamed {
  std::string query_name;
  std::uint32_t query_length;
  std::uint32_t query_start;
  std::uint32_t query_end;

  std::string target_name;
  std::uint32_t target_length;
  std::uint32_t target_start;
  std::uint32_t target_end;

  friend auto operator<=>(const OverlapNamed& lhs,
                          const OverlapNamed& rhs) = default;
};

auto ReverseOverlap(Overlap const& ovlp) -> Overlap;
auto OverlapLength(Overlap const& ovlp) -> std::uint32_t;
auto OverlapRatio(Overlap const& ovlp) -> double;
auto OverlapError(Overlap const& ovlp) -> double;

}  // namespace sniff
