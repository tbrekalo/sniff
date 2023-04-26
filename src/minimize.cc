#include "sniff/minimize.h"

#include <deque>

/* clang-format off */
constexpr static std::uint8_t kNucleotideCoder[] = {
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255,   0, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255,   0,   1 ,  1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255,
    255,   0,   1,   1,   0, 255, 255,   2,
      3, 255, 255,   2, 255,   1,   0, 255,
    255, 255,   0,   1,   3,   3,   2,   0,
    255,   3, 255, 255, 255, 255, 255, 255
};

constexpr static char kNucleotideDecoder[] = {
    'A', 'C', 'G', 'T'
};

/* clang-format on */

auto Hash(std::uint64_t val, std::uint64_t const kMask) -> std::uint64_t {
  val = ((~val) + (val << 21)) & kMask;
  val = val ^ (val >> 24);
  val = ((val + (val << 3)) + (val << 8)) & kMask;
  val = val ^ (val >> 14);
  val = ((val + (val << 2)) + (val << 4)) & kMask;
  val = val ^ (val >> 28);
  val = (val + (val << 31)) & kMask;
  return val;
}

namespace sniff {

auto Minimize(std::string_view sequence, std::uint32_t kmer_len,
              std::uint32_t window_len) -> std::vector<KMer> {
  auto dst = std::vector<KMer>();

  auto const mask =
      (1ULL << (static_cast<std::uint64_t>(kmer_len) * 2U)) - 1ULL;
  auto const shift_kmer = [mask](std::uint64_t kmer_val,
                                 char base) -> std::uint64_t {
    return ((kmer_val << 2ULL) | kNucleotideCoder[base]) & mask;
  };

  auto window = std::deque<std::pair<std::uint64_t, KMer>>();
  auto const window_push = [&window](std::uint64_t hash, KMer kmer) -> void {
    while (!window.empty() && window.back().first > hash) {
      window.pop_back();
    }

    window.emplace_back(hash, kmer);
  };

  auto const window_update = [&window](std::size_t window_start) {
    if (!window.empty() && window.front().second.position < window_start) {
      window.pop_front();
    }
  };

  auto kmer = std::uint64_t{};
  for (std::uint32_t i = 0; i < sequence.size(); ++i) {
    kmer = shift_kmer(kmer, sequence[i]);
    if (i >= kmer_len + window_len - 1) {
      window_update(i - (window_len + kmer_len - 1));
      if (dst.empty() || dst.back() != window.front().second) {
        dst.push_back(window.front().second);
      }
    }
    if (i >= kmer_len - 1) {
      window_push(Hash(kmer, mask),
                  KMer{.position = i - (kmer_len - 1), .value = kmer});
    }
  }

  return dst;
}

}  // namespace sniff
