#include "sniff/algo.h"

#include "biosoup/nucleic_acid.hpp"

namespace sniff {

static auto FindIndicesForLengthRelatedReads(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::uint32_t short_over_long_ratio_threshold) {
  auto related = std::vector<std::vector<std::size_t>>(reads.size());
  for (std::size_t i = 0; i < reads.size(); ++i) {
    for (std::size_t j = i + 1; j < reads.size(); ++j) {
      auto const len_short =
          std::min(reads[i]->inflated_len, reads[j]->inflated_len);
      auto const len_long =
          std::max(reads[i]->inflated_len, reads[j]->inflated_len);

      if (auto const ratio = static_cast<double>(len_short) / len_long;
          ratio > short_over_long_ratio_threshold) {
        related[i].push_back(j);
        related[j].push_back(i);
      }
    }
  }

  return related;
}

static auto ExtractReverseComplementSamples(
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> const& reads,
    std::vector<std::size_t> const& indices, std::uint32_t max_sample_len)
    -> std::vector<std::string> {
  auto dst = std::vector<std::string>(indices.size());
  for (std::size_t i = 0; i < indices.size(); ++i) {
    auto tmp = biosoup::NucleicAcid("", reads[indices[i]]->InflateData());
    tmp.ReverseAndComplement();

    dst.push_back(tmp.InflateData(std::min(tmp.inflated_len, max_sample_len)));
  }

  return dst;
}

// return indices for query samples that are similiar to target sample
static auto EvaluateSamples(std::string const& target,
                            std::vector<std::string> const& query_samples)
    -> std::vector<std::size_t> {
}

auto FindReverseComplementPairs(
    AlgoConfig cfg, std::vector<std::unique_ptr<biosoup::NucleicAcid>> reads)
    -> std::vector<std::pair<std::string, std::string>> {
  auto dst = std::vector<std::pair<std::string, std::string>>();
  auto length_related_reads =
      FindIndicesForLengthRelatedReads(reads, 1.0 - cfg.p);

  return dst;
}

}  // namespace sniff
