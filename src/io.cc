#include "sniff/io.h"

#include <algorithm>

// 3rd party
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "fmt/core.h"
#include "tbb/parallel_for.h"

// sniff
#include "sniff/minimize.h"

namespace sniff {

static constexpr auto kChunkSize = 1U << 26U;  // 64 MiB

static constexpr auto kFastaSuffxies =
    std::array<char const*, 4>{".fasta", "fasta.gz", ".fa", ".fa.gz"};

static constexpr auto kFastqSuffixes =
    std::array<char const*, 4>{".fastq", ".fastq.gz", ".fq", ".fq.gz"};

static auto IsSuffixFor(std::string_view const suffix,
                        std::string_view const query) -> bool {
  return suffix.length() <= query.length()
             ? suffix == query.substr(query.length() - suffix.length())
             : false;
}

static auto CreateParser(std::filesystem::path const& path)
    -> std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> {
  using namespace std::placeholders;
  if (std::filesystem::exists(path)) {
    if (std::any_of(kFastaSuffxies.cbegin(), kFastaSuffxies.cend(),
                    std::bind(IsSuffixFor, _1, path.c_str()))) {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastaParser>(path.c_str());
    } else if (std::any_of(kFastqSuffixes.cbegin(), kFastqSuffixes.cend(),
                           std::bind(IsSuffixFor, _1, path.c_str()))) {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<
          bioparser::FastqParser>(path.c_str());
    }
  }

  throw std::invalid_argument(
      "[camel::detail::CreateParser] invalid file path: " + path.string());
}

auto LoadReads(std::filesystem::path const& path)
    -> std::vector<std::unique_ptr<biosoup::NucleicAcid>> {
  auto timer = biosoup::Timer();

  timer.Start();
  auto parser = CreateParser(path);
  auto dst = std::vector<std::unique_ptr<biosoup::NucleicAcid>>();

  while (true) {
    auto reads = parser->Parse(kChunkSize);
    if (reads.empty()) {
      break;
    }
    dst.insert(dst.end(), std::make_move_iterator(reads.begin()),
               std::make_move_iterator(reads.end()));

    fmt::print(stderr,
               "\r[sniff::LoadSequences]({:12.3f}) loaded: {} sequences",
               timer.Lap(), dst.size());
  }

  std::sort(dst.begin(), dst.end(),
            [](std::unique_ptr<biosoup::NucleicAcid> const& lhs,
               std::unique_ptr<biosoup::NucleicAcid> const& rhs) -> bool {
              return lhs->id < rhs->id;
            });

  fmt::print(stderr,
             "\r[sniff::LoadSequences]({:12.3f}) loaded: {} sequences\n",
             timer.Stop(), dst.size());

  return dst;
}

}  // namespace sniff
