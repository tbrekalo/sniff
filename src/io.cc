#include "io.h"

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "fmt/core.h"

namespace sniff {

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
  auto dst = parser->Parse(std::numeric_limits<std::uint64_t>::max());
  fmt::print(stderr, "[sniff::LoadSequences]({:12.3f}) loaded: {} sequences\n",
             timer.Stop(), dst.size());

  return dst;
}

}  // namespace sniff
