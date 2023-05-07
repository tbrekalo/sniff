#include "sniff/io.h"

#include <iterator>

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

static constexpr auto kChunkSize = 1U << 30U;  // 1 GiB

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

auto LoadSketches(Config cfg, std::filesystem::path const& path)
    -> std::vector<Sketch> {
  auto timer = biosoup::Timer();

  timer.Start();
  auto parser = CreateParser(path);
  auto dst = std::vector<Sketch>();

  for (std::vector<Sketch> buff; true;) {
    auto reads = parser->Parse(kChunkSize);
    if (reads.empty()) {
      break;
    }

    buff.resize(reads.size());
    tbb::parallel_for(
        std::size_t(0), reads.size(),
        [&cfg, &reads, &buff](std::size_t idx) -> void {
          auto sample = reads[idx]->InflateData(0, cfg.sample_length);
          reads[idx]->ReverseAndComplement();
          auto rc_sample = reads[idx]->InflateData(0, cfg.sample_length);

          buff[idx] =
              Sketch{.read_identifier{.read_id = reads[idx]->id,
                                      .read_name = reads[idx]->name},

                     .read_len = reads[idx]->inflated_len,
                     .minimizers = Minimize(cfg.minimize_cfg, sample),
                     .rc_minimizers = Minimize(cfg.minimize_cfg, rc_sample)

              };
        });

    dst.insert(dst.end(), std::make_move_iterator(buff.begin()),
               std::make_move_iterator(buff.end()));

    fmt::print(stderr,
               "\r[sniff::LoadSequences]({:12.3f}) loaded: {} sequences",
               timer.Lap(), dst.size());
  }

  fmt::print(stderr,
             "\r[sniff::LoadSequences]({:12.3f}) loaded: {} sequences\n",
             timer.Stop(), dst.size());

  return dst;
}

}  // namespace sniff
