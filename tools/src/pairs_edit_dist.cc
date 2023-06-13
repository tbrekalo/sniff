#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <span>
#include <sstream>
#include <string>

#include "ankerl/unordered_dense.h"
#include "bindings/cpp/WFAligner.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "sniff/io.h"
#include "tbb/parallel_for.h"
#include "tbb/task_arena.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0;

using ReadMap =
    ankerl::unordered_dense::map<std::string,
                                 std::unique_ptr<biosoup::NucleicAcid>>;

struct ReadPair {
  std::string lhs;
  std::string rhs;
};

struct ReadPairEditRatio {
  ReadPair pair;
  double ratio;
};

static auto CommaSplit(std::string const& src) -> std::vector<std::string> {
  auto dst = std::vector<std::string>();
  auto istrm = std::istringstream(src);
  for (auto token = std::string(); std::getline(istrm, token, ',');) {
    dst.push_back(token);
  }

  return dst;
}

static auto LoadReads(std::filesystem::path const& reads_path) -> ReadMap {
  auto dst = ReadMap();
  auto reads = sniff::LoadReads(reads_path);

  dst.reserve(reads.size());
  for (auto& it : reads) {
    dst[it->name] = std::move(it);
  }

  return dst;
}

static auto LoadPairs(std::filesystem::path const& pairs_path)
    -> std::vector<ReadPair> {
  auto dst = std::vector<ReadPair>();
  auto ifstrm = std::ifstream(pairs_path);
  for (auto line = std::string(); std::getline(ifstrm, line);) {
    auto pair = CommaSplit(line);
    dst.push_back(
        ReadPair{.lhs = std::move(pair[0]), .rhs = std::move(pair[1])});
  }

  return dst;
}

static auto CreateRcString(std::unique_ptr<biosoup::NucleicAcid> const& read)
    -> std::string {
  auto dst = std::string(read->InflateData());
  for (auto i = 0U; i < dst.size(); ++i) {
    dst[i] = biosoup::kNucleotideDecoder[3 ^ read->Code(dst.size() - 1 - i)];
  }

  return dst;
}

static auto CreatePairsWithEditRatio(ReadMap const& reads,
                                     std::span<ReadPair const> pairs)
    -> std::vector<ReadPairEditRatio> {
  auto dst = std::vector<ReadPairEditRatio>(pairs.size());
  tbb::parallel_for(std::size_t(0), pairs.size(), [&](std::size_t idx) -> void {
    thread_local auto init = false;
    thread_local wfa::WFAlignerGapLinear aligner(
        -1, 1, 2, wfa::WFAligner::Score, wfa::WFAligner::MemoryUltralow);
    if (!init) {
      aligner.setHeuristicWFadaptive(10, 50, 10);
      aligner.setHeuristicZDrop(100, 100);
      aligner.setHeuristicBandedAdaptive(50, 50, 1);
      init = true;
    }

    auto const [lhs_name, rhs_name] = pairs[idx];
    auto lhs_str = reads.at(lhs_name)->InflateData();
    auto rhs_str = CreateRcString(reads.at(rhs_name));

    aligner.alignEnd2End(lhs_str, rhs_str);

    dst[idx] = ReadPairEditRatio{
        .pair = pairs[idx],
        .ratio = static_cast<double>(aligner.getAlignmentScore()) /
                 std::max(lhs_str.size(), rhs_str.size())};
  });

  return dst;
}

int main(int argc, char** argv) {
  try {
    auto options = cxxopts::Options(
        "sniff", "add edit distance to reverse complement pairs");
    /* clang-format off */
    options.add_options("general")
      ("h,help", "print help")
      ("t,threads", "number of threads to use",
        cxxopts::value<std::uint32_t>()->default_value("1"));
    options.add_options("input")
      ("reads", "input fasta/fastq reads", cxxopts::value<std::string>())
      ("pairs", "csv file containing reverse complements",
        cxxopts::value<std::string>());
    /* clang-format on */

    options.positional_help("<reads> <pairs>");
    options.parse_positional({"reads", "pairs"});
    options.show_positional_help();

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
      fmt::print(stderr, "{}\n", options.help());
      return EXIT_SUCCESS;
    }

    auto const reads_path =
        std::filesystem::path(result["reads"].as<std::string>());
    if (!std::filesystem::exists(reads_path)) {
      throw std::invalid_argument("invalid path: " + reads_path.string());
    }

    auto const pairs_path =
        std::filesystem::path(result["pairs"].as<std::string>());
    if (!std::filesystem::exists(pairs_path)) {
      throw std::invalid_argument("invalid path: " + pairs_path.string());
    }

    auto ta = tbb::task_arena(result["threads"].as<std::uint32_t>());
    ta.execute([&]() -> void {
      auto const reads = LoadReads(reads_path);
      auto const pairs = LoadPairs(pairs_path);

      fmt::print(stderr, "loaded {} reads and {} pairs\n", reads.size(),
                 pairs.size());

      auto const pairs_with_ratio = CreatePairsWithEditRatio(reads, pairs);
      for (auto [pair, ratio] : pairs_with_ratio) {
        fmt::print("{},{},{}\n", pair.lhs, pair.rhs, ratio);
      }
    });

  } catch (std::exception const& e) {
    fmt::print("{}\n", e.what());
  }

  return EXIT_SUCCESS;
}
