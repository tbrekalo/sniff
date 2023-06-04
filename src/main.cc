#include <sys/resource.h>

#include <cstdlib>

// 3rd party dependencies
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "tbb/task_arena.h"
#include "version.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0;

// sniff
#include "sniff/algo.h"
#include "sniff/io.h"

static auto GetPeakMemoryUsageKB() -> std::uint32_t {
  struct rusage rusage_info;
  getrusage(RUSAGE_SELF, &rusage_info);

  return rusage_info.ru_maxrss;
}

int main(int argc, char** argv) {
  try {
    auto options =
        cxxopts::Options("sniff", "pair up potential reverse complement reads");
    /* clang-format off */
    options.add_options("general")
      ("h,help", "print help")
      ("v,version", "print version")
      ("t,threads", "number of threads to use",
        cxxopts::value<std::uint32_t>()->default_value("1"));
    options.add_options("heuristic")
      ("p,percent",
        "maximum allowed difference in length as % of shorter read's length",
        cxxopts::value<double>()->default_value("0.01"));
    options.add_options("mapping")
      ("k,kmer-length", "kmer length used in mapping",
        cxxopts::value<std::uint32_t>()->default_value("15"))
      ("w,window-length", "window length used in mapping",
        cxxopts::value<std::uint32_t>()->default_value("5"))
      ("c,chain", "minimum chain length (in kmers)",
        cxxopts::value<std::uint32_t>()->default_value("4"))
      ("g,gap", "maximum gap between minimizers when chaining",
        cxxopts::value<std::uint32_t>()->default_value("500"));
    options.add_options("input")
      ("input", "input fasta/fastq file", cxxopts::value<std::string>());
    /* clang-format on */

    options.parse_positional({"input"});
    auto result = options.parse(argc, argv);

    auto early_quit = false;
    if (result.count("version")) {
      fmt::print(stderr, "{}.{}.{}\n", sniff_VERSION_MAJOR, sniff_VERSION_MINOR,
                 sniff_VERSION_PATCH);
      early_quit = true;
    }

    if (result.count("help")) {
      fmt::print(stderr, "{}\n", options.help());
      early_quit = true;
    }

    if (early_quit) {
      return EXIT_SUCCESS;
    }

    auto const n_threads = result["threads"].as<std::uint32_t>();
    auto const reads_path =
        std::filesystem::path(result["input"].as<std::string>());

    auto task_arena = tbb::task_arena(n_threads);
    auto timer = biosoup::Timer();
    timer.Start();

    task_arena.execute([&] {
      auto const cfg = sniff::Config{
          .p = result["percent"].as<double>(),
          .map_cfg =
              sniff::MapConfig{
                  .min_chain_length = result["chain"].as<std::uint32_t>(),
                  .max_chain_gap_length = result["gap"].as<std::uint32_t>(),
                  .kmer_len = result["kmer-length"].as<std::uint32_t>()},
          .minimize_cfg = sniff::MinimizeConfig{
              .kmer_len = result["kmer-length"].as<std::uint32_t>(),
              .window_len = result["window-length"].as<std::uint32_t>()}};

      /* clang-format off */
        fmt::print(stderr,
          "[sniff]\n"
          "\tp: {};\n"
          "\tk: {}; w: {}; chain: {}; gap: {};\n",
          cfg.p,
          cfg.minimize_cfg.kmer_len, cfg.minimize_cfg.window_len,
          cfg.map_cfg.min_chain_length, cfg.map_cfg.max_chain_gap_length);
      /* clang-format on */

      auto pairs =
          sniff::FindReverseComplementPairs(cfg, sniff::LoadReads(reads_path));
      for (auto const& [lhs, rhs] : pairs) {
        fmt::print("{},{}\n", lhs, rhs);
      }
    });

    fmt::print(stderr, "[sniff::main]({:12.3f}) peak rss {:0.3f} GB\n",
               timer.Stop(), static_cast<double>(GetPeakMemoryUsageKB()) / 1e6);
  } catch (std::exception const& e) {
    fmt::print(stderr, "{}\n", e.what());
  }

  return EXIT_SUCCESS;
}
